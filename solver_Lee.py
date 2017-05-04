import random
import sys

import numpy as np
import scipy.optimize, scipy.misc

def extract_chains(seq_data):
  alphas_per_well, betas_per_well = zip(*seq_data.well_data)
  return sorted(set(sum(alphas_per_well, []))), sorted(set(sum(betas_per_well, [])))

def solve(seq_data, iters=100, pair_threshold = 0.9):
  ## Computes a solution to the alpha-beta pairing problem, using the methods in Lee et al. (2017)
  def compute_well_pairings(alpha_idx, beta_idx, scores):
    from hungarian import solve_general_assignment
    
    # Reformulate problem as a general assignment problem
    # Then apply Hungarian algorithm
    ratings = [[int(scores[i][j]*10000) for j in beta_idx] for i in alpha_idx]
    assignment = solve_general_assignment(ratings)

    # Transform back into alpha- and beta-chain pairings
    pairings = []
    for i, row in enumerate(assignment):
      for j, v in enumerate(row):
        if v:  pairings.append((alpha_idx[i], beta_idx[j]))

    return pairings
  
  # Extract all distinct alpha- and beta-chains observed
  # TODO: might be better to extract the chains directly from the cells in the system
  all_alphas, all_betas = extract_chains(seq_data)
  
  # Create dictionary to look up alpha/beta index in constant time
  alpha_to_idx = {a: i for i,a in enumerate(all_alphas)}
  beta_to_idx = {b: i for i,b in enumerate(all_betas)}
  # Transform all well data to reference alpha- and beta-chains by index
  well_data = [[[alpha_to_idx[a] for a in data[0]], [beta_to_idx[b] for b in data[1]]] for data in seq_data.well_data]

  # Calculate association scores
  S = [[0 for j in range(len(all_betas))] for i in range(len(all_alphas))]
  for well in well_data:
    well_alpha_idx, well_beta_idx = well
    increment = 1./len(well_alpha_idx) + 1./len(well_beta_idx)
    for a_idx in well_alpha_idx:
      for b_idx in well_beta_idx:
        S[a_idx][b_idx] += increment

  overall_pairing_counts = {}
  well_pairings = [None]*len(well_data)
  percent_done = 0.0
  for i in range(iters):
    # Choose random subset of wells for this iter
    # Loop is to ensure that subset size is greater than 0 (can happen w/ small well count)
    wells_idx = []
    while len(wells_idx)==0:  wells_idx = [i for i in range(len(well_data)) if random.random()>0.5]
    

    # Compute well pairings for any well, if it hasn't been done already
    # Then accumulate the number of times each pair has been assigned in a well pairing
    pairing_counts = {}
    for well_idx in wells_idx:
      if well_pairings[well_idx] == None:
        well_pairings[well_idx] = compute_well_pairings(*well_data[well_idx], scores=S)
        percent_done += 100./len(well_pairings)
        print "Computing likely pairings... {0}%\r".format(int(percent_done)),
      for a,b in well_pairings[well_idx]:
        pairing_counts[(a,b)] = pairing_counts.get((a,b), 0) + 1

    # Compute filter cutoff (average of all nonzero pair counts)
    cutoff = float(sum(pairing_counts.values()))/len(pairing_counts.keys())

    # Extract all pairs with counts exceeding the cutoff
    good_pairs = [pair for pair in pairing_counts if pairing_counts[pair]>=cutoff]

    # For each pair exceeding the cutoff, increment the overall_pairing_counts number
    for pair in good_pairs:
      overall_pairing_counts[pair] = overall_pairing_counts.get(pair, 0) + 1

  overall_good_pairs = [pair for pair in overall_pairing_counts if overall_pairing_counts[pair]>=pair_threshold*iters]

  cells = [(all_alphas[a_idx], all_betas[b_idx]) for a_idx,b_idx in overall_good_pairs]
  cell_freqs, cell_freq_CIs = estimate_cell_frequencies(seq_data, cells)

  results = {
    'cells': cells,
    'cell_frequencies': cell_freqs,
    'cell_frequencies_CI': cell_freq_CIs
  }

  return results

def estimate_cell_frequencies(seq_data, cells):
  def extract_cells_per_well():
    cpw_distro = seq_data.metadata['cells_per_well_distribution']
    cpw_params = seq_data.metadata['cells_per_well_distribution_params']
    if cpw_distro == 'constant':
      cells_per_well = [cpw_params['cells_per_well']]*len(seq_data.well_data)
    elif cpw_distro == 'poisson':
      cells_per_well = [cpw_params['lam']]*len(seq_data.well_data) # This is an approx. Not sure if it'll affect results
    elif cpw_distro == 'explicit':
      cells_per_well = cpw_params['cells_per_well']
    else:
      print "Unknown cell/well distribution: {0} with parameters {1}".format(cpw_distro, cpw_params)
      return None, None, None
  
    # Gather distinct #s of cells per well (N) and count number of wells w/ each cell count
    N_dict = {}
    for cpw in cells_per_well:
      N_dict[cpw] = N_dict.get(cpw, 0) + 1
    N,W = zip(*sorted(N_dict.iteritems()))
    return cells_per_well, N, W

  def extract_cell_counts(cells_per_well, N, W):
    K = [[0]*len(N) for i in range(len(cells))]

    for well_size, well_data in zip(cells_per_well, seq_data.well_data):
      well_alphas = set(well_data[0])
      well_betas = set(well_data[1])
      N_idx = N.index(well_size)
      for i,(a,b) in enumerate(cells):
        if a in well_alphas and b in well_betas:
          K[i][N_idx] += 1

    return K

  def log_likelihood_func(f, a_idx, b_idx, N, W, K, Q_memo = {}, error_rate=10**-1):
    # Note: See Eqs (3) and (4) in Lee et al. for explanation of variables

    # Compute Q vector if not previously computed
    Q_key = (tuple(N), tuple(W), error_rate, f)
    if Q_key not in Q_memo:
      Q = []
      for n,w in zip(N,W):
        q = (1-f)**n + sum([
            (2*error_rate**m - error_rate**(2*m)) * scipy.misc.comb(n,m) * (f**m) * (1-f)**(n-m) 
        for m in range(1, n+1)])
        Q.append(q)
      Q_memo[Q_key] = Q

    # Retrieve Q from memoized dict of Q's
    # Note that Q only depends on N, W, error_rate, and f
    Q = Q_memo[Q_key]

    # Compute log likelihood as sum of Binomial probabilities
    # Note that the "combinations" in the binomial PDF is ignored as it does not affect
    # the location of the maximum
    return sum([np.log(scipy.misc.comb(w,k)) + k*np.log((1-q)) + (w-k)*np.log(q) for w,k,q in zip(W,K,Q)])

  cells_per_well, N, W = extract_cells_per_well()

  K = extract_cell_counts(cells_per_well, N, W)

  cell_freqs = []
  cell_freq_CIs = []
  for (a_idx, b_idx), k in zip(cells, K):
    L_func = lambda f: log_likelihood_func(f, a_idx, b_idx, N, W, k)

    # Find maximal likelihood
    f_opt = scipy.optimize.minimize_scalar(lambda f: -L_func(f), method='Bounded', bounds=(0,1)).x
    L_max = L_func(f_opt)

    # Find confidence interval, as specified in the paper
    f_min = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(0,f_opt)).x
    f_max = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(f_opt,1)).x

    cell_freqs.append(f_opt)
    cell_freq_CIs.append((f_min, f_max))

    print "Computing chain pair frequencies... {0}%\r".format(int(100.*len(cell_freqs)/len(cells))),
  
  return cell_freqs, cell_freq_CIs

