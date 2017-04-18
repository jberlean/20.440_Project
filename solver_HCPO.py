import math, random
import numpy as np

def extract_chains(seq_data):
  alphas_per_well, betas_per_well = zip(*seq_data.well_data)
  return sorted(set(sum(alphas_per_well, []))), sorted(set(sum(betas_per_well, [])))

def solve(seq_data, log_epsilon_prior = None, log_m_prior = None):
  def compute_epsilon(N, N_a, N_ax, n):
    if N_ax==0:  return float('inf')
    if N_ax==N:  return 0.0
    f_w_ax = float(N_ax)/N
    f_w_a = float(N_a)/(N - N_ax)
    f_c_ax = 1-(1-f_w_ax)**(1./n)
    f_c_a = 1-(1-f_w_a)**(1./n)
    return f_c_a/f_c_ax
  def process_alpha_chain(a_idx, pairs, done_alphas, done_betas, well_data, cells_per_well):
    # Find existing pairs involving a_idx
    prev_pairs = filter(lambda p: p[0]==a_idx, pairs)
    
    m = len(prev_pairs)

    N = len(well_data)
    wells_a, wells_ax = [],[]
    for i,well in enumerate(well_data):
      if a_idx not in well[0]:
        continue
      if any([b in well[1] for _,b in prev_pairs]):
        wells_ax.append(i)
      else:
        wells_a.append(i)
    N_a = len(wells_a)
    N_ax = len(wells_ax)

    cur_eps = compute_epsilon(N, N_a, N_ax, cells_per_well)
    cur_likelihood = log_epsilon_prior(cur_eps) + log_m_prior(m)

    new_pairs = []

    while True:
      # Compute beta chain frequency amongs wells_a
      b_freqs = {}
      for well_idx in wells_a:
        for b in well_data[well_idx][1]:
          if b not in done_betas:
            b_freqs[b] = b_freqs.get(b, 0) + 1
  
      # Get most common beta chain
      next_b = -1
      b_count_max = -1
      for b,count in b_freqs.items():
        if count > b_count_max:
          next_b = b
          b_count_max = count
  
      # If no beta chain options remained, we are done
      if next_b==-1:
        return new_pairs

      # Compute new likelihood contribution with this new pair
      new_eps = compute_epsilon(N, N_a - b_count_max, N_ax + b_count_max, cells_per_well)
      new_likelihood = log_epsilon_prior(new_eps) + log_m_prior(m+1)
  
      # Check if it's time to go home
      print m, cur_eps, new_eps
      print cur_likelihood, new_likelihood
      if new_likelihood <= cur_likelihood:
        return new_pairs

      # Otherwise, update all fields for next iter
      new_pairs.append((a_idx, next_b))

      cur_eps = new_eps
      m += 1
      cur_likelihood = new_likelihood

      for well_idx in wells_a:
        if next_b in well_data[well_idx]:
          wells_ax += [well_idx]
          wells_a.remove(well_idx)
  def process_beta_chain(b_idx, pairs, done_alphas, done_betas, well_data, cells_per_well):
    # Lazy, so I'm gonna just switch all alpha/beta chain data and run process_alpha_chain()
    pairs_rev = [(b,a) for a,b in pairs]
    well_data_rev = [[b_list, a_list] for a_list,b_list in well_data]
    new_pairs_rev = process_alpha_chain(b_idx, pairs_rev, done_betas, done_alphas, well_data_rev, cells_per_well)
    return [(a,b) for b,a in new_pairs_rev]
    
    

  if log_epsilon_prior is None:
    log_epsilon_prior = lambda x: -x*10**2
  if log_m_prior is None: # Note, this is not normalized
    log_m_prior = lambda x: -100*math.log(x) if x>=1 else float('-inf')

  # Extract # cells per well
  if seq_data.metadata['cells_per_well_distribution']=='constant':
    cpw = seq_data.metadata['cells_per_well_distribution_params']['cells_per_well']
  elif seq_data.metadata['cells_per_well_distribution']=='poisson':
    cpw = seq_data.metadata['cells_per_well_distribution_params']['lam']
  elif seq_data.metadata['cells_per_well_distribution']=='explicit':
    cpw = np.mean(seq_data.metadata['cells_per_well_distribution_params']['cells_per_well'])

  # Extract all distinct alpha- and beta-chains observed
  all_alphas, all_betas = extract_chains(seq_data)

  # Create dictionary to look up alpha/beta index in constant time
  alpha_to_idx = {a: i for i,a in enumerate(all_alphas)}
  beta_to_idx = {b: i for i,b in enumerate(all_betas)}
  # Transform all well data to reference alpha/beta chains by index
  well_data = [[set([alpha_to_idx[a] for a in data[0]]), set([beta_to_idx[b] for b in data[1]])] for data in seq_data.well_data]

  # Compute the frequency with which each alpha/beta chain appears in a well
  alpha_counts, beta_counts = [0]*len(all_alphas), [0]*len(all_betas)
  for well_alphas, well_betas in well_data:
    for a in well_alphas:  alpha_counts[a] += 1
    for b in well_betas:  beta_counts[b] += 1

  # Sort alpha and beta chains by frequency of occurrence
  alphas_sorted = sorted(range(len(all_alphas)), key=lambda i: -alpha_counts[i])
  betas_sorted = sorted(range(len(all_betas)), key=lambda i: -beta_counts[i])

  # Process the chains in order
  all_pairs = set()
  done_alphas = set()
  done_betas = set()

  i,j = 0,0
  while i < len(alphas_sorted) or j < len(betas_sorted):
    do_alpha = j>=len(betas_sorted) or (i<len(alphas_sorted) and alpha_counts[alphas_sorted[i]]>beta_counts[betas_sorted[j]])
    if do_alpha:
      a_idx = alphas_sorted[i]
      new_pairs = process_alpha_chain(a_idx, all_pairs, done_alphas, done_betas, well_data, cpw)
      done_alphas.add(a_idx)
      print "New pairs for alpha chain {0}:".format(a_idx), [(all_alphas[a_idx],all_betas[b_idx]) for a_idx,b_idx in new_pairs]
      i += 1
    else:
      b_idx = betas_sorted[j]
      new_pairs = process_beta_chain(b_idx, all_pairs, done_alphas, done_betas, well_data, cpw)
      done_betas.add(b_idx)
      print "New pairs for beta chain {0}:".format(b_idx), [(all_alphas[a_idx],all_betas[b_idx]) for a_idx,b_idx in new_pairs]
      j += 1
    all_pairs |= set(new_pairs)

  return [(all_alphas[a_idx], all_betas[b_idx]) for a_idx,b_idx in all_pairs]
