import random
import sys

def extract_chains(seq_data):
  alphas_per_well, betas_per_well = zip(*seq_data.well_data)
  return sorted(set(sum(alphas_per_well, []))), sorted(set(sum(betas_per_well, [])))

def solve(seq_data, iters=100, wells_per_iter=100, pair_threshold = 0.9):
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
  wells_per_iter = min(wells_per_iter, len(well_data))
  percent_done = 0.0
  for i in range(iters):
    # Choose random subset of wells for this iter
    wells_idx = [i for i in range(len(well_data)) if random.random()>0.5]

    # Compute well pairings for any well, if it hasn't been done already
    # Then accumulate the number of times each pair has been assigned in a well pairing
    pairing_counts = {}
    for well_idx in wells_idx:
      if well_pairings[well_idx] == None:
        well_pairings[well_idx] = compute_well_pairings(*well_data[well_idx], scores=S)
        percent_done += 100./len(well_pairings)
        print "{0}%\r".format(int(percent_done)),
        sys.stdout.flush()
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

  return [(all_alphas[a_idx], all_betas[b_idx]) for a_idx,b_idx in overall_good_pairs]
