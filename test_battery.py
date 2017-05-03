import sys

import numpy as np

from solver_Lee import solve
from seq_data import SequencingData as SD
from seq_generator import SequencingGenerator as SG

def run_tests(tests):
  ## tests is a list describing a set of runs to perform
  ## Each element of tests is of the following form:
  ##  (<num_reps>,
  ##   <seq_generator_func>,
  ##   <seq_generator_func_args>,
  ##   ((<solver1>, <solver1_args>, <solver1_stats_func>),
  ##    (<solver2>, <solver2_args>, <solver2_stats_func>),
  ##    ...)
  ##  )
  results = []
  for num_reps, seq_gen_func, seq_gen_args, solver_list in tests:
    test_results = []
    for i in range(num_reps):
      data = seq_gen_func(**seq_gen_args)
      rep_results = []
      for solver, args, stats_func in solver_list:
        res = solver(data, **args)
        stats = stats_func(data, res)
        rep_results.append(stats)
      test_results.append(rep_results)
    results.append(test_results)

  return results

def save_test_results(results, path):
  import json
  f = open(path, 'w')
  json.dump(results, f)
      

def run_Lee(data, **solver_kwargs):
  results = solve(data, **solver_kwargs)
  
  print "Ran Lee et al. solver with the following optional arguments:"
  for k,v in solver_kwargs.iteritems():
    print "  {0}: {1}".format(k,v)

  return results
def stats_Lee(data, results):
  pairs = results['cells']

  cells = data.metadata['cells']
  all_alphas, all_betas = zip(*cells)
  all_alphas,all_betas = set(all_alphas),set(all_betas)

  obs_alphas, obs_betas = zip(*data.well_data)
  obs_alphas, obs_betas = set(sum(obs_alphas, [])), set(sum(obs_betas, []))

  cells_set = set([(a,b) for a,b in cells])
  correct_pairs = [p for p in pairs if p in cells_set]
  incorrect_pairs = [p for p in pairs if p not in cells_set]

  pair_idxs = [cells.index(p) if p in cells else -1 for p in pairs]
  actual_freqs = [data.metadata['generated_data']['cell_frequencies'][i] if i!=-1 else 0.0 for i in pair_idxs]
  pred_freqs = results['cell_frequencies']
  pred_freqs_CI = results['cell_frequencies_CI']

  sorted_clones = sorted(zip(data.metadata['generated_data']['cell_frequencies'], data.metadata['cells']), reverse=True)
  p_s = 0
  for num_top_clones in range(len(sorted_clones)):
    p_s += sorted_clones[num_top_clones][0]
    num_top_clones += 1
    if p_s >= 0.5:  break
  _,top_clones = zip(*sorted_clones[:num_top_clones])
  _,tail_clones = zip(*sorted_clones[num_top_clones:])
  top_clones = set(top_clones)
  tail_clones = set(tail_clones)

  stats = {
    'num_cells': len(cells),
    'num_alphas': len(all_alphas),
    'num_betas': len(all_betas),
    'num_alphas_obs': len(obs_alphas),
    'num_betas_obs': len(obs_betas),
    'num_pairs': len(pairs),
    'num_pairs_correct': len(correct_pairs),
    'num_pairs_incorrect': len(incorrect_pairs),
    'false_negative': 1. - float(len(correct_pairs))/len(cells),
    'false_discovery': float(len(incorrect_pairs))/len(pairs),
    'freq_mse': np.mean([(f1-f2)**2 for f1,f2 in zip(actual_freqs, pred_freqs)]),
    'freq_ci_accuracy': np.mean([(f>=f_min and f<=f_max) for f,(f_min,f_max) in zip(actual_freqs, pred_freqs_CI)]),
    'depth_top': float(len([p for p in pairs if p in top_clones]))/len(top_clones),
    'depth_tail': float(len([p for p in pairs if p in tail_clones]))/len(tail_clones)
  }

  print "Solution statistics:"
  print "  Total cells (in system):", stats['num_cells']
  print "  Number of alpha chains (in system):", stats['num_alphas']
  print "  Number of beta chains (in system):", stats['num_betas']
  print "  Number of alpha chains (observed):", stats['num_alphas_obs']
  print "  Number of beta chains (observed):", stats['num_betas_obs']
  print "  Total pairs identified:", stats['num_pairs']
  print "  Correct pairs identified: {0} ({1}%)".format(stats['num_pairs_correct'], 100*(1 - stats['false_negative']))
  print "  Incorrect pairs identified: {0}".format(stats['num_pairs_incorrect'])
  print "  False discovery rate: {0}%".format(100*stats['false_discovery'])
  print "  Mean squared error of frequency guesses: {0}".format(stats['freq_mse'])
  print "  Percent of frequencies within confidence interval: {0}%".format(100.*stats['freq_ci_accuracy'])

  print "  Depth of top clones:", 100.*stats['depth_top']
  print "  Depth of tail:", 100.*stats['depth_tail']

  print

  return stats


def generate_cells(num_cells, max_alphas=None, max_betas=None):
  if max_alphas == None:  max_alphas = num_cells
  if max_betas == None:  max_betas = num_cells

  # Generate the degree for each alpha- and beta-chain from a given distribution
  a_sharing_probs=[0.816,0.085,0.021,0.007,0.033,0.005,0.033]
  b_sharing_probs=[0.859,0.076,0.037,0.019,0.009]
  #[0.8375, 0.0805, 0.029, 0.013, 0.021, 0.0025, 0.0165] # Averages from the Lee et al. paper
  adegs = np.random.choice(range(1,len(a_sharing_probs)+1), max_alphas, replace=True, p=a_sharing_probs)
  bdegs = np.random.choice(range(1,len(b_sharing_probs)+1), max_betas, replace=True, p=b_sharing_probs)
  
  # If you want to generate from a power law instead (not sure if this works as expected)
  #adegs = np.floor(np.random.pareto(2.1, size=max_alphas))+1
  #bdegs = np.floor(np.random.pareto(2.1, size=max_alphas))+1

  # Cut off at the desired number of cells
  alphas = sum([[i]*int(n) for i,n in enumerate(adegs)],[])[:num_cells] # this truncation will skew the distro a bit
  betas = sum([[i]*int(n) for i,n in enumerate(bdegs)], [])[:num_cells]

  # Randomly assign alpha- and beta-chains to each other
  np.random.shuffle(alphas)
  np.random.shuffle(betas)
  cells = list(set(zip(alphas, betas))) # Due to chance duplicates, there may be slightly less than num_cells cells
  # (A slightly more complex method could be used to ensure exactly num_cells cells)

  return cells

def generate_sequencing_data(num_cells, **seq_gen_args):
  gen = SG(**seq_gen_args)
  gen.cells = generate_cells(num_cells)
  #gen.set_cell_frequency_distribution(distro_type='explicit', frequencies=generate_cell_freqs(len(gen.cells),50))
  # TODO: Migrate gneerate_cell_freqs to seq_generator.py as "Lee" cell frequency distribution

  print "Generated data with the following parameters:"
  print "  Number of wells: {0}".format(gen.num_wells)
  print "  Cells/well distribution: {0} ({1})".format(gen.cells_per_well_distribution, gen.cells_per_well_distribution_params)
  
  alphas, betas = zip(*gen.cells)
  alphas,betas = set(alphas),set(betas)
  adeg = float(len(gen.cells))/len(alphas)
  bdeg = float(len(gen.cells))/len(betas)

  print "  Number of cells: {0}".format(len(gen.cells))
  print "  Number of unique alpha chains: {0}".format(len(alphas))
  print "  Number of unique beta chains: {0}".format(len(betas))
  print "  Average sharing of alpha chains: {0}".format(adeg)
  print "  Average sharing of beta chains: {0}".format(bdeg)

  print "  Chain deletion probability: {0}".format(gen.chain_deletion_prob)
  print "  Chain misplacement probability: {0}".format(gen.chain_misplacement_prob)
  print

  return gen.generate_data()
  

tests = [
  (5, 
   generate_sequencing_data,
   {'num_cells': 300,
    'chain_deletion_prob': 0.15,
    'num_wells': 96,
    'cells_per_well_distribution': 'constant',
    'cells_per_well_distribution_params': {'cells_per_well': 10},
    'cell_frequency_distribution': 'Lee',
    'cell_frequency_distribution_params': {'n_s': 50}
   },
   ((run_Lee, {'pair_threshold': 0.9}, stats_Lee),
    (run_Lee, {'pair_threshold': 0.6}, stats_Lee),
    (run_Lee, {'pair_threshold': 0.3}, stats_Lee))
  )
]
        
results = run_tests(tests)
