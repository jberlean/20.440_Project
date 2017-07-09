import sys
import itertools as it
import time

import random

import numpy as np
import scipy.stats as scipy_stats

from solver_Lee import solve as solve_Lee
from solver_440 import solve as solve_440
from seq_data import SequencingData as SD
from seq_generator import SequencingGenerator as SG

def roc(thresholds, correct, x_inc=None, y_inc=None):
  if x_inc is None: 1./(len(correct)-sum(correct))
  if y_inc is None: 1./sum(correct)

  thresholds = [th+np.random.normal(0, 1e-10) for th in thresholds]
  zipped = sorted(zip(thresholds, correct), reverse=True)

  x,y = [0],[0]
  for th, c in zipped:
    if c:
      x.append(x[-1])
      y.append(y[-1]+y_inc)
    else:
      x.append(x[-1]+x_inc)
      y.append(y[-1])

  x.append(1)
  y.append(1)
  auroc = sum([(x[i+1]-x[i])*(y[i]+y[i+1])/2. for i in range(len(x)-1)])

  return x,y,auroc

def coverage_vs_freq(data_list, res_list, num_bins = None, bins = None):
  freq_list = []
  coverage_list = []

  for data, res in zip(data_list, res_list):
    res_cells_set = set(res['cells'])
    freq_list.extend(data.metadata['generated_data']['cell_frequencies'])
    coverage_list.extend([c in res_cells_set for c in data.metadata['cells']])
  freq_list, coverage_list = zip(*sorted(zip(freq_list, coverage_list)))

  if num_bins is None:
    num_bins = 20
  if bins is None:
    log_min_freq, log_max_freq = np.log(freq_list[0]), np.log(freq_list[-1])
    bins = list(np.exp(np.arange(log_min_freq, log_max_freq, (log_max_freq-log_min_freq)/num_bins)[:(num_bins+1)]))
  bins.append(float('inf'))
  num_bins = len(bins)-1

  bin_counts, bin_successes = [0]*num_bins, [0]*num_bins
  bin_idx = 0
  for freq, covered in zip(freq_list, coverage_list):
    while freq >= bins[bin_idx+1]:
      bin_idx += 1 

    bin_counts[bin_idx] += 1
    if covered: bin_successes[bin_idx] += 1

  bin_mean = [float(succ)/(tot) if tot>0 else float('nan') for succ,tot in zip(bin_successes, bin_counts)]
  
  return bin_counts, bin_successes, bins[:-1]+[max(freq_list)]
      

def run_tests(tests, path):
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
  for testnum, (testname, num_reps, seq_gen_func, seq_gen_args, solver_list) in enumerate(tests):
    test_results = []
    print "*******"
    print "*******"
    print "STARTING TEST {} ({}/{})".format(testname, testnum, len(tests))
    print "*******"
    for i in range(num_reps):
      print "*******"
      print "TEST #{}, REP #{}".format(testnum, i)
      tstart = time.time()

      data = seq_gen_func(**seq_gen_args)
      rep_results = []
      for solver, args, stats_func in solver_list:
        res = solver(data, **args)
        stats = stats_func(data, res)
        rep_results.append(stats)
      test_results.append(rep_results)

      print "TEST #{}, REP #{} COMPLETED IN {} SECONDS".format(testnum, i, time.time()-tstart)
      print "*******"
    results.append(test_results)
    
    save_test_results(test_results, "{}_{}".format(testname, path))
  
  return results

def save_test_results(results, path):
  import json
  f = open(path, 'w')
  json.dump(results, f)
  f.close()
      

def run_Lee(data, **solver_kwargs):
  print "Running Lee et al. solver with the following optional arguments:"

  for k,v in solver_kwargs.iteritems():
    print "  {0}: {1}".format(k,v)

  return solve_Lee(data, **solver_kwargs)
def run_440(data, **solver_kwargs):
  print "Running 440 solver with the following optional arguments:"
  for k,v in solver_kwargs.iteritems():
    print "  {}: {}".format(k,v)

  return solve_440(data, **solver_kwargs)

def stats_Lee(data, results):
  def separate_top_tail_clones(fraction = 0.5):
    sorted_freqs, sorted_clones = zip(*sorted(zip(data.metadata['generated_data']['cell_frequencies'], data.metadata['cells']), reverse=True))
    idx_cutoff = len(filter(lambda f: f<fraction, np.cumsum(sorted_freqs)))
    top = sorted_clones[:idx_cutoff]
    tail = sorted_clones[idx_cutoff:]
    return top, tail

  def clones_to_pairings(clones):
    return set([p for alist,blist in clones for p in it.combinations([((a,),()) for a in alist]+[((),(b,)) for b in blist], 2)])
    

  true_clones = data.metadata['cells']
  true_pairings = clones_to_pairings(true_clones)
  true_top_clones, true_tail_clones = separate_top_tail_clones()

  true_dual_clones = [c for c in true_clones if len(c[0])>1 or len(c[1])>1]
  true_top_dual_clones = [c for c in true_top_clones if len(c[0])>1 or len(c[1])>1]
  true_tail_dual_clones = [c for c in true_tail_clones if len(c[0])>1 or len(c[1])>1]

  all_alphas, all_betas = zip(*true_clones)
  all_alphas,all_betas = set([a for alist in all_alphas for a in alist]), set([b for blist in all_betas for b in blist])
  obs_alphas, obs_betas = zip(*data.well_data)
  obs_alphas, obs_betas = set(sum(obs_alphas, [])), set(sum(obs_betas, []))

  clones = results['cells']
  dual_clones = [c for c in clones if len(c[0])>1 or len(c[1])>1]
  pairings = clones_to_pairings(clones)

  correct_clones = [c for c in clones if c in true_clones]
  correct_top_clones = [c for c in clones if c in true_top_clones]
  correct_tail_clones = [c for c in clones if c in true_tail_clones]

  candidate_dual_clones = [c for c in dual_clones if all([p in pairings for p in clones_to_pairings([c])]) and c in true_dual_clones]
  candidate_top_dual_clones = [c for c in candidate_dual_clones if c in true_top_dual_clones]
  candidate_tail_dual_clones = [c for c in candidate_dual_clones if c in true_tail_dual_clones]
  correct_dual_clones = [c for c in dual_clones if c in true_dual_clones]
  correct_top_dual_clones = [c for c in dual_clones if c in correct_top_clones]
  correct_tail_dual_clones = [c for c in dual_clones if c in correct_tail_clones]

  correct_pairings = [p for p in pairings if p in true_pairings]
  incorrect_pairings = [p for p in pairings if p not in true_pairings]

  clone_freqs = results['cell_frequencies']
  clone_freqs_CI = results['cell_frequencies_CI']
  clone_true_freqs = [data.metadata['generated_data']['cell_frequencies'][true_clones.index(c)] if c in true_clones else 0.0 for c in clones]

  #print "", correct_pairings
  #print "", [c for c in clones if c not in correct_clones]

  stats = {
    'num_cells': len(true_clones),
    'num_dual_alpha_cells': len(filter(lambda c: len(c[0])>1, true_clones)),
    'num_dual_beta_cells': len(filter(lambda c: len(c[1])>1, true_clones)),
    'num_alphas': len(all_alphas),
    'num_betas': len(all_betas),
    'num_alphas_obs': len(obs_alphas),
    'num_betas_obs': len(obs_betas),
    'num_top_clones': len(true_top_clones),
    'num_tail_clones': len(true_tail_clones),
    'num_top_dual_clones': len(true_top_dual_clones),
    'num_tail_dual_clones': len(true_tail_dual_clones),
    'num_pred_pairs': len(pairings),
    'num_pred_pairs_correct': len(correct_pairings),
    'num_pred_pairs_incorrect': len(incorrect_pairings),
    'num_pred_clones': len(clones),
    'num_pred_clones_correct': len(correct_clones),
    'num_pred_clones_incorrect': len(clones)-len(correct_clones),
    'num_pred_dual_clones': len(dual_clones),
    'num_pred_dual_clones_correct': len(correct_dual_clones),
    'num_pred_dual_clones_incorrect': len(dual_clones)-len(correct_dual_clones),
    'false_negative_pairs': 1. - float(len(correct_pairings))/len(true_pairings) if len(true_pairings)>0 else float('nan'),
    'false_negative_clones': 1. - float(len(correct_clones))/len(true_clones) if len(true_clones)>0 else float('nan'),
    'false_negative_top_clones': 1. - float(len(correct_top_clones))/len(true_top_clones) if len(true_top_clones)>0 else float('nan'),
    'false_negative_tail_clones': 1. - float(len(correct_tail_clones))/len(true_tail_clones) if len(true_tail_clones)>0 else float('nan'),
    'false_negative_dual_clones': 1. - float(len(correct_dual_clones))/len(true_dual_clones) if len(true_dual_clones)>0 else float('nan'),
    'false_negative_top_dual_clones': 1. - float(len(correct_top_dual_clones))/len(true_top_dual_clones) if len(true_top_dual_clones)>0 else float('nan'),
    'false_negative_tail_dual_clones': 1. - float(len(correct_tail_dual_clones))/len(true_tail_dual_clones) if len(true_tail_dual_clones)>0 else float('nan'),
    'false_negative_dual_clones_adj': 1. - float(len(correct_dual_clones))/len(candidate_dual_clones) if len(candidate_dual_clones)>0 else float('nan'),
    'false_negative_top_dual_clones_adj': 1. - float(len(correct_top_dual_clones))/len(candidate_top_dual_clones) if len(candidate_top_dual_clones)>0 else float('nan'),
    'false_negative_tail_dual_clones_adj': 1. - float(len(correct_tail_dual_clones))/len(candidate_tail_dual_clones) if len(candidate_tail_dual_clones)>0 else float('nan'),
    'false_discovery_pairs': float(len(incorrect_pairings))/len(pairings) if len(pairings)>0 else float('nan'),
    'false_discovery_clones': 1. - float(len(correct_clones))/len(clones) if len(clones)>0 else float('nan'),
    'false_discovery_dual_clones': 1. - float(len(correct_dual_clones))/len(dual_clones) if len(dual_clones)>0 else float('nan'),
    'freq_mse': np.mean([(f1-f2)**2 for f1,f2 in zip(clone_true_freqs, clone_freqs)]),
    'freq_ci_accuracy': np.mean([(f>=f_min and f<=f_max) for f,(f_min,f_max) in zip(clone_true_freqs, clone_freqs_CI)]),
  }

  print "Solution statistics:"
  print "  Total cells (in system):", stats['num_cells']
  print "  # clones in top 50%:", stats['num_top_clones']
  print "  # clones in bottom 50% (tail):", stats['num_tail_clones']
  print "  Total number of dual-alpha cells:", stats['num_dual_alpha_cells']
  print "  Total number of dual-beta cells:", stats['num_dual_beta_cells']
  print "  Number of alpha chains (in system):", stats['num_alphas']
  print "  Number of beta chains (in system):", stats['num_betas']
  print "  Number of alpha chains (observed):", stats['num_alphas_obs']
  print "  Number of beta chains (observed):", stats['num_betas_obs']
  print "  Total pairs identified:", stats['num_pred_pairs']
  print "  Correct pairs identified: {0} ({1}%)".format(stats['num_pred_pairs_correct'], 100*(1 - stats['false_negative_pairs']))
  print "  Incorrect pairs identified: {0}".format(stats['num_pred_pairs_incorrect'])
  print "  False discovery rate (over pairings): {0}%".format(100*stats['false_discovery_pairs'])
  print "  Total clones identified:", stats['num_pred_clones']
  print "  Correct clones identified: {} ({}%)".format(stats['num_pred_clones_correct'], 100*(1-stats['false_negative_clones']))
  print "  Incorrect clones identified:", stats['num_pred_clones_incorrect']
  print "  False discovery rate (over clonotypes): {0}%".format(100*stats['false_discovery_clones'])
  print "  Total dual clones identified:", stats['num_pred_dual_clones']
  print "  Correct dual clones identified: {} ({}%)".format(stats['num_pred_dual_clones_correct'], 100*(1-stats['false_negative_dual_clones']))
  print "  Incorrect dual clones identified:", stats['num_pred_dual_clones_incorrect']
  print "  False discovery rate (over dual-chain clonotypes): {0}%".format(100*stats['false_discovery_dual_clones'])
  print "  Overall depth (over pairs):", 100.*(1-stats['false_negative_pairs'])
  print "  Overall depth (over clones):", 100.*(1-stats['false_negative_clones'])
  print "  Depth of top clones:", 100.*(1-stats['false_negative_top_clones'])
  print "  Depth of tail clones:", 100.*(1-stats['false_negative_tail_clones'])
  print "  Depth of dual clones:", 100.*(1-stats['false_negative_dual_clones'])
  print "  Depth of top dual clones:", 100.*(1-stats['false_negative_top_dual_clones'])
  print "  Depth of tail dual clones:", 100.*(1-stats['false_negative_tail_dual_clones'])
  print "  Depth of dual clones (adj.):", 100.*(1-stats['false_negative_dual_clones_adj'])
  print "  Depth of top dual clones (adj.):", 100.*(1-stats['false_negative_top_dual_clones_adj'])
  print "  Depth of tail dual clones (adj.):", 100.*(1-stats['false_negative_tail_dual_clones_adj'])
  print "  False dual rate:", 100.*stats['false_discovery_dual_clones']

  print "  Mean squared error of frequency guesses: {0}".format(stats['freq_mse'])
  print "  Percent of frequencies within confidence interval: {0}%".format(100.*stats['freq_ci_accuracy'])

  print

  return stats
def stats_440(data, results):
  #results['cells'] = [((a,),(b,)) for a,b in results['cells']]
  return stats_Lee(data, results)



def generate_sequencing_data(num_cells, **seq_gen_args):
  gen = SG(**seq_gen_args)
  #gen.cells = SG.generate_cells(num_cells, alpha_dual_prob=0.3, beta_dual_prob=0.06)
  gen.cells = SG.generate_cells(num_cells, alpha_dual_prob=0.0, beta_dual_prob=0.00)
  #gen.set_cell_frequency_distribution(distro_type='explicit', frequencies=generate_cell_freqs(len(gen.cells),50))

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
  ("W={},CPW={}".format(w, cpw),
   5, 
   generate_sequencing_data,
   {'num_cells': 3000,
    'chain_deletion_prob': 0.15,
    'num_wells': w,
    'cells_per_well_distribution': 'constant',
    'cells_per_well_distribution_params': {'cells_per_well': cpw},
    'cell_frequency_distribution': 'power-law',
    'cell_frequency_distribution_params': {'alpha': -1}
   },
#   ((run_Lee, {'pair_threshold': 0.3, 'iters':25}, stats_Lee),
#    (run_Lee, {'pair_threshold': 0.6}, stats_Lee),
#    (run_Lee, {'pair_threshold': 0.3}, stats_Lee),
    ((run_440, {'pair_threshold': 0.75}, stats_440),)
  ) for w,cpw in  zip([96],[100]*5)+zip([500]*5,[5,10,20,50,100])#,1000,2000]#[5,10,50,100,300,500,1000,2000]
]
        

tstart = time.time()
#results = run_tests(tests, "TESTRESULTS_{}.txt".format(random.randint(100,999)))
data=tests[0][2](**tests[0][3])
#data = SD(path='test_dual2.txt')
res=tests[0][4][0][0](data, **tests[0][4][0][1])
tests[0][4][0][2](data, res)
print "*** Finished in {} seconds".format(time.time() - tstart)

obs_alphas = set([a for alist,_ in data.well_data for a in alist])
obs_betas = set([b for _,blist in data.well_data for b in blist])
num_correct = len(filter(lambda c: all([a in obs_alphas for a in c[0]]) and all([b in obs_betas for b in c[1]]), data.metadata['cells']))
num_incorrect = len(obs_alphas)*len(obs_betas)-num_correct

res2 = run_Lee(data,**tests[0][4][0][1])
tests[0][4][0][2](data,res2)

import matplotlib.pyplot as plt
plt.ion()

def compute_plot_coverage_vs_freq(data, res):
  bin_counts, bin_successes, bins = coverage_vs_freq([data],[res])
  bin_means = [float(s)/c if c>0 else float('nan') for s,c in zip(bin_successes,bin_counts)]
  bin_width = [bins[i+1]-bins[i] for i in range(len(bins)-1)]
  bin_left = [b+w/2. for b,w in zip(bins, bin_width)]
  return {'left': bin_left, 'height': bin_means, 'width': bin_width, 'color': (0,0,0,0), 'linewidth': 2}

plt.figure()
plt.xscale('log')
plt.bar(edgecolor='blue', **compute_plot_coverage_vs_freq(data, res))
plt.bar(edgecolor='black', **compute_plot_coverage_vs_freq(data,res2))

