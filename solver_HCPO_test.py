import sys

import numpy as np

from solver_HCPO import solve
from seq_data import SequencingData as SD
from seq_generator import SequencingGenerator as SG

def print_generator_args(gen):
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

def test_solver(data, **solver_kwargs):
  pairs = solve(data, **solver_kwargs)
  
  print "Solved with the following optional arguments:"
  for k,v in solver_kwargs.iteritems():
    print "  {0}: {1}".format(k,v)

  cells = data.metadata['cells']
  all_alphas, all_betas = zip(*cells)
  all_alphas,all_betas = set(all_alphas),set(all_betas)

  obs_alphas, obs_betas = zip(*data.well_data)
  obs_alphas, obs_betas = set(sum(obs_alphas, [])), set(sum(obs_betas, []))

  cells_set = set([(a,b) for a,b in cells])
  correct_pairs = [p for p in pairs if p in cells_set]
  incorrect_pairs = [p for p in pairs if p not in cells_set]

  fnr = 1.0 - float(len(correct_pairs))/len(cells)
  fdr = float(len(incorrect_pairs))/len(pairs)

  print "Solution statistics:"
  print "  Total cells (in system):", len(cells)
  print "  Number of alpha chains (in system):", len(all_alphas)
  print "  Number of beta chains (in system):", len(all_betas)
  print "  Number of alpha chains (observed):", len(obs_alphas)
  print "  Number of beta chains (observed):", len(obs_betas)
  print "  Total pairs identified:", len(pairs)
  print "  Correct pairs identified: {0} ({1}%)".format(len(correct_pairs), 100.*len(correct_pairs)/len(cells))
  print "  Incorrect pairs identified: {0}".format(len(incorrect_pairs))
  print "  False discovery rate: {0}%".format(100.*len(incorrect_pairs)/len(pairs))

  print ""

  return pairs, (fnr, fdr)

def get_stats(num_wells, cells_per_well):
  fnr_all = {}
  fdr_all = {}
  for nw, cpw in [(nw, cpw) for nw in num_wells for cpw in cells_per_well]:
    fnr = []
    fdr = []
    for i in range(10):
      #print "ITERATION", i
    
      gen = SG()
      gen.chain_deletion_prob=10**-1
      gen.num_wells = nw
      gen.set_cells_per_well('poisson', lam=cpw)
      cells = []
      a_idx, b_idx = 0,0
      while len(cells) < 1000:
        adeg = int(np.random.pareto(1)) + 1
        bdeg = int(np.random.pareto(1)) + 1
        cells += SG.generate_cells(int(np.random.pareto(adeg*bdeg)), adeg, bdeg, a_idx, b_idx)
        a_idx += bdeg
        b_idx += adeg
      gen.cells = cells[:1000]
      print_generator_args(gen)
      
      
      
      data = gen.generate_data()
      pairs, (fnr_temp, fdr_temp) = test_solver(data)
      fnr.append(fnr_temp)
      fdr.append(fdr_temp)
  
    print "# WELLS:", nw, "CPW:", cpw
    
    print "Mean fnr:", np.mean(fnr)
    print "Mean fdr:", np.mean(fdr)
    fnr_all[(nw, cpw)] = fnr
    fdr_all[(nw, cpw)] = fdr

  return fnr_all, fdr_all


num_wells = [1000]#[10, 50, 100, 200, 300, 500, 1000, 2000]
cells_per_well = [100]#[5,10,20, 50]#, 100]
fnr_all, fdr_all = get_stats(num_wells, cells_per_well)

import matplotlib.pyplot as plt

plt.ion()

plt.errorbar(num_wells, [np.mean(fnr_all[(nw, 100)]) for nw in num_wells], yerr = [np.std(fnr_all[(nw,100)]) for nw in num_wells], capsize=3, marker='*')
plt.errorbar(num_wells, [np.mean(fdr_all[(nw, 100)]) for nw in num_wells], yerr = [np.std(fdr_all[(nw,100)]) for nw in num_wells], capsize=3, marker='o')

