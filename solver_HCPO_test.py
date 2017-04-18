import sys

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

  return pairs

gen = SG(chain_deletion_prob=10**-1)
gen.num_wells = 1000
gen.set_cells_per_well('constant', cells_per_well=100)
gen.cells = SG.generate_cells(900, 1, 1) + SG.generate_cells(100, 2, 1, 900, 900)
data = gen.generate_data()
print_generator_args(gen)


pairs = test_solver(data)
