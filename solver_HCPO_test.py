import sys

from solver_HCPO import solve
from seq_data import SequencingData as SD

if len(sys.argv) < 2:
  print "Usage: python -i solver_HCPO_test.py <data_file>"
  print 
  exit()

def test_solver(data, **solver_kwargs):
  pairs = solve(data, **solver_kwargs)
  
  print "Solved with the following optional arguments:"
  for k,v in solver_kwargs.iteritems():
    print "  {0}: {1}".format(k,v)

  cells = data.metadata['cells']
  all_alphas, all_betas = zip(*cells)

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

data = SD(path=sys.argv[1])

pairs1 = test_solver(data)
#pairs2 = test_solver(data, wells_per_iter=50)
#pairs3 = test_solver(data, wells_per_iter=50, pair_threshold=.6)
#pairs4 = test_solver(data, wells_per_iter=90, pair_threshold=.1)
