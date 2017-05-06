
import numpy as np
from datetime import datetime

# This library consists of graphical tools that may be useful to visualize results

def assess_performance(results,data):
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
    print "  Mean squared error of frequency guesses: {0}".format(np.mean([(f1-f2)**2 for f1,f2 in zip(actual_freqs, pred_freqs)]))
    
def graphical_comparison(results1,results2):
    pass