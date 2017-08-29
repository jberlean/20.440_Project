
'''
Project PHASOR

Solver for unique sequences

'''


'''
Library Importation
'''

# standard libraries
import math
import operator as op
import time
import functools

# nonstandard libraries
import numpy as np
import scipy.misc
#import matplotlib.pyplot as plt
#from seq_generator import SequencingGenerator as SeqGen
#from seq_data import SequencingData


'''
Factory Methods
'''

# N choose K (OUTDATED, slower by ~25x)
def nCk_old(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

# N choose K
def nCk(n, r):
    return scipy.misc.comb(n,r)
    
'''
Frequency Estimators
'''

# non-match MLE estimator for f_a,f_b,f_ab
def nonmatch_frequency(w_ab,w_a,w_b,w_tot):
    return float(w_ab)/w_tot,float(w_a+w_ab)/w_tot,float(w_b+w_ab)/w_tot

# MLE estimator for f_a,f_b,f_ab
def match_frequency(w_ab,w_a,w_b,w_tot):
    if w_tot-w_ab == 0: f_a = 0
    else: f_a = float(w_a)/(w_tot-w_ab)
    if w_tot-w_ab == 0: f_b = 0
    else: f_b = float(w_b)/(w_tot-w_ab)
    f_ab = max((0,1. - (1.-(float(w_ab)/w_tot))/(1-f_a*f_b)))
    return f_ab,f_a,f_b

'''
Probability Calculators
'''

# prior probability
def nonmatch_probability(w_ab,w_a,w_b,w_tot):
    w_ab,w_a,w_b,w_tot = int(w_ab),int(w_a),int(w_b),int(w_tot)
    f_ab,f_a,f_b = nonmatch_frequency(w_ab,w_a,w_b,w_tot)
    prob = instance_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b)
    if prob == 0.: 
        return float('-inf')
    return math.log10(prob)

def match_probability(w_ab,w_a,w_b,w_tot):
    if w_ab == 0: return float('nan')
    w_ab,w_a,w_b,w_tot = int(w_ab),int(w_a),int(w_b),int(w_tot)
    f_ab,f_a,f_b = match_frequency(w_ab,w_a,w_b,w_tot)
    prob_total = 0
    for w in xrange(0,int(w_ab)+1):
        prob_total += (nCk(w_tot,w)*(f_ab**(w))*((1-f_ab)**(w_tot-(w))))*instance_probability(w_ab-w,w_a,w_b,w_tot-w,f_ab,f_a,f_b)        
    if prob_total == 0.: return float('-inf')
    return math.log10(prob_total)

def instance_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b):
    a = nCk(w_tot,w_a+w_ab)*(f_a**(w_a+w_ab))*((1-f_a)**(w_tot-(w_a+w_ab)))
    b = nCk(w_a+w_ab,w_ab)*(f_b**w_ab)*((1-f_b)**(w_a+w_ab-w_ab))
    c =  nCk(w_tot-(w_a+w_ab),w_b)*(f_b**w_b)*((1-f_b)**(w_tot-(w_a+w_b+w_ab)))
    return a*b*c

'''
Specialty Calculators
'''
    
def self_occurence(v):
    temp = np.matmul(v,np.transpose(v))
    np.fill_diagonal(temp,0)
    return temp

def match_score(w_ab,w_a,w_b,w_tot, match_prior = 0.5):
    if w_ab == 0: 
        return 0.,0.
    else:
        mp,nmp = match_probability(w_ab,w_a,w_b,w_tot),nonmatch_probability(w_ab,w_a,w_b,w_tot)
        f_ab,_,_ = match_frequency(w_ab,w_a,w_b,w_tot)
        return (10**mp * match_prior)/(10**mp * match_prior + 10**nmp * (1-match_prior)),f_ab

def possible_matches(real_matches,img_ab,a_uniques,b_uniques):
    w_ab = np.sum(img_ab[:,:,:],axis=2)
    real = []
    for rm in real_matches:
        try: 
            i,j = a_uniques.index(rm[0]),b_uniques.index(rm[1])
            if w_ab[i,j] > 0: real.append(rm)
        except ValueError: # if a real match never even occurs
            pass
    return real

'''
Important Functions
'''

def find_clonotypes(start_cell, img_a, img_b, threshold = 0.99, start_prior = 1.0, alpha_priors = None, beta_priors = None, results_memo = None):
  def compute_well_membership(clonotype, well_idxs = None):
    if well_idxs is None:
      chain_membership = [img_a[a] for a in clonotype[0]] + [img_b[b] for b in clonotype[1]]
    else:
      chain_membership = [img_a[a, well_idxs] for a in clonotype[0]] + [img_b[b, well_idxs] for b in clonotype[1]]
    return functools.reduce(np.multiply, chain_membership)
  def compute_well_count(clonotype, well_idxs = None):
    return np.sum(compute_well_membership(clonotype, well_idxs))
  def compute_score(chain_set_1, chain_set_2, prior, well_idxs = None, w_1=None, w_2=None, w_12=None):
    combined = ((chain_set_1[0]+chain_set_2[0]), (chain_set_1[1]+chain_set_2[1]))
    w_tot = len(well_idxs) if well_idxs is not None else img_a.shape[1]
    if w_1 is None:  w_1 = compute_well_count(chain_set_1, well_idxs = well_idxs) 
    if w_2 is None:  w_2 = compute_well_count(chain_set_2, well_idxs = well_idxs)
    if w_12 is None:  w_12 = compute_well_count((chain_set_1[0]+chain_set_2[0], chain_set_1[1]+chain_set_2[1]), well_idxs = well_idxs)

    return match_score(w_12, w_1-w_12, w_2-w_12, w_tot, prior)
  def get_prior(chain_set_1, chain_additions):
    if len(chain_additions[0])==1 and len(chain_additions[1])==0:
      return alpha_priors[len(chain_set_1[0])]
    elif len(chain_additions[0])==0 and len(chain_additions[1])==1:
      return beta_priors[len(chain_set_1[1])]
    else:
      print "Cannot compute prior to add {} to {}. Using prior of 0.5".format(chain_additions, chain_set_1)
      return 0.5
  def augment_chain_set(chain_set, addition):
    return (tuple(sorted(chain_set[0]+addition[0])), tuple(sorted(chain_set[1]+addition[1])))

  if results_memo is not None and start_cell in results_memo:
    #print "  memoized results used for {}".format(start_cell)
    return results_memo[start_cell]

  num_alphas = img_a.shape[0]
  num_betas = img_b.shape[0]

  if alpha_priors is None:  alpha_priors = [1./num_alphas, 0.3/num_alphas]
  if beta_priors is None:  beta_priors = [1./num_betas, 0.06/num_betas]

  max_alphas, max_betas = len(alpha_priors), len(beta_priors)

  w_1 = compute_well_count(start_cell)

  alist, blist = start_cell

  cells = []
  freqs = []
  scores = []

  if len(alist) < max_alphas:
    ##alpha_additions = [((a_idx,),()) for a_idx in range(alist[-1]+1 if len(alist)>0 else 0, num_alphas)]
    alpha_additions = [((a_idx,),()) for a_idx in range(num_alphas) if a_idx not in alist]
  else:
    alpha_additions = []

  if len(blist) < max_betas:
    ##beta_additions = [((),(b_idx,)) for b_idx in range(blist[-1]+1 if len(blist)>0 else 0, num_betas)]
    beta_additions = [((),(b_idx,)) for b_idx in range(num_betas) if b_idx not in blist]
  else:
    beta_additions = []

  for addition in alpha_additions+beta_additions:
      prior = get_prior(start_cell, addition)
      score, freq = compute_score(start_cell, addition, prior, w_1 = w_1)
      score *= start_prior

      new_cell = augment_chain_set(start_cell, addition)
      if score > threshold:
        ## Try to augment new cell with additional chains
        aug_clonotypes, aug_freqs, aug_scores = find_clonotypes(new_cell, img_a, img_b, threshold, start_prior = score, results_memo=results_memo)

        ## Check if there is still significant correlation within remaining cells
        unexplained_well_idxs_bool = functools.reduce(np.multiply, [1-compute_well_membership(c) for c in aug_clonotypes], np.ones(img_a.shape[1])).astype(int)
        unexplained_well_idxs = [i for i in range(len(unexplained_well_idxs_bool)) if unexplained_well_idxs_bool[i]==1]
        score_residual, running_cs = 1, ((new_cell[0][0],), ())
        additions = [((a,),()) for a in new_cell[0][1:]] + [((),(b,)) for b in new_cell[1]]
        for add in additions:
          score_temp, freq_temp = compute_score(running_cs, add, get_prior(running_cs, add), well_idxs=unexplained_well_idxs)
          #if len(unexplained_well_idxs) < img_a.shape[1]:
            #print running_cs, add, get_prior(running_cs, add), score_temp
            #print len(unexplained_well_idxs), compute_well_count(running_cs, well_idxs=unexplained_well_idxs), compute_well_count(add, well_idxs=unexplained_well_idxs), compute_well_count(augment_chain_set(running_cs, add), well_idxs=unexplained_well_idxs)
          score_residual *= score_temp
          running_cs = augment_chain_set(running_cs, add)

        #if len(unexplained_well_idxs) < img_a.shape[1]:
        #  print score, score_residual, aug_clonotypes

        #if new_cell[0][0] == 8:
        #  print '******', score_residual
        if score_residual > threshold:
          #if len(unexplained_well_idxs) < img_a.shape[1]:
          #  print "****"
          aug_clonotypes.append(new_cell)
          aug_freqs.append(freq_temp)
          aug_scores.append(score_residual)


        ## Add new clonotypes to list, if they are new
        for c, f, s in zip(aug_clonotypes, aug_freqs, aug_scores):
          c = (tuple(sorted(c[0])), tuple(sorted(c[1])))
          if len(c[0]) > 0 and len(c[1]) > 0:
            cells.append(c)
            freqs.append(f)
            scores.append(s)

  if results_memo is not None:
    results_memo[start_cell] = (cells, freqs, scores)

  return cells, freqs, scores

  

'''

Name: directional_matches
Function: Returns edges connecting matches that occur in group a -> b

Parameters:
    img_ab -
        description: co-occurence matrix for unique members of a and b, in each well
        size: (# unique a) x (# unique b) x (# wells)
    img_a -
        description: occurence matrix for unique members of a, in each well
        size: (# unique a) x (# wells)
    img_b -
        description: occurence matrix for unique members of b, in each well
        size: (# unique b) x (# wells)
    a_uniques -
        description: labels for unique members of a
        size: (# unique a) x (1)
    b_uniques -
        description: labels for unique members of a
        size: (# unique a) x (1)       
    threshold -
        description: crieria for Bayesian probability of match
        size: float in range (0 to 1)
    silent -
        description: whether function prints matches as it finds them
        size: boolean (True/False)   
Returns:
    predicted_ab - 
        description: list of tuples (a label,b label)
        
'''
# TODO: reduce image dimensions to 2D
# TODO: use sparse matrices
# TODO: add stochasicity to well dismissal

def directional_matches(
  #img_ab,
  img_a,img_b,a_uniques,b_uniques,threshold=0.99,silent=False):
    def compute_well_count(alist, blist):
      w_ab = np.sum(functools.reduce([img_a[a] for a in alist] + [img_b[b] for b in blist], np.product))
      return w_ab

    score = np.zeros((len(a_uniques),len(b_uniques)))
    frequency = np.zeros((len(a_uniques),len(b_uniques)))
    
    # important storage variable
    predicted_ab = []
    predicted_frequency = []
    predicted_score = []
    
    w_tot = img_a.shape[1]

    results_memo = {}
    for i in xrange(img_a.shape[0]):
        if not silent: print 'Starting analysis for alpha chain {}...\n'.format(a_uniques[i])

        cells, freqs, scores = find_clonotypes(((i,), ()), img_a, img_b, threshold, results_memo = results_memo)
        cells_adj = [(tuple(sorted([a_uniques[a] for a in alist])), tuple(sorted([b_uniques[b] for b in blist]))) for alist,blist in cells]
        predicted_ab.extend(cells_adj)
        predicted_frequency.extend(freqs)
        predicted_score.extend(scores)

        print 'Edge detection progress... {}%\r'.format(100*(i+1)/img_a.shape[0]),
        
    print ''


    predicted_ab, predicted_frequency, predicted_score = zip(*list(set(zip(predicted_ab, predicted_frequency, predicted_score))))

    print len(predicted_ab), len(set(predicted_ab))

    return predicted_ab,predicted_frequency,predicted_score # returns edges


def reduce_graph(all_edges,all_freqs,all_scores,all_uniques):
    # implicitly assumes order of sets in: [ab,ba,aa,bb]


    # initialize critical variables
    edges,freqs,scores = [],[],[]

    # TODO: add aa/bb solver ( this is non-trivial math )
    #a_clones,b_uniques = [(i) for i in all_uniques[0]],[(i) for i in all_uniques[1]]

    results = dict()
    # results['cells'] = [e[0] for edge in all_edges for e in edge] # eventually, when we start doing multiclones
    results['cells'] = [((e[0],),(e[1],)) for e in all_edges[0]]
    results['cell_frequencies'] = [f for freq in all_freqs for f in freq]
    results['threshold'] = [s for score in all_scores for s in score]
    results['cell_frequencies_CI'] = [(-s/2,s/2) for score in all_scores for s in score] # this is pseudo, and totally not right

    return results

    # implicitly assumes order of sets in: [ab,ba,aa,bb]
    

'''
The Meat-and-Potatoes Function
'''
# Verbose: on range 0 to 9
# TODO: verbose levels
def solve(data,pair_threshold = 0.99,verbose=0):
    
    if verbose >= 5: silent = False
    else: silent = True
    
    # find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # sort in ascending order
    sorted(a_uniques,key=int)
    sorted(b_uniques,key=int)

    # counts of each unique in wells
    w_tot = len(data.well_data)

    a,b = len(a_uniques),len(b_uniques)
    
    #img_ab = np.zeros((a,b,w_tot))
    img_a = np.zeros((a,w_tot))
    img_b = np.zeros((b,w_tot))
    

    if verbose >= 1: print 'Starting image creation...'
    
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        a_ind = [a_uniques.index(i) for i in well[0]]
        b_ind = [b_uniques.index(i) for i in well[1]]
        a_v,b_v = np.zeros((len(a_uniques),1)),np.zeros((len(b_uniques),1)) 
        np.put(a_v,a_ind,np.ones((len(a_ind))))
        np.put(b_v,b_ind,np.ones((len(b_ind))))
        # assign values to image layers
        #img_ab[:,:,w] = np.matmul(a_v,np.transpose(b_v))
        img_a[:,w] = np.squeeze(a_v)
        img_b[:,w] = np.squeeze(b_v)
            
        print 'Image generation progress... {}%\r'.format(100*(w+1)/w_tot),
    print ''

    if verbose >= 1: print 'Starting edge detection...'
            
    # Setup threshold values (TODO: better way to distinguish ab,ba from aa,bb
    t = pair_threshold
    t_shared = 1 - (1 - pair_threshold)**2
            
    
    # Find each type of available edge
    cells,freqs,scores = directional_matches(
        #img_ab,
        img_a,img_b,a_uniques,b_uniques,threshold=t,silent=silent)
    if verbose >= 2: print 'Finished AB edges!'
        
        
    if verbose >= 1: print 'Finished edge detection, analyzing graph...'
        
#    real_matches = data.metadata['cells']
#    
#    # checks to see these actually occur in data
#    potential_ab_matches = possible_matches(real_matches,img_ab,a_uniques,b_uniques)
#    
#    # solves for true edges
#    all_edges = [ab_edges]
#    all_freqs = [ab_freqs]
#    all_scores = [ab_scores]
#    all_uniques = [a_uniques,b_uniques]
#    
    # graph reduction
    #results = reduce_graph(all_edges,all_freqs,all_scores,all_uniques)
    results = {
      'cells': cells,
      'cell_frequencies': freqs,
      'cell_frequencies_CI': [(s,s) for s in scores],
      'threshold': scores
    }
    
    if verbose >= 1: print 'Finished!'
    
    return results
    

    '''
    Scanning for alpha -> beta graph edges
    '''




