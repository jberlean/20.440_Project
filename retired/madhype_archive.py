'''
THIS IS BEFORE REMOVING THE IMAGE METHOD FROM MAD-HYPE
'''

'''

Project: 
MAD-HYPE (Multicell Analytical Deconvolution for High Yield Paired-chain Evaluation)

Class(s): 
(none)

Function: 
Using analytics based in Bayesian systems, we attempt to identify T-cell clones from low-throughput well sequencing, with each well containing high-throughput numbers of cells

'''

'''
Library Importation
'''

# standard libraries
import math
import operator as op
import time

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

def directional_matches(img_ab,img_a,img_b,a_uniques,b_uniques,threshold=0.99,silent=False,data=None):
    score = np.zeros((len(a_uniques),len(b_uniques)))
    frequency = np.zeros((len(a_uniques),len(b_uniques)))
    
    # important storage variables
    predicted_ab = []
    predicted_frequency = []
    predicted_score = []
    
    w_tot = img_ab.shape[2]

    for i in xrange(img_ab.shape[0]):

        if not silent: print 'Starting analysis for alpha chain {}...\n'.format(a_uniques[i])

        unexplained_wells = list(xrange(w_tot))

        # create counts based on unexplained well data
        w_ab = np.sum(img_ab[i,:,unexplained_wells],axis=0)
        w_a = np.sum(img_a[i,unexplained_wells],axis=0)
        w_b = np.sum(img_b[:,unexplained_wells],axis=1)

        # assign scores
        for j in xrange(img_ab.shape[1]):
            n_ab = w_ab[j]
            n_a,n_b = w_a - n_ab,w_b[j] - n_ab
            n_tot = len(unexplained_wells)
            score[i,j],frequency[i,j] = match_score(n_ab,n_a,n_b,n_tot, 1./np.sqrt(img_ab.shape[0]*img_ab.shape[1]))
            #score[i,j],frequency[i,j] = match_score(n_ab,n_a,n_b,n_tot, 0.5) # effectively take out prior

            if score[i,j] > threshold:
              predicted_ab.append((a_uniques[i],b_uniques[j]))
              predicted_frequency.append(frequency[i,j])
              predicted_score.append(score[i,j])
              
              if data:
                  print 'Found match:'
                  print ' > Alpha - {}'.format(data.seqs_alpha[i])
                  print ' > Beta - {}'.format(data.seqs_beta[j])
                  print 'Number of wells with both chains: {}'.format(n_ab)
                  print 'Number of wells with alpha chain: {}'.format(n_a)
                  print 'Number of wells with beta chain: {}'.format(n_b)
                  print 'Number of wells total: {}'.format(n_tot)
                  print 'Confidence: {}'.format(score[i,j])
                  print 'Projected frequency: {}\n'.format(frequency[i,j])
              
        print 'Edge detection progress... {}%\r'.format(100*(i+1)/img_ab.shape[0]),
        
    print ''
    
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
def solve(data,pair_threshold = 0.99,verbose=0,real_data=False):
    
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
    
    img_ab,img_ba = np.zeros((a,b,w_tot)),np.zeros((b,a,w_tot))
    img_aa,img_bb = np.zeros((a,a,w_tot)),np.zeros((b,b,w_tot))
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
        img_ab[:,:,w] = np.matmul(a_v,np.transpose(b_v))
        img_ba[:,:,w] = np.matmul(b_v,np.transpose(a_v))
        img_aa[:,:,w] = self_occurence(a_v)
        img_bb[:,:,w] = self_occurence(b_v)
        img_a[:,w] = np.squeeze(a_v)
        img_b[:,w] = np.squeeze(b_v)
            
        print 'Image generation progress... {}%\r'.format(100*(w+1)/w_tot),
    print ''

    if verbose >= 1: print 'Starting edge detection...'
            
    # Setup threshold values (TODO: better way to distinguish ab,ba from aa,bb
    t = pair_threshold
    t_shared = 1 - (1 - pair_threshold)**2
            
    
    # Find each type of available edge
    ab_edges,ab_freqs,ab_scores = directional_matches(
        img_ab,img_a,img_b,a_uniques,b_uniques,threshold=t,silent=silent,data=data)
    if verbose >= 2: print 'Finished AB edges!'
        
    #aa_edges,aa_freqs,aa_scores = directional_matches(
    #    img_aa,img_a,img_a,a_uniques,a_uniques,threshold=t_shared,silent=silent)
    #if verbose >= 2: print 'Finished AA edges!'
        
    #bb_edges,bb_freqs,bb_scores = directional_matches(
    #    img_bb,img_b,img_b,b_uniques,b_uniques,threshold=t_shared,silent=silent)
    #if verbose >= 2: print 'Finished BB edges!'

        
    if verbose >= 1: print 'Finished edge detection, analyzing graph...'
    
    if real_data:
        for i,pair in enumerate(ab_edges):
            if False: #'CAASGD' in data.seqs_alpha[pair[0]] and 'CASSRDN' in data.seqs_beta[pair[1]]:
                print 'Pair {}:'.format(i+1)
                print '    Alpha - {}'.format(data.seqs_alpha[pair[0]])
                print '    Beta  - {}'.format(data.seqs_beta[pair[1]] )
                print '    Frequency - {}'.format(ab_freqs[i] )
                print '    Confidence: {}'.format(ab_scores[i])

    elif not real_data:
        real_matches = data.metadata['cells']
        # checks to see these actually occur in data
        potential_ab_matches = possible_matches(real_matches,img_ab,a_uniques,b_uniques)
        
        # solves for true edges
        all_edges = [ab_edges]#,aa_edges,bb_edges]
        all_freqs = [ab_freqs]#,aa_freqs,bb_freqs]
        all_scores = [ab_scores]#,aa_scores,bb_scores]
        all_uniques = [a_uniques,b_uniques]
        
        # graph reduction
        results = reduce_graph(all_edges,all_freqs,all_scores,all_uniques)
        
        if verbose >= 1: print 'Finished!'
        
        return results
    
    

    '''
    Scanning for alpha -> beta graph edges
    '''




