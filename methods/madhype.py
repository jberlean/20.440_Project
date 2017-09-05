
'''

Project: 
MAD-HYPE (Multicell Analytical Deconvolution for High Yield Paired-chain Evaluation)

Update (8/29/2017):
    Modified script to create dictionaries represeting each well, allows circumventation of checking many sequences

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
import os, sys
import pickle
import itertools

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

# TODO: what do I even do?
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

def precalculate_match_scores(w_tot,match_prior = 0.5):
    scores,freqs = np.zeros((w_tot+1,w_tot+1,w_tot+1)),np.zeros((w_tot+1,w_tot+1,w_tot+1))
    print ''
    for w_ab in xrange(w_tot+1):
        for w_a in xrange(w_tot+1):
            for w_b in xrange(w_tot+1):
                if w_ab + w_a + w_b > w_tot: 
                    continue
                else:
                    scores[w_ab][w_a][w_b],freqs[w_ab][w_a][w_b] = match_score(w_ab,w_a,w_b,w_tot,match_prior) # effectively take out prior
        #print 'Finished {}/{}...'.format(w_ab,w_tot)
    return scores,freqs

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
# TODO: add stochasicity to well dismissal

def find_all_clones(well_data, a_uniqs, b_uniqs, threshold=0.99, silent=False, distinct=False, chainset_wells_memo = dict()):
    def augment_chainsets(chainsets, well_data, threshold):
        # requires that chainset_wells_memo has the info on all elements of chainsets

        # important storage variables
        removed = set()
        added = set()

        # start iterating through well sets
        #for ind, well in enumerate(well_data):
        #    for chainset1 in chainsets:
        #        if ind not in chainset_wells_memo[chainset1]:  continue
        tstart = time.time()
        pairs_to_check = {c: set() for c in chainsets}
        for (well_ind, well),c in itertools.product(enumerate(well_data), chainsets):
                if well_ind in chainset_wells_memo[c]:
                    pairs_to_check[c] |= well
        t1 = time.time() - tstart
        #def abc(p):
        #    if p[0][0] in chainset_wells_memo[p[1]]:  pairs_to_check[p[1]] |= p[0][1]
        #map(abc, itertools.product(enumerate(well_data), chainsets))

        
        tstart = time.time()
        for ind, chainset1 in enumerate(chainsets):
#            for well_ind in chainset_wells_memo[chainset1]:
#                well = well_data[well_ind]
                for chainset2 in pairs_to_check[chainset1]:
                    merged = (tuple(sorted(chainset1[0]+chainset2[0])), tuple(sorted(chainset1[1]+chainset2[1])))#merge_chainsets(chainset1, chainset2)
                    #if merged == chainset1 or merged in added:  continue
                    #assert len(chainset1[0]+chainset1[1])+len(chainset2[0]+chainset2[1]) == len(merged[0]+merged[1])

                    # get score and frequency for parameter set
                    score,freq = score_chainset_merge(chainset1, chainset2, well_data, precalc_scores = scores, precalc_freqs = freqs)
                    #score *= scores_dict[chainset1] * scores_dict[chainset2]

                    if score > threshold:
                        #print len(chainset_wells_memo[chainset1]), len(chainset_wells_memo[chainset2]), len(chainset_wells_memo[merged])
                        #print chainset1, chainset2, merged, score
                        #print scores_dict[chainset1], scores_dict[chainset2], score
                        added.add(merged)
                        removed.add(chainset1)
                        scores_dict[merged] = score
                        freqs_dict[merged] = freq
    
                #print 'Finished {}/{}\r'.format(ind+1,len(chainsets)), 
            #sys.stdout.flush()
        t2 = time.time() - tstart            


        print ''
        print "{}% spent in prelims, {}% spent in work".format(t1/(t1+t2)*100,100*t2/(t1+t2))
        
        return removed, added
    def find_hidden_chainsets(potential_hidden, potential_hiders, well_data, threshold):
        ## TODO
        # for each clone removed in pass2:
        #   retrieve the clones that are overshadowing it
        #   determine the wells in which these clones are present (elim_wells)
        #   recalculate the chainset score with these wells removed
        #   if the score is good enough, add this to the set of hidden clones
        # return the set of hidden clones found
        return set()
    def score_chainset_merge(cs1, cs2, well_data, precalc_scores = None, precalc_freqs = None, elim_wells = set()):
        merged = (tuple(sorted(cs1[0]+cs2[0])), tuple(sorted(cs1[1]+cs2[1])))
    
        # TODO: replace with try (eafp)
#        if cs1 not in chainset_wells_memo:
#            chainset_wells_memo[cs1] = set([i for i in xrange(len(well_data)) if all([a in well_data[i][0] for a in cs1[0]]) and all([b in well_data[i][1] for b in cs1[1]])])
#        if cs2 not in chainset_wells_memo:
#            chainset_wells_memo[cs2] = set([i for i in xrange(len(well_data)) if all([a in well_data[i][0] for a in cs2[0]]) and all([b in well_data[i][1] for b in cs2[1]])])
#        if merged not in chainset_wells_memo:
#            chainset_wells_memo[merged] = chainset_wells_memo[cs1] & chainset_wells_memo[cs2]
#            
        try:
            w_a = len(chainset_wells_memo[cs1] - elim_wells)
        except:
            chainset_wells_memo[cs1] = set([i for i in xrange(len(well_data)) if all([a in well_data[i][0] for a in cs1[0]]) and all([b in well_data[i][1] for b in cs1[1]])])
            w_a = len(chainset_wells_memo[cs1] - elim_wells)
        try:
            w_b = len(chainset_wells_memo[cs2] - elim_wells)
        except:
            chainset_wells_memo[cs2] = set([i for i in xrange(len(well_data)) if all([a in well_data[i][0] for a in cs2[0]]) and all([b in well_data[i][1] for b in cs2[1]])])
            w_b = len(chainset_wells_memo[cs2] - elim_wells)
            
        try:
            w_ab = len(chainset_wells_memo[merged] - elim_wells)
        except:
            chainset_wells_memo[merged] = chainset_wells_memo[cs1] & chainset_wells_memo[cs2]
            w_ab = len(chainset_wells_memo[merged] - elim_wells)
            
        n_a = w_a - w_ab
        n_b = w_b - w_ab
        n_ab = w_ab
        n_tot = len(well_data) - len(elim_wells)

        if w_ab < 4:  return 0,0
        
        if precalc_scores is not None and precalc_freqs is not None:
            score = precalc_scores[n_ab,n_a,n_b]
            freq = precalc_freqs[n_ab,n_a,n_b]
        else:
            score, freq = match_score(n_ab,n_a,n_b,n_tot) # calc score manually

        return score, freq
    def merge_chainsets(cs1, cs2):
        return (tuple(sorted(set(cs1[0]+cs2[0]))), tuple(sorted(set(cs1[1]+cs2[1]))))

    # reformat well data
    well_data = [set([((a,),()) for a in well[0]] + [((),(b,)) for b in well[1]]) for well in well_data]
    w_tot = len(well_data)

    # retrieve pickled precalculated scores/freqs, or generate them anew if not found
    if not os.path.exists('./pickles'): os.makedirs('./pickles')
    if os.path.isfile('./pickles/val{}.p'.format(w_tot)): (scores,freqs) = pickle.load(open('./pickles/val{}.p'.format(w_tot),'r'))
    else: 
        scores,freqs = precalculate_match_scores(w_tot, match_prior=0.5)
        pickle.dump((scores,freqs),open('./pickles/val{}.p'.format(w_tot),'w'))

    # initialization
    scores_dict = {c: 1 for c in [((a,),()) for a in a_uniqs] + [((),(b,)) for b in b_uniqs]}
    freqs_dict = scores_dict.copy()
    for a in a_uniqs:  chainset_wells_memo[((a,),())] = set(filter(lambda i: ((a,),()) in well_data[i], xrange(w_tot)))
    for b in b_uniqs:  chainset_wells_memo[((),(b,))] = set(filter(lambda i: ((),(b,)) in well_data[i], xrange(w_tot)))

    # Do first pass (finding alpha-beta pairs)
    well_data_beta = [set(filter(lambda c: c[0]==(), well)) for well in well_data] # well_data only containing beta chains
    well_data_alpha = [set(filter(lambda c: c[1]==(), well)) for well in well_data] # well_data only containing alpha chains
    init_chainsets = [((a,),()) for a in a_uniqs]
    _, pass1 = augment_chainsets(init_chainsets, well_data_beta, threshold = threshold)
    _, pass1_2 = augment_chainsets(init_chainsets, well_data_alpha, threshold = threshold)
    _, pass1_3 = augment_chainsets([((),(b,)) for b in b_uniqs], well_data_beta, threshold=threshold) 
    pass1 |= pass1_2 | pass1_3
    #init_chainsets = [((a,),()) for a in a_uniqs] + [((),(b,)) for b in b_uniqs]
    #_, pass1 = augment_chainsets(init_chainsets, well_data, threshold = threshold)
    print '***', len(pass1)

    # Do second pass (adding alpha- or beta-chains to the a-b pairs)
    #pass1_removed, pass2 = augment_chainsets(pass1, well_data, threshold = threshold) # TODO: implement more stringent threshold
    #print len(pass1_removed), len(pass2)
    #for i,c in enumerate(pass2):
    #    if not all([merge_chainsets(*temp) in pass1 for temp in itertools.combinations([((a,),()) for a in c[0]] + [((),(b,)) for b in c[1]], 2)]):
    #        print '***',i, c
    #    for temp in itertools.combinations([((a,),()) for a in c[0]] + [((),(b,)) for b in c[1]], 2):
    #        if merge_chainsets(*temp) not in pass1:  print merge_chainsets(*temp)
    pass1_removed, pass2 = set(), set()

    # Do third pass (find clones with 4 chains, i.e. dual in both alpha and beta)
    #pass2_removed, pass3 = augment_chainsets(pass2, well_data, threshold = threshold) # TODO: implement more stringent threshold
    pass2_removed, pass3 = set(), set()

    # Compile preliminary list of clones to return
    clones_4chains = pass3
    clones_3chains = pass2 - pass2_removed
    clones_2chains = pass1 - pass1_removed

    # Recheck clones that were removed in prior passes due to being "overshadowed" by larger chainset
    pass2_hidden = find_hidden_chainsets(pass2_removed, clones_4chains, well_data, threshold = threshold)
    clones_3chains |= pass2_hidden

    pass1_hidden = find_hidden_chainsets(pass1_removed, clones_3chains, well_data, threshold = threshold)
    
    return clones_4chains|clones_3chains|clones_2chains, scores_dict, freqs_dict
           
       
# NOTE: This is still a relevant function, but I think it makes more sense to make a new fn for compile
# TODO: actually figure out graph solving feature

class CollectResults:
    def __init__(self):
        self.cells = []
        self.cell_frequencies = []
        self.threshold = []
        self.cell_frequencies_CI = []

    def add_results(self,all_edges,all_freqs,all_scores,id):

        # TODO: add aa/bb solver ( this is non-trivial math )

        # results['cells'] = [e[0] for edge in all_edges for e in edge] # eventually, when we start doing multiclones
        if id == 'ab' or id == 'AB':
            self.cells += [((e[0],),(e[1],)) for e in all_edges]
        elif id == 'aa' or id == 'AA':
            self.cells += [((e[0],e[1]),((),)) for e in all_edges]
        elif id == 'aa' or id == 'AA':
            self.cells += [(((),),(e[0],e[1])) for e in all_edges]

        self.cell_frequencies += [f for f in all_freqs]
        self.threshold += [s for s in all_scores]
        self.cell_frequencies_CI += [(-s/2,s/2) for s in all_scores] # this is pseudo, and totally not right

    def export(self):
        results = {
                    'cells':self.cells,
                    'cell_frequencies':self.cell_frequencies,
                    'threshold':self.threshold,
                    'cell_frequencies_CI':self.cell_frequencies_CI
                  }
        return results


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

    # create dictionaries for well presence
    a_wells = dict([(a_u,[]) for a_u in a_uniques])
    b_wells = dict([(b_u,[]) for b_u in b_uniques])
    
    if verbose >= 1: print 'Starting variable creation...'
    
    # TODO: Replace with well_dictionary method
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        for a in well[0]: a_wells[a].append(w) # drop well indices in each unique chains entry
        for b in well[1]: b_wells[b].append(w) # drop well indices in each unique chains entry
        print 'Image generation progress... {}%\r'.format(100*(w+1)/w_tot),

    print ''

    if verbose >= 1: print 'Starting edge detection...'
            
    # Setup threshold values (TODO: better way to distinguish ab,ba from aa,bb
    t = pair_threshold
    t_shared = 1 - (1 - pair_threshold)**2
            
    # Chain-set to well dictionary implementation; should obviate the use of a_wells and b_wells with some work
    chainset_wells = dict()
    for a, wells in a_wells.iteritems():  chainset_wells[((a,),())] = set(wells)
    for b, wells in b_wells.iteritems():  chainset_wells[((),(b,))] = set(wells)
    
    # Prep well data into sets (this is more efficient)
    well_data = [[set(alphas), set(betas)] for alphas,betas in data.well_data]
    
    # Find clones
    clones,scores_dict,freqs_dict = find_all_clones(
        well_data, a_uniques, b_uniques, threshold=t,silent=silent,distinct=True, chainset_wells_memo = chainset_wells)
    if verbose >= 2: print 'Finished clone-finding!'
    
    #aa_edges,aa_freqs,aa_scores = directional_matches(
    #    a_wells,a_wells,a_uniques,a_uniques,w_tot,threshold=t,silent=silent,distinct=False)
    #if verbose >= 2: print 'Finished AA edges!'
        
    #bb_edges,bb_freqs,bb_scores = directional_matches(
    #    b_wells,b_wells,b_uniques,b_uniques,w_tot,threshold=t,silent=silent,distinct=False)
    #if verbose >= 2: print 'Finished BB edges!'
    
    ab_edges = list(set(itertools.chain.from_iterable(itertools.product(c[0], c[1]) for c in clones)))
    ab_freqs = [freqs_dict[((a,),(b,))] for a,b in ab_edges]
    ab_scores = [scores_dict[((a,),(b,))] for a,b in ab_edges]
    
    aa_edges = list(set(itertools.chain.from_iterable(itertools.combinations(c[0], 2) for c in clones)))
    aa_freqs = [freqs_dict[((a1,a2),())] for a1,a2 in aa_edges]
    aa_scores = [scores_dict[((a1,a2),())] for a1,a2 in aa_edges]
    
    bb_edges = list(set(itertools.chain.from_iterable(itertools.combinations(c[1], 2) for c in clones)))
    bb_freqs = [freqs_dict[((),(b1,b2))] for b1,b2 in bb_edges]
    bb_scores = [scores_dict[((),(b1,b2))] for b1,b2 in bb_edges]
        
    if verbose >= 1: print 'Finished finding clones'
    
    if real_data:
        results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
        results_aa = {'edges':aa_edges,'freqs':aa_freqs,'scores':aa_scores}
        results_bb = {'edges':bb_edges,'freqs':bb_freqs,'scores':bb_scores}
        results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}
        return results

    elif not real_data:
        real_matches = data.metadata['cells']
        # checks to see these actually occur in data
        #potential_ab_matches = possible_matches(real_matches,img_ab,a_uniques,b_uniques)
        
        # solves for true edges
        #all_edges = [ab_edges]#,aa_edges,bb_edges]
        #all_freqs = [ab_freqs]#,aa_freqs,bb_freqs]
        #all_scores = [ab_scores]#,aa_scores,bb_scores]
        #all_uniques = [a_uniques,b_uniques]
       
        #results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
        #results_aa = {'edges':aa_edges,'freqs':aa_freqs,'scores':aa_scores}
        #results_bb = {'edges':bb_edges,'freqs':bb_freqs,'scores':bb_scores}

        #results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}
        
        # deposit results in the collect results function 
        compiler = CollectResults()
        compiler.add_results(ab_edges,ab_freqs,ab_scores,'AB')
        #compiler.add_results(aa_edges,aa_freqs,aa_scores,'AA')
        #compiler.add_results(bb_edges,bb_freqs,bb_scores,'BB')
        
        if verbose >= 1: print 'Finished!'
        
        return compiler.export() 
    
    

    '''
    Scanning for alpha -> beta graph edges
    '''




