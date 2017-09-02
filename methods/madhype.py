
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
import os
import pickle
import itertools

# nonstandard libraries
import numpy as np
import scipy.misc



'''
Specialty Calculators
'''


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

def directional_matches(a_wells,b_wells,a_uniqs,b_uniqs,w_tot,threshold=0.99,silent=False,distinct=False):

    # important storage variables
    predicted_ab = []
    predicted_frequency = []
    predicted_score = []

    print 'Making C++ data files...'
    print a_wells
    
    with open('chain_data_a.txt','w') as f:
        for w in a_uniqs: f.write('{}'.format(str(a_wells[w]))[1:-1]+'\n')
    with open('chain_data_b.txt','w') as f:
        for w in b_uniqs: f.write('{}'.format(str(b_wells[w]))[1:-1]+'\n')
    with open('uniques_a.txt','w') as f:
        for w in a_uniqs: f.write('{}\n'.format(w))
    with open('uniques_b.txt','w') as f:
        for w in b_uniqs: f.write('{}\n'.format(w))

    raw_input('Hold...')
        
    # create dictionaries that map each well dictionary to the reverse (well # -> chains)
    wells_a = dict([(i,set()) for i in xrange(w_tot)])
    wells_b = dict([(i,set()) for i in xrange(w_tot)])

    # iterate across both dictionaries
    print 'Starting reverse well dictionary formation...'
    for wells_x,x_wells in [(wells_a,a_wells),(wells_b,b_wells)]:
        for k,v in x_wells.items():
            for well in v:
                wells_x[well].add(k)
    print 'Finished!'

    # find a copy of precalulcated scores, build one if it doesn't exist 
    if not os.path.exists('./pickles'): os.makedirs('./pickles')
    if os.path.isfile('./pickles/val{}.p'.format(w_tot)): (scores,freqs) = pickle.load(open('./pickles/val{}.p'.format(w_tot),'r'))
    else: 
        scores,freqs = precalculate_match_scores(w_tot, match_prior=0.5)
        pickle.dump((scores,freqs),open('./pickles/val{}.p'.format(w_tot),'w'))


    # start iterating through well sets
    tested_pairs = set()
    for w in xrange(w_tot):
        print len(wells_a[w]),len(wells_b[w])
        p_tot = len(wells_a[w])*len(wells_b[w])
        pairs = itertools.product(wells_a[w],wells_b[w])
        for ind,p in enumerate(pairs):
            if not distinct and p[0] == p[1]: continue
            if p in tested_pairs:
                continue
            else:
                # calculate well occurances
                #w_ab = len([a for a in a_wells[p[0]] if a in b_wells[p[1]]]) 
                w_ab = len(set(a_wells[p[0]]).intersection(b_wells[p[1]])) # taps in to O(min(a,b)) for set intersection
                w_a,w_b = len(a_wells[p[0]]),len(b_wells[p[1]])

                # produce well counts (cross-over removed)
                n_ab = w_ab
                n_a,n_b = w_a - n_ab,w_b - n_ab
                n_tot = w_tot

                # create some catches for obviously bad cases
                if n_ab <= 3: continue # since atleast 4 occurances, better be higher than 3 matches
                f_ab = float(w_a*w_b)/(w_tot**2)

                # get score and frequency for parameter set
                score,frequency = scores[n_ab,n_a,n_b],freqs[n_ab,n_a,n_b] 
                 
                if score > threshold:
                    predicted_ab.append(p)
                    predicted_frequency.append(frequency)
                    predicted_score.append(score)

                tested_pairs.add(p)

            if ind%1000000 == 0:
                print 'Finished {}/{}'.format(ind+1,p_tot) 
        
    print ''
     
    return predicted_ab,predicted_frequency,predicted_score # returns edges


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
def solve(data,pair_threshold = 0.99,verbose=0,real_data=False,all_pairs=True):
    
    if verbose >= 5: silent = False
    else: silent = True

    w_tot = len(data.well_data)

    # Find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # Generate reverse dictionary
    a_wells = dict([(a_u,[]) for a_u in a_uniques])
    b_wells = dict([(b_u,[]) for b_u in b_uniques])
    
    if verbose >= 1: print 'Starting reverse dictionary creation...'
    
    # TODO: Replace with well_dictionary method
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        for a in well[0]: a_wells[a].append(w) # drop well indices in each unique chains entry
        for b in well[1]: b_wells[b].append(w) # drop well indices in each unique chains entry
        print 'Reverse dictionary progress... {}%\r'.format(100*(w+1)/w_tot),
    print ''

    # try to reduce dependencies on data processing here
    with open('./solver/chain_data_a.txt','w') as f:
        for w in a_uniques: f.write('{}'.format(str(a_wells[w]))[1:-1]+'\n')
    with open('./solver/chain_data_b.txt','w') as f:
        for w in b_uniques: f.write('{}'.format(str(b_wells[w]))[1:-1]+'\n')
    with open('./solver/uniques_a.txt','w') as f:
        for w in a_uniques: f.write('{}\n'.format(w))
    with open('./solver/uniques_b.txt','w') as f:
        for w in b_uniques: f.write('{}\n'.format(w))

    # C++ embedding  
    arg1,arg2,arg3 = str(w_tot),str(-math.log10(1.-pair_threshold)),str(1)

    os.system(os.getcwd() + '/solver/a.out {} {} {}'.format(arg1,arg2,arg3))
    
    # real data post-processing section
    if real_data:
        results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
        results_aa = {'edges':aa_edges,'freqs':aa_freqs,'scores':aa_scores}
        results_bb = {'edges':bb_edges,'freqs':bb_freqs,'scores':bb_scores}
        results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}
        return results

    # non-so-real data post-processing
    elif not real_data:
        # recalls real matches
        real_matches = data.metadata['cells']
        
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




