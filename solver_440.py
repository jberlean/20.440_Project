
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

# nonstandard libraries
import numpy as np
#import matplotlib.pyplot as plt
#from seq_generator import SequencingGenerator as SeqGen
#from seq_data import SequencingData


'''
Factory Methods
'''

# N choose K
def nCk(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom
    
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

def match_score(w_ab,w_a,w_b,w_tot):
    if w_ab == 0: return 0.
    else:
        mp,nmp = match_probability(w_ab,w_a,w_b,w_tot),nonmatch_probability(w_ab,w_a,w_b,w_tot)
        return (10**mp)/(10**mp+10**nmp)

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
        
Returns:
    predicted_ab - 
        description: list of tuples (a label,b label)
        
'''
def directional_matches(img_ab,img_a,img_b,a_uniques,b_uniques,threshold=0.99,silent=False):
    score = np.zeros((len(a_uniques),len(b_uniques)))
    w_tot = img_ab.shape[2]
    predicted_ab = []

    for i in xrange(img_ab.shape[0]):

        if not silent: print 'Starting analysis for alpha chain {}...\n'.format(a_uniques[i])

        unexplained_wells = list(xrange(w_tot))

        while True: # TODO: create exit condition
            # create counts based on unexplained well data
            w_ab = np.sum(img_ab[:,:,unexplained_wells],axis=2)
            w_a = np.sum(img_a[:,unexplained_wells],axis=1)
            w_b = np.sum(img_b[:,unexplained_wells],axis=1)

            # assign scores
            for j in xrange(img_ab.shape[1]):
                n_ab = w_ab[i,j]
                n_a,n_b = w_a[i] - n_ab,w_b[j] - n_ab
                n_tot = len(unexplained_wells)
                score[i,j] = match_score(n_ab,n_a,n_b,n_tot)

            # find best pair to add
            j = np.argmax(np.nan_to_num(score[i,:]))

            # check if the score is passable
            if score[i,j] > threshold:
                predicted_ab.append((a_uniques[i],b_uniques[j]))
                explained_wells = [k for k in xrange(w_tot) if img_ab[i,j,k] == 1]
                unexplained_wells = [item for item in unexplained_wells if item not in explained_wells]
                predicted_ab.append((a_uniques[i],b_uniques[j]))
                if not silent:
                    n_ab = w_ab[i,j]
                    n_a,n_b = w_a[i] - n_ab,w_b[j] - n_ab
                    n_tot = len(unexplained_wells)
                    print 'Accepted match ({},{}) with score {}.'.format(a_uniques[i],b_uniques[j],score[i,j])
                    print '{} newly explained cells, {} remaining.'.format(len(explained_wells),len(unexplained_wells))
                    print 'N_ab,N_a,N_b,N_tot:',n_ab,n_a,n_b,n_tot
            else:
                if not silent: print 'No more sufficient matches found for ({}), exiting...\n'.format(a_uniques[i])
                break

    return predicted_ab # returns edges

def generate_cells(num_cells, max_alphas=None, max_betas=None):
    if max_alphas == None:  max_alphas = num_cells
    if max_betas == None:  max_betas = num_cells

    # Generate the degree for each alpha- and beta-chain from a given distribution
    sharing_probs=[0.8375, 0.0805, 0.029, 0.013, 0.021, 0.0025, 0.0165] # Averages from the Lee et al. paper
    adegs = np.random.choice(range(1,8), max_alphas, replace=True, p=sharing_probs)
    bdegs = np.random.choice(range(1,8), max_betas, replace=True, p=sharing_probs)

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



'''
The Meat-and-Potatoes Function
'''
def solve(data,pair_threshold = 0.99,silent=False):
    # find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # sort in ascending order
    sorted(a_uniques,key=int)
    sorted(b_uniques,key=int)

    # counts of each unique in wells
    w_a = [float(sum([1 for well in data.well_data if i in well[0]])) for i in a_uniques] 
    w_b = [float(sum([1 for well in data.well_data if i in well[1]])) for i in b_uniques]

    img_ab,img_ba,img_aa,img_bb = None,None,None,None
    img_a = None
    img_b = None

    # creates all the necessary images of the data
    for well in data.well_data:
        a_ind = [a_uniques.index(i) for i in well[0]]
        b_ind = [b_uniques.index(i) for i in well[1]]
        a_v,b_v = np.zeros((len(a_uniques),1)),np.zeros((len(b_uniques),1)) 
        np.put(a_v,a_ind,np.ones((len(a_ind))))
        np.put(b_v,b_ind,np.ones((len(b_ind))))
        if img_ab == None: # if nothing has been assigned so far
            img_ab = np.matmul(a_v,np.transpose(b_v))
            img_ba = np.matmul(b_v,np.transpose(a_v))
            img_aa = self_occurence(a_v)
            img_bb = self_occurence(b_v)
            img_a = a_v
            img_b = b_v
        else: 
            img_ab = np.dstack((img_ab,np.matmul(a_v,np.transpose(b_v))))
            img_ba = np.dstack((img_ba,np.matmul(b_v,np.transpose(a_v))))
            img_aa = np.dstack((img_aa,self_occurence(a_v)))
            img_bb = np.dstack((img_bb,self_occurence(b_v)))
            img_a = np.hstack((img_a,a_v))
            img_b = np.hstack((img_b,b_v))



    real_matches = data.metadata['cells']

    t = pair_threshold
    t_shared = 1 - (1 - pair_threshold)**2

    ab_edges = directional_matches(img_ab,img_a,img_b,a_uniques,b_uniques,threshold=t,silent = silent)
    if not silent: print 'Finished AB edges!'
    ba_edges = directional_matches(img_ba,img_b,img_a,b_uniques,a_uniques,threshold=t,silent = silent)
    if not silent: print 'Finished BA edges!'
    aa_edges = directional_matches(img_aa,img_a,img_a,a_uniques,a_uniques,threshold=t_shared,silent = silent)
    if not silent: print 'Finished AA edges!'
    bb_edges = directional_matches(img_bb,img_b,img_b,b_uniques,b_uniques,threshold=t_shared,silent = silent)
    if not silent: print 'Finished BB edges!'


    # checks to see these actually occur in data
    potential_ab_matches = possible_matches(real_matches,img_ab,a_uniques,b_uniques)

    ab_edges_rev = [(a[1],a[0]) for a in ba_edges] # create symmetrical order
    predicted_ab_matches = list(set(ab_edges).intersection(ab_edges_rev))
    correct_ab_matches = [ab for ab in predicted_ab_matches if ab in potential_ab_matches]

    print 'Number of correct guesses: {}/{}'.format(len(correct_ab_matches),len(predicted_ab_matches))
    print 'Number of successes: {}/{}'.format(len(correct_ab_matches),len(potential_ab_matches))

    '''
    Scanning for alpha -> beta graph edges
    '''




