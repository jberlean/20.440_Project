
"""

MAD-HYPE Algorithm

This is the core scorer
Written in Python for testing...

"""


# Library importation
import sys

import math
import scipy.misc
from scipy.special import gammaln

######################
## DATA PROCESSING ###
##      METHODS    ###
######################

def intersection(lst1, lst2):
    return [x*y for x,y in zip(lst1,lst2)]

def intersection_sum(lst1, lst2):
    return sum(x*y for x,y in zip(lst1,lst2))

# Extract vector from string
def ParseCSString(string):
    return [int(v) for v in string.split(',')]


# Finds the unique indices in data
def LoadUniques(fname):
    file = open(fname)
    uniques = [int(v) for v in file]
    file.close()
    return uniques


# Well data load
def LoadChainData(fname, w_tot):
    chain_data = []

    file = open(fname)
    for line in file:
      data = [0]*w_tot
      line_vals = ParseCSString(line)
      for w in line_vals:  data[w] = 1
      chain_data.append(data)

    file.close()
    return chain_data

def LoadInitChainsets(fname):
    file = open(fname)

    init_chainsets = []
    for line in file:
      a1,a2,b1,b2  = ParseCSString(line)
      alist = () if a1==-1 else (a1,) if a2==-1 else (a1,a2)
      blist = () if b1==-1 else (b1,) if b2==-1 else (b1,b2)
      init_chainsets.append((alist, blist))

    file.close()
    return init_chainsets


####################
## COMPUTATIONAL ###
##    METHODS    ###
####################

# Returns N choose K integer

def nCk(n, k):
    return scipy.misc.comb(n,k)

def multinomial_prob(N, P):
    if any([n>0 and p==0 for n,p in zip(N,P)]):
      return 0

    coeff = gammaln(sum(N)+1) - sum([gammaln(n+1) for n in N])
    T = [n*math.log(p) if p!=0 else 0 for n,p in zip(N,P)]
    return math.exp(coeff + sum(T))

# Non-match MLE estimator for f_ab,f_a,f_b

def nonmatch_frequency(w_a, w_b, w_ab, w_tot):
    f_a = float(w_a+w_ab)/w_tot
    f_b = float(w_b+w_ab)/w_tot
    return f_a, f_b

def nonmatch_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot):
    w_o = w_tot - w_a - w_b - w_c - w_ab - w_ac - w_bc - w_abc;

    p_o = w_o/float(w_tot)
    p_nab = (w_o + w_c)/float(w_tot)
    p_nac = (w_o + w_b)/float(w_tot)
    p_nbc = (w_o + w_a)/float(w_tot)
    p_na = (w_o + w_b + w_c + w_bc)/float(w_tot)
    p_nb = (w_o + w_a + w_c + w_ac)/float(w_tot)
    p_nc = (w_o + w_a + w_b + w_ab)/float(w_tot)

    if p_o*p_nc > 0:  f_ab = max(0., 1 - p_nbc*p_nac/(p_o*p_nc))
    else:  f_ab = 0
    if p_o*p_nb > 0:  f_ac = max(0., 1 - p_nab*p_nbc/(p_o*p_nb))
    else:  f_ac = 0
    if p_o*p_na > 0:  f_bc = max(0., 1 - p_nab*p_nac/(p_o*p_na))
    else:  f_bc = 0

    if p_nbc > 0:  f_a = 1 - p_o/p_nbc
    else:  f_a = 0
    if p_nac > 0:  f_b = 1 - p_o/p_nac
    else:  f_b = 0
    if p_nab > 0:  f_c = 1 - p_o/p_nab
    else:  f_c = 0

    return f_a, f_b, f_c, f_ab, f_ac, f_bc


# Match MLE estimator for f_ab,f_a,f_b

def match_frequency(w_a, w_b, w_ab, w_tot):
    if (w_tot - w_ab - w_b == 0):  f_a = 0
    else:  f_a = float(w_a)/(w_tot-w_ab-w_b)

    if (w_tot - w_ab - w_a == 0):  f_b = 0
    else:  f_b = float(w_b)/(w_tot-w_ab-w_a)

    if f_a*f_b < 1:
      f_ab = max(0, (w_ab - f_a*f_b*w_tot) / ((1 - f_a*f_b)*w_tot))
    else:
      f_ab = 0

    return f_a, f_b, f_ab

def match_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot):
    w_o = w_tot - w_a - w_b - w_c - w_ab - w_ac - w_bc - w_abc

    if (w_a + w_o == 0):  f_a = 0
    else:  f_a = float(w_a) / (w_a+w_o)
    if (w_b + w_o == 0):  f_b = 0
    else:  f_b = float(w_b) / (w_b+w_o)
    if (w_c + w_o == 0):  f_c = 0
    else:  f_c = float(w_c) / (w_c+w_o)

    if (w_o+w_a+w_b+w_ab > 0 and f_a*f_b < 1):  f_ab = max(0., (float(w_ab) / (w_o+w_a+w_b+w_ab) - f_a*f_b) / (1 - f_a*f_b))
    else:  f_ab = 0
    if (w_o+w_a+w_c+w_ac > 0 and f_a*f_c < 1):  f_ac = max(0., (float(w_ac) / (w_o+w_a+w_c+w_ac) - f_a*f_c) / (1 - f_a*f_c))
    else:  f_ac = 0
    if (w_o+w_b+w_c+w_bc > 0 and f_b*f_c < 1):  f_bc = max(0., (float(w_bc) / (w_o+w_b+w_c+w_bc) - f_b*f_c) / (1 - f_b*f_c))
    else:  f_bc = 0

    f_false_abc = f_ab*(f_ac + (1-f_ac)*(f_bc + (1-f_bc)*f_c)) + (1-f_ab)*(f_ac*(f_bc + (1-f_bc)*f_b) + (1-f_ac)*(f_bc*f_a + (1-f_bc)*f_a*f_b*f_c))
    if f_false_abc < 1:
      f_abc = max(0., (float(w_abc) / w_tot - f_false_abc) / (1 - f_false_abc))
    else:
      f_abc = 0.

    return f_a, f_b, f_c, f_ab, f_ac, f_bc, f_abc

# Instantaneous probability for nonmatch instance

def nonmatch_instant_probability(w_a, w_b, w_ab, w_tot, f_a, f_b):
    val =  multinomial_prob([w_tot-w_a-w_b-w_ab, w_a, w_b, w_ab], [(1-f_a)*(1-f_b), f_a*(1-f_b), f_b*(1-f_a), f_a*f_b])
    return val;


# Instantaneous probability for match instance

def match_instant_probability(w_a, w_b, w_ab, w_tot, f_a, f_b, f_ab):
    val =  multinomial_prob([w_tot-w_a-w_b-w_ab, w_a, w_b, w_ab], [(1-f_a)*(1-f_b)*(1-f_ab), f_a*(1-f_b)*(1-f_ab), f_b*(1-f_a)*(1-f_ab), 1 - (1-f_a*f_b)*(1-f_ab)])
    return val


# Non-match probability calculation

def nonmatch_probability(w_a, w_b, w_ab, w_tot):
    f_a, f_b = nonmatch_frequency(w_a, w_b, w_ab, w_tot)
    prob = nonmatch_instant_probability(w_a, w_b, w_ab, w_tot, f_a, f_b)
    return prob
def nonmatch_probability_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot):
    f_a, f_b, f_c, f_ab, f_ac, f_bc = nonmatch_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)

    n_a, n_b, n_c, n_ab, n_ac, n_bc = 1-f_a, 1-f_b, 1-f_c, 1-f_ab, 1-f_ac, 1-f_bc

    p_abc = f_ab*(f_ac + n_ac*f_bc + n_ac*n_bc*f_c) + n_ab*(f_ac*(f_bc + n_bc*f_b) + n_ac*(f_bc*f_a + n_bc*f_a*f_b*f_c))
    p_ab = n_ac*n_bc*n_c*(f_ab + n_ab*f_a*f_b)
    p_ac = n_ab*n_bc*n_b*(f_ac + n_ac*f_a*f_c)
    p_bc = n_ab*n_ac*n_a*(f_bc + n_bc*f_b*f_c)
    p_a = n_ab*n_ac*n_bc*f_a*n_b*n_c
    p_b = n_ab*n_ac*n_bc*n_a*f_b*n_c
    p_c = n_ab*n_ac*n_bc*n_a*n_b*f_c
    p_o = n_ab*n_ac*n_bc*n_a*n_b*n_c

    prob = multinomial_prob([w_tot-w_a-w_b-w_c-w_ab-w_ac-w_bc-w_abc, w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc], [p_o, p_a, p_b, p_c, p_ab, p_ac, p_bc, p_abc])
    return prob


# Match probability calculation

def match_probability(w_a, w_b, w_ab, w_tot):
    f_a, f_b, f_ab = match_frequency(w_a, w_b, w_ab, w_tot)
    if (f_ab==0):  return 0.0
    prob = match_instant_probability(w_a, w_b, w_ab, w_tot, f_a, f_b, f_ab)
    return prob
def match_probability_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot):
    f_a, f_b, f_c, f_ab, f_ac, f_bc, f_abc = match_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)

    n_a, n_b, n_c, n_ab, n_ac, n_bc, n_abc = 1-f_a, 1-f_b, 1-f_c, 1-f_ab, 1-f_ac, 1-f_bc, 1-f_abc

    p_abc = f_abc + n_abc*(f_ab*(f_ac + n_ac*f_bc + n_ac*n_bc*f_c) + n_ab*(f_ac*(f_bc + n_bc*f_b) + n_ac*(f_bc*f_a + n_bc*f_a*f_b*f_c)))
    p_ab = n_abc*n_ac*n_bc*n_c*(f_ab + f_a*f_b*n_ab)
    p_ac = n_abc*n_ab*n_bc*n_b*(f_ac + f_a*f_c*n_ac)
    p_bc = n_abc*n_ab*n_ac*n_a*(f_bc + f_b*f_c*n_bc)
    p_a = n_abc*n_ab*n_ac*n_bc*f_a*n_b*n_c
    p_b = n_abc*n_ab*n_ac*n_bc*n_a*f_b*n_c
    p_c = n_abc*n_ab*n_ac*n_bc*n_a*n_b*f_c
    p_o = n_abc*n_ab*n_ac*n_bc*n_a*n_b*n_c

    prob = multinomial_prob([w_tot-w_a-w_b-w_c-w_ab-w_ac-w_bc-w_abc, w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc], [p_o, p_a, p_b, p_c, p_ab, p_ac, p_bc, p_abc])

    return prob

# Match score calculator

def match_score(w_a, w_b, w_ab, w_tot):
    # If there are two or fewer matches, its unlikely to be real
    # JB - changed to 3 or fewer matches to be same as backup
    if ( w_ab <= 3 ):
        score = -float('inf')
        freq = 0.
    else:
        mp = match_probability(w_a, w_b, w_ab, w_tot)
        nmp = nonmatch_probability(w_a, w_b, w_ab, w_tot)
        _, _, freq = match_frequency(w_a, w_b, w_ab, w_tot)
        if mp==0 or nmp==0 or freq==0:
          score = -float('inf')
          freq = 0
        else:
          score = math.log10(mp) - math.log10(nmp)
        #print w_a, w_b, w_ab, mp, nmp, score
    return score, freq

def match_score_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot):
    # Calculates the match score for the chains a, b, and c
    # If there are three or fewer matches, its unlikely to be real
    if ( w_abc <= 3 ):
        score = -float('inf')
        freq = 0.
    else:
#        mp = match_probability_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)
#        nmp = nonmatch_probability_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)
#        _, _, _, _, _, _, freq = match_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)
#        if mp==0 or nmp==0 or freq==0 or w_tot-w_abc-w_ab-w_ac-w_bc-w_a-w_b-w_c==0:
#          score = -float('inf')
#          freq = 0
#        else:
#          score = math.log10(mp) - math.log10(nmp)
        _, _, _, _, _, _, freq = match_frequency_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)
        if freq > 0:
          mp1 = match_probability(w_a+w_ac, w_b+w_bc, w_ab+w_abc, w_tot)
          nmp1 = nonmatch_probability(w_a+w_ac, w_b+w_bc, w_ab+w_abc, w_tot)
          mp2 = match_probability(w_ab, w_c + w_ac + w_bc, w_abc, w_tot)
          nmp2 = nonmatch_probability(w_ab, w_c + w_ac + w_bc, w_abc, w_tot)
          #mp2 = match_probability(w_a+w_ab, w_c + w_bc, w_ac+w_abc, w_tot)
          #nmp2 = nonmatch_probability(w_a+w_ab, w_c + w_bc, w_ac + w_abc, w_tot)
          score1_exp = mp1/nmp1
          score2_exp = mp2/nmp2
          score = -float('inf') if score1_exp==0 or score2_exp==0 else -math.log10(1/score1_exp + 1/score2_exp + 1/(score1_exp*score2_exp))
        else:
          score = -float('inf')

        if match_score(w_a+w_ac, w_b+w_bc, w_ab+w_abc, w_tot)[0]<2.9 or match_score(w_a+w_ab, w_c+w_bc, w_ac+w_abc, w_tot)[0]<2.9 or match_score(w_b+w_ab, w_c+w_ac, w_bc+w_abc, w_tot)[0]<2.9:
          score = -float('inf')
        
#        if (score>=7):
#          if match_score(w_a+w_ac, w_b+w_bc, w_ab+w_abc, w_tot)[0]<3 or match_score(w_a+w_ab, w_c+w_bc, w_ac+w_abc, w_tot)[0]<3 or match_score(w_b+w_ab, w_c+w_ac, w_bc+w_abc, w_tot)[0]<3:
#            print w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, mp, nmp, score

    return score, freq

def chainset_num_alphas(c):
    return len(c[0])
def chainset_num_betas(c):
    return len(c[1])
def chainset_size(c):
    return len(c[0]) + len(c[1])

def add_chain_to_chainset(cs, c):
    return (cs[0]+c[0], cs[1]+c[1])

def make_result_string(cs, freq, score, uniques_a, uniques_b):
    a1 = -1 if chainset_num_alphas(cs)==0 else uniques_a[cs[0][0]]
    a2 = -1 if chainset_num_alphas(cs)<2 else uniques_a[cs[0][1]]
    b1 = -1 if chainset_num_betas(cs)==0 else uniques_b[cs[1][0]]
    b2 = -1 if chainset_num_betas(cs)<2 else uniques_b[cs[1][1]]
    return "{}\t{}\t{}\t{}\t{}\t{}".format(score, freq, a1, a2, b1, b2)

# Attempt to add chains to an existing chainset
def try_chain_additions_nondual(chainset, additions, threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot):
    # Calculate chainset well data
    if (chainset_num_alphas(chainset) > 0):
        chainset_welldata = chain_data_a[chainset[0][0]]
    else:
        chainset_welldata = chain_data_b[chainset[1][0]]

    results = []
    for chain in additions:
        # Retrieve data on chain being added
        chain_welldata = chain_data_a[chain[0][0]] if chainset_num_alphas(chain)>0 else chain_data_b[chain[1][0]]

        # Score probability of chainset+chain
        w_ab = intersection_sum(chainset_welldata, chain_welldata)
        w_a = sum(chainset_welldata) - w_ab
        w_b = sum(chain_welldata) - w_ab
        score, freq = match_score(w_a, w_b, w_ab, w_tot)
    
        # Store result if it's good enough
        if (score > threshold):
            new_chainset = add_chain_to_chainset(chainset, chain)
            results.append(make_result_string(new_chainset, freq, score, uniques_a, uniques_b));

    return results



def try_chain_additions_dual(chainset, additions, threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot):
    # Calculate chainset well data

    # Assumes the chainset has exactly 1 alpha and 1 beta. We could generalize this without too much work.
    chainset_a_welldata = chain_data_a[chainset[0][0]]
    chainset_b_welldata = chain_data_b[chainset[1][0]]

    chainset_ab_welldata = intersection(chainset_a_welldata, chainset_b_welldata)

    results = []
    for chain in additions:
        # Retrieve data on chain being added
        chain_welldata = chain_data_a[chain[0][0]] if chainset_num_alphas(chain)>0 else chain_data_b[chain[1][0]]

        # Score probability of chainset+chain
        w_abc = intersection_sum(chainset_ab_welldata, chain_welldata)
        w_ab = sum(chainset_ab_welldata) - w_abc
        w_ac = intersection_sum(chainset_a_welldata, chain_welldata) - w_abc
        w_bc = intersection_sum(chainset_b_welldata, chain_welldata) - w_abc
        w_a = sum(chainset_a_welldata) - w_ac - w_ab - w_abc
        w_b = sum(chainset_b_welldata) - w_bc - w_ab - w_abc
        w_c = sum(chain_welldata) - w_ac - w_bc - w_abc
        score, freq = match_score_dual(w_a, w_b, w_c, w_ab, w_ac, w_bc, w_abc, w_tot)
    
        # Store result if it's good enough
        if (score > threshold):
            new_chainset = add_chain_to_chainset(chainset, chain)
            results.append(make_result_string(new_chainset, freq, score, uniques_a, uniques_b))

    return results

def try_chain_additions(chainset, additions, threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot):
    if chainset_size(chainset)==1:
      return try_chain_additions_nondual(chainset, additions,threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot)
    elif chainset_size(chainset)==2:
      return try_chain_additions_dual(chainset, additions,threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot)
    elif chainset_size(chainset)==3:
      print chainset
      assert False, "***"
      pass


def main(args):
    print "starting... "
    w_tot = int(args[0])
    threshold = float(args[1])
    try_alphas = bool(int(args[2]))
    try_betas = bool(int(args[3]))
    index = int(args[4])
   
    # Input parameters
    print "Starting Process-{} with parameters:".format(index)
    print "  W_tot:", w_tot
    print "  Threshold:", threshold
    print "  try_alphas:", try_alphas
    print "  try_betas:", try_betas
    
   
    # Declare data variables
    fname_init_chainsets = "./solver/initial_" + str(index) + ".txt"
    fname_data_a = "./solver/chain_data_a.txt"
    fname_data_b = "./solver/chain_data_b.txt"
    fname_uniques_a = "./solver/uniques_a.txt"
    fname_uniques_b = "./solver/uniques_b.txt"

    # Load unique chains for a/b
#    print "Loading unique chains..."
    uniques_a = LoadUniques(fname_uniques_a)
    uniques_b = LoadUniques(fname_uniques_b)
#    print "{}/{} unique a/b chains loaded!".format(len(uniques_a), len(uniques_b))
    
    # Load well data for a/b 
#    print "Loading well djata..."
    chain_data_a = LoadChainData(fname_data_a, w_tot)
#    print "Finished loading chain data A!"
    chain_data_b = LoadChainData(fname_data_b, w_tot)
#    print "Finished loading chain data B!"

#    print "Loading starting chainsets from {}:".format(fname_init_chainsets)
    init_chainsets = LoadInitChainsets(fname_init_chainsets)
#    print len(init_chainsets), "initial chain sets loaded!"

    results = []

    # Create vector of chain additions to be tested
    additions = []
    if (try_alphas):
        additions.extend([((i,),()) for i in xrange(len(uniques_a))])
    if (try_betas):
        additions.extend([((),(i,)) for i in xrange(len(uniques_b))])
#    print "trying", len(additions), "additions"
    
    for i, init_chainset in enumerate(init_chainsets):
      results.extend(try_chain_additions(init_chainset, additions, threshold, uniques_a, chain_data_a, uniques_b, chain_data_b, w_tot))
      print "Finished {}/{}\r".format(i+1, len(init_chainsets)),
      sys.stdout.flush()

    # Output results to txt file
    outfile = open("./solver/results_"+str(index) + ".txt", "w")
    for res in results:  outfile.write(res + "\n")

main(sys.argv[1:])
