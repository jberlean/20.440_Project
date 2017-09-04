
# Overall results

core_count = 8
lines = []

# get all the results into one list
for i in xrange(core_count):
    with open('results_{}.txt'.format(i+1),'r') as f:
        lines += [l.strip().split('\t') for l in f.readlines()]

lines = [[(int(l[2]),int(l[3])),float(l[0]),float(l[1])] for l in lines] # fix data types

# safe in clean form 
ab_edges = [(int(l[2]),int(l[3])) for l in lines]
ab_freqs = [float(l[1]) for l in lines]
ab_scores = [float(l[0]) for l in lines]

results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}


    
