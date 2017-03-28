import operator as op
import matplotlib.pyplot as plt
import numpy as np


def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


w = 1000  # number of wells
n = 35  # number of cells/well
fs = [10**(-f/5.) for f in xrange(2,31)]
#fs = [float(f)/2000 for f in xrange(1,21)]
#fs = [.05]
alpha = 1

# generate distribution
occurs = xrange(1,20001)
p = [o**-alpha for o in occurs]
p = [i/sum(p) for i in p]
p_present = [1 - (1-i)**n for i in p]
fs_present = [1 - (1-i)**n for i in fs]

# results

#misses2 = []
#misses3 = []

results = []

for f,f_orig in zip(fs_present,fs):
    misses = []

    for i in xrange(0,w+1): # cut out not appearing case and one appearance case (noise)
        misses.append((1 - reduce(op.mul,[1 - ((f1**i)*((1-f1)**(w-i))) for f1 in p_present],1))*ncr(w,i)*(f**i)*((1-f)**(w-i)))
        #misses2.append((1 - reduce(op.mul,[1 - ((f1**i)*((1-f1)**(w-i))) for f1 in p_present],1)))
        #misses3.append(ncr(w,i)*(f**i)*((1-f)**(w-i)))
        
    #misses.append(ncr(w,0)*((1-f)**(w)))

    print 'Added clone with frequency {}%: Failure rate of {}%'.format(round(100*f_orig,4),round(100*sum(misses),4))     
    results.append(sum(misses))
    
'''
plt.plot(xrange(2,2+len(misses2)), misses2)
plt.xlabel('Number of well appearances (k)')
plt.ylabel('Probability of match occurring')
plt.title('Match distribution for f = {} as a function of k'.format(f))
print 'Probability of failure: {}'.format(sum(misses))
'''
plt.semilogx(fs,results)
plt.xlabel('Frequency of target clone')
plt.ylabel('Probability of convoluted signals')
plt.title('Probability of convoluted signals for n = {}'.format(n))
#'''
plt.show()
