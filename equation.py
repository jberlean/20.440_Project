# This program is intended to solve for probability two chains match

import random
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

# true population characteristics
# per well basis
f_ax_true = 0.05
f_xb_true = 0.05
f_ab_true = 0.05

w = 1000

ab_wells = []

for i in xrange(w):
     rands = [random.random() for r in xrange(3)]
     a,b = rands[0] < f_ax_true or rands[2] < f_ab_true,rands[1] < f_xb_true or rands[2] < f_ab_true
     ab_wells.append((a,b))

w_a = sum([1 for ab in ab_wells if ab[0] == True])
w_b = sum([1 for ab in ab_wells if ab[1] == True])
w_ab = sum([1 for ab in ab_wells if ab == (True,True)])

f_a = float(w_a)/w
f_b = float(w_b)/w

p_ab = []
f_ab = [(float(w_ab)/w)*float(i)/100 for i in xrange(101)]

for f in f_ab:
    f_ax,f_xb = (f_a-f)/(1-f),(f_b-f)/(1-f)
    a,b,c = np.arange(0,w_ab+1),np.arange(w_ab,-1,-1),np.arange(w,w-(w_ab+1),-1)
    p = [stats.binom.pmf(a[i],w,f)*stats.binom.pmf(b[i],c[i],f_ax*f_xb) for i in xrange(w_ab+1)]
    p_ab.append(sum(p))

#for f,p in zip(f_ab,p_ab): print f,p

print 'Wells:',w
print 'Wells a:',w_a
print 'Wells b:',w_b
print 'Wells ab:',w_ab
print 'Max(P(f_ab = {})) = {}'.format(f_ab[p_ab.index(max(p_ab))],max(p_ab))

plt.plot(f_ab,p_ab)


plt.show(block=False)
raw_input('Hit enter to close...')
plt.close()
