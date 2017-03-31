
'''
Project (working title): 
    Multicell sequence decoupling
Function: 
    Simulates the ability to decouple information about a population of paired chain sequences using combinatorics
Purpose: 
    Joseph and Patrick 440 project
'''

'''
Library imports and initialization
'''

# standard libraries
import operator as op

# nonstandard libraries
import matplotlib.pyplot as plt
import numpy as np

'''
Main function execution (if file not imported, this will run)
'''

def main():
    mcs = Multicell_Sequencing(w=1000,n=100)
    mcs.build_repertoire(alpha=1.0,count=20000)
    mcs.fs_logspace(-9,-1,40)
    fs,results = mcs.start_test()

    print 'Clone repertoires :',mcs.p[:10]
    print 'Clone repertoires:',mcs.p[-10:]

    visualize_1D(fs,results)
    visualize_repertoire(mcs.p) 

'''
Factory methods
'''

# N choose K evaluation, with approximated results for N >> 1
def ncr_exact(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

# approximation of binomial distribution
def ncr_approx(n, r, f):
    avg,var = n*f,n*f*(1-f)
    return (1/((2*var*np.pi)**0.5))*np.exp(-((r-avg)**2)/(2*var))

# creates a series of numbers representing binomial distribution (approx. and exact)
def binomial_list(w,f):
    if False: #w*f > 50 and w*(1-f) > 50: # condition for inaccurate approximation
        temp = [ncr_approx(n,i,f) for i in xrange(0,w+1)]
        return [t/sum(temp) for t in temp]
    else:
        return [ncr_exact(w,i)*(f**i)*((1-f)**(w-i)) for i in xrange(0,w+1)]

# probability a permutation of i choose w elements is met (frequency f)
def perm_match(w,f,i):
    return (f**i)*((1-f)**(w-i))


'''
Class: Multicell Sequencing Analysis
Function: Test for well properties of particular parameter sets
'''

class Multicell_Sequencing:
    def __init__(self,w=1000,n=1000):
        self.w,self.n = w,n
        self.p,self.fs = None,None
    
    def build_repertoire(self,alpha=1,count=20000):
        p = [(o+1)**-alpha for o in xrange(count)]
        self.p = [i/sum(p) for i in p]

    def fs_linspace(self,start,stop,num=5):
        self.fs = np.linspace(start,stop,num)

    def fs_logspace(self,start,stop,num=5):
        self.fs = np.logspace(start,stop,num) 

    def start_test(self,silent=True):
        # check to make sure important sets are define
        if self.p is None:
            print 'Repertoire not defined!'
            return None
        elif self.fs is None:
            print 'Clone frequencies not defined!'
            return None 
        
        self.p_present = [1 - (1-i)**self.n for i in self.p]
        self.fs_present = [1 - (1-i)**self.n for i in self.fs]

        # prep for simulation
        results = []
        repertoire_perm_match = [(1 - reduce(op.mul,[1 - perm_match(self.w,f1,i) for f1 in self.p_present],1)) for i in xrange(0,self.w+1)] # p_present modification

        # main simulation loop
        for f,f_orig in zip(self.fs_present,self.fs):
            binomial = binomial_list(self.w,f)
            results.append(sum([a*b for a,b in zip(binomial,repertoire_perm_match)]))
            if not silent: print 'Added clone with frequency {}%: Failure rate of {}%'.format(round(100*f_orig,4),round(100*results[-1],4))     

        return self.fs,results 

'''
Visualization methods (maybe put this into a class eventually
'''

def visualize_1D(fs,results):
    plt.semilogx(fs,results)
    plt.xlabel('Frequency of target clone')
    plt.ylabel('Probability of convoluted signals')
    plt.show(block=False)
    raw_input('Hit enter to close...')
    plt.close()

def visualize_repertoire(p):
    plt.loglog(xrange(1,len(p)+1),p)
    plt.title('Repertoire clonal frequencies')
    plt.xlabel('Clone Index (#)')
    plt.ylabel('Frequency of Clone (% total)')
    plt.show(block=False)
    raw_input('Hit enter to close...')
    plt.close()

'''
Namespace main execution
'''

if __name__ == '__main__':
    main()


