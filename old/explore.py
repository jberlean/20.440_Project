
'''
This file is primarily for exploring data in the Howie et al. study
'''

import os,sys
import time
import cPickle as pickle
import distance

def main():
    params = {'threshold': 25,
              'percent_distance': 0.05,
              'silent': False,
              'dirname': 'subjectX'}
    explore_wells(params)

def numeric_hamming4(num1, num2):
    assert len(num1) == len(num2)
    return numarray.sum(num1 != num2)

def filter_well_sequences(well,params):
    with open(well) as f: lines = f.readlines() #read each line in file

    if len(lines) < 14: return [] # if theres no data in this file

    header = lines[12]
    lines = [x.strip().split('\t') for x in lines[13:]]
    if not params['silent']: print '{} prefiltered unique sequences'.format(len(lines))

    threshold_lines = [x for x in lines if all((int(x[2]) > params['threshold'], 
                                                x[1].startswith('C'),
                                                x[1].endswith(('F','W'))))]
    threshold_lines = [x[0:2]+[int(x[2])]+[x[3:]] for x in threshold_lines]

    if not params['silent']: print '{} unique sequences'.format(len(threshold_lines))
    if len(lines) < 14: return [] # if theres no data in this file

    # fold up low frequency clones
    saved_lines = []
    limit = params['percent_distance']*len(threshold_lines[0][0])
    total_lines = len(threshold_lines)
    for l in xrange(len(threshold_lines)-1,-1,-1):
        print 'Clustering progress: {}%\r'.format((100*(total_lines-l))/total_lines),
        match = False
        for l2 in xrange(l-1,-1,-1):
            if distance.hamming(threshold_lines[l][0],threshold_lines[l2][0]) <= limit:
                threshold_lines[l2][2] += threshold_lines[l][2]
                match = True
                break
        if match == False: saved_lines.append(threshold_lines[l])

    if not params['silent']: 
        print '\n{} clustered unique sequences'.format(len(saved_lines))
        print 'Finished!\n'
    return saved_lines
    

def explore_wells(params):
    ### Acquire paired chains discovered by Howie ###
    # save original directory
    start_dirname = os.getcwd()
    # Get it the right place
    try: os.chdir('./{}/cdna_data'.format(params['dirname']))
    except OSError: raise Exception("Directory doesn't exist!")
    # collect well data files
    all_well_data = [[f for f in os.listdir('./') if f.endswith('.tsv')],
                 [f for f in os.listdir('./') if f.endswith('.tsv (2)')]]
    
    ### Unpack and store sequences (with filters) ###
    sequences_by_chain = []
    
    # threshold parameters
    labels = ['alpha','beta']
    data = {}
    
    for l,well_data in zip(labels,all_well_data):
        sequences_by_well = []
        for i,well in enumerate(well_data):
            print 'Start anaylsis on {}...'.format(i+1)
            sequences_by_well.append(filter_well_sequences(well,params))
        print 'Finished {} chain wells.'
        data[l] = sequences_by_well
        
    storage = {'data':data,'params':params}
    os.chdir(start_dirname)
    
    for i in xrange(1,101):
        fname = '{}_compact_{}.p'.format(params['dirname'],i)
        if not os.path.exists(fname):
            pickle.dump(storage,open(fname,"wb"))
            break
                  
    print 'Saved!'
                
                    
                        
            

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
def explore_experiment1():
    # get into the correct directory
    os.chdir('./experiment1')

    ### Acquire paired chains discovered by Howie ###
    # find lines of real sequences
    with open('tcr_pairseq_fdr1pct.pairs') as f:
        lines = f.readlines()
    lines = [x.strip().split('\t') for x in lines]
    header,lines = lines[0],lines[1:]

    # print first few
    print '{} pairs predicted.'.format(len(lines))

    ### Start pulling out data
    # find all sequences in each well 
    os.chdir('./cdna_data')
    well_data_a = [f for f in os.listdir('./') if f.endswith('.tsv') and
                    f.startswith('TCRA')]
    well_data_b = [f for f in os.listdir('./') if f.endswith('.tsv') and
                    f.startswith('TCRB')]
    print '{} alpha chain well files detected.'.format(len(well_data_a))
    print '{} beta chain well files detected.'.format(len(well_data_b))
    # iterate through files, to find sequence densities
    with open(well_data_a[0]) as f:
        lines = f.readlines()
    lines = [x.strip() for x in lines]
    for i in xrange(5): print lines[i]
    print len(lines)

if __name__ == '__main__':
    main()
