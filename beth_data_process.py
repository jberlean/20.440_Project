
'''
Data formatting and interpretation
MAD-HYPE Trial data
'''



# standard libraries
import os
from collections import Counter

# nonstandard libraries
from datetime import datetime

# homegrown libraries
from methods import * # libraries: madhype,alphabetr,pairseq



# data characteristics
def main():
    dir_name = './plate-data'
    data = Data(dir_name)
    data.process_files()

    # actual simulation
    startTime = datetime.now()
    new_results = madhype.solve(data,pair_threshold=0.995,verbose=2,real_data=True) # not stringent
    print 'MAD-HYPE took {} seconds.\n'.format(datetime.now()-startTime)

class Data:
    # TODO: make a dictionary passing feature to assign multiple parameters
    def __init__(self,dir_name,params={}):
        self.dir_name = dir_name 
        # default parameters
        self.empty_wells = ['A3','B4','C5','D6','E7','F8','G9','H10']
        self.except_wells = ['D4','G10']
        self.count_threshold = 4 # each sequence must have atleast this number of read counts
        self.well_threshold = 3 # each sequence must appear in this number of wells to be analyzed

    def process_files(self):
        # modify well names to slip in a C
        empty_wells,except_wells = [e[0]+'C'+e[1]+'.txt' for e in self.empty_wells],[e[0]+'C'+e[1]+'.txt' for e in self.except_wells]

        # access all text files in directory
        text = {'bg_wells':[],'wells':[]}
        for file in os.listdir(self.dir_name):
            if file.endswith('.txt'):
                if any([e in file for e in empty_wells]):
                    print 'Considering this well as background',file
                    with open(os.path.join(self.dir_name,file),'rb') as f:
                        data = [line.split('\t') for line in f.read().split('\r\n')]
                        text['bg_wells'].append(data)
                elif any([e in file for e in except_wells]):
                    print 'Skipping file',file
                    continue
                else:
                    with open(os.path.join(self.dir_name,file),'rb') as f:
                        data = [line.split('\t') for line in f.read().split('\r\n')]
                        text['wells'].append(data)

        # process data to acquire counts (real wells)
        data_alpha,data_beta = {},{}
        for i,d in enumerate(text['wells']): # check for each well file in text
            data_alpha[i],data_beta[i]  = {},{} # set up dictionary for each well
            for l in d[1:]: 

                # filtering and thresholding
                if len(l) != 11: continue
                if int(l[0]) < self.count_threshold: continue

                # add to alpha dictionary
                if any(['A' in s for s in l[4:7]]) and all([not 'B' in s for s in l[4:7]]): # had no beta IDs and atleast one alpha ID
                    if l[3] in data_alpha[i].keys():
                        data_alpha[i][l[3]] += int(l[0]) # adds to sequences containing identical CDR3 regions 
                    else:
                        data_alpha[i][l[3]] = int(l[0]) # creates new entry corresponding to a CDR3 region

                # add to beta dictionary
                if any(['B' in s for s in l[4:7]]) and all([not 'A' in s for s in l[4:7]]): # had no alpha IDs and atleast one beta ID
                    if l[3] in data_beta[i].keys():
                        data_beta[i][l[3]] += int(l[0]) # adds to sequences containing identical CDR3 regions 
                    else:
                        data_beta[i][l[3]] = int(l[0]) # creates new entry corresponding to a CDR3 region

        # generate quick histogram of frequency of clones across wells
        all_seqs_alpha = [a for b in [d[1].keys() for d in data_alpha.items()] for a in b]
        all_seqs_beta  = [a for b in [d[1].keys() for d in data_beta.items() ] for a in b]
        seqs_alpha = [a[0] for a in Counter(all_seqs_alpha).items() if a[1] >= self.well_threshold]
        seqs_beta= [a[0] for a in Counter(all_seqs_beta).items()  if a[1] >= self.well_threshold] 

        # generate lists of object presense
        well_data_alpha = [[seqs_alpha.index(item[0]) for item in well[1].items() if item[0] in seqs_alpha] for well in data_alpha.items()]
        well_data_beta =  [[seqs_beta.index(item[0])  for item in well[1].items() if item[0] in seqs_beta ] for well in data_beta.items() ]
        
        # save important attributes to class
        self.well_data = [[a,b] for a,b in zip(well_data_alpha,well_data_beta)]
        self.seqs_alpha = seqs_alpha
        self.seqs_beta  = seqs_beta

# namespace catch
if __name__ == '__main__':
    main()
