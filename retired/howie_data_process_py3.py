 
'''
Testing
Function: Figuring out how to process Howie data
'''

""" LIBRARY IMPORTATION """
# standard libraries
import os
import sys
import gzip 
import csv
import pickle as pickle
from collections import Counter

# nonstandard libraries
from datetime import datetime
import matplotlib.pyplot as plt

# homegrown libraries
from methods import * # libraries: madhype,alphabetr,pairseq


""" FACTORY METHODS """

def flatten(l):
    """ Flatten a list into a 1D list """
    return [item for sublist in l for item in sublist]


""" NOTES """
# According to Howie data notes:
# > cDNA refers to sample of origin assignments
# > gDNA refers to repertoire frequency measurements

def main(mode='coast2coast'):

    if mode == 'coast2coast':
        ### directory assignments
        dirnameX,dirnameY = './data/howie/subjectX','./data/howie/subjectY'
        dirname_exp = './data/howie/experiment2' 
        origin_dict = []

        ### analysis on presented data
        filesX,filesY = subjectXYdata(dirnameX,dirnameY) # returns dictionaries
        origin_dict = catalog_repertoire(filesX,filesY,overwrite=False) 
        data = data_assignment(dirname_exp,origin_dict,overwrite=False) # no save due to memory

        ### run analyis (mad-hype)
        #startTime = datetime.now()
        results_madhype = madhype.solve(data,pair_threshold=0.995,verbose=7,real_data=True)
        pickle.dump(results_madhype,open('./pickles/results_{}.p'.format(dirname_exp[-1]),'wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        print('MAD-HYPE took {} seconds.\n'.format(datetime.now()-startTime))
        
        ### process results
        results = pickle.load(open('./pickles/results_{}.p'.format(dirname_exp[-1]),'r'))
        interpret_results(data,results,dirname_exp)
            
    elif mode == 'analysis':
        ### directory assignments
        ### EDIT BELOW ###
        dirnameX,dirnameY = './data/howie/subjectX','./data/howie/subjectY'
        dirname_exp = './data/howie/experiment1' 
        origin_dict = []
        ### EDIT ABOVE ###

        ### process results
        origin_dict = []
        data = data_assignment(dirname_exp,origin_dict,overwrite=False)
        results = pickle.load(open('./pickles/results_{}.p'.format(dirname_exp[-1]),'r'))
        all_results = interpret_results(data,results,dirname_exp,graph=False) # only returns edges above threshold
         
        ### graph features
        all_results 
        
        


def interpret_results(data,results,dirname,graph=True):
    
    # prepare howie results
    with open(dirname + '/tcr_pairseq_fdr1pct.pairs','r') as f: 
        howie_results = [((a[0],a[2]),0.0,-float(a[5])) for i,a in enumerate(csv.reader(f, dialect="excel-tab")) if i != 0]

    # prepare madhype results
    results_ab = [(((a[0]),(a[1])),b,c,'AB') for a,b,c in zip(
        results['AB']['edges'],results['AB']['freqs'],results['AB']['scores'])]
    results_aa = [(((a[0],a[1]),()),b,c,'AA') for a,b,c in zip(
        results['AA']['edges'],results['AA']['freqs'],results['AA']['scores'])]
    results_bb = [(((),(a[0],a[1])),a,b,c,'BB') for a,b,c in zip(
        results['BB']['edges'],results['BB']['freqs'],results['BB']['scores'])]
    all_results = results_ab  + results_aa + results_bb

    # graphing feature, if requested
    if graph:
        # map results using data dictionaries
        xH,yH,tH = map_results(howie_results,data)
        xM,yM,tM = map_results(all_results,data)

        # plot results
        plt.plot(xH['X'],yH['X'],label='Howie (X)')
        plt.plot(xH['Y'],yH['Y'],label='Howie (Y)')
        plt.plot(xM['X'],yM['X'],label='MAD-HYPE (X)')
        plt.plot(xM['Y'],yM['Y'],label='MAD-HYPE (Y)')
        
        plt.xlabel('Number of false positives')
        plt.ylabel('Number of true positives')
        plt.legend()
        plt.show(block=False)
        input('Press enter to continue...')

    return [edge for edge in all_results if edge[2] >= tH] 

def map_results(results,data):
    results.sort(key=lambda x:-float(x[2]))
    print('Total matches:',len(results))
    x,y = {'X':[0],'Y':[0]},{'X':[0],'Y':[0]}
    
    # create reverse dictionary in case of need
    seq_chain_A = dict([(v,k) for k,v in list((data.chain_seq['A']).items())])
    seq_chain_B = dict([(v,k) for k,v in list((data.chain_seq['B']).items())])

    for edge in results:
        if type(edge[0][0]) == int:
            if len(edge) == 4:
                x_origin = data.chain_origin[edge[3][0]][edge[0][0]]
                y_origin = data.chain_origin[edge[3][1]][edge[0][1]]
            else:
                x_origin,y_origin = data.chain_origin['A'][edge[0][0]],data.chain_origin['B'][edge[0][1]]
        elif type(edge[0][0]) == str:
            try:
                #print seq_chain_A.keys()
                x_ind,y_ind = seq_chain_A[edge[0][0]],seq_chain_B[edge[0][1]]
                x_origin,y_origin = data.chain_origin['A'][x_ind],data.chain_origin['B'][y_ind]
            except KeyError:
                continue
        # check if same origin patient (exclusive)
        if x_origin == ['X'] and y_origin == ['X']:
            x['X'].append(x['X'][-1])
            y['X'].append(y['X'][-1]+1)
        if x_origin == ['Y'] and y_origin == ['Y']:
            x['Y'].append(x['Y'][-1])
            y['Y'].append(y['Y'][-1]+1)
        elif x_origin == ['X'] and y_origin == ['Y'] or x_origin == ['Y'] and y_origin == ['X']: 
            x['X'].append(x['X'][-1]+1)
            x['Y'].append(x['Y'][-1]+1)
            y['X'].append(y['X'][-1])
            y['Y'].append(y['Y'][-1])
            if x['X'][-1] == 50: break

    return x,y,edge[2]

def subjectXYdata(dirnameX,dirnameY):
    """ Pulls out file names corresponding the directories submitted, subjects X/Y """
    # initialize dictionaries
    subject_x_files,subject_y_files = {'gdna':{},'cdna':{}},{'gdna':{},'cdna':{}}
    
    # iterate across dirnames and file dictionaries
    for dirname,files in zip([dirnameX,dirnameY],[subject_x_files,subject_y_files]):
        for data_type in ['gdna','cdna']:
            # story directory contents
            dirfiles = os.listdir(dirname+'/{}_data'.format(data_type))
            # check across files for file in dirfiles:
            for file in dirfiles:
                # check if starting file is gzipped
                if file.endswith('.gz') and not any([file[:-3] in d for d in dirfiles if not d.endswith('.gz')]):
                    with gzip.open(dirname+'/{}_data/'.format(data_type)+file,'rb') as f:
                        with open(dirname+'/{}_data/'.format(data_type)+file[:-3],'wb') as f_new:
                            f_new.write(f.read(file[:-3]))
                    print('Unzipped {}.'.format(dirname+'/{}_data/'.format(data_type)+file))
                    files[data_type][file[file.index('TCR')+3] + file[file.index('.well')+5:file.index('.results')]] = \
                            dirname+'/{}_data/'.format(data_type)+file[:-3]

                # otherwise store file location
                elif file.endswith('.tsv'):
                    files[data_type][file[file.index('TCR')+3] + file[file.index('.well')+5:file.index('.results')]] = \
                            dirname+'/{}_data/'.format(data_type)+file
            
    return subject_x_files,subject_y_files

def catalog_repertoire(filesX,filesY,overwrite=False):
    """ Creates a dictionary that takes each unique sequence and assigns ownership between subject X/Y """
    """ This analysis is very slow, so I'm pulling as many tricks out as I can """
    if overwrite == True or not os.path.isfile('./pickles/origin_dict.p'):
        origin_dict = {} 
        # iterate across repertoires to make independent dictionaries
        for patient_id,files in zip(['X','Y'],[filesX,filesY]):
            origin_dict[patient_id] = {}
            for chain_id in ['A','B']: 
                found,keep = set(),[] # create some storage variables
                add,app = found.add,keep.append
                # return 
                for file_id,file in list(files['cdna'].items()): 
                    if not chain_id in file_id: continue # catch to weed out nonmatching chains
                    print('Processing sample {} for patient {}'.format(file_id,patient_id)) 
                    # go through file and pull sequences
                    with open(file,'rb') as f:
                        for line in csv.reader(f, dialect="excel-tab"):
                            if line[0] not in found:
                                add(line[0]) # circumvent attribute lookup
                                app(line[0]) # circumvent attribute lookup
                origin_dict[patient_id][chain_id] = keep 
                pickle.dump(origin_dict,open('./pickles/origin_dict.p','wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        # merge dictionaries
        final_dict = {'A':{},'B':{}}
        for patient_id,chain_dict in list(origin_dict.items()):
            for chain_id,seqs in list(chain_dict.items()):
                print('{} sequences found for chain {}, patient {}'.format(len(seqs),chain_id,patient_id))
                for seq in seqs:
                    try:
                        final_dict[chain_id][seq].append(patient_id)
                    except KeyError:
                        final_dict[chain_id][seq] = [patient_id]
        pickle.dump(final_dict,open('./pickles/origin_dict.p','wb'),protocol=pickle.HIGHEST_PROTOCOL) 
    
    # if there already exists a dictionary that isn't to be overwritten 
    else:
        print('Loading existing origin dictionary...')
        final_dict = pickle.load(open('./pickles/origin_dict.p','r'))
        print('Finished loading!')

    # final counts on chain sequence popularity 
    for chain_id,seq_dict in list(final_dict.items()):
        print('Total sequences for chain {}: {}'.format(chain_id,len(list(seq_dict.keys()))))

    return final_dict

def data_assignment(dirname,origin_dict={},threshold=4,overwrite=True,silent=False):
    """ Pulls out data from experiment directory and assigns TCR wells """
    """ Only valid for two subject testing """
    """ Note: threshold mimicked from Howie """
    # initialize dictionaries
    well_dict = {'A':[],'B':[]} # will hold lists of lists containing chain indices
    chain_2_origin = {'A':{},'B':{}}  # will hold dictionary of chain index mapping to subject origin
    chain_2_sequence = {'A':{},'B':{}} # will hold dictionary of chain index mapping to sequence
    sequence_2_chain = {'A':{},'B':{}} # will hold dictionary of sequence mapping to chain index 
    files = {'A':{},'B':{}}

    # save files
    

    """ Start unpackaging data, create the four dictionaries """ 
    if overwrite == True or not os.path.isfile('./pickles/well_dict_{}.p'.format(dirname[-1])):
        dirfiles = os.listdir(dirname+'/cdna_data/')
        """ Find the file names of all the data, unzip where needed """
        # check across all files, unzip as needed
        for file in dirfiles:
            # check if starting file is gzipped
            if file.endswith('.gz') and not any([file[:-3] in d for d in dirfiles if not d.endswith('.gz')]):
                with gzip.open(dirname+'/cdna_data/'+file,'rb') as f:
                    with open(dirname+'/cdna_data/'+file[:-3],'wb') as f_new:
                        f_new.write(f.read(file[:-3]))
                print('Unzipped {}.'.format(dirname+'/cdna_data/'+file))
                files[file[file.index('TCR')+3]][int(file[file.index('.well')+5: \
                        file.index('.results')])] = dirname+'/cdna_data/'+file[:-3]

            # otherwise store file location
            elif file.endswith('.tsv'):
                files[file[file.index('TCR')+3]][int(file[file.index('.well')+5: \
                        file.index('.results')])] = dirname+'/cdna_data/'+file
        ###
        if not origin_dict:
            print('No origin dictionary provided, exiting...')
            return None
        for chain_id,chain_files in list(files.items()): # iterate across file locations 
            # create a list of all sequences in origin dictionary
            origin_seqs = origin_dict[chain_id] 
            # iterate across repertoires to make independent dictionaries
            for well_id in sorted(list(chain_files.keys()),key=int):
                # go through file and pull sequences
                with open(files[chain_id][well_id],'rb') as f:
                    if not silent: print('Analyzing well {} for chain {}...'.format(well_id,chain_id))
                    well_dict[chain_id].append([])
                    for line in csv.reader(f, dialect="excel-tab"):
                        if line[0] in origin_seqs:    
                            try: # try to assign the chain to a sequence index
                                well_dict[chain_id][-1].append(sequence_2_chain[chain_id][line[0]])
                            except KeyError: # if there isn't an existing index
                                # TODO: change to use origin dict
                                well_dict[chain_id][-1].append(len(chain_2_origin[chain_id]))
                                chain_2_origin[chain_id][len(chain_2_origin[chain_id])] = \
                                        origin_dict[chain_id][line[0]] 
                                chain_2_sequence[chain_id][len(chain_2_sequence[chain_id])] = line[0]
                                sequence_2_chain[chain_id][line[0]] = len(sequence_2_chain[chain_id])


        # remove chains that occur less than threshold 
        print('Adjusting well data for occurance threshold...')
        chain_keys = {'A':[k for k,v in list(Counter(flatten(well_dict['A'])).items()) if v >= threshold],
                      'B':[k for k,v in list(Counter(flatten(well_dict['B'])).items()) if v >= threshold]}
        well_dict['A'] = [[i for i in j if i in chain_keys['A']] for j in well_dict['A']]
        well_dict['B'] = [[i for i in j if i in chain_keys['B']] for j in well_dict['B']]

        if not silent: print('Finished data for chain {}, creating savepoint...'.format(chain_id))
        pickle.dump(well_dict,open('./pickles/well_dict_{}.p'.format(dirname[-1]),'wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        pickle.dump(chain_2_origin,open('./pickles/chain_2_origin_{}.p'.format(dirname[-1]),'wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        pickle.dump(chain_2_sequence,open('./pickles/chain_2_sequence_{}.p'.format(dirname[-1]),'wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        print('Finished saving!')
    
    # if there already exists a dictionary that isn't to be overwritten 
    else:
        print('Loading existing variable names dictionary...')
        well_dict = pickle.load(open('./pickles/well_dict_{}.p'.format(dirname[-1]),'r'))
        chain_2_origin = pickle.load(open('./pickles/chain_2_origin_{}.p'.format(dirname[-1]),'r'))
        chain_2_sequence = pickle.load(open('./pickles/chain_2_sequence_{}.p'.format(dirname[-1]),'r'))
        print('Finished loading!')

    return Results(well_dict,chain_2_origin,chain_2_sequence)
    
class Results:
    """ Quick class to package output data """
    def __init__(self,well_dict,chain_origin,chain_seq):
        self.well_data = [[a,b] for a,b in zip(well_dict['A'],well_dict['B'])]
        self.chain_origin = chain_origin
        self.chain_seq = chain_seq

# script call catch
if __name__ == '__main__':
    if len(sys.argv) > 1:
        print('Starting script under command: {}'.format(sys.argv[1]))
        main(sys.argv[1])
    else:
        main()

