

'''

Project: MAD-HYPE
Class(s): (none) 
Function: Organizes main pipeline execution and testing for paired-sequence analysis via 3 methods:

    - MAD-HYPE (Holec, Berleant)
    - ALPHABETR (Lee, et al.)
    - pairSEQ (Howie, et al.)

Author: Patrick V. Holec
Date Created: 6/26/2017
Date Updated: 6/26/2017

'''

'''
Library Importation
'''

print 'Importing libraries...'

# standard libraries
import time
from datetime import datetime

# nonstandard libraries
import numpy as np

# personal libraries
from methods import * # libraries: madhype,alphabetr,pairseq,backup_madhype,madhype_python
from graphers import * # libraries: graphical_tools
from datasim import * # libaries: seq_data,seq_generator

print 'Finished!'


'''
Testing Fucntion
'''

def main():
    test = Testing()
    test.solve()
    test.probe()

    return test


'''
Evaluation class
'''

class Performance:
    def __init__(self,results,data):

        # 
        # Function: use the standardized results/data classes to assess overall performance
        # Results (attributes called):
        #   cell_frequencies,cell_frequencies_CI
        # Data (attributes called):
        #   metadata.cells,well_data,metadata.generated_data.cell_frequencies
        #

        self.pred_cells = results['cells']
        self.cells = data.metadata['cells']
        
        all_alphas, all_betas = zip(*self.cells)
        self.all_alphas,self.all_betas = set([a for alist in all_alphas for a in alist]),set([b for blist in all_betas for b in blist])

        obs_alphas, obs_betas = zip(*data.well_data)
        self.obs_alphas, self.obs_betas = set(sum(obs_alphas, [])), set(sum(obs_betas, []))

        cell_idx_dict = {c: i for i,c in enumerate(self.cells)}

        # Pair data
        self.pairs = [(a,b) for c in self.cells for a in c[0] for b in c[1]]
        pairs_set = set(self.pairs)
        self.pred_pairs = set([(a,b) for c in self.pred_cells for a in c[0] for b in c[1]])
        self.correct_pairs = [p for p in self.pred_pairs if p in pairs_set]
        self.incorrect_pairs = [p for p in self.pred_pairs if p not in pairs_set]

        # Overall cell data
        self.correct_cells = [c for c in self.pred_cells if c in cell_idx_dict]
        self.incorrect_cells = [c for c in self.pred_cells if c not in cell_idx_dict]

        # Non-dual cell data
        is_nondual = lambda c:  len(c[0])==1 and len(c[1])==1
        self.nondual_cells = filter(is_nondual, self.cells)
        self.pred_nondual_cells = filter(is_nondual, self.pred_cells)
        self.correct_nondual_cells = [c for c in self.pred_nondual_cells if c in cell_idx_dict]
        self.incorrect_nondual_cells = [c for c in self.pred_nondual_cells if c not in cell_idx_dict]

        # Dual cell data
        is_dual = lambda c:  not(is_nondual(c))
        self.dual_cells = filter(is_dual, self.cells)
        self.pred_dual_cells = filter(is_dual, self.pred_cells)
        self.correct_dual_cells = [c for c in self.pred_dual_cells if c in cell_idx_dict]
        self.incorrect_dual_cells = [c for c in self.pred_dual_cells if c not in cell_idx_dict]
        
        # Compute stats on frequency predictions
        cell_idxs = [cell_idx_dict[c] if c in cell_idx_dict else -1 for c in self.pred_cells]
        freqs = [data.metadata['generated_data']['cell_frequencies'][i] if i!=-1 else 0.0 for i in cell_idxs]
        pred_freqs = results['cell_frequencies']
        pred_freqs_CI = results['cell_frequencies_CI']

        #self.correct_pairs_percent =  100.*len(self.correct_pairs)/len(self.cells)
        #self.fdr = 100.*len(self.incorrect_pairs)/len(self.pairs)
        self.rmse = np.sqrt(np.mean([(f1-f2)**2 for f1,f2 in zip(freqs, pred_freqs)]))

    def __call__(self):
        print "Solution statistics:"
        print "  Total cells (in system): {} ({} non-dual, {} dual)".format(len(self.cells), len(self.nondual_cells), len(self.dual_cells))
        print "  Number of alpha chains (in system):", len(self.all_alphas)
        print "  Number of beta chains (in system):", len(self.all_betas)
        print "  Number of alpha chains (observed):", len(self.obs_alphas)
        print "  Number of beta chains (observed):", len(self.obs_betas)
        print "  Pair statistics:"
        print "    Total in system:", len(self.pairs)
        print "    Total identified:", len(self.pred_pairs)
        print "    Correctly identified: {} ({}%)".format(len(self.correct_pairs), 100.*len(self.correct_pairs)/len(self.pairs))
        print "    Incorrectly identified:", len(self.incorrect_pairs)
        print "    FDR: {}%".format(100.*len(self.incorrect_pairs)/len(self.pred_pairs))
        print "  Overall cell statistics:"
        print "    Total in system:", len(self.cells)
        print "    Total identified:", len(self.pred_cells)
        print "    Correctly identified: {} ({}%)".format(len(self.correct_cells), 100.*len(self.correct_cells)/len(self.cells))
        print "    Incorrect identified:", len(self.incorrect_cells)
        print "    FDR: {}%".format(100.*len(self.incorrect_cells)/len(self.pred_cells))
        print "    RMSE of frequency guesses:", self.rmse
        if len(self.nondual_cells)>0 and len(self.dual_cells)>0:
            print "  Non-dual cell statistics:"
            print "    Total in system:", len(self.nondual_cells)
            print "    Total identified:", len(self.pred_nondual_cells)
            print "    Correctly identified: {} ({}%)".format(len(self.correct_nondual_cells), 100.*len(self.correct_nondual_cells)/len(self.nondual_cells))
            print "    Incorrectly identified:", len(self.incorrect_nondual_cells)
            if len(self.pred_nondual_cells)>0:  print "    FDR: {}%".format(100.*len(self.incorrect_nondual_cells)/len(self.pred_nondual_cells))
            print "  Dual cell statistics:"
            print "    Total in system:", len(self.dual_cells)
            print "    Total identified:", len(self.pred_dual_cells)
            print "    Correctly identified: {} ({}%)".format(len(self.correct_dual_cells), 100.*len(self.correct_dual_cells)/len(self.dual_cells))
            print "    Incorrectly identified:", len(self.incorrect_dual_cells)
            if len(self.pred_dual_cells)>0:  print "    FDR: {}%".format(100.*len(self.incorrect_dual_cells)/len(self.pred_dual_cells))




# TODO: Convert to class

class Testing:

    def default_model(self):
        # basically every parameter defined in one dictionary
        default_params = {
                         ### Solver Parameters ###
                         'solver_methods':['madhype','madhype_python', 'backup_madhype','alphabetr'],
                         ### Experimental Parameters ###
                         'well_total':96,
                         'cell_per_well_distribution':'constant', # distribution of cells in each well
                         'cell_per_well_total':500, # number of cells per well
                         ### Experimental Noise Parameters ###
                         'chain_misplacement_prob':0.00, # migration to a different well
                         'chain_deletion_prob':0.01, # disappearance of a chain in well
                         ### Repertoire Parameters ###
                         'cell_frequency_distro':'power-law', # can be constant,power-law,Lee,explicit
                         'cell_frequency_alpha':-1, # defining coefficient in the cell frequency distribution
                         'repertoire_cell_total':[2000], # number of cells simulated in repertoire
                         'alpha_chain_total':[1],
                         'beta_chain_total':[1],
                         ### Metadata Parameters
                         'save_filename':None # name of file that instructions are saved under
                          }

        # apply all changes
        self.update_model(default_params)

    # use a dictionary to update class attributes
    def update_model(self,params={}):
        # makes params dictionary onto class attributes
        for key, value in params.items():
            setattr(self, key, value)
            
        # check for inconsistencies in the parameter selection
        assert all(len(x) == len(self.repertoire_cell_total) for x in [self.alpha_chain_total,self.beta_chain_total]), \
                'Not all distribution parameters are the same length'
    
    def __init__(self,params = None):
        
        ### Settings management ###
        
        # set all default parameters
        self.default_model()
        
        # check to see if there is an update
        if params:
            print 'Updating model with new parameters:'
            for key,value in params.iteritems(): print '  - {}: {}'.format(key,value)
            self.update_model(params)
            
            
        ### Generate landscapes ###
        
        gen = seq_generator.SequencingGenerator()

        # Set distribution for cell frequencies
        gen.set_cell_frequency_distribution(self.cell_frequency_distro, alpha=self.cell_frequency_alpha)
        
        gen.chain_misplacement_prob = self.chain_misplacement_prob # Prob of a chain migrating to another well
        gen.chain_deletion_prob = self.chain_deletion_prob # Prob of a chain failing to be amplified
        gen.num_wells = self.well_total # total number of wells

        gen.set_cells_per_well(self.cell_per_well_distribution, cells_per_well=self.cell_per_well_total)

        # iterate through unique declarations of cell packages
        #cells,a_ind,b_ind = [],0,0
        #for r,a,b in zip(self.repertoire_cell_total,self.alpha_chain_total,self.beta_chain_total):
        #    cells += gen.generate_cells(r,alpha_start_idx=a_ind, beta_start_idx=b_ind) # removed sharing #'s (a/b)
        #    a_ind,b_ind = a_ind + a*r,b_ind + b*r
        #gen.cells = cells
        # Use SequencingGenerator's built-in cell generator (for slightly more realistic repertoire)
        gen.cells = seq_generator.SequencingGenerator.generate_cells(sum(self.repertoire_cell_total), alpha_dual_prob=0.3, beta_dual_prob=0.0)

        ## Save data to a file
        self.data = gen.generate_data()
        if self.save_filename: data.save_data(self.save_filename)
            
        
    def solve(self):
       
        results = {}
   
    
        if 'madhype' in self.solver_methods:
            startTime = datetime.now()
            #results['madhype'] = madhype.solve(self.data,pair_threshold=0.999,verbose=0) # not stringent
            results['madhype'] = madhype.solve(self.data, {'pair_threshold':0.}) # not stringent
            print 'MAD-HYPE took {} seconds.\n'.format(datetime.now()-startTime)

        if 'madhype_python' in self.solver_methods:
            startTime = datetime.now()
            results['madhype_python'] = madhype_python.solve(self.data, {'pair_threshold':0.}) # not stringent
            print 'MAD-HYPE (python) took {} seconds.\n'.format(datetime.now()-startTime)

        #if 'backup_madhype' in self.solver_methods:
        #    startTime = datetime.now()
        #    results['backup_madhype'] = backup_madhype.solve(self.data,pair_threshold=0.999,verbose=0) # not stringent
        #    print 'MAD-HYPE took {} seconds.\n'.format(datetime.now()-startTime)

        #if 'alphabetr' in self.solver_methods:
        #    startTime = datetime.now()
        #    results['alphabetr'] = alphabetr.solve(self.data,pair_threshold=0.1) # not stringent
        #    print 'ALPHABETR took {} seconds.\n'.format(datetime.now()-startTime)


        for k,v in results.items():
            print '\n{} performance:'.format(k)
            performance = Performance(v,self.data)
            performance()

        self.results = list(results.values())


    def probe(self):
         
        #graphical_tools.graphical_data_summary(self.data,self.results) 

        graphical_tools.graphical_auroc(self.results[0], self.results[1], self.data)

            

if __name__ == '__main__':

    test = main()


