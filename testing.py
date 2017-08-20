
from solver_Lee import solve as lee_solve
from solver_440 import solve as new_solve
from seq_generator import SequencingGenerator as SeqGen
from seq_data import SequencingData
import numpy as np
from datetime import datetime
import graphical_tools

print 'Finished importation!'

'''
Generate landscapes
'''

gen = SeqGen()

'''
Set distributions
'''

w_tot = 480

#gen.set_cell_frequency_distribution('Lee')
gen.set_cell_frequency_distribution('power-law', alpha=-1)
gen.set_options(chain_misplacement_prob=0,chain_deletion_prob=10**-1)
#gen.chain_misplacement_prob = 10**5 # Prob of a chain migrating to another well
#gen.chain_deletion_prob = 10**5 # Prob of a chain failing to be amplified

gen.num_wells = w_tot

# Make SequencingGenerator object
gen.set_cells_per_well('constant', cells_per_well=50)
#gen.set_cells_per_well('poisson', lam=100)

#gen.cells = SeqGen.generate_cells_lee(600)
gen.cells = SeqGen.generate_cells(500)
#cells += SeqGen.generate_cells(250, 2, 1, alpha_start_idx=250, beta_start_idx=251) 
#cells += SeqGen.generate_cells(250, 1, 2, alpha_start_idx=2000, beta_start_idx=2000) 
#gen.cells = cells

## Save data to a file
data = gen.generate_data()
data.save_data('patrick_testing.txt')

print 'Finished!'


startTime = datetime.now()
new_results = new_solve(data,pair_threshold=0.00001,verbose=2) # not stringent
print 'PHASOR took {} seconds.\n'.format(datetime.now()-startTime)

startTime = datetime.now()
lee_results = lee_solve(data,pair_threshold=0.000001) # not stringent
print 'Lee took {} seconds.\n'.format(datetime.now()-startTime)

print '\nPHASOR Performance:'
graphical_tools.assess_performance(new_results,data)

print '\nLEE Performance:'
graphical_tools.assess_performance(lee_results,data)


graphical_tools.graphical_frequencies(new_results,data,True)
graphical_tools.graphical_frequencies(lee_results,data,False)

#graphical_tools.graphical_auroc(new_results,lee_results,data.metadata['cells'])
#graphical_tools.graphical_network(data,new_results,lee_results)


