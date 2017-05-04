from seq_generator import SequencingGenerator
from seq_data import SequencingData

## Create a SequencingGenerator object with all the default values:
#    -num_wells = 96
#    -cells_per_well = 35 (constant)
#    -1000 cells, no shared alpha- or beta-chains
#    -cell frequencies distributed by power law
#    -no sequencing error
gen = SequencingGenerator()
data = gen.generate_data()

## Save data to a file
data.save_data("test.txt")
# Use the following to load existing data:
data2 = SequencingData(path = "test.txt") # Make sure to specify the path using the keyword arg 'path'


## Set the number of wells:
gen.num_wells = 1000

## Set the distribution of cells per well
# Available distributions are:
#  'constant' [parameters: cells_per_well (integer)]
#  'poisson'  [parameters: lam (positive real)]
#  'explicit' [parameters: cells_per_well (list of ints)]
gen.set_cells_per_well(distro_type = 'poisson', lam=5) # poisson distribution with avg of 5 cells/well

## Set the cells in the system
# Creates a list of cells with 100 cells with unique alpha- and beta-chains along with 200 cells with alpha-chains each associated with 2 beta-chains
cells = SequencingGenerator.generate_cells(100) # Generates 100 cells with default distribution of alpha/beta chain sharing
cells = SequencingGenerator.generate_cells(100, alpha_start_idx=100, beta_start_idx=100) # Generates 100 cells with alpha- and beta-chains indexed starting at 100
cells = SequencingGenerator.generate_cells(100, alpha_dual_prob = 0.2) # Generates 100 cells, each with a 20% of having dual alpha-chains
# Set the cells used by the generator
gen.cells = cells

#  also try:
# cells = SequencingGenerator.generate_cells(100, 1, 2) # Generates 200 cells with each beta chain associated with 2 alpha chains
# cells = SequencingGenerator.generate_cells(100, 2, 2) # Generates 400 cells with pairs of alpha chains associated with pairs of beta chains (e.g. a1+b1, a1+b2, a2+b1, a2+b2)


## Set the cell frequency distribution
# Available distributions are:
#  'constant'
#  'power-law' [parameters: alpha] **PVH, I need to talk to you about this, b/c numpy is acting weird --JB
#  'explicit'  [parameters: cells_per_well (list of ints)]
gen.set_cell_frequency_distribution('power-law', alpha=-1)


## Set the error parameters
gen.chain_misplacement_prob = 10**-5 # Prob of a chain migrating to another well
gen.chain_deletion_prob = 10**-3 # Prob of a chain failing to be amplified


## Set multiple options at once:
gen.set_options(
  num_wells = 96,
  cells_per_well_distribution = 'constant',
  cells_per_well_distribution_params = {'cells_per_well': 10}
)


## Generate new data based on new parameters and save to file
data3 = gen.generate_data()

## Access generated well data
# well_data is a list of length <num_wells>, where each element is two lists, the first being the list of alpha-chains and the second being the list of beta-chains
print "Generated data in well 0:", data3.well_data[0]

