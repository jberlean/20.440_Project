import itertools as it
import numpy as np

from seq_data import SequencingData

class SequencingGenerator(object):
  """ Encapsulates parameters for generating sequencing data from a simulated experiment.
  Accepts the following parameters:
    Number of wells
      - num_wells [READ/WRITE]
    Distribution of cells/well
      - cells_per_well_distribution [READ ONLY]
      - cells_per_well_distribution_params [READ ONLY]
      - set_cells_per_well()
    Cell information
      - cells [READ/WRITE]
    Cell frequency distribution
      - cell_frequency_distribution [READ ONLY]
      - cell_frequency_distribution_params [READ ONLY]
      - set_cell_frequency_distribution()
    Noise stuff
      - chain_misplacement_prob [READ/WRITE]
  """ 

  def __init__(self, **kwargs):

    ## Set default parameter values
    self.num_wells = 96
    self.set_cells_per_well(distro_type = 'constant', cells_per_well = 35)
    self.cells = SequencingGenerator.generate_cells(1000, 1, 1)
    self.set_cell_frequency_distribution(distro_type = 'power-law', alpha = -1)
    self.chain_misplacement_prob = 0
    self.chain_deletion_prob = 0

    ## If any parameters are specified in kwargs, modify them from the defaults
    self.set_options(**kwargs)


  ## Well parameters
  
  # Number of wells
  @property
  def num_wells(self):
    return self._num_wells
  @num_wells.setter
  def num_wells(self, val):
    self._num_wells = val

  # Distribution of cells/well
  @property
  def cells_per_well_distribution(self):
    return self._cpw_distro
  @property
  def cells_per_well_distribution_params(self):
    return self._cpw_params
  @property
  def cells_per_well(self):
    return self._cpw
  def set_cells_per_well(self, distro_type, **distro_params):
    self._cpw_distro = distro_type
    self._cpw_params = distro_params


  ## Cell parameters
  
  # Cell info - list of (alpha_id, beta_id) representing distinct cell types
  @property
  def cells(self):
    return self._cells
  @cells.setter
  def cells(self, val):
    self._cells = val[:]

  # Cell frequencies
  @property
  def cell_frequency_distribution(self):
    return self._cell_freq_distro
  @property
  def cell_frequency_distribution_params(self):
    return self._cell_freq_params
  def set_cell_frequency_distribution(self, distro_type, **distro_params):
    self._cell_freq_distro = distro_type
    self._cell_freq_params = distro_params


  ## Noise parameters
  # Probability that a chain will relocate to a random other well
  @property
  def chain_misplacement_prob(self):
    return self._cmp
  @chain_misplacement_prob.setter
  def chain_misplacement_prob(self, prob):
    self._cmp = prob
  # Probability that a chain will not be amplified properly (i.e. deleted)
  @property
  def chain_deletion_prob(self):
    return self._cdp
  @chain_deletion_prob.setter
  def chain_deletion_prob(self, prob):
    self._cdp = prob
  

  def set_options(self, **kwargs):
    if 'num_wells' in kwargs:
      self.num_wells = kwargs['num_wells']
    if 'cells_per_well_distribution' in kwargs and 'cells_per_well_distribution_params' in kwargs:
      self.set_cells_per_well(kwargs['cells_per_well_distribution'], **kwargs['cells_per_well_distribution_params'])
    if 'cells' in kwargs:
      self.cells = kwargs['cells']
    if 'cell_frequency_distribution' in kwargs and 'cell_frequency_distribution_params' in kwargs:
      self.set_cell_frequency_distribution(kwargs['cell_frequency_distribution'], **kwargs['cell_frequency_distribution_params'])
    if 'chain_misplacement_prob' in kwargs:
      self.chain_misplacement_prob = kwargs['chain_misplacement_prob']
    if 'chain_deletion_prob' in kwargs:
      self.chain_deletion_prob = kwargs['chain_deletion_prob']


  @staticmethod
  def generate_cells(reps, alpha_degree, beta_degree, alpha_start_idx=0, beta_start_idx=0):
    # alpha_degree: number of beta chains associated with each alpha chain
    # beta_degree: number of alpha chains associated with each beta chain
    # reps: number of iterations to run (total cells created = alpha_degree*beta_degree*reps)

    # TODO: explain what is going on
    cells = []
    for i in range(reps):
      alpha_indices = xrange(alpha_start_idx, alpha_start_idx+beta_degree)
      beta_indices = xrange(beta_start_idx, beta_start_idx+alpha_degree)
      cells.extend(it.product(alpha_indices, beta_indices))

      alpha_start_idx += beta_degree
      beta_start_idx += alpha_degree
    return cells

  def _sample_cells_per_well(self):
    # Generate a list of number of cells for each well, based on the specified distribution
    # This is called to generate a new list each time generate_sequencing_data() is called.
    distro = self.cells_per_well_distribution
    params = self.cells_per_well_distribution_params
    if distro == 'constant':
      return [params['cells_per_well']] * self.num_wells
    elif distro == 'poisson':
      return list(np.random.poisson(params['lam'], self.num_wells))
    elif distro == 'explicit':
      return params['cells_per_well']
    else:
      assert False, "Unknown distribution of cells/well: {0}".format(distro)
  def _sample_cell_freqs(self):
    # Generate a list of cell frequencies based on the specified distribution
    # This is called each time generate_sequencing_data() is called.
    distro = self.cell_frequency_distribution
    params = self.cell_frequency_distribution_params
    if distro == 'constant':
      freqs = np.array([1]*len(self.cells))
    elif distro == 'power-law':
      freqs = np.random.pareto(-params['alpha'], len(self.cells)) ## TODO: there's something screwy abt this distro, talk to holec abt it
    elif distro == 'explicit':
      freqs = np.array(params['frequencies'])
    else:
      assert False, "Unknown distribution of cell frequencies: {0}".format(distro)

    freqs = freqs / np.sum(freqs) # Normalize freqs so it is a probability distro
    return list(freqs)

  def generate_data(self):
    # Generates sequencing data based on this SequencingGenerator object's parameter values.
    # Results are returned in a SequencingData object, which can be saved to a file with seq_data.save()

    cells_per_well = self._sample_cells_per_well()
    cell_freqs = self._sample_cell_freqs()

    misplaced_alphas = []
    misplaced_betas = []
    cells_per_well_idx = []
    well_data = []
    for cpw in cells_per_well:
      # Pick cpw cells based on distro in cell_freqs
      # TODO: use trees to do this in O(log n) time?
      # (tbh: i actually don't know the efficiency of numpy's algorithm)
      cells_idx = np.random.choice(range(len(self.cells)), cpw, replace=True, p=cell_freqs)
      cells = [self.cells[idx] for idx in cells_idx]

      # Extract alpha and beta chains in the well
      alphas, betas = zip(*cells)
      alphas, betas = list(alphas), list(betas)

      # Determine if any chains are misplaced
      # If so, move them to the appropriate list of misplaced chains
      for a in alphas:
        if np.random.uniform() < self.chain_deletion_prob:
          alphas.remove(a)
        elif np.random.uniform() < self.chain_misplacement_prob:
          misplaced_alphas.append(a)
          alphas.remove(a)
      for b in betas:
        if np.random.uniform() < self.chain_deletion_prob:
          betas.remove(b)
        elif np.random.uniform() < self.chain_misplacement_prob:
          misplaced_betas.append(b)
          betas.remove(b)

      # Remove duplicate chains and add to well_data
      well_data.append([sorted(set(alphas)), sorted(set(betas))])

      # Store actual cells in well for metadata
      cells_per_well_idx.append(list(cells_idx))

    # Put misplaced chains in random wells
    for a in misplaced_alphas:  well_data[np.random.randint(0, len(well_data))][0].append(a)
    for b in misplaced_betas:  well_data[np.random.randint(0, len(well_data))][1].append(b)
      
    metadata = {
      'num_wells': self.num_wells,
      'cells_per_well_distribution': self.cells_per_well_distribution,
      'cells_per_well_distribution_params': self.cells_per_well_distribution_params,
      'cells': self.cells,
      'cell_frequency_distribution': self.cell_frequency_distribution,
      'cell_frequency_distribution_params': self.cell_frequency_distribution_params,
      'generated_data': {
        'cells_per_well': cells_per_well,
        'cells_per_well_idx': cells_per_well_idx,
        'cell_frequencies': cell_freqs
      }
    }
    seq_data = SequencingData(well_data = well_data, metadata = metadata)
    return seq_data
