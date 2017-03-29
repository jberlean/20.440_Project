from numpy import random as nprand

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
    TODO: add noise parameters
  """ 

  def __init__(self, **kwargs):
    ## Set default parameter values
    self.num_wells = 96
    self.set_cells_per_well(distro_type = 'constant', cells_per_well = 1)
    self.cells = SequencingGenerator.generate_cells(1000, 1, 1)
    self.set_cell_frequency_distribution(distro_type = 'power-law', alpha = -1)


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

  ## TODO
  

  def set_options(self, **kwargs):
    if 'num_wells' in kwargs:
      self.num_wells = kwargs['num_wells']
    if 'cells_per_well_distribution' in kwargs and 'cells_per_well_distribution_params' in kwargs:
      self.set_cells_per_well_distribution(kwargs['cells_per_well_distribution'], **kwargs['cells_per_well_distribution_params'])
    if 'cells' in kwargs:
      self.cells = kwargs['cells']
    if 'cell_frequency_distribution' in kwargs and 'cell_frequency_distribution_params' in kwargs:
      self.set_cell_frequency_distribution(kwargs['cell_frequency_distribution'], **kwargs['cell_frequency_distribution_params'])
    ## TODO: allow setting of noise parameters


  @staticmethod
  def generate_cells(self, reps, alpha_degree, beta_degree, alpha_start_idx=0, beta_start_idx=0):
    # alpha_degree: number of beta chains associated with each alpha chain
    # beta_degree: number of alpha chains associated with each beta chain
    # reps: number of iterations to run (total cells created = alpha_degree*beta_degree*reps)

    # TODO: explain what is going on
    cells = []
    for i in range(reps):
      alpha_indices = xrange(alpha_start_idx, alpha_start_idx+beta_degree)
      beta_indices = xrange(beta_start_idx, beta_start_idx+alpha_degree)
      cells.extend(it.product(alpha_indices, beta_indices))
    return cells

  def _sample_cells_per_well(self):
    # Generate a list of number of cells for each well, based on the specified distribution
    # This is called to generate a new list each time generate_sequencing_data() is called.
    distro = self.cells_per_well_distribution
    params = self.cells_per_well_distribution_params
    if distro == 'constant':
      return [params['cells_per_well']] * self.num_wells
    elif distro == 'poisson':
      return list(nprand.poisson(params['lambda'], self.num_wells))
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
      freqs = nprand.power(params['alpha']+1, len(self.cells))
    elif distro == 'explicit':
      freqs = np.array(params['frequencies'])
    else:
      assert False, "Unknown distribution of cell frequencies: {0}".format(distro)

    freqs = freqs / np.sum(freqs) # Normalize freqs so it is a probability distro
    return freqs

  def generate_sequencing_data(self, path):
    cells_per_well = self._sample_cells_per_well()
    cell_freqs = self._sample_cell_freqs()
    
    well_data = []
    for cpw in cells_per_well:
      # Pick cpw cells based on distro in cell_freqs
      # TODO: use trees to do this in O(log n) time?
      # (fyi: i actually don't know the efficiency of numpy's algorithm)
      cells = nprand.choice(self.cells, cpw, replace=True, p=cell_freqs)

      # Extract alpha and beta chains in the well
      alphas, betas = zip(*cells)
      # Remove duplicate chains and add to well_data
      well_data.append([sorted(set(alphas)), sorted(set(betas))])

    metadata = {
      'num_wells': self.num_wells,
      'cells_per_well_distribution': self.cells_per_well_distribution,
      'cells_per_well_distribution_params': self.cells_per_well_distribution_params,
      'cells': self.cells,
      'cell_frequency_distribution': self.cell_frequency_distribution,
      'cell_frequency_distribution_params': self.cell_frequency_distribution_params,
      'generated_data': {
        'cells_per_well': cells_per_well,
        'cell_frequencies': cell_freqs
      }
    }
    seq_data = SequencingData(well_data = well_data, metadata = metadata)
    seq_data.save(path)
