from seq_data import SequencingData

class SequencingGenerator(object):
  """ Encapsulates parameters for generating sequencing data from a simulated experiment.
  Accepts the following parameters:
    Number of wells
      - num_cells [READ/WRITE]
    Distribution of cells/well
      - cells_per_well_distribution [READ ONLY]
      - cells_per_well_distribution_params [READ ONLY]
      - set_cells_per_well()
    Cell information
      - cell_info [READ/WRITE]
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
    self.cell_info = SequencingGenerator.generate_cell_info(1000, 1, 1)
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
  def cell_info(self):
    return self._cell_info
  @cell_info.setter
  def cell_info(self, val):
    self._cell_info = val[:]

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
    if 'cell_info' in kwargs:
      self.cell_info = kwargs['cell_info']
    if 'cell_frequency_distribution' in kwargs and 'cell_frequency_distribution_params' in kwargs:
      self.set_cell_frequency_distribution(kwargs['cell_frequency_distribution'], **kwargs['cell_frequency_distribution_params'])
    ## TODO: allow setting of noise parameters


  def _update_well_info(self):
    ## TODO
    pass
  def _update_cell_info(self):
    ## TODO
    pass

  @staticmethod
  def generate_cell_info(self, reps, alpha_degree, beta_degree, alpha_start_idx=0, beta_start_idx=0):
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

  def generate_sequencing_data(self):
    ## TODO
    pass


