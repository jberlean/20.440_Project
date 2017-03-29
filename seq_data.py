import json

class SequencingData(object):
  def __init__(self, well_data=None, metadata=None):
    self.well_data = well_data
    self.metadata = metadata

  def load_data(self, path):
    data = json.load(open(path, 'r'))
    self.well_data = data['well_data']
    self.metadata = data['metadata']
  def save_data(self, path):
    data = {'well_data': self.well_data, 'metadata': self.metadata}
    json.dump(data, open(path, 'w'))

  def get_well_data(self, well_id = None):
    if well_id == None:
      return self.well_data
    else:
      return self.well_data[well_id]

  def get_metadata(self):
    return self.metadata
