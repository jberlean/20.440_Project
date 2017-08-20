

"""

MAD-HYPE : Runtime Interface

Creator: Patrick V. Holec
Created: August 10th, 2017
Updated: August 10th, 2017

"""

# standard libraries
from time import sleep
from csv import reader

# nonstandard libraries
from Tkinter import *


class RuntimeInterface:
    def __init__(self, master, params={}):
        self.master = master
        update = Toplevel(Frame(master))
        update.wm_title("Progress Update")

        # pass params dictionary into class attributes
        for k,v in params.items():
            setattr(self,k,v)

        self.progress_text = ''
        self.progress_label = StringVar() # string variable for progress
        self.progress = Label(update, textvariable=self.progress_label) # assign label 
        self.progress.pack(side="top",fill="both",expand=True,padx=50,pady=50)
        
        self.add_output('Importing barcodes...')
        bc_dict = load_barcodes(self.barcode_fname)
        self.add_output(bc_dict['warnings'])
        self.add_output('Finished importing!')
        print 'Results:',bc_dict
        self.add_output('Something else again...')
        


    def add_output(self,new_text):
        """ Adds a str/list of text to the output window """ 
        if not new_text: pass 
        elif type(new_text) == str:
            self.progress_text += '\n' + new_text
        elif type(new_text) == list:
            self.progress_text += '\n'+'\n'.join(new_text)
        self.progress_label.set(self.progress_text)
        self.master.update_idletasks() 

    def clear_output(self,new_text):
        """ Clears the output window """
        self.progress_text = ''
        self.progress_label.set(self.progress_text)
        self.master.update_idletasks() 



def check_dna(my_str):
    if not my_str: return False
    return all([a in 'autcgAUTCG' for a in my_str])

def load_barcodes(fname):

    """ Loads barcodes in template form into MAD-HYPE """

    bc_dict = {'A':{},'B':{},'warnings':[],'crash':False}

    with open(fname,'rb') as f:

        r = reader(f,delimiter='\t')
        fixed_regions = ["A5 Fixed", "A3 Fixed", "B5 Fixed", "B3 Fixed"] # fixed regions for id
        switch = False # switches from fixed region declaration to barcodes

        # assign barcodes to appropriate well/chain dictionary slots
        for i,line in enumerate(r):
            
            # reading the conservation section
            if not switch: # skip the first line

                # look for conserved region identifiers 
                print line[0]
                print fixed_regions
                print line[0] in fixed_regions
                if line[0] in fixed_regions:
                    bc_dict[line[0]] = line[1] 
                    fixed_regions.remove(line[0])
               
                if line[0] == 'BARCODING SECTION':
                    switch = True
                    continue

            # reading the barcoding section
            else:
                if line[0] == 'Well': continue # skip the first line
                elif line[1] == '':
                    bc_dict['A'][line[0]] = [l.upper() for l in line[2:] if check_dna(l)]
                    bc_dict['B'][line[0]] = [l.upper() for l in line[2:] if check_dna(l)]
                elif line[1].upper() in 'AB': # check for chain identifiers (case-insensitive)
                    bc_dict[line[1].upper()][line[0]] = [l.upper() for l in line[2:] if check_dna(l)]
                else:
                    bc_dict['warnings'].append(
                            'WARNING: Well {} chain name not found ({})...'.format(line[0],line[1]))

        # delete wells that do not have A/B components
        if not set(bc_dict['A'].keys()) == set(bc_dict['B'].keys()):
            bc_dict['warnings'].append('Some wells lack both A/B chains, deleting...')
            set_a,set_b = set(bc_dict['A']),set(bc_dict['B']) 
            del_a,del_b = [a for a in set_a if not a in set_b],[b for b in set_b if b not in set_a]
            for d in del_a: del d['A'][d]
            for d in del_b: del d['B'][d]
        
        # raise warnings if certain fixed regions are missed
        if len(fixed_regions) != 0:
            for region in fixed_regions:
                bc_dict['warnings'].append(
                        'ERROR: Conserved region ({}) not found in file.'.format(region))
                bc_dict['crash'] = True 
    
        # raise warnings if no barcoding section was found
        if not switch:
            bc_dict['warnings'].append(
                    'ERROR: No BARCODING SECTION line was found, barcodes not added.')
            bc_dict['crash'] = True
    
    # return resulting dictionary
    return bc_dict
            
            
def load_sequences(fname,bc_dict): 
    # identify file format
    if fname.endswith('.fasta'): marker = '>'
    if fname.endswith('.fastq'): marker = '@'
       
    # identify barcodes
    
    for seq in fasta_iter(fname):
        print 


""" Helper Functions """

def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq






