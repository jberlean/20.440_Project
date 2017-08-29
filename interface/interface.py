
"""

MAD-HYPE : User Interface

Creator: Patrick V. Holec
Created: August 8th, 2017
Updated: August 10th, 2017

"""

# standard libraries
import os
import time
import random

# nonstandard libraries
from Tkinter import Tk, Frame, Label, Button, Entry, StringVar, DISABLED, NORMAL, END, N, E, W, S 
from Tkinter import LEFT, IntVar, Checkbutton 
from Tkinter import *
import tkMessageBox

# homegrown libraries
from madhype import *

class UserInterface:
    def __init__(self, master):

        self.master = master
        self.error_message = ""
        #Frame.__init__(Frame(master))
        #self.frame = Frame(master)

        sgr1,sgr2 = 4,9 # stagger for menu

        master.title("MAD-HYPE")

        for i in xrange(1,3): master.grid_columnconfigure(i, weight=1, uniform='foo')

        self.solver = IntVar()
        self.silent= IntVar()
        self.label_text = StringVar()
        vcmd = master.register(self.validate)

        self.label_text.set(self.error_message)
        self.label = Label(master, textvariable=self.label_text, fg='red')

        self.locations_text = Label(master,text='File Locations')
        self.fasta_text = Label(master,text='Sequence file:')
        self.barcode_text = Label(master,text='Barcode file:')
        self.settings_text = Label(master,text='Settings')
        self.fdr_text = Label(master,text='FDR (%):')
        self.error_text= Label(master,text='Errors:')

        self.barcode_entry = Entry(master)
        self.fasta_entry = Entry(master)
        self.fdr_entry = Entry(master, validate="key", validatecommand=(vcmd, '%P'))

        self.barcode_entry.insert(END, os.getcwd() + '/my_barcode.csv')
        self.fasta_entry.insert(END, os.getcwd() + '/my_sequences.fasta')
        self.fdr_entry.insert(END, '1.0')

        self.solver_entry = Checkbutton(master, text='Clone Solver', variable=self.solver,
                onvalue=1, offvalue=0, anchor=CENTER)
        self.silent_entry = Checkbutton(master, text='Silent Mode', variable=self.silent,
                onvalue=1, offvalue=0, anchor=CENTER)

        self.help_button = Button(master, text="Help", command=self.help, state=NORMAL)
        self.info_button = Button(master, text="Info", command=self.info, state=NORMAL)
        self.reset_button = Button(master, text="Reset", command=self.reset, state=NORMAL)
        self.run_button = Button(master, text="Run", command=self.run, state=NORMAL,bg='green')
        
        self.locations_text.grid(row=0, column=1, columnspan=2, sticky='news')
        self.fasta_text.grid(row=1, column=0, columnspan=1, sticky='nes')
        self.barcode_text.grid(row=2, column=0, columnspan=1, sticky='nes')

        for i in xrange(sgr1-3): Label(master,text='').grid(row=3+i,column=0,columnspan=3) 
        self.settings_text.grid(row=sgr1, column=1, columnspan=2, sticky='news')
        self.fdr_text.grid(row=sgr1+1, column=0, columnspan=1, sticky='nes')

        for i in xrange(sgr2-sgr1-4): Label(master,text='').grid(row=sgr1+4+i,column=0,columnspan=3) 
        self.error_text.grid(row=sgr2, column=1, columnspan=2, sticky='news')
        self.label.grid(row=sgr2+1, column=1, rowspan=1, columnspan=2, sticky='s')

        self.barcode_entry.grid(row=1, column=1, columnspan=2, sticky='news')
        self.fasta_entry.grid(row=2, column=1, columnspan=2, sticky='news')
        self.fdr_entry.grid(row=sgr1+1, column=1, columnspan=1, sticky='news')
        self.solver_entry.grid(row=sgr1+1, column=2, columnspan=1, sticky='news')
        self.silent_entry.grid(row=sgr1+2, column=2, columnspan=1, sticky='news')

        # Shift buttons over
        self.help_button.grid(row=sgr2-2, column=4)
        self.info_button.grid(row=sgr2-1, column=4)
        self.reset_button.grid(row=sgr2, column=4)
        self.run_button.grid(row=sgr2+2, column=4)


    def validate(self, new_text):
        """ Validates key choice for FDR """
        if not new_text: # the field is being cleared
            return True
        try:
            if new_text == '-' or new_text == '.':
                return True
            self.entered_number = float(new_text)
            return True
        except ValueError:
            return False


    def help(self):
        tkMessageBox.showinfo('About MAD-HYPE','This is information about MAD-HYPE.')

    def info(self):
        tkMessageBox.showinfo('Info','GUI Created by:\nPatrick V. Holec')

    def reset(self):
        self.barcode_entry.delete(0,'end')
        self.fasta_entry.delete(0,'end')
        self.fdr_entry.delete(0,'end')

        self.error_message = ""
        self.label_text.set(self.error_message)

    def run(self):

        new_errors,new_errors_full = [],[]
        
        # check for a bunch of errors
        if not os.path.isfile(self.fasta_entry.get()):
            new_errors_full.append('Fasta file not located at /{}'.format(self.fasta_entry.get()))
            new_errors.append('Fasta file not located')
        if not os.path.isfile(self.barcode_entry.get()):
            new_errors_full.append('Barcode file not located at /{}'.format(self.barcode_entry.get()))
            new_errors.append('Barcode file not located')
        if not self.fasta_entry.get().endswith('.fasta') or self.fasta_entry.get().endswith('.fastq'): 
            new_errors_full.append('Sequence file does not contain correct extension (.fasta/.fastq)')
            new_errors.append('Wrong extension (.fasta/.fastq)')
        if not self.barcode_entry.get().endswith('.csv'): 
            new_errors_full.append('Barcode file does not contain correct extension (.csv)')
            new_errors.append('Wrong extension (.csv)')
        try:
            if len(self.fdr_entry.get()) == 0:
                new_errors_full.append('No FDR given!')
                new_errors.append('No FDR given')
            elif (0. > float(self.fdr_entry.get())) or (float(self.fdr_entry.get()) > 100.):
                new_errors_full.append('Invalid FDR')
                new_errors.append('FDR out of valid range ({}%)!'.format(self.fdr_entry.get()))
        except ValueError:
            new_errors_full.append('Invalid FDR')
            new_errors.append('Invalid FDR')
        
        # check if there are any errors
        if not len(new_errors) > 0:
            tkMessageBox.showinfo('Execution Error','\n'.join(new_errors_full))
            self.label_text.set('\n'.join(new_errors))
        else:
            self.label_text.set('')
            params = {
                      'fasta_fname':self.fasta_entry.get(),
                      'barcode_fname':self.barcode_entry.get(),
                      'fdr':self.fdr_entry.get(),
                      'use_solver':self.solver,
                      'silent':self.silent
                     }

        RuntimeInterface(self.master,params)

root = Tk()
interface = UserInterface(root)
root.mainloop()


