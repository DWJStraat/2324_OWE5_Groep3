"""
A pipeline for building and searching HMMs.
Written by: David Straat
Date: 25-sept-2023
"""
import os
import argparse

import Bio
from Bio import SCOP, SeqIO
from hydrophobicity import Hydrophobicity
from conservation import Conservation
import pandas as pd


class Pipeline:
    """
    This class is a pipeline for building and searching HMMs.
    """

    def __init__(self, file, nr_database, delta_database, iterations=1):
        self.frame = None
        self.sequences = None
        self.cons = None
        self.hydro = None
        self.plot = None
        self.domains = {}
        self.blast = None
        self.file = file
        self.nrDatabase = nr_database
        self.DeltaBase = delta_database
        self.iter = iterations
        self.msa_file = f'{self.file}-msa.fna'

    def readFasta(self):
        self.sequences = []
        self.headers = []
        self.sequences.extend(record.seq for record in SeqIO.parse(
            self.file, "fasta"))
        self.headers.extend(record.id for record in SeqIO.parse(
            self.file, "fasta"))

    def align(self):
        """
        Aligns the sequences in the file using MAFFT.
        """
        os.system(f'mafft {self.file} > {self.msa_file}')

    def hmmBuild(self):
        """
        Builds a HMM from the aligned sequences.
        """
        os.system(f'hmmbuild {self.file}.hmm {self.msa_file}')

    def hmmSearch(self):
        """
        Searches the database for sequences that match the HMM.
        """
        os.system(f'hmmsearch --tblout {self.file}-output.txt '
                  f'{self.file}.hmm {self.nrDatabase}')

    def run(self, iterations: int = None, get_plot: bool = False):
        """
        Runs the pipeline.
        :param iterations: The number of iterations to run. Defaults to the
        parameter passed when initiating the class.
        :param get_plot: If True, a plot of the hydrophobicity and
        conservation of the sequences will be made. Defaults to False.
        :return: None
        """
        self.readFasta()
        if iterations is not None:
            self.iter = iterations
        for _ in range(self.iter):
            self.loop()
        if get_plot:
            self.getPlot()

    def loop(self):
        """
        Runs the pipeline once.
        """
        self.align()
        self.hmmBuild()
        self.hmmSearch()

    def getDomains(self):
        """
        Gets the domains of the sequences.
        """
        with open(f'{self.file}-output.txt') as file:
            content = file.readlines()
        for line in content:
            print(line)
        sequences = []
        for sequence in sequences:
            scop_handle = SCOP.search(sequence)
            self.domains[sequence] = scop_handle

    def getPlot(self, hydrophobicity: bool = True, conservation: bool = True,
                graphType: str = 'line', name: str = None, ref_seq: int =
                0):
        """
        Gets a plot of the hydrophobicity and/or conservation of the sequences.
        """
        ax = None
        height = 10
        width = 50
        descriptor = ""
        if hydrophobicity:
            self.hydro = Hydrophobicity(self.sequences[ref_seq])
            hydroFrame = self.hydro.calculateAllHydrophobicity()
            ax = hydroFrame.plot(kind=graphType, figsize=(width, height))
            descriptor += "Hydrophobicity"
        if conservation:
            self.cons = Conservation(self.file, ref=ref_seq)
            consPanda = self.cons.getConservation()
            if hydrophobicity:
                descriptor += " and "
            if ax is None:
                ax = consPanda.plot(kind = graphType, figsize=(width, height))
            else:
                consPanda.plot(ax=ax, secondary_y=True)
            descriptor += "Conservation"
        self.plot = ax
        self.plot.set_xlabel("Position")
        title = f"Plot of {descriptor}"
        if name is not None:
            title += f" for {name}"
        self.plot.set_title(title)

        self.plot = self.plot.get_figure()


    def showPlot(self):
        if self.plot is None:
            self.getPlot()
        self.plot.show()

    def savePlot(self, file_name: str = None):
        if self.plot is None:
            self.getPlot()
        if file_name is None:
            file_name = f"{self.file}-plot.png"
        self.plot.savefig(file_name)


if "__name__" == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('file', help='The name of the file to use.')
    argparse.add_argument('database', help='The name of the database to use.')
    argparse.add_argument('iterations',
                          help='The number of iterations to run.')
    argparse.add_argument('-p', '--plot', help='Show the plot of the '
                                               'hydrophobicity and '
                                               'conservation of the sequences.',
                          action='store_true')
    args = argparse.parse_args()
    pipe = Pipeline(args.file, args.database, args.iterations)
    pipe.run()
    if args.plot:
        pipe.savePlot()
