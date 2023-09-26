"""
A pipeline for building and searching HMMs.
Written by: David Straat
Date: 25-sept-2023
"""
import os
import argparse


class Pipeline:
    """
    This class is a pipeline for building and searching HMMs.
    """

    def __init__(self, file, database, iterations=1):
        self.blast = None
        self.file = file
        self.database = database
        self.iter = iterations
        self.msa_file = f'{self.file}-msa.fna'

    def align(self):
        """
        Aligns the sequences in the file using MAFFT.
        """
        os.system(f'mafft {self.file} > {self.msa_file}')

    def hmmbuild(self):
        """
        Builds a HMM from the aligned sequences.
        """
        os.system(f'hmmbuild {self.file}.hmm {self.msa_file}')

    def hmmsearch(self):
        """
        Searches the database for sequences that match the HMM.
        """
        os.system(f'hmmsearch --tblout {self.file}-output.txt '
                  f'{self.file}.hmm {self.database}')

    def run(self, iterations=None):
        """
        Runs the pipeline.
        """
        if iterations is not None:
            self.iter = iterations
        for _ in range(self.iter):
            self.loop()

    def loop(self):
        """
        Runs the pipeline once.
        """
        self.align()
        self.hmmbuild()
        self.hmmsearch()

if "__name__" == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('file', help='The name of the file to use.')
    argparse.add_argument('database', help='The name of the database to use.')
    argparse.add_argument('iterations', help='The number of iterations to run.')
    args = argparse.parse_args()
    pipe = Pipeline(args.file, args.database, args.iterations)
    pipe.run()