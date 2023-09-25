import os

class Pipeline:
    def __init__(self, file,database):
        self.blast = None
        self.file = file
        self.database = database
        self.msa_file = f'{self.file}-msa.fna'

    def align(self):
        os.system(f'mafft {self.file} > {self.msa_file}')

    def hmmbuild(self):
        os.system(f'hmmbuild {self.file}.hmm {self.msa_file}')

    def hmmsearch(self):
        os.system(f'hmmsearch --tblout {self.file}-output.txt '
                  f'{self.file}.hmm {self.database}')
