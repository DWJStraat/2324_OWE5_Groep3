import sys
import os
from Bio.Align.Applications import MafftCommandline
import pyHMMER as hmmer


class main():
    def __init__(self, file, database):
        self.hmm = None
        self.hits = None
        self.file = file
        self.database = database
        self.msa_file = f'{self.file}-output.fna'
        self.hmm_file = f'{self.file}.hmm'
        self.alphabet = hmmer.easel.Alphabet.amino()
        self.background = hmmer.plan7.Background(self.alphabet)

    def mafft(self):
        MafftCommandline(input=self.file, auto=True, fastaout=True,
                         out=self.msa_file)

    def hmmbuild(self):

        with hmmer.easel.SequenceFile(self.msa_file, self.alphabet) as \
                seqfile:
            msa = seqfile.read()
        builder = hmmer.plan7.Builder()
        self.hmm, __, __ = builder.build(msa, self.background)
        with hmmer.easel.HMMFile(self.hmm_file, 'w') as hmmfile:
            hmmfile.write(self.hmm)

    def hmmsearch(self):
        pipeline = hmmer.plan7.Pipeline(self.alphabet, self.background)
        with hmmer.easel.SequenceFile(self.database, self.alphabet) as \
                seqfile:
            self.hits = pipeline.search_hmm(self.hmm, seqfile)
        with open(f'{self.file}-output.txt', 'w') as output:
            for hit in self.hits:
                output.write(f'{hit}\n')




target_file = sys.argv[1]
database = sys.argv[2]





