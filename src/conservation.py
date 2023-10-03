from src import pycanal


class conservation():
    def __init__(self, file_path):
        self.file_path = file_path
        self.canal = pycanal.Canal(self.file_path,
                                   ref = 0)

    def getConservation(self):
        return self.canal.analysis(include=None, method='relative')

    def getDistribution(self, site, file_name= None):
        if file_name is None:
            return self.canal.plotSiteDistribution(site=site)
        return self.canal.plotSiteDistribution(site=site, saveplot=file_name)

    def getConsensusSequence(self, file_name = None):
        return self.canal.getConsensusSequence(savefasta=file_name)

a = conservation('outputaaaa.fna')
b = a.getConservation()