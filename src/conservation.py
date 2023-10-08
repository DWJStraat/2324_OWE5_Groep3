"""
This module is used to get the conservation of a protein sequence and plot it.
"""

from src import pycanal


class Conservation:
    """
    A wrapper for the PyCanal package. This class is used to get the
    conservation of a protein sequence and plot it.
    """

    def __init__(self, file_path, ref=0):
        self.file_path = file_path
        self.canal = pycanal.Canal(self.file_path, ref=ref)

    def getConservation(self):
        """
        Get the conservation of the protein sequence.
        :return: A pandas dataframe with the conservation of each position.
        """
        try:
            analysis = self.canal.analysis()
            return analysis.rename(columns={"relative": "Relative Entropy"})
        except FileNotFoundError:
            print("The file was not found. Please run the pipeline first.")

    def getConservationPlot(self, file_name=None, show=False, name=None,
                            graphType='line'):
        """
        Plot the conservation of the protein sequence.
        :param file_name:  The name of the file to save the plot to.
        defaults to None.
        :param show: The plot will be shown if True. Defaults to False.
        :param name: The name of the protein sequence to be added to the
        title.
        Defaults to None.
        :param graphType: The type of graph to be plotted. Defaults to 'line'.
        :return plot: The plot of the conservation of the protein sequence.
        """
        title = "Conservation Plot"
        conservation_frame = self.getConservation()
        plot = conservation_frame.plot(kind=graphType)
        plot.set_xlabel("Position")
        plot.set_ylabel("Relative Entropy")
        if name is not None:
            title += f" for {name}"
        plot.set_title(title)
        plot = plot.get_figure()
        if show:
            plot.show()
        if file_name is not None:
            plot.savefig(file_name)
        return plot

    def getDistribution(self, site, file_name=None):
        """
        Get the distribution of the amino acids at a given site.
        :param site: The site to get the distribution of.
        :param file_name: The name of the file to save the plot to.
        defaults to None.
        :return: Returns a pandas plot with the distribution of the
        amino acids at the given site.
        """
        if file_name is None:
            return self.canal.plotSiteDistribution(site=site)
        return self.canal.plotSiteDistribution(site=site, saveplot=file_name)

    def getConsensusSequence(self, file_name=None):
        """
        Get the consensus sequence of the protein sequence.
        :param file_name: The name of the file to save the consensus sequence
        to. Defaults to None.
        :return: Returns the consensus sequence of the protein sequence.
        """
        return self.canal.getConsensusSequence(savefasta=file_name)

    def getAverageConservation(self):
        """
        Get the average conservation of the protein sequence.
        :return: Returns the average conservation of the protein sequence.
        """
        conservation = self.getConservation()
        return conservation["Relative Entropy"].mean()


if __name__ == "__main__":
    from tkinter import filedialog as fd

    file = fd.askopenfilename()
    cons = Conservation(file)
