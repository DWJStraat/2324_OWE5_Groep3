"""
A module to calculate the hydrophobicity of a protein sequence utilizing the
Kyte-Doolittle scale.
"""

import pandas as pd


class Hydrophobicity:
    """
    A class to calculate the hydrophobicity of a protein sequence utilizing
    the Kyte-Doolittle scale.
    """

    def __init__(self, sequence):
        self.aminozuur_dict = {
            "A": 1.800,
            "R": -4.500,
            "N": -3.500,
            "D": -3.500,
            "C": 2.500,
            "Q": -3.500,
            "E": -3.500,
            "G": -0.400,
            "H": -3.200,
            "I": 4.500,
            "L": 3.800,
            "K": -3.900,
            "M": 1.900,
            "F": 2.800,
            "P": -1.600,
            "S": -0.800,
            "T": -0.700,
            "W": -0.900,
            "Y": -1.300,
            "V": 4.200
        }
        self.sequence = sequence.replace("-", "")

    def calculateHydrophobicity(self, site: int):
        """
        Calculate the hydrophobicity of a site.
        :param site: The site to calculate the hydrophobicity of.
        :return: The hydrophobicity of the site.
        """
        return self.aminozuur_dict[self.sequence[site]]

    def calculateAverageHydrophobicity(self):
        """
        Calculate the average hydrophobicity of the sequence.
        :return: The average hydrophobicity of the sequence.
        """
        return sum(
            self.aminozuur_dict[aa] for aa in self.sequence
        ) / len(self.sequence)

    def calculateAllHydrophobicity(self):
        """
        Calculate the hydrophobicity of all sites in the sequence.
        :return: A pandas dataframe of the hydrophobicity of all sites in the
        sequence.
        """
        return pd.DataFrame(
            [self.aminozuur_dict[aa] for aa in self.sequence],
            columns=["Hydrophobicity"])

    def plotHydrophobicity(self, file_name: str = None, show: bool = False,
                           name: str = None):
        """
        Plot the hydrophobicity of the sequence.
        :param file_name: The name of the file to save the plot to. Will not
        be saved by default.
        :param show: The plot will be shown if True. Defaults to False.
        :param name: The name of the protein sequence to be added to the
        title. Defaults to None.
        :return: Returns a pandas plot of the hydrophobicity of the sequence.
        """
        df = self.calculateAllHydrophobicity()
        plot = df.plot()
        plot.set_xlabel("Position")
        plot.set_ylabel("Hydrophobicity")
        title = "Hydrophobicity Plot"
        if name is not None:
            title += f" for {name}"
        plot.set_title(title)
        plot = plot.get_figure()
        if show:
            plot.show()
        if file_name is not None:
            plot.savefig(file_name)
        return plot

