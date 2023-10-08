"""
This file is used to run the pipeline in a conventional form, rather than a
commandline form. It takes in the name of the file, the name of the
database, and the number of iterations as input.
"""

import pipeline

file_name = input('Enter the name of the file: ')
database_name = input('Enter the name of the database: ')
iterations = int(input('Enter the number of iterations: '))
plot = input('Do you want to plot the hydrophobicity and conservation? (y/n) ')

pipe = pipeline.Pipeline(file_name, database_name, iterations)
pipe.run()
if plot == 'y':
    pipe.savePlot()
