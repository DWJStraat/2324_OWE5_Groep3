import pipeline

file_name = input('Enter the name of the file: ')
database_name = input('Enter the name of the database: ')
iterations = int(input('Enter the number of iterations: '))

pipe = pipeline.Pipeline(file_name, database_name, iterations)
pipe.run()