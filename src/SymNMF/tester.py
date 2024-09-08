import pandas as pd
from scipy.spatial.distance import euclidean, pdist, squareform
import numpy as np
import pprint
import math
# np.random.seed(1234)
import logging

FILE_NAME = "input.txt"


def set_log_config():
	logging.basicConfig(format="[%(levelname)s] %(asctime)s %(message)s",
						datefmt="%Y-%m-%d %H:%M:%S",
						encoding='utf-8',
						handlers=[logging.FileHandler("my_logs.log"),
								  logging.StreamHandler()],
						level=logging.DEBUG)


def create_data_vectors(n, bound):
	data = [[np.random.uniform(-bound, bound) for j in range(n)] for i in range(n)]
	data = np.array(data)
	with open(FILE_NAME, 'w') as file:
		for line in data:
			file.write(f"{','.join(str('%.4f' % element) for element in line)}\n")
	print(data)
	logging.info("Created data vectors")
	return data


def print_matrix(matrix):
	for line in matrix:
		print(",".join(str("%.4f" % element) for element in line))
def read_input_file(file_name):
	data_points = pd.read_csv(file_name, header=None)
	data_points = data_points.to_numpy()
	return data_points

def similarity_func(vector1, vector2):
	return math.exp(-(math.pow(euclidean(vector1, vector2), 2) / 2))

def sym(data_points):
	dists = pdist(data_points, similarity_func)
	sym_matrix = squareform(dists)
	return sym_matrix


def main():
	flag = True
	set_log_config()
	n = 10
	range = 2
	if flag:
		create_data_vectors(n, range)
	data_points = read_input_file(FILE_NAME)
	print("Data points:")
	print_matrix(data_points)
	sym_matrix = sym(data_points)
	print("Similarity matrix:")
	print_matrix(sym_matrix)



if __name__ == '__main__':
	main()
