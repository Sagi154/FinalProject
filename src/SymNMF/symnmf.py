import sys
import pandas as pd
import numpy as np
import math
import symnmfmodule as c
np.random.seed(1234)


def avg_W_entries(normalized_similarity_matrix, vectors_count):
	total_matrix_sum = 0
	for line in normalized_similarity_matrix:
		for element in line:
			total_matrix_sum += element
	return total_matrix_sum / (vectors_count * vectors_count)


def initialize_decomposition_matrix_H(vectors_count, m, K):
	# Upper bound not tight
	decomposition_matrix = [[np.random.uniform(low=0, high=(2 * math.sqrt(m/K) + 1e-10)) for j in range(K)] for i in range(vectors_count)]
	return decomposition_matrix


def read_input_file(file_name):
	data_points = pd.read_csv(file_name, header=None)
	data_points = data_points.iloc[:, 1:].to_numpy()
	# data_points = data_points.to_numpy() # if first column isn't a key
	return data_points


def parse_arguments():
	K = int(float(sys.argv[1]))
	goal = sys.argv[2]
	file_name = sys.argv[3]
	return K, goal, file_name


def main():
	K, goal, file_name = parse_arguments()
	try:
		data_points = read_input_file(file_name)
		vectors_count = len(data_points)
		result_matrix = None
		if goal == "symnmf":
			normalized_similarity_matrix = c.norm(data_points.tolist(), vectors_count)
			if normalized_similarity_matrix is None:
				print("An Error Has Occurred")
				return
			m = avg_W_entries(normalized_similarity_matrix, vectors_count)
			result_matrix = c.symnmf(K, initialize_decomposition_matrix_H(vectors_count, m, K), normalized_similarity_matrix, vectors_count)
		elif goal == "sym":
			result_matrix = c.sym(data_points.tolist(), vectors_count)
		elif goal == "ddg":
			result_matrix = c.ddg(data_points.tolist(), vectors_count)
		elif goal == "norm":
			result_matrix = c.norm(data_points.tolist(), vectors_count)
		if result_matrix is None:
			print("An Error Has Occurred")
			return
		else:
			for line in result_matrix:
				print(",".join(str("%.4f" % element) for element in line))
	except Exception as e:
		print("An Error Has Occurred")