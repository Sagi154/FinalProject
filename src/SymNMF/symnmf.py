import sys
import pandas as pd
import numpy as np
import math
import symnmfmodule as c
from tester import *
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
	data_points = data_points.to_numpy()
	return data_points.tolist()


def parse_arguments():
	K = int(float(sys.argv[1]))
	goal = sys.argv[2]
	file_name = sys.argv[3]
	return K, goal, file_name


def main():
	K, goal, file_name = parse_arguments()
	# try:
	# create_data_vectors(file_name, 4, 4, 2)
	data_points = read_input_file(file_name)
	if data_points is None:
		print("An Error Has Occurred")
		return
	vectors_count = len(data_points)
	vector_length = len(data_points[0])
	result_matrix = None
	expected = None
	if goal == "symnmf":
		normalized_similarity_matrix = c.norm(data_points, vectors_count, vector_length)
		if normalized_similarity_matrix is None:
			return
		m = avg_W_entries(normalized_similarity_matrix, vectors_count)
		initial_H = initialize_decomposition_matrix_H(vectors_count, m, K)
		result_matrix = c.symnmf(K, initial_H, normalized_similarity_matrix, vectors_count)
	elif goal == "sym":
		result_matrix = c.sym(data_points, vectors_count, vector_length)
		expected = sym(data_points)
	elif goal == "ddg":
		result_matrix = c.ddg(data_points, vectors_count, vector_length)
		expected = ddg(data_points)
	elif goal == "norm":
		result_matrix = c.norm(data_points, vectors_count, vector_length)
		expected = norm(data_points)
	if result_matrix is None:
		return
	else:
		print("printing result matrix received from C")
		for line in result_matrix:
			print(",".join(str("%.4f" % element) for element in line))
		# print("printing expected matrix from tester")
		# print_matrix(expected)
		# comparison = compare_results(expected, result_matrix)
		# print(f"Comparison: {comparison}")

	# except Exception as e:
	# 	print("An Error Has Occurred")


if __name__ == "__main__":
	main()
