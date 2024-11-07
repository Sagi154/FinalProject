import pandas as pd
from scipy.spatial.distance import euclidean, pdist, squareform
from scipy.sparse.csgraph import laplacian
import numpy as np
import pprint
import math
# np.random.seed(1234)
import logging
import os

import Prev_final_100.symnmf as prev
import Prev_final_100.analysis as prev_anal
import symnmf as our

FILE_NAME = "input.txt"

current_dir = "."
prev_dir = "Prev_final_100"
current_symnmf_path = "${current_dir}/symnmf.py"
prev_symnmf_path = "${prev_dir}/symnmf.py"


def set_log_config():
	logging.basicConfig(format="[%(levelname)s] %(asctime)s %(message)s",
						datefmt="%Y-%m-%d %H:%M:%S",
						encoding='utf-8',
						handlers=[logging.FileHandler("my_logs.log"),
								  logging.StreamHandler()],
						level=logging.DEBUG)


def create_data_vectors(file_name, vectors_count, vector_length, bound):
	data = [[np.random.uniform(-bound, bound) for j in range(vector_length)] for i in range(vectors_count)]
	data = np.array(data)
	with open(file_name, 'w') as file:
		for line in data:
			file.write(f"{','.join(str('%.4f' % element) for element in line)}\n")
	logging.info("Created data vectors")
	return data


def transform_matrix_elements_to_4_digits(matrix):
	for line in matrix:
		for element in line:
			element = float("{:.4f}".format(element))


def compare_results(expected, actual):
	transform_matrix_elements_to_4_digits(expected)
	transform_matrix_elements_to_4_digits(actual)
	for i, line in enumerate(actual):
		for j, element in enumerate(line):
			if abs(element - expected[i][j]) > 0.0001:
				print(f"Error at ({i},{j}) element")
				return False
	return True



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


def ddg(data_points):
	n = len(data_points)
	sym_matrix = sym(data_points)
	# print("printing sym matrix from tester")
	# print_matrix(sym_matrix)
	degrees = np.sum(sym_matrix, axis=1)
	ddg_matrix = np.diag(degrees)
	return ddg_matrix


def norm(data_points):
	n = len(data_points)
	sym_matrix = sym(data_points)
	ddg_matrix = ddg(data_points)
	ddg_inv_sqrt = [[1.0 / math.sqrt(ddg_matrix[i][j]) if ddg_matrix[i][j] != 0.0 else 0.0 for j in range(n)] for i in range(n)]
	ddg_inv_sqrt = np.array(ddg_inv_sqrt)
	norm_matrix = ddg_inv_sqrt @ sym_matrix
	norm_matrix = norm_matrix @ ddg_inv_sqrt
	return norm_matrix

def reached_error(expected, actual, goal):
	print(f"Comparison in goal:{goal} failed")
	print("Expected matrix:")
	print_matrix(expected)
	print("actual matrix:")
	print_matrix(actual)

def run_test(vectors_count_limit, vector_length_limit, cord_value_limit):
	print("Running test...")
	count = 0
	for i in range(2, vectors_count_limit):
		for j in range(2, vector_length_limit):
			vectors_count = i
			vector_length = j
			file_name = f"testing/input_{i}_{j}.txt"
			data_points = create_data_vectors(file_name, vectors_count, vector_length, cord_value_limit).tolist()
			for K in range(1, min(vectors_count, 5) + 1):
				print(f"Running test for matrix ({i},{j})")
				sym_matrix_our, sym_matrix_prev = our.perform_goal(data_points, i, j, K, "sym")
				if not compare_results(sym_matrix_prev, sym_matrix_our):
					reached_error(sym_matrix_prev, sym_matrix_our, "sym")
					count += 1
					return
				ddg_matrix_our, ddg_matrix_prev = our.perform_goal(data_points, i, j, K, "ddg")
				if not compare_results(ddg_matrix_prev, ddg_matrix_our):
					reached_error(ddg_matrix_prev, ddg_matrix_our, "ddg")
					count += 1
					return
				norm_matrix_our, norm_matrix_prev = our.perform_goal(data_points, i, j, K, "norm")
				if not compare_results(norm_matrix_prev, norm_matrix_our):
					reached_error(norm_matrix_prev, norm_matrix_our, "norm")
					count += 1
					return
				symnmf_matrix_our, symnmf_matrix_prev = our.perform_goal(data_points, i, j, K, "symnmf")
				if not compare_results(symnmf_matrix_prev, symnmf_matrix_our):
					reached_error(symnmf_matrix_prev, symnmf_matrix_our, "symnmf")
					count += 1
			# return
				# os.remove(f"{file_name}")
				# if successful and reached here, delete previous test file
	print(f"Errors count: {count}")
	print("Test successful")


def main():
	vectors_count_limit = 5
	vector_length_limit = 4
	cord_value_limit = 2
	run_test(vectors_count_limit, vector_length_limit, cord_value_limit)


if __name__ == '__main__':
	main()
