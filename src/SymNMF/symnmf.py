import sys
import pandas as pd
import numpy as np
import math
import symnmfmodule as c
np.random.seed(1234)


def initialize_decomposition_matrix_H(vectors_count: int, m: float, K: int) -> np.ndarray:
	"""
	Performing step 1.4.1 of the algorithm.
	:param vectors_count: Number of vectors in the data set.
	:param m: The mean of all entries of W.
	:param K: Number of clusters.
	:return: Initial decomposition matrix.
	"""
	# Upper bound not tight
	decomposition_matrix = [[2 * math.sqrt(m / K) * np.random.uniform(0, 1) for j in range(K)] for i in range(vectors_count)]
	# decomposition_matrix = [[np.random.uniform(low=0, high=(2 * math.sqrt(m / K) + 1e-10)) for j in range(K)] for i in range(vectors_count)]
	return decomposition_matrix


def read_input_file(file_name: str) -> list:
	"""
	Load dataset from a .txt file
	:param file_name: The name of the .txt file.
	:return: A numpy array representing the data.
	"""
	data_points = pd.read_csv(file_name, header=None)
	data_points = data_points.to_numpy()
	return data_points.tolist()


def parse_arguments():
	"""
	Parsing the arguments from the command line.
	:return: K, goal and file name variables.
	"""
	K = int(float(sys.argv[1]))
	goal = sys.argv[2]
	file_name = sys.argv[3]
	return K, goal, file_name


def perform_goal(data_points: list, vectors_count: int, vector_length: int, K: int, goal: str) -> list:
	"""
	Performs the required goal
	:param data_points: A numpy array representing the data vectors.
	:param vectors_count: Number of vectors in the data.
	:param vector_length: Length of the vectors.
	:param K: Number of clusters.
	:param goal: Goal client wants to perform.
	:return: result matrix of said goal.
	"""
	result_matrix = None
	if goal == "symnmf":
		normalized_similarity_matrix = c.norm(data_points)
		if normalized_similarity_matrix is None:
			return
		m = np.mean(np.array(normalized_similarity_matrix))
		initial_H = initialize_decomposition_matrix_H(vectors_count, m, K)
		result_matrix = c.symnmf(K, initial_H, normalized_similarity_matrix, vectors_count)
	elif goal == "sym":
		result_matrix = c.sym(data_points)
	elif goal == "ddg":
		result_matrix = c.ddg(data_points)
	elif goal == "norm":
		result_matrix = c.norm(data_points)
	return result_matrix


def main():
	"""
	Main function
	"""
	K, goal, file_name = parse_arguments()
	try:
		data_points = read_input_file(file_name)
		if data_points is None:
			print("An Error Has Occurred")
			return
		vectors_count, vector_length = len(data_points), len(data_points[0])
		result_matrix = perform_goal(data_points, vectors_count, vector_length, K, goal)
		if result_matrix is None:
			print("An Error Has Occurred")
			return
		else:
			for line in result_matrix:
				print(",".join(str("%.4f" % element) for element in line))
	except Exception as e:
		print("An Error Has Occurred")


if __name__ == "__main__":
	main()
