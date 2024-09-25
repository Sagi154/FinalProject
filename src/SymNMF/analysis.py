import sys

from src.SymNMF.symnmf import parse_arguments

def parse_arguments():
    K = int(float(sys.argv[1]))
    file_name = sys.argv[2]
    return K, file_name

def symnmf_cluster_assignment():


def main():
    K, file_name = parse_arguments()