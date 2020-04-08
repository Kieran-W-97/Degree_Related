"""A series of functions used in the hydrophobic analysis"""
# function to load the csv gene sequence into a list:
import csv
import numpy as np
def csv_reader(filename):
    gene_seq = []
    with open(filename+'.csv') as csvfile:
        readcsv = csv.reader(csvfile)
        for row in readcsv:
            gene_seq = row
    return gene_seq
# a function to create lists of bins:
def bin_lister(B,hydro_list):
    # B is the number of elements in each bin. At the moment only odd values of B are taken.
    L = len(hydro_list)
    N = int(float((B - 1)) / 2)
    x_bin = range(N+1, L - N + 1)
    bin_dict = {}
    for i in x_bin:
        bin_list = []
        for j in range(B):
            bin_list.append(hydro_list[i + j - N - 1])
        bin_dict[i] = bin_list # create a dict of lists of elements of ith bin.
    y_bin = []
    # create array of weights (triangular weighting):
    e_w = float(input('Enter edge weight relative to centre (as decimal): '))
    w = np.zeros(B)
    # calculate incremental increase of weighting:
    d = float(1-e_w)/(float(B-1)/2)
    for i in range(int(float(B-1)/2)):
        w[i] = e_w
        e_w = e_w + d
    w[N] = 1
    for i in range(int(float(B+1)/2),B):
        e_w = e_w - d
        w[i] = e_w
    for i in x_bin:
        y_bin.append(np.average(bin_dict[i], weights = w))
    return np.array([x_bin,y_bin])