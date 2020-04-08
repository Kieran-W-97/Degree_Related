"""AUTHOR: Kieran Wachsmuth, 19th Jan 2020
A code to plot the Hydrophobicity of Amino Acid Sequence"""
from hydro_fns import csv_reader
from hydro_fns import bin_lister
import matplotlib.pyplot as plt
gene_seq = csv_reader('gene_sequence')
L = len(gene_seq)
# create dict of hydrophobicity values for each Amino Acid:
hydro_dict = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'W':-0.9,'S':-0.8,'Y':-1.3,
              'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5}
hydro_list = []
for i in range(L):
    hydro_list.append(hydro_dict[gene_seq[i]])
# create bins for averaging hydrophobicity:
B = int(input('Enter size of bins (odd number only): '))
bins_array_1 = bin_lister(B,hydro_list)
# why not have a second bin:
B_2 = int(input('Enter size of bins (2nd bining), (Odd number only):'))
bins_array_2 = bin_lister(B_2,hydro_list)
# plot:
x_bin = bins_array_1[0]
y_bin = bins_array_1[1]
x_bin_2 = bins_array_2[0]
y_bin_2 = bins_array_2[1]
x = range(1,L+1)
y = hydro_list
#plt.xticks(x, gene_seq)
plt.xticks(range(0,500,25))
plt.plot(x,y,'bo', x_bin,y_bin,'g', x_bin_2,y_bin_2,'r')
#plt.plot(x_bin,y_bin, 'r')
plt.xlabel('Amino Acid')
plt.ylabel('Relative Hydrphobicity')
plt.axhline(color = 'k')
plt.grid(True)
plt.show()