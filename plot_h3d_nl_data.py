from __future__ import print_function
import ctypes
import pickle
from functools import reduce
from scipy.integrate import simps
import shelve
from random import choice
from numpy import exp, power
from numpy.random import randint, rand, seed
import logging
from matplotlib import pyplot as plt, colors
from scipy.optimize import curve_fit
import time
import csv

data_file = "data.txt"

def plot_MQ_v_B(filename=None):
    B_vals = []
    M_vals = []
    Q_vals = []

    if not filename:
        csvfile = open(data_file, newline='')
    else:
        csvfile = open(filename, newline='')

    datareader = csv.reader(csvfile, delimiter=',')

    for row in datareader:
        B_vals.append(float(row[0]))
        M_vals.append(float(row[1]))
        Q_vals.append(float(row[2])*100)

    plt.plot(B_vals, M_vals, 'o-b')
    plt.plot(B_vals, Q_vals, "og")
    plt.show()


if __name__ == "__main__":
    plot_MQ_v_B()
