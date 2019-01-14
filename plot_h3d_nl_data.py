from __future__ import print_function
import ctypes
import pickle
from functools import reduce
from scipy.integrate import simps
import shelve
from random import choice
import numpy as np
from numpy import exp, power
from numpy.random import randint, rand, seed
import logging
from matplotlib import pyplot as plt, colors
from scipy.optimize import curve_fit
import time
import csv

data_file = "out67_data.txt"

def plot_MQ_v_B(filename=None):
    B_vals = []
    M_vals = []
    Q_vals = []


    if not filename:
        csvfile = open(data_file, newline='')

        datareader = csv.reader(csvfile, delimiter=',')

        for row in datareader:
            B_vals.append(float(row[0]))
            M_vals.append(float(row[1]))
            Q_vals.append(float(row[2])*40000/900)
    else:
        for f in filename:
            csvfile = open(f, newline='')
            datareader = csv.reader(csvfile, delimiter=',')

            if len(Q_vals) == 0:
                for i,row in enumerate(datareader):
                    Q_vals.append(float(row[2])*40000/900);
                    B_vals.append(float(row[0]))
                    M_vals.append(float(row[1]))
            else:
                for i,row in enumerate(datareader):
                    Q_vals[i] += (float(row[2])*40000/900);

        Q_vals = [x/len(filename) for x in Q_vals]





    plt.plot(B_vals[:int(len(B_vals)/2)], M_vals[:int(len(B_vals)/2)], '.-r', label="$M ->$")
    plt.plot(B_vals[int(len(B_vals)/2):], M_vals[int(len(B_vals)/2):], '.-b', label="$M <-$")

    plt.plot(B_vals[:int(len(B_vals)/2)], Q_vals[:int(len(B_vals)/2)], ".m", label="Q/4e4 spins ->")
    plt.plot(B_vals[int(len(B_vals)/2):], Q_vals[int(len(B_vals)/2):], ".g", label="Q/4e4 spins <-")

    plt.xlabel("$B/|J|$", fontdict = {"fontsize": 16})
    plt.legend()
    plt.show()

def plot_X_inv_v_T(file_list=None, label_list = None, style_list=None):
    for i,f in enumerate(file_list):
        T_vals = []
        X_vals = []

        if not f:
            csvfile = open(data_file, newline='')
        else:
            csvfile = open(f, newline='')

        datareader = csv.reader(csvfile, delimiter=',')

        for row in datareader:
            if not float(row[0]) > 2.0:
                T_vals.append(float(row[0]))
                X_vals.append(float(row[1]))


        #plt.plot(T_vals[:-1], X_vals[:-1], '.-b', label="$J=-1, N=40$")
        plt.plot(T_vals, X_vals, style_list[i], label=label_list[i])

    plt.ylabel("$c$", fontdict = {"fontsize": 16})
    plt.xlabel("$T/|J|$", fontdict = {"fontsize": 16})
    plt.legend(title="B")
    plt.show()

def plot_heat_map():

    star_points = [
    [0.38001,0.3],
    [0.4,0.45],
    [0.4,0.6],
    [0.36,0.85],
    [0.36,1],
    [0.36,1.15],
    [0.38,1.3],
    [0.28,1.45]]

    filename = "heatmap.csv"
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        data = []
        #print(reader.fieldnames)
        #for fn in reader.fieldnames:
        #    data.append([])


        for i,row in enumerate(reader):
            if i == 0:
                continue
            if float(row[0]) > .5:
                break
            #for i,n in enumerate(reader.fieldnames):
            #    data[i].append(row[n])
            #if not i%2 == 0:
            data.append(list(reversed(row[1:10])))

        print("data: ", data[:][2])

        fig, ax = plt.subplots(1, 1, figsize=(10,10))   # Magnetization per spin, all lattices

        im = ax.imshow([[float(item)*(1000.0/400) for item in np.array(data)[:,i]] for i in range(len(data[0]))], cmap=plt.get_cmap('hot_r'), interpolation='catrom' ,origin='lower', extent = [.025, .5, .15, 1.45], aspect="auto")
        ax.set_xlabel("$T / |J|$", fontdict = {"fontsize": 16})
        ax.set_ylabel("$B / |J|$", fontdict = {"fontsize": 16})
        c = fig.colorbar(im)
        c.set_label("Q/1000 spins", fontdict = {"fontsize": 16})

        ax.plot([x[0] for x in star_points], [x[1] for x in star_points], "*w")
        plt.show()


if __name__ == "__main__":
    #plot_MQ_v_B(["out52_data.txt", "out53_data.txt", "out54_data.txt", "out55_data.txt", "out56_data.txt", "out57_data.txt", "out58_data.txt", "out59_data.txt", "out60_data.txt"])
    plot_MQ_v_B()
    #plot_X_inv_v_T([ "af_c_45.txt","af_c_60.txt", "af_c_85.txt",  "af_c_100.txt", "af_c_115.txt", "af_c_130.txt", "af_c_145.txt"], [".45",".60",".85", "1.0", "1.15", "1.3", "1.45"],[".-r",".-y", ".-k",".-c", ".-m",".-b", ".-g"])
    #plot_X_inv_v_T([ "af_c_45_micro.txt", "af_c_85_micro.txt",   "af_c_130_micro.txt", "af_c_145_micro.txt"], [".45",".85", "1.3", "1.45"],[".-r",".-k",".-b", ".-g"])
    #plot_X_inv_v_T([ "af_c_45_mini.txt","af_c_60_mini.txt", "af_c_85_mini.txt",  "af_c_100_mini.txt", "af_c_115_mini.txt", "af_c_130_mini.txt", "af_c_145_mini.txt"], [".45",".60",".85", "1.0", "1.15", "1.3", "1.45"],[".-r",".-y", ".-k",".-c", ".-m",".-b", ".-g"])

    #plot_X_inv_v_T([]"curie_chi_data_40.txt"])
    #plot_heat_map()
