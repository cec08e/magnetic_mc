from __future__ import print_function
import ctypes
import pickle, json
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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from matplotlib.colors import Normalize

data_file = "out167.txt"

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
                    Q_vals.append(float(row[2])*1000/900);
                    B_vals.append(float(row[0]))
                    M_vals.append(float(row[1]))
            else:
                for i,row in enumerate(datareader):
                    Q_vals[i] += (float(row[2])*1000/900);
                    B_vals[i]+= (float(row[0]))
                    M_vals[i]+= (float(row[1]))

        Q_vals = [x/len(filename) for x in Q_vals]
        M_vals = [x/len(filename) for x in M_vals]
        B_vals = [x/len(filename) for x in B_vals]






    plt.plot(B_vals[:int(len(B_vals)/2)], M_vals[:int(len(B_vals)/2)], '.-r', label="$M ->$")
    #plt.plot(B_vals[int(len(B_vals)/2):], M_vals[int(len(B_vals)/2):], '.-b', label="$M <-$")

    plt.plot(B_vals[:int(len(B_vals)/2)], Q_vals[:int(len(B_vals)/2)], ".m", label="Q/1000 spins ->")
    #plt.plot(B_vals[int(len(B_vals)/2):], Q_vals[int(len(B_vals)/2):], ".g", label="Q/1000 spins <-")

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


def visualize_impurities(sim_num = 0):

    filename = "impurity_sites_182.txt"

    f = open(filename, "r")
    data = eval(f.read())
    print(data)
    ################# FIXXXXX #################
    grid = [[0 for item in range(30)] for row in range(30)]
    for tup in data:
        grid[tup[1]][tup[2]] = 1
        grid[(tup[1]+1)%len(grid)][tup[2]] = 1
        grid[tup[1]][(tup[2]+1)%len(grid[0])] = 1
        grid[(tup[1]+1)%len(grid)][(tup[2]+1)%len(grid[0])] = 1

        #grid[0][0] = 10
    plt.imshow(grid, cmap=plt.get_cmap('Spectral'))

    plt.axis("off")
    plt.show()

def visualize_spins(sim_num = 0):

    filename = "lattice_record.txt"

    data = json.load(open(filename, 'rb'))
    print(len(data[0]))
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make the grid
    for i, row in enumerate(data[0]):
        for j, spin in enumerate(row):
            print("Spin: ", spin)
            ax.quiver(i*2, j*2, 0, spin[0], spin[1], spin[2], length=1, pivot = 'tail', normalize=True)
    ax.set_zlim(-30,30)
    plt.show()

def visualize_lattice(sim_num=0):
    #filename = "sim_results/sim_" + str(sim_num) + ".pickle"
    filename = "lattice_record.txt"

    data = json.load(open(filename, 'rb'))
    print(len(data[0]))
    grid = [[item[2] for item in row] for row in data[0]]
    #grid[0][0] = 10
    plt.imshow(grid, cmap=plt.get_cmap('Spectral'))
    plt.axis("off")
    plt.show()

def visualize_TC_density(sim_num = 0):
    filename = "lattice_record.txt"

    data = json.load(open(filename, 'rb'))
    print(data)
    # triangulate lattice and calculate the TC vals and plot
    TC_points = []
    x = []
    y = []
    t = []
    for i, row in enumerate(data[0]):
        for j,spin in enumerate(row):
            print("i,j: ", i, ", ", j)
            #index of point given by i*(len(row)) + j
            x.append(i)
            y.append(j)
            # Calc TC for (i,j) --> (i, j-1) ---> (i+1, j) and plot at (i+1/3, j-1/3)
            if j-1 >= 0 and i+1 < len(data[0]):
                #TC_points.append((i+(1/3.0), j-(1/3.0), calc_SA(data[0][i][j], data[0][i][(j-1)%len(row)], data[0][(i+1)%len(data[0])][j])))
                TC_points.append(calc_SA(data[0][i][j], data[0][i][(j-1)%len(row)], data[0][(i+1)%len(data[0])][j]))
                t.append((i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))
                print("Appending triangle: ", (i*len(row) + j, i*len(row) + j - 1, (i+1)*len(row) + j))

            # Calc TC for (i,j) --> (i, j+1) ---> (i-1, j) and plot at (i-1/3, j+1/3)
            if i-1 >= 0 and j+1 < len(row):
                #TC_points.append((i-(1/3.0), j+(1/3.0), calc_SA(data[0][i][j], data[0][i][(j+1)%len(row)], data[0][(i-1)%len(data[0])][j])))
                TC_points.append(calc_SA(data[0][i][j], data[0][i][(j+1)%len(row)], data[0][(i-1)%len(data[0])][j]))
                t.append((i*len(row) + j, i*len(row) + j + 1, (i-1)*len(row) + j))

    print(len(t))
    print(len(TC_points))
    print(reduce(lambda x,y: x+y, TC_points))

    #print(t[31])
    #TC_points[31] = 10

    triang = Triangulation(x,y,t)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    tpc = ax.tripcolor(triang, TC_points, cmap=plt.get_cmap('PiYG'))
    c = fig.colorbar(tpc)
    c.set_label("Q", fontdict = {"fontsize": 16})

    plt.gca().invert_yaxis()



    plt.show()

    '''
    cmap = matplotlib.cm.get_cmap('cool')
    for elem in TC_points:
        plt.plot(elem[0], elem[1], cmap=plt.get_cmap('Spectral'))
    plt.show()
    '''
def calc_SA(n1, n2, n3):
    n2_cross_n3 = np.cross(n2, n3)
    n1_dot_n2 = np.dot(n1, n2)
    n2_dot_n3 = np.dot(n2, n3)
    n3_dot_n1 = np.dot(n3, n1)


    rho = np.sqrt(2*(1+n1_dot_n2)*(1+n2_dot_n3)*(1+n3_dot_n1));
    temp = complex((1.0/rho)*(1 + n1_dot_n2 + n2_dot_n3 + n3_dot_n1), (1/rho)*np.dot(n1, n2_cross_n3))

    Omega = (2*np.log(temp).imag)/(4*np.pi)
    #print("Omega: ", Omega)

    return Omega

def coercive_plot():
    B = np.linspace(-1,1,10000)
    M1 = [np.cos(np.pi  - .1/(x+(5/2))) for x in B]
    M2 = [np.cos(np.pi - .25/(x+(5/2))) for x in B]
    M3 = [np.cos(np.pi - .4/(x+(5/2))) for x in B]
    plt.plot( B, M1, 'g', label="D = .1")
    plt.plot( B, M2, 'b', label="D = .25")
    plt.plot( B, M3, 'r', label="D = .4")
    plt.legend()
    plt.xlabel("$B$", fontdict = {"fontsize": 16})
    plt.ylabel(r"$cos(\theta )$", fontdict = {"fontsize": 16})



    plt.show()



if __name__ == "__main__":
    #plot_MQ_v_B(["out115.txt", "out116.txt", "out117.txt", "out118.txt","out119.txt", "out120.txt", "out21.txt","out22.txt", "out23.txt","out24.txt","out25.txt", "out26.txt", "out27.txt", "out28.txt", "out29.txt", "out30.txt", "out131.txt", "out132.txt", "out133.txt", "out134.txt",
    #             "out140.txt", "out141.txt", "out142.txt", "out143.txt", "out144.txt", "out145.txt", "out146.txt", "out147.txt", "out148.txt",
    #             "out149.txt", "out150.txt", "out151.txt", "out152.txt", "out153.txt", "out154.txt", "out155.txt", "out166.txt", "out167.txt", "out168.txt", "out169.txt", "out170.txt",
    #             "out171.txt", "out172.txt", "out173.txt", "out174.txt", "out175.txt"])
    #plot_MQ_v_B(["out176.txt", "out177.txt", "out178.txt", "out179.txt","out180.txt"])
    #plot_MQ_v_B(["out181.txt", "out182.txt", "out183.txt", "out184.txt","out185.txt"])

    #plot_MQ_v_B()
    #coercive_plot()
    #plot_X_inv_v_T([ "af_c_45.txt","af_c_60.txt", "af_c_85.txt",  "af_c_100.txt", "af_c_115.txt", "af_c_130.txt", "af_c_145.txt"], [".45",".60",".85", "1.0", "1.15", "1.3", "1.45"],[".-r",".-y", ".-k",".-c", ".-m",".-b", ".-g"])
    #plot_X_inv_v_T([ "af_c_45_micro.txt", "af_c_85_micro.txt",   "af_c_130_micro.txt", "af_c_145_micro.txt"], [".45",".85", "1.3", "1.45"],[".-r",".-k",".-b", ".-g"])
    #plot_X_inv_v_T([ "af_c_45_mini.txt","af_c_60_mini.txt", "af_c_85_mini.txt",  "af_c_100_mini.txt", "af_c_115_mini.txt", "af_c_130_mini.txt", "af_c_145_mini.txt"], [".45",".60",".85", "1.0", "1.15", "1.3", "1.45"],[".-r",".-y", ".-k",".-c", ".-m",".-b", ".-g"])

    #plot_X_inv_v_T([]"curie_chi_data_40.txt"])
    #plot_heat_map()
    #visualize_spins()
    visualize_lattice()
    visualize_TC_density()
    visualize_impurities()
