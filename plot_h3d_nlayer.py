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
from matplotlib import pyplot, colors
from scipy.optimize import curve_fit
import time

SIM_NUM = 1
NUM_L = 1
ROWS = 20
COLS = 20

PDOUBLE = ctypes.POINTER(ctypes.c_double)
PPDOUBLE = ctypes.POINTER(PDOUBLE)

_h3d = ctypes.CDLL('./heisenberg3d_nlayer.so')

class Spin_t(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double),
                ("y", ctypes.c_double),
                ("z", ctypes.c_double)]

PSPIN_T = ctypes.POINTER(Spin_t)
PPSPIN_T = ctypes.POINTER(PSPIN_T)
PPPSPIN_T = ctypes.POINTER(PPSPIN_T)



_h3d.initialize_lattice.argtypes = ()
_h3d.M_v_B.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.M_v_B.restype = ctypes.c_int
_h3d.record_lattice.argtypes = (ctypes.POINTER(PPSPIN_T),)
_h3d.record_lattice.restype = None
_h3d.sample_mag.argtypes = (PDOUBLE, ctypes.c_int)
_h3d.sample_mag.restype = None

_h3d.mag_v_temp.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.mag_v_temp.restype = ctypes.c_int
_h3d.suscept_v_temp.argtypes = (ctypes.POINTER(PDOUBLE),)
_h3d.suscept_v_temp.restype = ctypes.c_int


#_h3d.staggered_mag_v_temp.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.staggered_mag_v_temp.restype = ctypes.c_int
#_h3d.staggered_suscept_v_temp.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.staggered_suscept_v_temp.restype = ctypes.c_int
#_h3d.M_v_T.argtypes = (ctypes.POINTER(PDOUBLE), ctypes.c_double, ctypes.c_double, ctypes.c_double)
#_h3d.M_v_T.restype = ctypes.c_int
#_h3d.M_v_K.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.M_v_K.restype = ctypes.c_int
#_h3d.M_v_J.argtypes = (ctypes.POINTER(PDOUBLE),)
#_h3d.M_v_J.restype = ctypes.c_int

def autocorrelate_mag(sweeps, selection):
    # Calc and plot magnetization autocorrelation function as a function of sweeps
    # Average magnetization per site

    # self.mag_vals should hold all the magnetization per spin values calculated at the end of
    # each sweep.
    global _h3d
    mag_vals_arr = ctypes.c_double*sweeps
    mag_vals = mag_vals_arr()

    _h3d.initialize_lattice()
    _h3d.sample_mag(mag_vals, sweeps)

    avg_mag = reduce((lambda x,y: x+y), mag_vals)/len(mag_vals)
    avg_mag_sq = power(avg_mag, 2)

    chi_0 = autocorrelate_int(avg_mag_sq, 0, sweeps, mag_vals) # Initial autocorrelation function chi(0)
    x = [i for i in range(selection)]
    y = [autocorrelate_int(avg_mag_sq, i, sweeps, mag_vals, norm=chi_0) for i in range(selection)]
    pyplot.plot(x, y, 'b')
    pyplot.ylabel('Magnetization autocorrelation $\chi(t)$')
    pyplot.xlabel('Sweeps')

    pyplot.show()

    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump([x,y], fp)


def autocorrelate_int(avg_mag_sq, t, sweep_num, mag_vals, norm = 1):
    # Performs autocorrelation integral at sweep t
    # Can only integrate over t'=0 to t'=(self.step_num/self.lat_size)-t ?
    #t_prime = [i for i in range((int(self.step_num/self.lat_size) - t))] # x samples
    print("Autocorrelating with t = ", t)
    t_prime = [i for i in range(sweep_num - t)] # x samples
    y_samples = [((mag_vals[i]*mag_vals[i+t]) - avg_mag_sq)/norm for i in t_prime]
    return simps(y_samples, t_prime)

'''
def plot_staggered_mag():
    global _h3d

    mag_vals_entry = ctypes.c_double*2
    mag_vals_arr = PDOUBLE*10000

    mag_vals = mag_vals_arr()

    for i in range(10000):
        mag_vals[i] = mag_vals_entry()

    _h3d.initialize_lattice()
    num_samples = _h3d.staggered_mag_v_temp(mag_vals)

    data = [[] for i in range(2)]
    for i in range(num_samples):
        data[0].append(mag_vals[i][0])
        data[1].append(mag_vals[i][1])


    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)
'''
def plot_mag():
    global _h3d

    mag_vals_entry = ctypes.c_double*2
    mag_vals_arr = PDOUBLE*10000

    mag_vals = mag_vals_arr()

    for i in range(10000):
        mag_vals[i] = mag_vals_entry()

    _h3d.initialize_lattice()
    num_samples = _h3d.mag_v_temp(mag_vals)

    data = [[] for i in range(2)]
    for i in range(num_samples):
        data[0].append(mag_vals[i][0])
        data[1].append(mag_vals[i][1])


    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)

'''
def plot_staggered_suscept():
    global _h3d

    s_vals_entry = ctypes.c_double*2
    s_vals_arr = PDOUBLE*10000

    s_vals = s_vals_arr()

    for i in range(10000):
        s_vals[i] = s_vals_entry()

    _h3d.initialize_lattice()
    num_samples = _h3d.staggered_suscept_v_temp(s_vals)

    data = [[] for i in range(2)]
    for i in range(num_samples):
        data[0].append(s_vals[i][0])
        data[1].append(s_vals[i][1])

    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)
'''
def plot_suscept():
    global _h3d

    s_vals_entry = ctypes.c_double*2
    s_vals_arr = PDOUBLE*10000

    s_vals = s_vals_arr()

    for i in range(10000):
        s_vals[i] = s_vals_entry()

    _h3d.initialize_lattice()
    num_samples = _h3d.suscept_v_temp(s_vals)

    data = [[] for i in range(2)]
    for i in range(num_samples):
        data[0].append(s_vals[i][0])
        data[1].append(s_vals[i][1])

    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)



def plot_lattice():
    global _h3d

    col_arr = Spin_t*COLS
    row_arr = PSPIN_T*ROWS
    lattice_arr = PPSPIN_T*NUM_L

    record = lattice_arr()
    for i in range(NUM_L):
        record[i] = row_arr()
        for j in range(ROWS):
            record[i][j] = col_arr()

    _h3d.initialize_lattice()
    _h3d.record_lattice(record)

    data = [[] for i in range(NUM_L)]
    for i in range(NUM_L):
        for j in range(ROWS):
            data[i].append([])
            for k in range(COLS):
                data[i][j].append((record[i][j][k].x, record[i][j][k].y, record[i][j][k].z))

    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)


def plot_M_v_B(max_samples = 10000):
    """ Plots the mean magnetization per spin vs. external field strength.
        Data is acquired by calling the M_v_B() function in the C implementation.

        Plots will include separate sub-figures for total magnetization M, and
        individual layer magnetizations (identified as M1, M2, M3, M4).

        The number of data points returned by the C function is available as
        num_samples. The data will populate the results list, where results[i][0]
        is the external field strength of the data point and results[i][1] is the
        total magnetization. The data results[i][2:5] represent the magnetization
        values of the individual layers.
    """
    global _h3d

    # Create results list to pass by reference
    sample_arr = ctypes.c_double*(NUM_L + 3) # B, M, M1, M2, M3, M4, ... array for specific B value
    data_arr = PDOUBLE*max_samples

    results = data_arr()
    for i in range(max_samples):
        results[i] = sample_arr()

    _h3d.initialize_lattice()
    num_samples = _h3d.M_v_B(results) # pass results list by ref to be populated by C function
    ''''
    B_vals = []
    M_vals = []
    M1_vals = []
    M2_vals = []
    M3_vals = []
    M4_vals = []
    for i in range(num_samples):
        B_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])
    '''
    data = [[] for i in range((NUM_L + 3))]
    for i in range(num_samples):
        for j in range(NUM_L + 3):
            data[j].append(results[i][j])
    '''
        B_vals.append(results[i][0])
        M_vals.append(results[i][1])
        M1_vals.append(results[i][2])
        M2_vals.append(results[i][3])
        M3_vals.append(results[i][4])
        M4_vals.append(results[i][5])
    '''
    '''
    plt.subplot(511)   # Magnetization per spin, all lattices
    plt.plot(B_vals[0:int(num_samples/2)], M_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)

    plt.ylabel("M")
    #plt.xlabel("B")
    plt.xticks([])

    plt.subplot(512)   # Magnetization per spin, first lattice
    plt.plot(B_vals[0:int(num_samples/2)], M1_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M1_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{1}$")
    plt.xticks([])
    #plt.xlabel("B")

    plt.subplot(513)   # Magnetization per spin, second lattice
    plt.plot(B_vals[0:int(num_samples/2)], M2_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M2_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{2}$")
    plt.xticks([])

    #plt.xlabel("B")

    plt.subplot(514)   # Magnetization per spin, third lattice
    plt.plot(B_vals[0:int(num_samples/2)], M3_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M3_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{3}$")
    #plt.xlabel("B")
    plt.xticks([])


    plt.subplot(515)   # Magnetization per spin, fourth lattice
    plt.plot(B_vals[0:int(num_samples/2)], M4_vals[0:int(num_samples/2)], 'r')
    plt.plot(B_vals[int(num_samples/2):], M4_vals[int(num_samples/2):], 'b')
    plt.ylim(-1.0, 1.0)


    plt.ylabel("$M_{4}$")
    plt.xlabel("B (J_inter = {.05,.05,.05,0}, K = {.05,.05,.05,.05}, J_intra = {1.0,1.0,1.0,1.0})")
    #plt.xlabel("$H$ ($J_{FM} = 1.0$, $J_{AFM} = .05$, $K_{1} = .06$, $K_{234} = .05$)")
    plt.savefig("sim_results/sim_"+str(SIM_NUM)+".png")
    '''
    with open("sim_results/sim_"+str(SIM_NUM)+".pickle", 'wb') as fp:
        pickle.dump(data, fp)

    #plt.show()

if __name__ == "__main__":
    #plot_lattice()
    #autocorrelate_mag(300000, 3000)
    #plot_mag()
    plot_M_v_B()
    #plot_suscept()
    #plot_staggered_mag()
