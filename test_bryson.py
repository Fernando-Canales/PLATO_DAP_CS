# Test of various ways of estimating the error in the centroid coordinates
# I use Monte-Carlo simulations to establish the "ground-truth"
# Then compare with the Bryson estimate (Eq. 2 from Bryson et al.)
# and with what I think the correct version of Bryson et al. should be
#
# some parts of this code are adapted from the imagette.py from
# Fernando's repository
#
# Aaron Birch, June 15, 2024

import numpy as np

def make_2d_grid():
    M = 6
    x,y = np.meshgrid(np.arange(0, M), np.arange(0, M))
    return x,y

def measure_centroid_from_an_imagette(I):
    x, y = make_2d_grid()
    c_x = np.sum(x * I) / np.sum(I)
    c_y = np.sum(y * I) / np.sum(I)
    return c_x, c_y

def compute_bryson_variance_in_c_x(I):
    x, y = make_2d_grid()
    c_x, c_y = measure_centroid_from_an_imagette(I)
    f_tot = np.sum(I)
    return np.sum( (x**2+c_x**2) * I) / np.sum(I)**2

def compute_aaron_variance_in_c_x(I):
    x, y = make_2d_grid()
    c_x, c_y = measure_centroid_from_an_imagette(I)
    return (np.sum((x-c_x)** 2 * I)) / np.sum(I)**2

def generate_a_single_realization_of_an_imagette(A):
    x, y = make_2d_grid()
    Ibar = A*np.exp(-( (x-2)**2+(y-2)**2)/2**2)
    return np.random.poisson(lam=Ibar, size=Ibar.shape)

def run_the_MC_simulation_to_estimate_std_c_x(N_trials):
    A = 1000 # expectation value of max value of imagette
    c_x = np.zeros(N_trials)
    for j in range(N_trials):
        I = generate_a_single_realization_of_an_imagette(A)
        if j == 0:
            M = I.shape[0]
            I_save = np.zeros((N_trials,I.shape[0],I.shape[1]))
        I_save[j,:,:] = I
        c_x[j], not_used = measure_centroid_from_an_imagette(I)
    print('mean c_x from the MC experiment = %.4f' % np.mean(c_x))
    return np.std(c_x), I_save
#===============

N_trials = 10000

MC_std_c_x, all_I = run_the_MC_simulation_to_estimate_std_c_x(N_trials)

average_I = np.mean(all_I, axis=0)
var_I = np.var(all_I, axis=0)

c_x, not_used = measure_centroid_from_an_imagette(average_I)
print('c_x measured from mean imagette = %.4f' % c_x)

np.set_printoptions(precision=1)
x, y = make_2d_grid()
print('x grid')
print(x)
print('mean imagette:')
print(average_I)
print('variance of imagette:')
print(var_I)

var_c_x_bryson = compute_bryson_variance_in_c_x(average_I)
var_c_x_aaron = compute_aaron_variance_in_c_x(average_I)
#---
print('Monte-Carlo std c_x = %.4f' % MC_std_c_x)
print('Bryson std c_v = %.4f' % np.sqrt(var_c_x_bryson))
print('Aaron std c_v = %.4f ' % np.sqrt(var_c_x_aaron))




