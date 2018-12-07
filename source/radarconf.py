import os
import matplotlib.pyplot as plt
import numpy as np
from functions import *
import itertools
from scipy.constants import pi as pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from time import time, gmtime, strftime
import h5py

# put latex in figures
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

t = time()

# Parameters
mp_tol = 1e-1

radar = 'JONES'

if radar == 'JONES':
    freq = 31
    lambda0 = lambda0(frequency=freq)
    xycoords = np.array([[0, 2],
                         [0, -2.5],
                         [-2, 0],
                         [2.5, 0],
                         [0, 0]])
elif radar == 'symmetric1':
    freq = 31
    d = 3
    lambda0 = lambda0(frequency=freq)
    xycoords = np.array([[d, 0],
                         [-d, 0],
                         [0, d],
                         [0, -d],
                         [d / np.sqrt(2), d/np.sqrt(2)],
                         [-d / np.sqrt(2), d/np.sqrt(2)],
                         [d / np.sqrt(2), -d/np.sqrt(2)],
                         [-d / np.sqrt(2), -d/np.sqrt(2)],
                         [0, 0]])
elif radar == 'Ydist':
    freq = 31
    lambda0 = lambda0(frequency=freq)
    d = 3
    xycoords = np.array([[d * np.cos(np.radians(67.5)), d * np.sin(np.radians(67.5))],
                        [d * np.cos(np.radians(112.5)), d * np.sin(np.radians(112.5))],
                        [0, -d],
                        [0, 0]])

if not os.path.exists('../results/'+radar):
    os.makedirs('../results/'+radar)

# generate and write .log file
log = open('../results/' + radar + '/' + radar + '.log', mode='w', newline='\r\n')
log.write('Solution for ' + radar + ' radar configuration is calculated. \r\n')
log.write('.................................................... \r\n')
log.write('frequency (f) = ' + str(freq) + 'MHz \r\n')
log.write('array positions: ' + '\r\n')



# CODE STARTS HERE

# Do not count the last group as it is defined as the origin
# Sn : sensor_groups
Sn = np.shape(xycoords)[0] - 1

##
R = R_cal(sensor_groups=Sn, xycoords=xycoords)

# Calculate linear coefficients
K = linCoeff_cal(R=R)

# Calculate base numbers
n0 = np.floor(2 * K)

# K's from 1 to *no + 1
k_length = 2 * n0 + 1

# Create all possible permutations from 3 first sets of planes
PERMS_number = np.prod(k_length[0:3])

iterables3 = [range(1, int(k_length[0])+1), range(1, int(k_length[1])+1), range(1, int(k_length[2])+1)]
PERMS_J = list(itertools.product(*iterables3))  # list of tuples, each tuples containing 3 int.


print('First intersection calculation: %1d permutations of 3 planes \n' % PERMS_number)
log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Started ' + str(PERMS_number)
          + ' permutations of 3 planes. \r\n')

# 3 first sets of planes

I = ([0, 1, 2])

W_matrix = np.transpose(R[:, 0:3]/np.tile(np.linalg.norm(R[:, 0:3], axis=0), (3, 1)))

pinv_norm = np.zeros((int(PERMS_number), 1))
intersection_line = np.zeros((3, int(PERMS_number)))

# Start looping over all combinations

for aux, j in enumerate(PERMS_J):

    b_vector = np.zeros((len(I), 1))

    for i in I:
        b_vector[i] = -np.dot(nvec_j(I[i], R), p0_jk(I[i], j[i], R, n0, K))

    intersection_line[:, aux], pinv_norm[aux] = mooore_penrose_solution(W=W_matrix,
                                                                        b=b_vector)


# Choose out from all possible combinations the one that are valid
intersection_line_norm = np.transpose(np.sqrt([np.sum(intersection_line ** 2, axis=0)]))

intersections = intersections_cal(pinv_norm=pinv_norm,
                                  mp_tol=mp_tol,
                                  PERMS_J=PERMS_J,
                                  intersection_line=intersection_line,
                                  R=R, norm=intersection_line_norm)

SURVIVORS = np.zeros((Sn - 2, 1))
SURVIVORS[0] = intersections['number']

log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Finished ' + str(PERMS_number)
          + ' permutations of 3 planes. \r\n')
log.write('                      There were ' + str(intersections['number']) + ' survivors \r\n')


print('Done with first three permutations')

for ii in range(3, Sn):

    PERMS_J_base = PERMS_J

    # Create all possible permutations from ii first sets of planes with only the surviving set + all new and
    # recursively iterate

    # PERMS_number = int(intersections_n * k_length[ii])
    PERMS_number = int(intersections['number'] * k_length[ii])

    print(
        'Starting plane intersections for new sensor %1d of %1d with %1d permutations on %1d remaining solutions \n' % (
            ii+1, Sn, PERMS_number, intersections['number']))

    log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Starting plane intersections for new sensor ' + str(ii+1)
              + ' of ' + str(Sn) + ' with ' + str(PERMS_number) + ' permutations on ' + str(intersections['number'])
              + ' remaining solutions \r\n')

    PERMS_J = permutations_create(permutations_base=PERMS_J_base,
                                  intersections_ind=intersections['indexes'],
                                  k_length=k_length,
                                  permutation_index=ii)

    I = range(0, ii+1)

    W_matrix = np.transpose(R[:, 0:ii+1] / np.tile(np.linalg.norm(R[:, 0:ii+1], axis=0), (3, 1)))

    pinv_norm = np.zeros((PERMS_number, 1))
    intersection_line = np.zeros((3, PERMS_number))

    PERMS_J_prime = []
    b_vector = np.zeros((len(I), PERMS_number))

    for aux, j in enumerate(PERMS_J):
        PERMS_J_prime.append(np.hstack((j[0], j[1])))
        for i in I:
            b_vector[i, aux] = -np.dot(nvec_j(I[i], R), p0_jk(I[i], np.hstack((j[0], j[1]))[i], R, n0, K))

    mooore_penrose_solution_par(W=W_matrix,
                                b_set=b_vector,
                                pnum=10,
                                niter=PERMS_number,
                                intersection_line_set=intersection_line,
                                pinv_norm_set=pinv_norm)

    PERMS_J = PERMS_J_prime

    intersections = intersections_cal(pinv_norm, mp_tol, PERMS_J, intersection_line, R)

    SURVIVORS[ii-2] = intersections['number']

    log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Finished plane intersections for new sensor ' + str(ii+1)
              + ' of ' + str(Sn) + ' with ' + str(PERMS_number) + ' permutations on ' + str(intersections['number'])
              + ' remaining solutions \r\n')
    if not ii+1 == Sn:
        log.write('                      There were ' + str(intersections['number']) + ' survivors \r\n')
    else:
        log.write('                      At the end there are ' + str(intersections['number'])
                  + ' brave survivors \r\n')


intersections_last = intersections['number']

print('Last number of valid intersections ' + str(intersections_last))

intersections_integers_complete = np.vstack((np.zeros(
    (1, np.shape(intersections['integers'])[1])), intersections['integers']))

AmbiguityDistances = dict()
AmbiguityDistances['int_form_mat'] = np.abs(intersections_integers_complete - np.round(intersections_integers_complete))
AmbiguityDistances['int_form_mean'] = np.mean(AmbiguityDistances['int_form_mat'], axis=0)
AmbiguityDistances['wave_form_mat'] = np.exp(1j * 2 * pi * intersections_integers_complete) - \
                                      np.exp(1j * 2 * pi * np.round(intersections_integers_complete))
AmbiguityDistances['wave_form'] = np.sqrt(
    np.sum(AmbiguityDistances['wave_form_mat'] * np.conjugate(AmbiguityDistances['wave_form_mat']), axis=0)).real

log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Time elapsed for main calculations was ' + str(time()-t)
          + ' seconds \r\n')

with h5py.File('../processed_data/' + radar + '.h5', 'w') as hdf:
    hdf.create_dataset('intersections_integers_complete', data=intersections_integers_complete)
    hdf.create_dataset('ambiguity_distances_INT_FORM_MAT', data=AmbiguityDistances['int_form_mat'])
    hdf.create_dataset('ambiguity_distances_INT_FORM_mean', data=AmbiguityDistances['int_form_mean'])
    hdf.create_dataset('ambiguity_distances_WAVE_FORM_MAT', data=AmbiguityDistances['wave_form_mat'])
    hdf.create_dataset('ambiguity_distances_WAVE_FORM', data=AmbiguityDistances['wave_form'])
    hdf.create_dataset('intersections_line', data=intersection_line[:, intersections['indexes'][0]])

log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Results exported to /processed_data/' + radar + '.h5 \r\n')

log.close()

# TODO create a separate function from here till the end.
k0 = k0(el0=50, az0=270)

cutoff_ph_ang = pi/2

# Find all s-lines that intersect with the cap by range check

cap_intersections_of_slines = slines_intersections(k0=k0,
                                                   intersections_ind=intersections['indexes'][0],
                                                   intersection_line=intersection_line,
                                                   cutoff_ph_ang=cutoff_ph_ang)

# From knowing what lines intercept with cap, find al possible DOA ambiguities that are part of this
ambiguity_distances_explicit, ambiguity_distances_normal, k_finds \
    = explicit(intersection_line=intersection_line,
               intersections_ind=intersections['indexes'][0],
               cap_intersections_of_slines=cap_intersections_of_slines,
               xy=xycoords,
               k0=k0)

print('Done with all other permutations')

# PLOTS

fig1, ax1 = plt.subplots()
ax1.scatter(xycoords[:, 0] * lambda0, xycoords[:, 1] * lambda0, s=85, alpha=0.85, marker='o', label='Sensor position')
ax1.grid(which='both')
ax1.set_xlabel('x [m]', fontsize=14)
ax1.set_ylabel('y [m]', fontsize=14)
ax1.set_title(r'\textbf{' + radar + ' radar sensor configuration}', fontsize=14)
fig1.savefig('../results/'+radar+'/figure1', format='eps')

fig2, ax2 = plt.subplots()
ax2.plot(range(3, Sn+1), SURVIVORS, 'bs')
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
ax2.grid(which='both')
ax2.set_xlabel('Number of sensors included', fontsize=14)
ax2.set_ylabel('Number of common plane intersections', fontsize=14)
ax2.set_title(r'\textbf{Intersection plane counts}', fontsize=14)
fig2.savefig('../results/'+radar+'/figure2', format='eps')

fig3, ax3 = plt.subplots()
ax3.scatter(intersection_line[0, intersections['indexes'][0]], intersection_line[1, intersections['indexes'][0]], s=40)
ax3.grid(which='both')
ax3.set_xlabel(r'$s_{x}$ \ [1]', fontsize=14)
ax3.set_ylabel(r'$s_{y}$ \ [1]', fontsize=14)
ax3.set_title(r'\textbf{Intersection lines}', fontsize=14)
ax3.set_aspect('equal')
fig3.savefig('../results/'+radar+'/figure3', format='eps')

fig4 = plt.figure()
ax4 = fig4.gca(projection='3d')
for S_ind in range(0, len(intersections['indexes'][0])):

    s_point = intersection_line[:, intersections['indexes'][0][S_ind]]
    s_line = np.repeat(np.transpose([s_point]), repeats=100, axis=1)
    s_line[2, :] = np.linspace(start=-np.sqrt(4 - np.dot(s_point, s_point)),
                               stop=np.sqrt(4 - np.dot(s_point, s_point)),
                               num=100,
                               endpoint=True)
    if S_ind == 25:
        ax4.plot(s_line[0, :], s_line[1, :], s_line[2, :], '.r')
    else:
        ax4.plot(s_line[0, :], s_line[1, :], s_line[2, :], '.b')
ax4.set_xlabel(r'$s_{x} \ [1]$', fontsize=14)
ax4.set_ylabel(r'$s_{y} \ [1]$', fontsize=14)
ax4.set_zlabel(r'$s_{z} \ [1]$', fontsize=14)
ax4.set_title(r'\textbf{Solution set $\Omega$}', fontsize=14)
fig4.savefig('../results/'+radar+'/figure4', format='eps')

circ_cutoff_ph_ang_x = np.sin(cutoff_ph_ang) * np.cos(np.linspace(start=0, stop=2*pi, num=100))
circ_cutoff_ph_ang_y = np.sin(cutoff_ph_ang) * np.sin(np.linspace(start=0, stop=2*pi, num=100))

fig6, ax6 = plt.subplots()
for S_ind in range(0, len(intersections['indexes'][0])):
    s_point = intersection_line[:, intersections['indexes'][0][S_ind]]
    ax6.scatter(s_point[0], s_point[1], s=20, c=(0, 0, 1), alpha=0.6)
    plt.text(s_point[0]+0.03, s_point[1]+0.03, "%0.2f" % AmbiguityDistances['wave_form'][S_ind])
ax6.plot(circ_cutoff_ph_ang_x, circ_cutoff_ph_ang_y)
ax6.scatter(k0[0], k0[1], facecolors='none', edgecolors='r', s=20)
ax6.grid(which='both')
ax6.set_xlabel(r'$s_{x}$ \ [1]', fontsize=14)
ax6.set_ylabel(r'$s_{y}$ \ [1]', fontsize=14)
ax6.set_title(r'\textbf{Solution set $\Omega$}', fontsize=14)
ax6.set_aspect('equal')
fig6.savefig('../results/'+radar+'/figure6', format='eps')

fig7, ax7 = plt.subplots()
ax7.scatter(k0[0], k0[1], facecolors='none', edgecolors='r', s=20)
for i in range(0, np.shape(k_finds)[1]):
    ax7.scatter(k_finds[0, i], k_finds[1, i], s=20, c=(0, 0, 1), alpha=0.6)
    plt.text(k_finds[0, i] + 0.1, k_finds[1, i] + 0.1, "%0.2f" % ambiguity_distances_explicit[i])
ax7.grid(which='both')
ax7.set_xlabel(r'$k_{x}$ \ [1]', fontsize=14)
ax7.set_ylabel(r'$sk_{y}$ \ [1]', fontsize=14)
ax7.set_title(r'\textbf{Solution set $\Omega (k)$}', fontsize=14)
ax7.set_aspect('equal')
fig7.savefig('../results/'+radar+'/figure7', format='eps')

print(time()-t)
plt.show()
