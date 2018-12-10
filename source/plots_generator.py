import os
import matplotlib.pyplot as plt
import numpy as np
from functions import *
from radarconf import radar_conf
from scipy.constants import pi as pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
import h5py
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')


def generate_plots(radar_name, frequency, elevation, azimuth):

    """
    Import results summary from ambiguity_calculator algorithm and plot results for a DOA given by azimuth and elevation
    angles.

    :param radar_name: choose radar configuration
    :param elevation: DOA elevation angle [º]
    :param azimuth: DOA azimuth angle [º]
    """

    if not os.path.exists('../results/'+radar_name):
        os.makedirs('../results/'+radar_name)

    radar = radar_name
    k0 = k0_cal(el0=elevation, az0=azimuth)

    lambda0, xycoords = radar_conf(radar_name=radar, frequency=frequency)

    AmbiguityDistances = dict()
    intersections = dict()

    with h5py.File('../processed_data/' + radar + '/' + radar + '.h5', 'r') as hdf:
        G1 = hdf.get('trivial_calculations')
        Sn = np.array(G1.get('sensor_groups'))
        G2 = hdf.get('results_permutations')
        AmbiguityDistances['int_form_mat'] = np.array(G2.get('ambiguity_distances_INT_FORM_MAT'))
        AmbiguityDistances['wave_form'] = np.array(G2.get('ambiguity_distances_WAVE_FORM'))
        SURVIVORS = np.array(G2.get('survivors'))
        intersection_line = np.array(G2.get('intersection_line'))
        intersections['indexes'] = np.array(G2.get('intersection_indexes'))


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



    fig1, ax1 = plt.subplots()
    ax1.scatter(xycoords[:, 0] * lambda0, xycoords[:, 1] * lambda0, s=85, alpha=0.85, marker='o',
                label='Sensor position')
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
    ax3.scatter(intersection_line[0, intersections['indexes'][0]], intersection_line[1, intersections['indexes'][0]],
                s=40)
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

    plt.show()


# generate_plots(radar_name='Ydist', frequency=31, elevation=50, azimuth=270)