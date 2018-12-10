import os
import numpy as np
from functions import *
import itertools
from scipy.constants import pi as pi
from time import time, gmtime, strftime
import h5py


def ambiguities_calculate(radar_name, frequency):
    """
    Calculates the solution for the radar ambiguity problem by implementing developed by Daniel Kastinen in his paper
    "Determining all ambiguities in direction of arrival measured by radar systems". The output of the calculations is
    summarized in a .h5 file with HDF5 format and saved into the /processed_data folder

    :param radar_name: Name of the radar to be studied out of a list of defined configurations.
    :param frequency: Operating frequency [MHz]
    """

    radar_name = radar_name

    lambda0, xycoords = radar_conf(radar_name=radar_name, frequency=frequency)

    if not os.path.exists('../processed_data/'+radar_name):
        os.makedirs('../processed_data/'+radar_name)

    # generate and write .log file
    log = open('../processed_data/' + radar_name + '/' + radar_name + '.log', mode='w', newline='\r\n')
    log.write('Solution for ' + radar_name + ' radar configuration is calculated. \r\n')
    log.write('.................................................... \r\n')
    log.write('frequency (f) = ' + str(frequency) + 'MHz \r\n')
    log.write('array positions: ' + '\r\n')

    # CODE STARTS HERE
    mp_tol = 1e-1

    # Do not count the last group as it is defined as the origin
    # Sn : sensor_groups
    Sn = np.shape(xycoords)[0] - 1

    ##
    R = R_cal(sensor_groups=Sn, xycoords=xycoords)

    # Calculate linear coefficients
    K = linCoeff_cal(R=R)

    # Calculate base numbers
    n0 = np.floor(2 * K)

    # K's from 1 to no + 1
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
            'Starting plane intersections for new sensor %1d of %1d with %1d permutations on %1d remaining solutions \n'
            % (ii+1, Sn, PERMS_number, intersections['number']))

        log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())
                  + ' : Starting plane intersections for new sensor ' + str(ii+1)
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

        log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())
                  + ' : Finished plane intersections for new sensor ' + str(ii+1)
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
    AmbiguityDistances['int_form_mat'] = np.abs(intersections_integers_complete
                                                - np.round(intersections_integers_complete))
    AmbiguityDistances['int_form_mean'] = np.mean(AmbiguityDistances['int_form_mat'], axis=0)
    AmbiguityDistances['wave_form_mat'] = np.exp(1j * 2 * pi * intersections_integers_complete) - \
                                          np.exp(1j * 2 * pi * np.round(intersections_integers_complete))
    AmbiguityDistances['wave_form'] = np.sqrt(
        np.sum(AmbiguityDistances['wave_form_mat'] * np.conjugate(AmbiguityDistances['wave_form_mat']), axis=0)).real

    log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ' : Time elapsed for main calculations was ' + str(time()-t)
              + ' seconds \r\n')

    with h5py.File('../processed_data/' + radar_name + '/' + radar_name + '.h5', 'w') as hdf:

        G1 = hdf.create_group('trivial_calculations')
        G1.create_dataset('sensor_groups', data=Sn)
        G1.create_dataset('subgroup_phase_center', data=R)
        G1.create_dataset('linear_coefficients', data=K)
        G1.create_dataset('base_numbers', data=n0)
        G1.create_dataset('k_length', data=k_length)

        G2 = hdf.create_group('results_permutations')
        G2.create_dataset('intersections_integers_complete', data=intersections_integers_complete)
        G2.create_dataset('ambiguity_distances_INT_FORM_MAT', data=AmbiguityDistances['int_form_mat'])
        G2.create_dataset('ambiguity_distances_INT_FORM_mean', data=AmbiguityDistances['int_form_mean'])
        G2.create_dataset('ambiguity_distances_WAVE_FORM_MAT', data=AmbiguityDistances['wave_form_mat'])
        G2.create_dataset('ambiguity_distances_WAVE_FORM', data=AmbiguityDistances['wave_form'])
        G2.create_dataset('intersection_line', data=intersection_line)
        G2.create_dataset('intersection_indexes', data=intersections['indexes'])
        G2.create_dataset('survivors', data=SURVIVORS)

    log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())
              + ' : Results exported to /processed_data/' + radar_name + '.h5 \r\n')

    log.close()

    print(time()-t)


if __name__ == "__main__":

    t = time()
    ambiguities_calculate(radar_name='JONES', frequency=31)
