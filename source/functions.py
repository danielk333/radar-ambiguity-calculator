from scipy.constants import speed_of_light as c # [m/s]
from scipy.constants import pi as pi
import numpy as np
import itertools


def lambda0(frequency):

    return c * 1e-6 / frequency


def rRrho_cal(sensor_groups, subgroup_size, xycoords, xpos, ypos, zpos):

    r = np.zeros((3, subgroup_size, sensor_groups))

    r[0, :] = xpos[0, :sensor_groups]
    r[1, :] = ypos[1, :sensor_groups]
    r[2, :] = zpos[2, :sensor_groups]

    R = np.zeros((3, sensor_groups))
    R[0, :] = xycoords[:sensor_groups, 0]
    R[1, :] = xycoords[:sensor_groups, 1]
    R[2, :] = np.multiply(0, xycoords[:sensor_groups, 0])

    rho = np.zeros((3, subgroup_size, sensor_groups))

    for i in range(0, sensor_groups):

        rho[:, :, i] = r[:, :, i] - np.tile(np.reshape(R[:, i], (3, 1)), (1, subgroup_size))

    return r, R, rho


#TODO Am I really using this function somewhere?
def spherical_angle(theta, phi):

    return [np.cos(theta) * np.sin(phi)], [np.sin(theta) * np.sin(phi)], [np.cos(phi)]


def linCoeff_cal(R):

    return np.sum(R ** 2, axis=0) / np.linalg.norm(R, axis= 0)

#TODO test function p0_jk
def p0_jk(j, k, R, n0, K):

    return (n0[j] + 1 - k) / K[j] / np.linalg.norm(R[:, j], axis=0) * R[:, j]

#TODO test function nvec_j
def nvec_j(j, R):

    return R[:, j] / np.linalg.norm(R[:, j], axis=0)

#TODO test function moore_penrose_solution
def mooore_penrose_solution (W, b):

    Moore_Penrose_solution_check = np.linalg.multi_dot([W, np.linalg.pinv(W), b]) - b
    intersection_line = np.dot(np.linalg.pinv(W), b)[:, 0]
    pinv_norm= np.linalg.norm(Moore_Penrose_solution_check)

    return intersection_line, pinv_norm


def intersections_cal (pinv_norm, mp_tol, PERMS_J, intersection_line, R, **kwargs):

    intersections = dict()

    if 'norm' in kwargs:
        intersections['indexes'] = np.where(
            (pinv_norm < mp_tol) & (kwargs['norm'] <= 2) & (kwargs['norm'] != 0))
    else:
        intersections['indexes'] = np.where(pinv_norm < mp_tol)

    intersections['number'] = len(intersections['indexes'][0])
    intersections['permutations'] = np.array(PERMS_J)[intersections['indexes'][0]]

    integers = (
        np.repeat(np.transpose(R)[:, :, np.newaxis], repeats=intersections['number'], axis=2)[:, :, None] * np.reshape(
            intersection_line[:, intersections['indexes'][0]], (3, 1, intersections['number']))).sum(axis=1)

    intersections['integers'] = np.reshape(integers, (np.shape(integers)[0], intersections['number']))

    return intersections


def permutations_create(permutations_base, intersections_ind, k_length, permutation_index):
    iterables = [np.array(permutations_base)[intersections_ind[0], :], range(1, int(k_length[permutation_index]) + 1)]

    return list(itertools.product(*iterables))

#TODO test function k0
def k0(el0, az0):

    return [np.sin(np.radians(az0)) * np.cos(np.radians(el0)), np.cos(np.radians(az0)) * np.cos(
        np.radians(el0)), np.sin(np.radians(el0))]

#TODO test function slines_intersections
def slines_intersections(k0, intersections_ind, intersection_line, cutoff_ph_ang):

    cap = np.repeat([[k0[0]], [k0[1]]], repeats=len(intersections_ind), axis=1) - \
          intersection_line[0:2, intersections_ind]
    cap = np.sqrt(np.sum(cap ** 2, axis=0))
    cap = np.where(cap <= np.sin(cutoff_ph_ang))

    return cap

# TODO test function explicit
def explicit(intersection_line, intersections_ind, cap_intersections_of_slines, xy, k0):
    s_sel = intersection_line[:, intersections_ind[cap_intersections_of_slines]]

    aux1 = np.repeat([[k0[0]], [k0[1]]], repeats=np.shape(s_sel)[1], axis=1) - s_sel[0:2, :]
    aux2 = np.sqrt(1 - aux1[0, :] ** 2 - aux1[1, :] ** 2)

    k_finds = np.vstack((aux1, aux2))

    subgroup_signal_k0 = np.exp(-1j * 2 * pi * (xy[:, 0] * k0[0] + xy[:, 1] * k0[1]))
    subgroup_signal = np.exp(-1j * 2 * pi * (
            np.transpose(np.repeat([xy[:, 0]], repeats=np.shape(k_finds)[1], axis=0)) * np.repeat(
        [k_finds[0, :]], repeats=np.shape(xy)[0], axis=0) + np.transpose(np.repeat([xy[:, 1]], repeats=np.shape(
        k_finds)[1], axis=0)) * np.repeat([k_finds[1, :]], repeats=np.shape(xy)[0], axis=0)))

    ambiguity_distances_explicit = np.linalg.norm(np.transpose(np.repeat([subgroup_signal_k0], repeats=np.shape(
        k_finds)[1], axis=0)) - subgroup_signal, axis=0)

    ambiguity_normal_explicit = (np.transpose(np.repeat([subgroup_signal_k0], repeats=np.shape(
        k_finds)[1], axis=0)) - subgroup_signal) / np.repeat([ambiguity_distances_explicit], repeats=np.shape(
        subgroup_signal)[0], axis=0)

    return ambiguity_distances_explicit, ambiguity_normal_explicit
