from scipy.constants import speed_of_light as c # [m/s]
from scipy.constants import pi as pi
import numpy as np


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


def spherical_angle(theta, phi):

    return [np.cos(theta) * np.sin(phi)], [np.sin(theta) * np.sin(phi)], [np.cos(phi)]


def linCoeff_cal(R):

    return np.sum(R ** 2, axis=0) / np.linalg.norm(R, axis= 0)


def p0_jk(j, k, R, n0, K):

    return (n0[j] + 1 - k) / K[j] / np.linalg.norm(R[:, j], axis=0) * R[:, j]


def nvec_j(j, R):
    return R[:, j] / np.linalg.norm(R[:, j], axis=0)


def k0(el0, az0):

    return [np.sin(np.radians(az0)) * np.cos(np.radians(el0)), np.cos(np.radians(az0)) * np.cos(
        np.radians(el0)), np.sin(np.radians(el0))]


def slines_intersections(k0, intersections_ind, intersection_line, cutoff_ph_ang):

    cap = np.repeat([[k0[0]], [k0[1]]], repeats=len(intersections_ind), axis=1) - intersection_line[0:2, intersections_ind]
    cap = np.sqrt(np.sum(cap ** 2, axis=0))
    cap = np.where(cap <= np.sin(cutoff_ph_ang))

    return cap


# FIXME check at symmetries and improve for performance

class AmbiguityDistances:
    def __init__(self, intersections_integers_complete):
        self.intersections_integers_complete = intersections_integers_complete

    def int_form_mat(self):

        return np.abs(self.intersections_integers_complete - np.round(self.intersections_integers_complete))

    def int_form_mean(self):

        return np.mean(self.int_form_mat(), axis=0)

    def wave_form_mat(self):

        return np.exp(1j * 2 * pi * self.intersections_integers_complete) - np.exp(1j * 2 * pi * np.round(
            self.intersections_integers_complete))

    def wave_form(self):

        return np.sqrt(np.sum(self.wave_form_mat() * np.conjugate(self.wave_form_mat()), axis=0)).real

    def explicit(self, intersection_line, intersections_ind, cap_intersections_of_slines, xy, k0):

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