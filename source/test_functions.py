import unittest
from functions import *
import numpy as np


class TestFunctions(unittest.TestCase):

    def test_lambda0(self):

        self.assertEqual(lambda_cal(frequency=31), 9.670724451612903)


    def test_intersections_cal(self):

        pinv_norm = np.array([0.06, 0.01, 0.03, 0.15, 0.04])
        intersection_line_norm = np.array([0.75, 0, 3, 1, 0.75])
        PERMS_J = [([0, 0], 1), ([1, 1], 1), ([2, 2], 2), ([3, 3], 3), ([4, 4], 4)]
        R = [[0, 0, -2, 2.5], [2, -2.5, 0, 0], [0, 0, 0, 0]]
        intersection_line = np.array([[2, 1.5, 1, 0.5, 0],
                                      [-2.22e-16, -2.22e-16, -2.22e-16, -2.22e-16, -2.22e-16],
                                      [0, 0, 0, 0, 0]])

        intersections = intersections_cal(pinv_norm=pinv_norm, mp_tol=0.1, PERMS_J=PERMS_J,
                                          intersection_line=intersection_line, R=R, norm=intersection_line_norm)

        np.testing.assert_array_equal(intersections['indexes'], [[0, 4]])
        np.testing.assert_array_equal(intersections['permutations'], [PERMS_J[0], PERMS_J[4]])
        self.assertEqual(intersections['number'], 2)
        self.assertEqual(np.shape(intersections['integers'])[1], 2)

    def test_permutations_create(self):

        PERMS_J_base = [[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 1, 4]]
        indexes = [[1, 3], [0, 0]]
        k_length = [2, 3]
        ii = 1

        PERMS_J = permutations_create(permutations_base=PERMS_J_base, intersections_ind=indexes, k_length=k_length,
                                      permutation_index=ii)

        self.assertEqual(len(PERMS_J), len(np.array(PERMS_J_base)[indexes[0], :]) * k_length[1])

    def test_R_cal(self):

        Sn = 4

        xycoords = np.array([[0, 2], [0, -2.5], [-2, 0], [2.5, 0], [0, 0]])

        R = R_cal(sensor_groups=Sn,  xycoords=xycoords)

        np.testing.assert_array_equal(np.array(np.shape(R)), np.array([3, Sn]))

    def test_linCoeff_cal(self):

        R = np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3]])

        K = linCoeff_cal(R)

        self.assertEqual(len(K), 3)
        np.testing.assert_array_almost_equal(K, np.array([np.sqrt(3*1), np.sqrt(3*4), np.sqrt(3*9)]))

    def test_k0_cal(self):

        a = k0_cal(el0=50, az0=270)

        np.testing.assert_array_almost_equal(np.array(a), np.array([-0.642788, 0, 0.76604]), decimal=4)

    def test_p0_jk(self):

        j = 1
        k = 1
        R = np.array([[1, 2], [3, 4]])
        n0 = np.array([1, 2, 3])
        K = np.array([4, 5, 6])

        a = p0_jk(j, k, R, n0, K)

        np.testing.assert_array_almost_equal(a, np.array([0.178885438199983, 0.357770876399966]), decimal=8)

    def test_nvec_j(self):

        j = 1
        R = np.array([[1, 2], [3, 4]])

        a = nvec_j(j, R)

        np.testing.assert_array_almost_equal(a, np.array([0.447213595499958, 0.894427190999916]), decimal=8)

    def test_explicit(self):

        intersection_line = np.array([[.1, .2, .3, 1], [.4, .5, .6, .2], [.7, .8, .9, .13]])
        intersections_ind = np.array([0, 2, 3])
        cap_intersections_of_slines = [0, 2]
        xy = np.array([[-0.25, 0.25], [0.25, -0.25]])
        k0 = [.8, .15]

        distances, normal, k_finds = explicit(intersection_line, intersections_ind, cap_intersections_of_slines, xy, k0)

        np.testing.assert_array_almost_equal(distances, np.array([0.6602832, 1.66250775]), decimal=8)
        np.testing.assert_array_almost_equal(normal, np.array([[0.67249851-0.21850801j, -0.27059805+0.65328148j],
                                                               [0.67249851+0.21850801j, -0.27059805-0.65328148j]]),
                                             decimal=8)

    def test_slines_intersections(self):

        k0 = [0.15, 0.20]
        intersections_ind = [0, 2, 4]
        intersection_line = np.array([[.1, .2, .3, .4, .5],
                                      [.6, .7, .8, .9, .8],
                                      [.7, .6, .5, .4, .3],
                                      [.2, .1, .2, .3, .4]])
        cutoff_ph_ang = np.radians(40)

        cap = slines_intersections(k0, intersections_ind, intersection_line, cutoff_ph_ang)

        np.testing.assert_array_equal(cap, np.array([[0, 1]]))


if __name__ == '__main__':
    unittest.main()
