import numpy as np
from functions import lambda_cal


def radar_conf(radar_name, frequency):

    """
    Select radar radar configuration after name.

    :param radar_name: input by user
    :param frequency: frequency at which the radar array is being operated [Mhz]
    :return lambda0: wave lenght corresponding to radar frequency [m]
    :return xycoords: coordinates of subarray centers w.r.t. center of radar configuration in wavelengths.
    """

    if radar_name == 'JONES':
        lambda0 = lambda_cal(frequency=frequency)
        xycoords = np.array([[0, 2],
                             [0, -2.5],
                             [-2, 0],
                             [2.5, 0],
                             [0, 0]])
    elif radar_name == 'symmetric1':
        d = 3
        lambda0 = lambda_cal(frequency=frequency)
        xycoords = np.array([[d, 0],
                             [-d, 0],
                             [0, d],
                             [0, -d],
                             [d / np.sqrt(2), d/np.sqrt(2)],
                             [-d / np.sqrt(2), d/np.sqrt(2)],
                             [d / np.sqrt(2), -d/np.sqrt(2)],
                             [-d / np.sqrt(2), -d/np.sqrt(2)],
                             [0, 0]])
    elif radar_name == 'Ydist':
        lambda0 = lambda_cal(frequency=frequency)
        d = 3
        xycoords = np.array([[d * np.cos(np.radians(67.5)), d * np.sin(np.radians(67.5))],
                            [d * np.cos(np.radians(112.5)), d * np.sin(np.radians(112.5))],
                            [0, -d],
                            [0, 0]])
    else:
        raise Exception('The name you have given is not recognize by the defined radar configurations. ')

    return lambda0, xycoords
