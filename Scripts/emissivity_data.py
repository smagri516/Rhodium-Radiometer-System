import scipy.io
import numpy as np

# Load initial MATLAB file
y_vs_e = np.genfromtxt('WavelengthVsEmissivity.csv', delimiter=',', skip_header=1, missing_values=.2)

y = y_vs_e[:,0] * 1e-6 # All wavelengths from optical to lwir spectrum in
e = y_vs_e[:, 1] # All emissivities corresponding to y
y_vs_e = np.column_stack((y, e))


def resonsivityVsWavelength(yWorking):
    """
    Return a 2d array of wavelength vs Responsivity. This is an approximation based
    on the data sheet
    :param yWorking: array of wavelengths in working region (aka what the det can see)
                    the 2nd column represents WL, while the first is original index
    :return: 3column array of with R on third column
    """

    # THE R VALUES ARE HARDCODED IN FOR THIS SPECIFIC DETECTOR LOOKING AT THE SPECIFIC WLs

    #Region 1
    slope1 = (20 - 8) / (500E-9 - 380E-9)
    R1 = slope1 * (yWorking[:, 1] - 380E-9) + 8;

    # # Region 2
    # R2 = np.ones((48-24))
    # R2.fill(24.5)
    #
    # # Region 3
    # slope3 = (2 - 21) / (1000E-9 - 700E-9);
    # R3 = slope3*(yWorking[49:, 1] - 700E-9) + 21;

    R = R1
    yWorking = np.column_stack([yWorking, R])

    return yWorking







