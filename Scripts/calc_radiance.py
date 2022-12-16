import numpy as np
import math
from scipy.integrate import simps
pi = np.pi

def calc_radiance(T, y_e_band, broadening_width=5e-9, return_yvsBand=False):
    """
    Calculate the radiance of a given light band
    :param T: Temperature in K of material
    :param y_e_band: np array of y vs e for the band
    :return: float of total radiance for a given band
    """

    y = y_e_band[:, 0] # Isolate wavelengths
    e = y_e_band[:, 1] # Isolate emissivities

    # Planck's constant (in J*s)
    h = 6.62607015e-34

    # Speed of light (in m/s)
    c = 299792458

    # Boltzmann constant (in J/K)
    k = 1.380649e-23

    num = 2 * h * c**2 / (y**5)
    denom = np.exp(h*c / (y * k * T)) - 1
    Ly_bb = num / denom


    Ly_e = Ly_bb * e

    # This is for calculation of line emission radiance
    if int(len(y)) == 1:
        L_total = Ly_e.item() * broadening_width
    else:
        L_total = simps(Ly_e, y) # otherwise, do numerical intergration over band

    # This if statement is you want to return the total radiance and the yvsLy_e array
    if return_yvsBand:
        y_vs_Lye = np.column_stack((y, Ly_e))
        return L_total, y_vs_Lye

    return L_total


def calc_line_emission_L(peakLybb, lineWidth, initialLy):
    """
    Calculate the total radiance of line emission with given width. This just uses a constant initial L
    and line width. THIS IS NOT FOR USE IN THE NON-CONSTANT EMISSION DATA
    :param peakLybb: Peak spectral radiance of a black body
    :param lineWidth: Band width to be integrated over
    :return: Total spectral radiance of line emission + thermal emission
    """

    Lytotal = peakLybb + initialLy
    return Lytotal * lineWidth

def update_yVsLye(yvsLye, peakLybb, index):
    """
    Take the yVsLye data and update the specific spectral radiances that haave line emission
    :param yvsLye: numpy array with y vs. Lye
    :param peakLybb: peak black body spectral radiance
    :param index: index to be updated
    :return: None
    """

    yvsLye[index][1] += peakLybb

def calc_idc(isig, idark, ibkg):

    idc = 0
    idc+= ibkg
    idc+= idark
    idc += isig
    return idc



"""def peak_bb_emission(y, T):
    
    # Calculate peak bb emission,
    # :param y: total emission wavelengths
    # :param T: Temperature in K
    # :return: Ly_bb_peak Peak spectral radiance for a black body
    #          y_bb_peak  Peak emission line for bb
    #

    # Planck's constant (in J*s)
    h = 6.62607015e-34

    # Speed of light (in m/s)
    c = 299792458

    # Boltzmann constant (in J/K)
    k = 1.380649e-23

    num = 2 * h * c ** 2 / (y ** 5)
    denom = np.exp(h * c / (y * k * T)) - 1
    Lybb = num / denom

    y_vs_Lybb = np.column_stack((y, Lybb))
    max_Lybb = np.max(Lybb)
    max_index = np.where((Lybb[:] == max_Lybb))

    return max_Lybb, y_vs_Lybb[max_index,0], max_index"""






