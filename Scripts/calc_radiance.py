import numpy as np
import math
from scipy.integrate import simps
pi = np.pi

def calc_radiance(T, y_e_band, broadening_width=5e-9):
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

    if int(len(y)) == 1:
        L_total = Ly_e.item() * broadening_width
    else:
        L_total = simps(Ly_e, y)


    return L_total




