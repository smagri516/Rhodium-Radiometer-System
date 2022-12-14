# Sebastian Magri

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

from emissivity_data import y, e, y_vs_e
from calc_radiance import calc_radiance, calc_line_emission_L, update_yVsLye #, peak_bb_emission
from power_to_det import system_specs, radianceToIncidentPower
from plotting import plot_y_vs_Ly


T = 1677  # Temperature in K of metal sample
def main(y):

    # TODO change peak_ybb and peak_Lybb to not manual entry but to API or calculated
    peak_y = 1.727e-6 # in m -- peak spectral wavelength
    peak_Lybb = 54325700000 # in W/m^2/sr/m


    # Get bands to be working with. The outputs are y_vs_e of band
    opti_y_vs_e, opti_y, opti_e = get_band(.38e-6, .78e-6) # Optical band
    swir_y_vs_e, swir_y, swir_e = get_band(1.4e-6, 3e-6) # SWIR band
    mwir_y_vs_e, mwir_y, mwir_e = get_band(3e-6, 8e-6) # MWIR band
    lwir_y_vs_e, lwir_y, lwir_e = get_band(8e-6,14.99999e-6) # LWIR band

    # Calculate thermal emission radiances
    Lopti = calc_radiance(T, opti_y_vs_e)
    Lswir = calc_radiance(T, swir_y_vs_e)
    Lmwir = calc_radiance(T, mwir_y_vs_e)
    Llwir = calc_radiance(T, lwir_y_vs_e)
    Ltotal, y_vs_Lye = calc_radiance(T, y_vs_e, return_yvsBand=True)




    # Now we need to add spectral line radiances to the thermal
    # Spectral Line radiances

    spectral_lines = [.385e-6, .395e-6, .413e-6, .421e-6, .437e-6] # LIST DESIRED SPECTRAL LINES HERE
    broad_width = 5e-9 # SPECIFY BROADENING WIDTH
    #
    sline_indicies = get_center_emission_band_indicies(spectral_lines, broad_width)

    lines_emission_L = [] # Array of each spectral line and its total radiance
    for sline in sline_indicies:
        L_line = calc_line_emission_L(peak_Lybb, broad_width, sline[1])
        update_yVsLye(y_vs_Lye, peak_Lybb, sline[1])
        lines_emission_L.append([sline[0], L_line])



    plot_y_vs_Ly(y_vs_Lye) # This saves y_vs_Ly plot on MATLAB
    # plt.plot(y_vs_Lye) # This plots in matplotlib ,but I like MATLAB better, so I'd rather plot there



    # Calculate power onto detector system First calculate the system specifications given your detector radius,
    # spot, Fno, and reduced distance from spot to lens
    Rdet = .0001 # Radius of detector in meters
    Rspot = 0.3 # m
    Fno = 4
    zML_red = 1.9833333
    Rlens, f1, zLD_red = system_specs(Rdet,Rspot, Fno, zML_red) # outputs the radius of lens, focal length, and reduced distance from lens to det
    rhoa = .04
    rhob = .04
    aPrime = .02 # mm^-1
    tWindow = 2

    # TIME TO CALCULATE SIGNAL
    # yvsR = resonsivityVsWavelength(y)
    Llines = [row[1] for row in lines_emission_L]
    Lsignal = sum(Llines)
    incidentDetectorSignal = radianceToIncidentPower(Lsignal,rhoa, rhob, aPrime, tWindow,
                                                         Rdet, Rlens, zLD_red)





def get_band(ymin, ymax, y_vs_e=y_vs_e):
    band_indicies = np.where((y_vs_e[:, 0] >= ymin) & (y_vs_e[:, 0] <= ymax))
    band = y_vs_e[band_indicies]
    y_band = band[:, 0]
    e_band = band[:, 1]
    return band, y_band, e_band

def get_center_emission_band_indicies(slines, width, y_vs_e=y_vs_e):
    """
    Get band from a center wavelength
    :param slines: array of slines of band [wavelength1, wavelength2, etc]
    :param width: width of band
    :return:
    """

    slines_indicies = []
    for sline in slines:
        ymin = sline - (width/2)
        ymax = sline + (width/2)
        band_indicies = np.where((y_vs_e[:, 0] >= ymin) & (y_vs_e[:, 0] <= ymax))

        # Rerun if there is no index that is inside the given range
        if int(len(band_indicies[0])) == 0:

            s_index_scalar = get_center_emission_band_indicies([sline], width*2)[0]
            s_index = np.ndarray((1,))
            s_index.fill(s_index_scalar)
            print("Help")
        # If there are multiple indexes, take the middle one
        elif int(len(band_indicies[0])) > 1:
            ns = int(len(band_indicies[0]))
            mid_index_val = ns//2
            index = band_indicies[0][mid_index_val]
            s_index = np.ndarray((1,))
            s_index.fill(index)

            return s_index
        # else take the only index
        else:
            s_index = band_indicies[0]

        slines_indicies.append([sline, int(s_index[0].item())])


    return slines_indicies



main(y)
