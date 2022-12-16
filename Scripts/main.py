# Sebastian Magri
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

from emissivity_data import y, e, y_vs_e, resonsivityVsWavelength
from slice_band import defineWorkingRegion
from calc_radiance import calc_radiance, calc_line_emission_L, update_yVsLye, calc_idc #, peak_bb_emission
from power_to_det import system_specs, radianceToIncidentPower, systemCoeff, calculateSignal, calc_totalNoise, calculate_1f, calc_ijohnson
from plotting import plot_y_vs_Ly


c = 3e8 # speed of light
q = 1.6e-19 # charge of electron
kb = 1.38e-23 # Boltzmann's constant

def main(y):

    # TODO change peak_ybb and peak_Lybb to not manual entry but to API or calculated
    T = 2000  # Temperature in K of metal sample
    peak_y = 1.727e-6 # in m -- peak spectral wavelength
    peak_Lybb = 131066e6 # in W/m^2/sr/m

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
    # plt.xlabel("Wavelength")
    # plt.ylabel("Spectral Radiance (W/m^2/sr/m)")
    # plt.show()

    # Calculate power onto detector system First calculate the system specifications given your detector radius,
    # spot, Fno, and reduced distance from spot to lens

    Rdet = .0001 # Radius of detector in meters
    Rspot = 0.5 # m
    Fno = 1.8
    zML_red = 1
    Rlens, f1, zLD_red = system_specs(Rdet,Rspot, Fno, zML_red) # outputs the radius of lens, focal length, and reduced distance from lens to det
    rhoa = .04
    rhob = .04
    aPrime = .02 # mm^-1
    tWindow = 2
    Coeff = systemCoeff(rhoa, rhob, aPrime, tWindow, Rdet, Rlens, zLD_red)
    # TIME TO CALCULATE SIGNAL
    bandBottom = 380e-9
    bandTop =440e-9
    f_top = c/bandBottom
    f_bottom = c/bandTop

    y_working = defineWorkingRegion(y, bandBottom, bandTop)
    arrLength = int(len(y_working[:, 0]))  # Define array length to slice later
    y_working = resonsivityVsWavelength(y_working) # Add responsivity to matrix
    y_working = np.column_stack([y_working, y_vs_Lye[:arrLength,1]])

    # The following loop calcluates the detector signal from the emission lines
    # The detector signal is in amps
    i_Signal = 0
    for j in range(len(sline_indicies)):
        Lsignal = lines_emission_L[j][1]
        incidentDetectorPhi = radianceToIncidentPower(Lsignal,rhoa, rhob, aPrime, tWindow,
                                                         Rdet, Rlens, zLD_red)
        line_emission_phi = incidentDetectorPhi * y_working[j,2] # Take our power and then multiply
                        # Responsivity at that wavelength
        i_Signal += line_emission_phi

    # NOISE CALCULATIONS

    noises = []
    # Calculate background noise as a result of thermal emission
    working_band_y_vs_e = get_band(bandBottom, bandTop)
    # Take Ly and integrate with R(y) and multiply by Coeff
    i_bkgnd_noise = 1e-1*calculateSignal(y_working , Coeff)
    noises.append(i_bkgnd_noise)
    i_dark = 50e-12 # Dark noise arising from 150 V Reverse Bias
    noises.append(i_dark)
    bandwidth = c/380e-9 - c/450e-9# Bandpass filter
    bandwidth = 90e6-50e3
    i_shot = math.sqrt(2 * q *(i_Signal + i_bkgnd_noise + i_dark) * bandwidth) # to agree with Lye units
    noises.append(i_shot)
    Rload = 50 # ohms
    T_det = 253 # in K at detector
    i_johnson = calc_ijohnson(T_det, bandwidth, Rload)
    noises.append(i_johnson)
    idc = calc_idc(i_dark, i_Signal, i_bkgnd_noise)
    i_1f = calculate_1f(idc, f_top, f_bottom)
    noises.append(i_1f)
    i_noise = calc_totalNoise(noises)

    print("Background Noise: ", i_bkgnd_noise)
    print("Dark Noise: ", i_dark)
    print("Shot Noise", i_shot)
    print("Johnson Noise: ", i_johnson)
    print("1/f Noise: ", i_1f)
    print("Total Noise: ", i_noise)
    print("Total Signal: ", i_Signal)


    SNR = (i_Signal / i_noise)

    print("Calculated SNR: ", SNR)




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
