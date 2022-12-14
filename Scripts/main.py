# Sebastian Magri
from emissivity_data import y, e, y_vs_e
import numpy as np
from calc_radiance import calc_radiance #, peak_bb_emission
from power_to_det import system_specs

T = 1677  # Temperature in K of metal sample
def main(y):

    # TODO change peak_y and peak_Ly to not manual entry but to API or calculated
    peak_y = 1.727e-6 # in m -- peak spectral wavelength
    peak_Ly = 54325700000 # in W/m^2/sr/m


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
    Ltotal = calc_radiance(T, y_vs_e)

    # Now we need to add spectral line radiances to the thermal
    # Spectral Line radiances

    spectral_lines = [.385e-6, .395e-6, .413e-6, .421e-6, .437e-6] # LIST DESIRED SPECTRAL LINES HERE
    broad_width = 5e-9 # SPECIFY BROADENING WIDTH

    L_slines = []


    for s_line in spectral_lines: # Iterate through each spectral lines and append name and value to L_slines
        # Each spectral line has 100% of the peak black body spectral radiance
        name = str(s_line)
        ye, y, e = get_center_emission_band(s_line, broad_width)
        Lline = calc_radiance(T,ye)
        L_slines.append([name, Lline])


    # Calculate power onto detector system First calculate the system specifications given your detector radius,
    # spot, Fno, and reduced distance from spot to lens
    Rdet = .0001 # Radius of detector in meters
    Rspot = 0.3 # m
    Fno = 4
    zML_red = 1.9833333
    Rlens, f1, zLD_red = system_specs(Rdet,Rspot, Fno, zML_red) # outputs the radius of lens, focal length,
    # and reduced distance from lens to det





def get_band(ymin, ymax, y_vs_e=y_vs_e):
    band_indicies = np.where((y_vs_e[:, 0] >= ymin) & (y_vs_e[:, 0] <= ymax))
    band = y_vs_e[band_indicies]
    y_band = band[:, 0]
    e_band = band[:, 1]
    return band, y_band, e_band

def get_center_emission_band(ycenter, width, y_vs_e=y_vs_e):
    """
    Get band from a center wavelength
    :param y_center: center of band
    :param width: width of band
    :return:
    """
    ymin = ycenter - (width/2)
    ymax = ycenter + (width/2)
    band_indicies = np.where((y_vs_e[:, 0] >= ymin) & (y_vs_e[:, 0] <= ymax))
    if int(len(band_indicies[0])) == 0:

        return get_center_emission_band(ycenter, width*2)

    band = y_vs_e[band_indicies]
    y_band = band[:, 0]
    e_band = band[:, 1]
    return band, y_band, e_band



main(y)
