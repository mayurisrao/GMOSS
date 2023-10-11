import healpy as hp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.special import kve
from scipy.integrate import quad

GSPAN = 100
speed_of_light  = 2.99792458e08
mass_of_electron  = 9.10938356e-31
charge_of_electron = 1.6e-19
sine_alpha  = 1.0
magnetic_field = 1e-9

scale_gam_nu = (
            3.0 * charge_of_electron * magnetic_field * sine_alpha
        ) / (4.0 * np.pi * mass_of_electron * speed_of_light)

def combined_concave_convex(df1, df2):
    combined_df = pd.concat([df1, df2], ignore_index=True)
    combined_df_sorted = combined_df.sort_values(by="PIXEL", ascending=True, ignore_index = True)
    return combined_df_sorted

def concave_func(nu, fnorm_1, fnorm_2, alpha_1, alpha_2, T_x, T_e, nu_t):
    """
    Function Name: 
    -------------
    "concave_func()"
    
    Description:
    -------------
    Calculates the value of the concave spectral model at the given frequency.

    Parameters:
    -------------
    - nu (float): The frequency/frequencies at which to calculate the model.
    - fnorm_1 (float): The first scaling/normalization parameter.
    - fnorm_2 (float): The second scaling/normalization parameter.
    - alpha_1 (float): The low frequency spectral index.
    - alpha_2 (float): The high frequency spectral index.
    - T_x (float): parameter representing optically thin free free emission index.
    - T_e (float): The electron temparature.
    - nu_t (float): The frequecny of thermal absorption turnover.

    Returns:
    -------------
    The value of the concave spectral model at the given frequency.
    """
    expo = np.exp(-1.0 * (nu_t / nu) ** 2.1)
    temp = nu**-alpha_1 + (fnorm_2 / fnorm_1) * nu**-alpha_2
    FFIT = fnorm_1 * (temp + T_x * nu**-2.1) * expo + T_e * (1.0 - expo)
    return FFIT

def modbessik2(u):  
    xnu = 5.0 / 3.0
    ORDER = xnu
    N = 1
    ARG = np.divide(1, u)
    BK = kve(ORDER, ARG)
    xrk = BK / np.exp(ARG)
    return np.divide(xrk, (u * u))

def fofx1(gamafloat, nu, C1):
    gama = float(gamafloat)
    nu_c = (gama * gama * scale_gam_nu) / 1.0e9
    x = nu / nu_c
    xl = 0.0
    xu = 1 / x
    rint, _ = quad(modbessik2, xl, xu)
    p1 = (2 * C1) - 3.0
    integ = rint * (gama ** (-1.0 * p1)) * x
    return integ

def convex_func(nu, fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t):
    """
    Function Name:
    -------------
    "convex_func()"

    Description:
    -------------
    Calculates the value of the convex spectral model at the given frequency.

    Parameters:
    -------------
    - nu (float): The frequency at which to calculate the model.
    - fnorm (float): The scaling/normalization parameter.
    - alpha1 (float): The low frequency spectral index.
    - alpha2 (float): The high frequency spectral index.
    - nu_break (float): The break frequency.
    - Tx (float): parameter representing optically thin free free emission index.
    - Te (float): The electron temparature.
    - nu_t (float): The frequecny of thermal absorption turnover.

    Returns:
    -------------
    The value of the convex spectral model at the given frequency.
    """
    nu_min = nu * 1e9 / GSPAN
    nu_max = nu * 1e9 * GSPAN
    gama_min = np.sqrt(nu_min / scale_gam_nu)
    gama_max = np.sqrt(nu_max / scale_gam_nu)
    gama_break = np.sqrt(nu_break / scale_gam_nu)
    xl = gama_min
    xu = gama_max
    xb = gama_break

    if xl > xb:
        C1 = alpha2
        rint, _ = sp.integrate.quad(fofx1, xl, xu, args=(nu, C1))
        rint *= gama_break ** (2 * C1 - 3)
    elif xu < xb:
        C1 = alpha1
        rint, _ = sp.integrate.quad(fofx1, xl, xu, args=(nu, C1))
        rint *= gama_break ** (2 * C1 - 3)
    else:
        xu = xb
        C1 = alpha1
        rint1, _ = sp.integrate.quad(fofx1, xl, xu, args=(nu, C1))
        rint1 *= gama_break ** (2 * C1 - 3)
        xl = xb
        xu = gama_max
        C1 = alpha2
        rint2, _ = sp.integrate.quad(fofx1, xl, xu, args=(nu, C1))
        rint2 *= gama_break ** (2 * C1 - 3)
        rint = rint1 + rint2

    extn = np.exp(-1.0 * (nu_t / nu) ** 2.1)
    FFIT = fnorm * ((nu**-2.0 * rint) + Tx * nu**-2.1) * extn + Te * (
        1.0 - extn
    )

    return FFIT

def all_sky_map(path_to_concave_fits, path_to_convex_fits, frequency, number_of_pixels = 3072, plot_all_sky = False):
    """
    Plots an all-sky map of a given frequency using the concave pixel fits data provided and saves the allsky map as a png file.

    Parameters:
    -------------
        DATA (str): Path to the directory where the concave pixel fits data is stored.
        frequency (int): The frequency of the data to be plotted.

    Returns:
    -------------
        None.
    """
    concave_pixel_fits_df = pd.read_csv(path_to_concave_fits)
    convex_pixel_fits_df = pd.read_csv(path_to_convex_fits)
    PIXELS = np.arange(1,number_of_pixels+1)

    # PIXELS = concave_pixel_fits_df.loc[:, "PIXEL"].values
    concave_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)
    convex_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)

    b_temp = np.array([])
    for i in range(1, number_of_pixels+1):
        if i in concave_pixel_fits_df.loc[:,"PIXEL"].values:
            b_temp = np.append(
                b_temp,
                concave_func(
                    frequency,
                    concave_pixel_fits_df.loc[i, "FNORM1"],
                    concave_pixel_fits_df.loc[i, "FNORM2"],
                    concave_pixel_fits_df.loc[i, "ALPHA_1"],
                    concave_pixel_fits_df.loc[i, "ALPHA_2"],
                    concave_pixel_fits_df.loc[i, "T_X"],
                    concave_pixel_fits_df.loc[i, "T_E"],
                    concave_pixel_fits_df.loc[i, "NU_T"],
                ),
            )
        elif i in convex_pixel_fits_df.loc[:,"PIXEL"].values:
            b_temp = np.append(
                b_temp,
                convex_func(
                    frequency,
                    convex_pixel_fits_df.loc[i, "FNORM1"],
                    convex_pixel_fits_df.loc[i, "FNORM2"],
                    convex_pixel_fits_df.loc[i, "ALPHA_1"],
                    convex_pixel_fits_df.loc[i, "ALPHA_2"],
                    convex_pixel_fits_df.loc[i, "T_X"],
                    convex_pixel_fits_df.loc[i, "T_E"],
                    convex_pixel_fits_df.loc[i, "NU_T"],
                ),
            )
        
        else:
            print(f"model value for pixel {i} was not added, hence will be taken as 0")
            b_temp = np.append(b_temp, 0)

       
    if plot_all_sky == True: 
        NSIDE = 16
        print(
            "Approximate resolution at NSIDE {} is {:.2} deg".format(
                NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
            )
        )
        NPIX = hp.nside2npix(NSIDE)
        m = np.arange(NPIX)
        from matplotlib import cm

        hp.mollview(
            b_temp,
            title=f"{frequency} GHz - Concave Fits",
            nest=True,
            norm="hist",
            cmap=cm.jet,
        )
        hp.graticule()
        plt.savefig(f"concave_fit_mollview_{frequency}_GHz.png")


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os
    load_dotenv()
    DATA = os.environ.get("DATA")
    frequency = 1200 * 1e-3  # GHz
    all_sky_map(f"{DATA}concave_pixel_fits.csv",f"{DATA}convex_pixel_fits.log", frequency)

