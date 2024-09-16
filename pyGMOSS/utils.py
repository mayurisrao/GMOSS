import healpy as hp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.special import kve
from scipy.integrate import quad
import pygmoss.pygmoss_consts as pygc

GSPAN = pygc.GSPAN
speed_of_light = pygc.SPEED_OF_LIGHT
mass_of_electron = pygc.MASS_Of_ELECTRON
charge_of_electron = pygc.CHARGE_OF_ELECTRON
sine_alpha = pygc.SINE_ALPHA
magnetic_field = pygc.MAGNETIC_FIELD
scale_gam_nu = pygc.SCALE_GAM_NU


def combine_df(df1, df2):
    """
    Function Name:
    -------------
    "combine_df"

    Description:
    -------------
    Combines two data frames by it's pixel value. Make sure the Pixel values are in the 0th column.

    Parameters:
    -------------
    df1(pd.DataFrame): df1 to be combined.
    df2(pd.DataFrame): df2 to be combined.

    Returns:
    -------------
    The combined data frame in the form of a pandas.DataFrame.
    """
    df1.columns = range(df1.shape[1])
    df2.columns = range(df2.shape[1])
    combined_df = pd.concat([df1, df2], ignore_index=True)
    combined_df_sorted = combined_df.sort_values(
        by="PIXEL", ascending=True, ignore_index=True
    )
    return combined_df_sorted


def concave_func(nu, fnorm_1, fnorm_2, alpha_1, alpha_2, T_x, T_e, nu_t) -> float:
    """
    Function Name:
    -------------
    "concave_func()"

    Description:
    -------------
    Calculates the value of the concave spectral model at the given frequency.

    Parameters:
    -------------
    - nu (float): The frequency/frequencies at which to calculate the model(in GHz).
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


def convex_func(nu, fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t) -> float:
    """
    Function Name:
    -------------
    "convex_func()"

    Description:
    -------------
    Calculates the value of the convex spectral model at the given frequency.

    Parameters:
    -------------
    - nu (float): The frequency at which to calculate the model(in GHz).
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
    FFIT = fnorm * ((nu**-2.0 * rint) + Tx * nu**-2.1) * extn + Te * (1.0 - extn)

    return FFIT


def all_sky_map(
    path_to_concave_fits,
    path_to_convex_fits,
    frequency,
    number_of_pixels=3072,
    plot_all_sky=False,
) -> pd.Series:
    """
    Function name
    -------------
    "all_sky_map()"

    Description
    -------------
    Generates the all sky map for a given frequency and writes it to a pandas.Series using
    the concave pixel fits and convex pixel fits data provided. The function can also
    plot the all sky map(in .png format) of the input frequency using Healpy.

    Parameters:
    -------------
    path_to_concave_fits (str): Path to the directory where the concave pixel fits file
    path_to_convex_fits (str): Path to the directory where the concave pixel fits file
    frequency (int): The frequency (in GHz) to plot the all sky map.
    number_of_pixels(int): number of pixels and it is defaulted to 3072.
    plot_all_sky(bool): defaulted to False, set the parameter to True if all sky plot using Healpy is desired.

    Returns:
    -------------
        A pandas.Series containing the all sky map at the desired frequency.

    Example Usage:
    -------------
    ```
    load_dotenv()
    DATA = os.environ.get("DATA")
    frequency = 1200 * 1e-3  # GHz
    series = all_sky_map(
        f"{DATA}concave_pixel_fits.csv",
        f"{DATA}convex_pixel_fits.log",
        frequency,
        number_of_pixels=3072,
        plot_all_sky=True,
    )
    series.to_csv(f"{DATA}all_sky_map_{frequency/1e-3}.csv", index=False, header=False)
    ```
    """
    concave_pixel_fits_df = pd.read_csv(path_to_concave_fits)
    convex_pixel_fits_df = pd.read_csv(path_to_convex_fits)
    PIXELS = np.arange(1, number_of_pixels + 1)

    concave_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)

    convex_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)

    b_temp = np.array([])
    print(concave_pixel_fits_df.index.values)
    print(convex_pixel_fits_df.index.values)
    for i in range(1, number_of_pixels + 1):
        print(f"Pixel number: {i}")
        if i in concave_pixel_fits_df.index.values:
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
        elif i in convex_pixel_fits_df.index.values:
            b_temp = np.append(
                b_temp,
                convex_func(
                    frequency,
                    convex_pixel_fits_df.loc[i, "FNORM"],
                    convex_pixel_fits_df.loc[i, "ALPHA1"],
                    convex_pixel_fits_df.loc[i, "ALPHA2"],
                    convex_pixel_fits_df.loc[i, "NU_BREAK"],
                    convex_pixel_fits_df.loc[i, "TX"],
                    convex_pixel_fits_df.loc[i, "TE"],
                    convex_pixel_fits_df.loc[i, "NU_T"],
                ),
            )

        else:
            print(f"model value for pixel {i} was not added, hence will be taken as 0")
            b_temp = np.append(b_temp, 0)

    series_plot = pd.Series(b_temp)

    if plot_all_sky == True:
        NSIDE = 16
        # print(
        #     "Approximate resolution at NSIDE {} is {:.2} deg".format(
        #         NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
        #     )
        # )
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
        plt.savefig(f"mollview_{frequency}_GHz.png")

    return series_plot


def btemp_vs_frequency(
    path_to_concave_fits,
    path_to_convex_fits,
    path_to_the_brightness_temperature_file,
    pixel_val,
):
    """
    Function name
    -------------
    "btemp_vs_frequency"

    Description
    -------------
    Plots brightness temperature vs frequency for a specified pixel val using the concave or convex model depending on the pixel and saves it as a .png file.

    Parameters:
    -------------
    path_to_concave_fits (str): Path to the the concave pixel fits file(csv file)
    path_to_convex_fits (str): Path to the concave pixel fits file(csv file)
    path_to_the_brightness_temperature_file: Path to the brightness temperature file.
    number_of_pixels(int): number of pixels(it is defaulted to 3072). Recommended not change the default value.
    pixel_val(int): desired pixel value to be plotted.

    Returns:
    -------------
        none

    Example Usage:
    -------------
    ```
    load_dotenv()
    DATA = os.environ.get("DATA")
    pixel = 2060
    btemp_vs_freqency(
        f"{DATA}concave_pixel_fits.csv",
        f"{DATA}convex_pixel_fits.log",
        f"{DATA}brightness_temp_per_pixel.csv",
        pixel_val=pixel,
    )
    ```
    """
    concave_pixel_fits_df = pd.read_csv(path_to_concave_fits)
    convex_pixel_fits_df = pd.read_csv(path_to_convex_fits)
    brightness_temperature_df = pd.read_csv(path_to_the_brightness_temperature_file)
    concave_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)
    convex_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)
    brightness_temperature_df.set_index("PIXEL", inplace=True, drop=True)

    brightness_temperature_cols = brightness_temperature_df.columns.values
    print(brightness_temperature_cols)

    frequencies = np.array([])
    frequencies_string = np.array([])
    for i in range(len(brightness_temperature_cols)):
        values = brightness_temperature_cols[i][:-3]
        frequencies_string = np.append(frequencies_string, values)
        frequencies = np.append(frequencies, int(values))

    frequencies = frequencies * 1e-3  # GHz
    print(frequencies)

    b_temps = brightness_temperature_df.loc[pixel_val].values
    print(b_temps)

    x_ax = np.linspace(2.2e-2, 23, 100)
    y_ax = np.array([])
    convexity = "null"

    if pixel_val in concave_pixel_fits_df.index.values:
        for fre in x_ax:
            y_ax = np.append(
                y_ax,
                concave_func(
                    fre,
                    concave_pixel_fits_df.loc[pixel_val, "FNORM1"],
                    concave_pixel_fits_df.loc[pixel_val, "FNORM2"],
                    concave_pixel_fits_df.loc[pixel_val, "ALPHA_1"],
                    concave_pixel_fits_df.loc[pixel_val, "ALPHA_2"],
                    concave_pixel_fits_df.loc[pixel_val, "T_X"],
                    concave_pixel_fits_df.loc[pixel_val, "T_E"],
                    concave_pixel_fits_df.loc[pixel_val, "NU_T"],
                ),
            )
        convexity = "concave"

    if pixel_val in convex_pixel_fits_df.index.values:
        for fre in x_ax:
            y_ax = np.append(
                y_ax,
                convex_func(
                    fre,
                    convex_pixel_fits_df.loc[pixel_val, "FNORM"],
                    convex_pixel_fits_df.loc[pixel_val, "ALPHA1"],
                    convex_pixel_fits_df.loc[pixel_val, "ALPHA2"],
                    convex_pixel_fits_df.loc[pixel_val, "NU_BREAK"],
                    convex_pixel_fits_df.loc[pixel_val, "TX"],
                    convex_pixel_fits_df.loc[pixel_val, "TE"],
                    convex_pixel_fits_df.loc[pixel_val, "NU_T"],
                ),
            )
        convexity = "convex"

    plt.plot(x_ax, y_ax, label="optimized model")

    plt.plot(frequencies, b_temps, "r*")
    plt.grid()
    plt.legend()
    plt.xlabel("log b_temp")
    plt.ylabel("log frequency")
    plt.title(f"{convexity} plots for pixel {pixel_val}")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(f"{convexity}_pixel_{pixel_val}.png")


def concave_convex_fits_generator(
    path_to_input_param_file_from_GMOSS,
    convert_to_csv=True,
) -> pd.DataFrame:
    """
    Function name
    -------------
    "concave_convex_fits_generator"

    Description
    -------------
    Generates Concave and Convex fits from GMOSS(Dr. Mayuri et. al.).

    Parameters:
    -------------
    - path_to_input_param_file_from_GMOSS (str): Path to the the concave and convex fits file generated from GMOSS(.txt format)
    - convert_to_csv (bool): The default value of the parameter is set to True, when True it generates two files convex_pixel_fits
      and concave_pixel_fits. If the parmeter is set to False no files are generated, the function returns the convex_fits and concave_fits
      as two data frames

    Returns:
    -------------
        if convert_to_csv = True then the function generates concave and convex fits from the input_parameter_file
        else if convert_to_csf = False then the function returns two data frames containing the convex fits and concave fits.


    Example Usage:
    -------------
    ```
    from pygmoss.utils import concave_convex_fits_generator

    all_sky_map_generator_from_optimized_input_parm_file("GMOSS_params_allpix.txt")

    ```
    """
    df_1_for_convex = pd.read_csv(
        path_to_input_param_file_from_GMOSS, delimiter=" ", header=None
    )
    df_2_for_concave = df_1_for_convex

    convex_column_names = [
        "PIXEL",
        "flag",
        "norm_param",
        "alpha1_param",
        "dalpha_param",
        "nubrk_param",
        "Tx_param",
        "Te_param",
        "nut_param",
        "chisquared",
    ]
    df_1_for_convex.columns = convex_column_names
    convex_df = df_1_for_convex[df_1_for_convex.flag == 0]

    concave_column_names = [
        "PIXEL",
        "flag",
        "fnorm1_param",
        "alpha1_param",
        "dalpha_param",
        "fnorm2_param",
        "Tx_param",
        "Te_param",
        "nut_param",
        "chisquared",
    ]

    df_2_for_concave.columns = concave_column_names
    concave_df = df_2_for_concave[df_2_for_concave.flag == 1]

    # Creating df for convex

    PIXEL = pd.Series(convex_df["PIXEL"].values, name="PIXEL")

    FNORM = pd.Series(10 ** convex_df["norm_param"].values, name="FNORM")
    ALPHA1 = pd.Series(10 ** convex_df["alpha1_param"].values, name="ALPHA1")
    d_alph = pd.Series(10 ** convex_df["dalpha_param"].values)
    ALPHA2 = ALPHA1 + d_alph
    ALPHA2.name = "ALPHA2"
    NU_BREAK = pd.Series(10 ** convex_df["nubrk_param"].values, name="NU_BREAK")
    TX = pd.Series(10 ** convex_df["Tx_param"].values, name="TX")
    TE = pd.Series(10 ** convex_df["Te_param"].values, name="TE")
    NU_T = pd.Series(10 ** convex_df["nut_param"].values, name="NU_T")
    CHISQUARED = pd.Series(convex_df["chisquared"].values, name="CHISQUARED")

    convex_refined_df = pd.concat(
        [PIXEL, FNORM, ALPHA1, ALPHA2, NU_BREAK, TX, TE, NU_T, CHISQUARED], axis=1
    )

    # Creating df for concave

    PIXEL = pd.Series(concave_df["PIXEL"].values, name="PIXEL")
    FNORM1 = pd.Series(concave_df["fnorm1_param"].values, name="FNORM1")
    FNORM2 = pd.Series(concave_df["fnorm2_param"].values, name="FNORM2")
    ALPHA_1 = pd.Series(10 ** concave_df["alpha1_param"].values, name="ALPHA_1")
    d_alph = pd.Series(10 ** concave_df["dalpha_param"].values)
    ALPHA_2 = ALPHA_1 - d_alph
    ALPHA_2.name = "ALPHA_2"
    T_X = pd.Series(10 ** concave_df["Tx_param"].values, name="T_X")
    T_E = pd.Series(10 ** concave_df["Te_param"].values, name="T_E")
    NU_T = pd.Series(10 ** concave_df["nut_param"].values, name="NU_T")
    CHISQUARED = pd.Series(concave_df["chisquared"].values, name="CHISQUARED")

    concave_refined_df = pd.concat(
        [PIXEL, FNORM1, FNORM2, ALPHA_1, ALPHA_2, T_X, T_E, NU_T, CHISQUARED], axis=1
    )
    if convert_to_csv == True:
        concave_refined_df.to_csv("concave_pixel_fits.csv", index=False)
        convex_refined_df.to_csv("convex_pixel_fits.csv", index=False)
    else:
        return (convex_refined_df, concave_refined_df)


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.environ.get("DATA")
    frequency = 1200 * 1e-3  # GHz
    series = all_sky_map(
        f"{DATA}concave_pixel_fits.csv",
        f"{DATA}convex_pixel_fits.log",
        frequency,
        number_of_pixels=3072,
        plot_all_sky=True,
    )
    series.to_csv(f"{DATA}all_sky_map_{frequency/1e-3}.csv", index=False, header=False)

    pixel = 2060
    btemp_vs_frequency(
        f"{DATA}concave_pixel_fits.csv",
        f"{DATA}convex_pixel_fits.log",
        f"{DATA}brightness_temp_per_pixel.csv",
        pixel_val=pixel,
    )

    # to just plot all sky maps
    concave_convex_fits_generator("GMOSS_params_allpix.txt")
    frequency = 1200 * 1e-3  # GHz
    sky_map_1200 = all_sky_map(
        "concave_pixel_fits.csv",
        "convex_pixel_fits.csv",
        0.05,
        number_of_pixels=3072,
        plot_all_sky=False,
    )
    sky_map_1200.to_csv(f"all_sky_map_{frequency/1e-3}.csv", index=False, header=False)
