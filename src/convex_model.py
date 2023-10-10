import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.integrate as integrate
import scipy as sp
from scipy.special import kve
from scipy.integrate import quad
import time
import logging


class ConvexModel:
    """
    Class Name:
    -------------
    ConvexModel

    Desciption:
    -------------
    This class represents the convex spectral model.

    Parameters:
    -------------
    - path_to_brightness_temperature_file : str
        The path to the brightness temperature file.
    - path_to_convexity_file : str
        The path to the convexity file.
    - path_to_initial_parameters_file : str
        The path to the initial parameters file.
    - path_to_convex_fits_file : str
        The path to the convex_fits file to log data. The input file should be in .log format
    - GSPAN : int, optional
        The GSPAN value (default is 100).
    - NHPIX : int, optional
        The number of pixels (default is 3072).
    - speed_of_light : float, optional
        The speed of light in meters per second (default is 2.99792458e08).
    - mass_of_electron : float, optional
        The mass of the electron in kilograms (default is 9.10938356e-31).
    - charge_of_electron : float, optional
        The charge of the electron in coulombs (default is 1.6e-19).
    - sine_alpha : float, optional
        The sine of alpha (default is 1.0).
    - magnetic_field : float, optional
        The magnetic field in tesla (default is 1e-9).
    - verbose : bool, optional
        Whether to output verbose information (default is False).
    """

    def __init__(
        self,
        path_to_brightness_temperature_file: str,
        path_to_convexity_file: str,
        path_to_initial_parameters_file: str,
        path_to_convex_fits_file: str,
        GSPAN: float = 100,
        speed_of_light: float = 2.99792458e08,
        mass_of_electron: float = 9.10938356e-31,
        charge_of_electron: float = 1.6e-19,
        sine_alpha: float = 1.0,
        magnetic_field: float = 1e-9,
        verbose=False,
        start_pix: int = 0,
        stop_pix: int = 0,
    ):
        self.brightness_temperature_df = pd.read_csv(
            path_to_brightness_temperature_file
        )
        self.convexity_df = pd.read_csv(path_to_convexity_file)
        self.initial_parameters_df = pd.read_csv(path_to_initial_parameters_file)
        self.convex_fits_file = path_to_convex_fits_file
        self.GSPAN = GSPAN
        self.speed_of_light = speed_of_light
        self.mass_of_electron = mass_of_electron
        self.charge_of_electron = charge_of_electron
        self.sine_alpha = sine_alpha
        self.magnetic_field = magnetic_field
        self.verbose = verbose
        self.start_pix = start_pix
        self.stop_pix = stop_pix

        self.frequencies = np.array([22.0, 45.0, 150.0, 408.0, 1420.0, 23000.0])
        self.frequencies *= 1e-3
        self.frequencies_string = np.array(
            ["22", "45", "150", "408", "1420", "23000"]
        )  # MHz
        self.df = pd.merge(
            self.brightness_temperature_df, self.convexity_df, on="PIXEL"
        )
        self.df_convex = self.df.loc[self.df.loc[:, "Concave/Convex"] == "Convex", :]
        self.df_convex.reset_index(inplace=True, drop=True)
        self.pixels = self.df_convex.loc[:, "PIXEL"].values
        self.df_convex.set_index("PIXEL", inplace=True)

        self.scale_gam_nu = (
            3.0 * self.charge_of_electron * self.magnetic_field * self.sine_alpha
        ) / (4.0 * np.pi * self.mass_of_electron * self.speed_of_light)

    def modbessik2(self, u):  
        xnu = 5.0 / 3.0
        ORDER = xnu
        N = 1
        ARG = np.divide(1, u)
        BK = kve(ORDER, ARG)
        xrk = BK / np.exp(ARG)
        return np.divide(xrk, (u * u))

    def fofx1(self, gamafloat, nu, C1):
        gama = float(gamafloat)
        nu_c = (gama * gama * self.scale_gam_nu) / 1.0e9
        x = nu / nu_c
        xl = 0.0
        xu = 1 / x
        rint, _ = quad(self.modbessik2, xl, xu)
        p1 = (2 * C1) - 3.0
        integ = rint * (gama ** (-1.0 * p1)) * x
        return integ

    def convex_func(self, nu, fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t):
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
        nu_min = nu * 1e9 / self.GSPAN
        nu_max = nu * 1e9 * self.GSPAN
        gama_min = np.sqrt(nu_min / self.scale_gam_nu)
        gama_max = np.sqrt(nu_max / self.scale_gam_nu)
        gama_break = np.sqrt(nu_break / self.scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            rint, _ = sp.integrate.quad(self.fofx1, xl, xu, args=(nu, C1))
            rint *= gama_break ** (2 * C1 - 3)
        elif xu < xb:
            C1 = alpha1
            rint, _ = sp.integrate.quad(self.fofx1, xl, xu, args=(nu, C1))
            rint *= gama_break ** (2 * C1 - 3)
        else:
            xu = xb
            C1 = alpha1
            rint1, _ = sp.integrate.quad(self.fofx1, xl, xu, args=(nu, C1))
            rint1 *= gama_break ** (2 * C1 - 3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            rint2, _ = sp.integrate.quad(self.fofx1, xl, xu, args=(nu, C1))
            rint2 *= gama_break ** (2 * C1 - 3)
            rint = rint1 + rint2

        extn = np.exp(-1.0 * (nu_t / nu) ** 2.1)
        FFIT = fnorm * ((nu**-2.0 * rint) + Tx * nu**-2.1) * extn + Te * (
            1.0 - extn
        )

        return FFIT

    def chisquared(self, params, xobs, yobs):
        chisq = 0
        yobs_length = len(yobs)
        for i in range(yobs_length):
            chisq = (
                chisq + ((yobs[i] - self.convex_func(xobs[i], *params)) / yobs[i]) ** 2
            )
        return chisq / 6

    def fit(self):
        """
        Function Name:
        -------------
        "fit()"

        Description:
        -------------
        This method fits the convex model to the brightness temperature data for each pixel in the dataset and write the fitted/optimized parameters into a log file.

        Parameters:
        -------------
        This method takes no parameters.

        Returns:
        -------------
        - This method writes the fitted/optimized parameters parameters into the user specified .log file(If the log file does not exist the program creates the .log file).
          The .log file has the following columns:
        - "PIXEL": The pixel number for which the parameters were fitted.
        - "FNORM": The scaling/normalization parameter.
        - "ALPHA1": The low frequency spectral index.
        - "ALPHA2": The high frequency spectral index.
        - "NU_BREAK": The break frequency.
        - "TX": parameter representing optically thin free free emission index.
        - "TE": The electron temparature.
        - "NU_T": The frequecny of thermal absorption turnover.

        Example:
        -------------
        ```python
        model = ConvexModel("brightness_temperature.csv", "convexity.csv", "convex_model_initial_parameters.csv", "convex_pixel_fits.log")
        model.fit()
        ```
        In this example, the "ConvexModel" class is instantiated with the "brightness_temperature.csv",  "convexity.csv", "convex_model_initial_parameters.csv" and "convex_pixel_fits.log" files. Then, the `fit()` method is called to fit the convex model to the data
        and write the fitted parameters into "convex_pixel_fits.log" file. The resulting DataFrame is assigned to the `result_df` variable.
        """
        logging.basicConfig(
            filename=self.convex_fits_file,
            level=logging.INFO,
            format="",
        )

        if os.stat(self.convex_fits_file).st_size == 0:
            logging.info("PIXEL,FNORM,ALPHA1,ALPHA2,NU_BREAK,TX,TE,NU_T,CHISQUARE")
            start_pixel_index = 0

        else:
            frame = pd.read_csv(self.convex_fits_file)
            last_end_pixel = frame.loc[:, "PIXEL"].values[-1]
            (start_pixel_index,) = np.where(self.pixels == last_end_pixel)[0] + 1

        stop_pixel = self.pixels[-1]
        (stop_pixel_index,) = np.where(self.pixels == stop_pixel)[0]

        if self.start_pix != 0 and self.stop_pix != 0:
            start_pixel_index = np.where(self.pixels == self.start_pix)[0]
            start_pixel_index = start_pixel_index[0]
            stop_pixel_index = np.where(self.pixels == self.stop_pix)[0]
            stop_pixel_index = stop_pixel_index[0]

        for pixel in self.pixels[start_pixel_index : stop_pixel_index + 1]:
            start_time = time.time()
            if self.verbose:
                print(f"Pixel: {pixel}")

            fnorm = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "FNORM"
            ].values[0]
            alpha1 = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA1"
            ].values[0]
            alpha2 = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA2"
            ].values[0]
            nu_break = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "NU_BREAK"
            ].values[0]
            T_X = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "TX"
            ].values[0]
            T_E = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "TE"
            ].values[0]
            nu_t = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "NU_T"
            ].values[0]

            x0 = [fnorm, alpha1, alpha2, nu_break, T_X, T_E, nu_t]
            print(f"initial params: {x0}")

            self.b_temp = np.array([])
            for f in self.frequencies_string:
                self.b_temp = np.append(
                    self.b_temp, self.df_convex.loc[pixel, f"{f}MHz"]
                )
            print(self.b_temp)
            bounds = (
                [-np.inf, np.inf],
                [2, 3],
                [2, 3],
                [-np.inf, np.inf],
                [-np.inf, np.inf],
                [0, 10000],
                [-np.inf, np.inf],
            )
            result = minimize(
                self.chisquared,
                args=(self.frequencies, self.b_temp),
                x0=x0,
                method="Nelder-Mead",
                bounds=bounds,
                options={"maxiter": 10000},
            )

            print(
                f"optimized params:- Pixel: {pixel}, Fnorm: {result.x[0]}, Alpha1: {result.x[1]}, Alpha2: {result.x[2]}, Nu_break: {result.x[3]}, Tx: {result.x[4]}, Te: {result.x[5]}, Nu_T: {result.x[6]}, Chisquared: {result.fun}\n"
            )
            logging.info(
                f"{pixel},{result.x[0]},{result.x[1]},{result.x[2]},{result.x[3]},{result.x[4]},{result.x[5]},{result.x[6]},{result.fun}"
            )
            stop_time = time.time()
            print(f"Iteration time for the pixel: {(start_time - stop_time)/60} min")


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.environ.get("DATA")

    path_to_brightness_temperature_file = f"{DATA}brightness_temp_per_pixel.csv"
    path_to_convexity_file = f"{DATA}convexity.csv"
    path_to_initial_parameters_file = f"{DATA}convex_model_initial_parameters.csv"
    path_to_convex_fits_file = f"{DATA}convex_pixel_fits.log"
    convex_model = ConvexModel(
        path_to_brightness_temperature_file,
        path_to_convexity_file,
        path_to_initial_parameters_file,
        path_to_convex_fits_file,
        verbose=True,
        start_pix=0,
        stop_pix=0,
    )
    convex_model.fit()


