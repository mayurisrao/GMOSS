import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.integrate as integrate
import scipy as sp


class ConvexModel:
    """
    This class represents the convex spectral model.

    Parameters
    ----------
    - path_to_brightness_temperature_file : str
        The path to the brightness temperature file.
    - path_to_convexity_file : str
        The path to the convexity file.
    - path_to_initial_parameters_file : str
        The path to the initial parameters file.
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
        GSPAN=100,
        NHPIX=3072,
        speed_of_light=2.99792458e08,
        mass_of_electron=9.10938356e-31,
        charge_of_electron=1.6e-19,
        sine_alpha=1.0,
        magnetic_field=1e-9,
        verbose=False,
    ):
        self.brightness_temperature_df = pd.read_csv(
            path_to_brightness_temperature_file
        )
        self.convexity_df = pd.read_csv(path_to_convexity_file)
        self.initial_parameters_df = pd.read_csv(path_to_initial_parameters_file)
        self.GSPAN = GSPAN
        self.NHPIX = NHPIX
        self.speed_of_light = speed_of_light
        self.mass_of_electron = mass_of_electron
        self.charge_of_electron = charge_of_electron
        self.sine_alpha = sine_alpha
        self.magnetic_field = magnetic_field
        self.verbose = verbose

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
        self.pixels = self.df_convex.loc[:, "PIXEL"]
        self.df_convex.set_index("PIXEL", inplace=True)

        self.scale_gam_nu = (
            3.0 * self.charge_of_electron * self.magnetic_field * self.sine_alpha
        ) / (4.0 * np.pi * self.mass_of_electron * self.speed_of_light)

    def F(self, x):
        if x < 3:
            one = (np.pi * x) / np.sqrt(3)
            two = (9 * (x ** (11 / 3)) * gamma(-2 / 3)) / (160 * 2 ** (2 / 3))
            three = ((x ** (1 / 3)) * (16 + (3 * x**2)) * gamma(-1 / 3)) / (
                24 * 2 ** (1 / 3)
            )
            return -one + two - three
        else:
            exponential_term = np.exp(-x) / (967458816 * np.sqrt(2) * x ** (5 / 2))
            const = 13 * np.sqrt(np.pi)
            quad_term = 2429625 + 2 * x * (
                -1922325 + (5418382 * x) + 83221732 * (x**2)
            )
            error_function_term = (
                1196306216
                * np.exp(x)
                * np.pi
                * (x ** (7 / 2))
                * sp.special.erfc(np.sqrt(x))
            )
            return exponential_term * ((const * quad_term) - error_function_term)

    def integrand_for_convex(self, gama, alpha, nu):
        nu_c = (self.scale_gam_nu * (gama**2)) / 1e9
        x = nu / nu_c
        integrand_ = self.F(x) * x * np.power(gama, -1 * (2 * alpha - 3))
        return integrand_

    def func(self, parameters) -> float:  # pp is an array
        # TO DO: To extract initial parameters from file
        fnorm, alpha1, alpha2, nu_break, T_X, T_E, nu_t = parameters

        if alpha1 < 2.0 or alpha1 > 3.0:
            alpha1 = 1000000000
        if alpha2 < 2.0 or alpha2 > 3.0:
            alpha2 = 1000000000
        if T_E < 0.0 or T_E > 10000:
            T_E = 1000000000

        chisq = 0.0
        for i in range(len(self.frequencies)):
            nu = self.frequencies[i]
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
                I, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha2, nu)
                )
                I *= np.power(gama_break, 2 * C1 - 3)

            elif xu < xb:
                C1 = alpha1
                I, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha1, nu)
                )
                I *= np.power(gama_break, 2 * C1 - 3)
            else:
                xu = xb
                C1 = alpha1
                I1, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha1, nu)
                )
                I1 *= np.power(gama_break, 2 * C1 - 3)
                xl = xb
                xu = gama_max
                C1 = alpha2
                I2, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha2, nu)
                )
                I2 *= np.power(gama_break, 2 * C1 - 3)
                I = I1 + I2
            extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
            FFIT = fnorm * (
                (np.power(nu, -2) * I) + T_X * np.power(nu, -2.1)
            ) * extn + T_E * (1.0 - extn)
            DDIF = (self.b_temp[i] - FFIT) / (self.b_temp[i])
            if i <= 5:
                chisq += DDIF * DDIF
        chisq /= 6.0

        if self.verbose:
            print(f"-- Chi Square: {chisq}")
        return chisq

    def fit(self):
        to_save_list = []
        for pixel in self.pixels[35:36]:
            if self.verbose:
                print(f"Pixel: {pixel}")
            to_save = {}
            to_save["PIXEL"] = pixel

            # PIXEL,FNORM,ALPHA1,ALPHA2,NU_BREAK,TX,TE,NU_T
            fnorm = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "FNORM"
            ]
            alpha1 = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA1"
            ]
            alpha2 = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA2"
            ]
            nu_break = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "NU_BREAK"
            ]
            T_X = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "TX"
            ]
            T_E = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "TE"
            ]
            nu_t = self.initial_parameters_df.loc[
                self.initial_parameters_df.loc[:, "PIXEL"] == pixel, "NU_T"
            ]

            x0 = [fnorm, alpha1, alpha2, nu_break, T_X, T_E, nu_t]

            self.b_temp = np.array([])
            for f in self.frequencies_string:
                self.b_temp = np.append(
                    self.b_temp, self.df_convex.loc[pixel, f"{f}MHz"]
                )

            bounds = (
                [-np.inf, np.inf],
                [2, 3],
                [2, 3],
                [0.001, np.inf],
                [0.001, np.inf],
                [0.001, 10000],
                [0.001, np.inf],
            )
            result = minimize(
                self.func,
                x0=x0,
                method="Nelder-Mead",
                bounds=bounds,
                options={"verbose": 1, "maxiter": 10000},
            )

            to_save["FNORM"] = result.x[0]
            to_save["ALPHA1"] = result.x[1]
            to_save["ALPHA2"] = result.x[2]
            to_save["NU_BREAK"] = result.x[3]
            to_save["TX"] = result.x[4]
            to_save["TE"] = result.x[5]
            to_save["NU_T"] = result.x[6]

            to_save_list.append(to_save)

        return pd.DataFrame(to_save_list)

    def convex_func(self, parameters, frequencies=None):
        fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = parameters
        b_temps = []

        if frequencies == None:
            frequencies = self.frequencies

        for i in range(len(frequencies)):
            nu = self.frequencies[i]
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
                I, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha2, nu)
                )
                I *= np.power(gama_break, 2 * C1 - 3)

            elif xu < xb:
                C1 = alpha1
                I, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha1, nu)
                )
                I *= np.power(gama_break, 2 * C1 - 3)
            else:
                xu = xb
                C1 = alpha1
                I1, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha1, nu)
                )
                I1 *= np.power(gama_break, 2 * C1 - 3)
                xl = xb
                xu = gama_max
                C1 = alpha2
                I2, _ = integrate.quad(
                    self.integrand_for_convex, xl, xu, args=(alpha2, nu)
                )
                I2 *= np.power(gama_break, 2 * C1 - 3)
                I = I1 + I2
            extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
            FFIT = fnorm * (
                (np.power(nu, -2) * I) + Tx * np.power(nu, -2.1)
            ) * extn + Te * (1.0 - extn)
            b_temps.append(FFIT)
        return b_temps


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.environ.get("DATA")

    path_to_brightness_temperature_file = f"{DATA}brightness_temp_per_pixel.csv"
    path_to_convexity_file = f"{DATA}convexity.csv"
    path_to_initial_parameters_file = f"{DATA}convex_model_initial_parameters.csv"
    convex_model = ConvexModel(
        path_to_brightness_temperature_file,
        path_to_convexity_file,
        path_to_initial_parameters_file,
        verbose=True,
    )
    df = convex_model.fit()

    df.to_csv(f"{DATA}convex_pixel_fits.csv", index=False)
