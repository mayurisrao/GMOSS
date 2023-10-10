import numpy as np
import pandas as pd
from scipy.special import gamma
import scipy.integrate as integrate
import scipy as sp
from scipy.special import kve
from scipy.integrate import quad


class ConvexModelInitialParameters:
    """
    Class Name:
    -------------
    ConvexModelInitialParameters

    Description:
    -------------
    Takes the brightness temparature and convex pixels files and outputs initial parameters for the convex model.

    Parameters:
    ---------------
    - path_to_brightness_temparature_file : str
        Path to folder where the brightness temparature (CSV) file is located with column names [PIXEL,1420MHz,150MHz,22MHz,23000MHz,408MHz,45MHz]
    - path_to_convexity_file : str
        Path to the folder where the convexity file (CSV) is located with column names [PIXEL,ALPHA_1,ALPHA_2,Concave/Convex]
    - speed_of_light: float, optional
    - mass_of_electron: float, optional
    - charge_of_electron: float, optional
    - sine_alpha: float, optional
    - magnetic_field: float, optional
    - GSPAN: float, optional
    - verbose: bool, optional

    Returns:
    ---------------
    df: pandas.DataFrame with initial parameters [PIXEL,FNORM,ALPHA1,ALPHA2,NU_BREAK,TX,TE,NU_T] for the convex model.
    """

    def __init__(
        self,
        path_to_brightness_temparature_file: str,
        path_to_convexity_file: str,
        speed_of_light: float = 2.99792458e08,
        mass_of_electron: float = 9.1e-31,
        charge_of_electron: float = 1.6e-19,
        sine_alpha: float = 1.0,
        magnetic_field: float = 1e-9,
        GSPAN: float = 100.0,
        verbose: bool = False,
    ) -> pd.DataFrame:
        self.path_to_brightness_temparature_file = path_to_brightness_temparature_file
        self.path_to_convexity_file = path_to_convexity_file
        self.speed_of_light = speed_of_light
        self.mass_of_electron = mass_of_electron
        self.charge_of_electron = charge_of_electron
        self.sine_alpha = sine_alpha
        self.magnetic_field = magnetic_field
        self.GSPAN = GSPAN
        self.verbose = verbose

        self.scale_gam_nu = (
            3.0 * self.charge_of_electron * self.magnetic_field * self.sine_alpha
        ) / (4.0 * np.pi * self.mass_of_electron * self.speed_of_light)

    def modbessik2(self, u):  # can take in an array or an integer
        xnu = 5.0 / 3.0
        ORDER = xnu
        N = 1
        ARG = np.divide(1, u)
        BK = kve(ORDER, ARG)
        xrk = BK / np.exp(ARG)
        return np.divide(xrk, (u * u))

    def integrand_for_param(self, gama, alpha, nu):
        nu_c = (self.scale_gam_nu * (gama**2)) / 1e9
        x = nu / nu_c
        integrand_ = self.F(x) * x * np.power(gama, -1 * (2 * alpha - 3))
        return integrand_

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

    def convex_model_initial_parameter_generator(self):
        """
        Function Name:
        -------------
        "convex_model_initial_parameter_generator()"

        Description:
        -------------
        This method generates the initial parameters for the convex model.

        Parameters:
        -------------
        This method takes no parameters.

        Returns:
        -------------
        - This method returns a DataFrame with the initial(to be optimized) parameters.
          The DataFrame has the following columns:
        - "PIXEL": The pixel number for which the parameters will be fitted.
        - "FNORM": initial scaling/normalization parameter.
        - "ALPHA1": initial low frequency spectral index parameter.
        - "ALPHA2": initial high frequency spectral index parameter.
        - "NU_BREAK": initial break frequency parameter.
        - "TX": initial parameter representing optically thin free free emission index.
        - "TE": electron temparature initial parameter.
        - "NU_T": frequecny of thermal absorption turnover initial paramter.

        Example:
        -------------
        ```python
        convex_model_initial_parameters = ConvexModelInitialParameters("brightness_temp.csv", "convexity.csv", verbose=True)
        result_df = convex_model_initial_parameters.convex_model_initial_parameter_generator()
        df.to_csv("convex_model_initial_parameters.csv")
        ```
        In this example, the "ConvexModelInitialParameters" class is instantiated with the "brightness_temperature.csv" and "convexity.csv". Then, the 'convex_model_initial_parameter_generator'
        method is called to generate the initial parameters. The resulting DataFrame is assigned to the 'result_df' variable which is later converted to a csv file using
        'to_csv' function.
        """
        df = pd.read_csv(self.path_to_brightness_temparature_file)
        convexity_df = pd.read_csv(self.path_to_convexity_file)

        convexity_df = convexity_df.loc[
            convexity_df.loc[:, "Concave/Convex"] == "Convex", :
        ]
        df = pd.merge(df, convexity_df, on="PIXEL", how="right")
        pixels = df.loc[:, "PIXEL"].values
        alpha_1_list = df.loc[:, "ALPHA_1"].values
        alpha_2_list = df.loc[:, "ALPHA_2"].values

        # check the units for frequency for calculation of initial_parameters
        frequency = np.array([22.0, 45.0, 150.0, 408.0, 1420.0, 23000.0])  # x_values
        frequency *= 1e-3
        frequency_string = np.array(["22", "45", "150", "408", "1420", "23000"])

        b_temp = np.zeros(len(frequency_string))
        to_save_list = []
        for idx, pixel in enumerate(pixels):
            if self.verbose:
                print(f"pixel number: {pixel}")

            to_save = {}

            to_save["PIXEL"] = pixel
            for i, f in enumerate(frequency_string):
                b_temp[i] = df.loc[df.loc[:, "PIXEL"] == pixel, f"{f}MHz"]

            Te = 8000.0
            nu_t = 0.001
            nu_break = np.sqrt(0.150 * 0.408) * 1e9
            extn = np.exp(-1.0 * np.power((nu_t / frequency[4]), 2.1))
            alpha1, alpha2 = alpha_1_list[idx], alpha_2_list[idx]

            nu = 1.420
            nu_min = nu * 1e9 / self.GSPAN
            nu_max = nu * 1e9 * self.GSPAN
            gama_min = np.sqrt((nu_min) / self.scale_gam_nu)
            gama_max = np.sqrt((nu_max) / self.scale_gam_nu)
            gama_break = np.sqrt((nu_break) / self.scale_gam_nu)

            xb = gama_break
            xl = gama_min
            xu = gama_max

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

            fnorm = (b_temp[4] - (Te * (1.0 - extn))) / (
                (np.power(nu, -2.0) * float(rint)) * extn
            )

            nu = 22.690
            nu_min = nu * 1e9 / self.GSPAN
            nu_max = nu * 1e9 * self.GSPAN
            gama_min = np.sqrt((nu_min) / self.scale_gam_nu)
            gama_max = np.sqrt((nu_max) / self.scale_gam_nu)
            gama_break = np.sqrt((nu_break) / self.scale_gam_nu)

            xb = gama_break
            xl = gama_min
            xu = gama_max

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

            extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
            temp1 = fnorm * extn

            Tx = (
                ((b_temp[5] - Te * (1.0 - extn)) / temp1)
                - (np.power(nu, -2.0) * float(rint))
            ) / np.power(frequency[4], -2.1)
            if Tx <= 0:
                Tx = 1.0e-10

            to_save["FNORM"] = fnorm
            to_save["ALPHA1"] = alpha1
            to_save["ALPHA2"] = alpha2
            to_save["NU_BREAK"] = nu_break
            to_save["TX"] = Tx
            to_save["TE"] = Te
            to_save["NU_T"] = nu_t
            to_save_list.append(to_save)
        return pd.DataFrame(to_save_list)


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.environ.get("DATA")
    convex_model_initial_parameters = ConvexModelInitialParameters(
        f"{DATA}brightness_temp_per_pixel.csv", f"{DATA}convexity.csv", verbose=True
    )
    df = convex_model_initial_parameters.convex_model_initial_parameter_generator()
    df.to_csv(f"{DATA}convex_model_initial_parameters.csv")
