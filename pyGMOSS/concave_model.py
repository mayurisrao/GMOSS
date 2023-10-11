import numpy as np
import pandas as pd
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt


class ConcaveModel:
    """
    Class Name:
    -------------
    ConcaveModel

    Description:
    -------------
    This class represents the concave spectal model

    Parameters:
    -------------
        - path_to_brightness_temperature_file (str): 
            The path to the file containing brightness temperature data.
        - path_to_convexity_file (str): 
            The path to the file containing convexity data.
        - T_e (float): 
            The electron temperature initial param to use in the model (default 8000).
        - nu_t (float): 
            The frequency of thermal absorbtion turnover initial param (default 0.001).
    """

    def __init__(
        self,
        path_to_brightness_temperature_file: str,
        path_to_convexity_file: str,
        T_e: float = 8000,
        nu_t: float = 0.001,
    ):
        self.T_e = T_e
        self.nu_t = nu_t
        self.brightness_temperature_df = pd.read_csv(
            path_to_brightness_temperature_file
        )
        self.convexity_df = pd.read_csv(path_to_convexity_file)
        self.frequencies_string = np.array(
            ["22", "45", "150", "408", "1420", "23000"]
        )  # MHz
        self.frequencies = np.array([22, 45, 150, 408, 1420, 23000])  # MHz
        self.frequencies = self.frequencies * 1e-3  # GHz
        self.df = pd.merge(
            self.brightness_temperature_df, self.convexity_df, on="PIXEL"
        )
        self.df_concave = self.df.loc[self.df.loc[:, "Concave/Convex"] == "Concave", :]
        self.df_concave.reset_index(inplace=True, drop=True)
        self.pixels = self.df_concave.loc[:, "PIXEL"]

    def concave_func(self, nu, fnorm_1, fnorm_2, alpha_1, alpha_2, T_x, T_e, nu_t):
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

    def chisquared(self, params, xobs, yobs):
        chisq = 0
        yobs_length = len(yobs)
        for i in range(yobs_length):
            chisq = (
                chisq + ((yobs[i] - self.concave_func(xobs[i], *params)) / yobs[i]) ** 2
            )
        return chisq / 6

    def fit(self):
        """
        Function Name:
        -------------
        `fit()`

        Description:
        -------------
        This method fits a concave model to the brightness temperature data for each pixel in the dataset and returns the resulting fitted parameters.

        Parameters:
        -------------
        This method takes no parameters.

        Returns:
        -------------
        - This method returns a Pandas DataFrame object containing the fitted parameters/optimized for each pixel in the dataset. 
          The DataFrame has the following columns:
        - "PIXEL": The pixel number for which the parameters were fitted.
        - "FNORM1": The first scaling constant.
        - "FNORM2": The second scaling constant.
        - "ALPHA_1": The low frequency spectral index.
        - "ALPHA_2": The high frequency spectral index.
        - "T_X": parameter representing optically thin free free emission index.
        - "T_E": The electron temperature.
        - "NU_T": The frequecny of thermal absorption turnover.
        

        Example:
        -------------
        ```python
        model = ConcaveModel("brightness_temperature.csv", "convexity.csv")
        result_df = model.fit()
        ```
        In this example, the "ConcaveModel" is instantiated with the "brightness_temperature.csv" and "convexity.csv" files. Then, the `fit()` method is called to fit the concave model to the data and return a DataFrame containing the fitted parameters. 
        The resulting DataFrame is assigned to the `result_df` variable.
        """
        to_save_list = []
        for pixel in self.pixels:
            to_save = {}
            print(f"Pixel Number: {pixel}")

            self.b_temp = np.array([])
            for f in self.frequencies_string:
                self.b_temp = np.append(
                    self.b_temp,
                    self.df_concave.loc[
                        self.df_concave.loc[:, "PIXEL"] == pixel, f"{f}MHz"
                    ],
                )

            alpha_1 = self.df_concave.loc[
                self.df_concave.loc[:, "PIXEL"] == pixel, "ALPHA_1"
            ].values[0]
            alpha_2 = self.df_concave.loc[
                self.df_concave.loc[:, "PIXEL"] == pixel, "ALPHA_2"
            ].values[0]

            fnorm_1 = self.b_temp[2] / np.power(self.frequencies[2], -1 * alpha_1)
            fnorm_2 = (
                self.b_temp[3] / np.power(self.frequencies[3], -1 * alpha_2)
            ) / fnorm_1
            expo = np.exp(-1 * np.power(self.nu_t / self.frequencies[5], 2.1))
            T_x = (1.0 / np.power(self.frequencies[5], -2.1)) * (
                (self.b_temp[5] - self.T_e * (1.0 - expo)) / (fnorm_1 * expo)
            ) - (
                (np.power(self.frequencies[5], -1.0 * alpha_1))
                + (fnorm_2 * np.power(self.frequencies[5], -1.0 * alpha_2))
            )

            if T_x <= 0:
                T_x = 1e-10

            x0 = [fnorm_1, fnorm_2, alpha_1, alpha_2, T_x, self.T_e, self.nu_t]
            bounds = (
                [-np.inf, np.inf],
                [-np.inf, np.inf],
                [2, 3],
                [2, 3],
                [-np.inf, np.inf],
                [0, 10000],
                [-np.inf, np.inf],
                
            )

            self.result = minimize(
                self.chisquared,
                args=(self.frequencies, self.b_temp),
                x0=x0,
                method="Nelder-Mead",
                bounds=bounds,
            )

            to_save["PIXEL"] = pixel
            to_save["FNORM1"] = self.result.x[0]
            to_save["FNORM2"] = self.result.x[1]
            to_save["ALPHA_1"] = self.result.x[2]
            to_save["ALPHA_2"] = self.result.x[3]
            to_save["T_X"] = self.result.x[4]
            to_save["T_E"] = self.result.x[5]
            to_save["NU_T"] = self.result.x[6]
            
            to_save["CHISQUARED"] = self.chisquared(
                self.result.x, self.frequencies, self.b_temp
            )

            to_save_list.append(to_save)

        return pd.DataFrame(to_save_list)


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os
    load_dotenv()
    DATA = os.environ.get("DATA")

    model = ConcaveModel(f"{DATA}brightness_temp_per_pixel.csv", "{DATA}convexity.csv")
    df = model.fit()
    df.to_csv("{DATA}concave_pixel_fits_testing.csv", index=False)
