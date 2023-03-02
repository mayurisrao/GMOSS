import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class PlotBrightnessTemperatureModel:
    def __init__(
        self,
        path_to_brightness_temperature_file,
        path_to_concave_fit_parameters_file,
        path_to_convex_fit_parameters_file,
    ):
        self.path_to_brightness_temperature_file = path_to_brightness_temperature_file
        self.path_to_concave_fit_parameters_file = path_to_concave_fit_parameters_file
        # self.path_to_convex_fit_parameters_file = path_to_convex_fit_parameters_file

        self.frequencies = np.array([22, 45, 150, 408, 1420, 23000]) * 1e-3
        self.frequencies_string = np.array(["22", "45", "150", "408", "1420", "23000"])
        self.frequencies_for_plot = np.arange(22, 25e3, 10) * 1e-3

        self.brightness_temperature_df = pd.read_csv(
            path_to_brightness_temperature_file
        )
        self.concave_fit_parameters_df = pd.read_csv(
            path_to_concave_fit_parameters_file
        )
        # self.convex_fit_parameters_df = pd.read_csv(path_to_convex_fit_parameters_file)

    def concave_func(self, x, C_1, C_2, alpha_1, alpha_2, I_x, nu_t, T_e):
        one = np.power(x, -alpha_1)
        two = (C_2 / C_1) * np.power(x, -alpha_2)
        three = I_x * np.power(x, -2.1)
        expo = np.exp(-1 * np.power((nu_t / x), 2.1))
        eqn_one = C_1 * (one + two + three) * expo
        eqn_two = T_e * (1 - expo)
        return eqn_one + eqn_two

    def plot_concave_pixel_fit(self, pixel):
        FNORM1 = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "FNORM1"
        ]
        FNORM2 = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "FNORM2"
        ]
        ALPHA_1 = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA_1"
        ]
        ALPHA_2 = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "ALPHA_2"
        ]
        T_X = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "T_X"
        ]
        NU_T = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "NU_T"
        ]
        T_E = self.concave_fit_parameters_df.loc[
            self.concave_fit_parameters_df.loc[:, "PIXEL"] == pixel, "T_E"
        ]

        brightness_temperature_data_points = np.zeros(len(self.frequencies))
        for i, f in enumerate(self.frequencies_string):
            brightness_temperature_data_points[i] = self.brightness_temperature_df.loc[
                self.brightness_temperature_df.loc[:, "PIXEL"] == pixel, f"{f}MHz"
            ]

        brightness_temperature_model_points = np.zeros(len(self.frequencies_for_plot))
        for i, f in enumerate(self.frequencies_for_plot):
            brightness_temperature_model_points[i] = self.concave_func(
                f, FNORM1, FNORM2, ALPHA_1, ALPHA_2, T_X, NU_T, T_E
            )

        fig, ax = plt.subplots()
        ax.scatter(self.frequencies, brightness_temperature_data_points, c="red")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.plot(self.frequencies_for_plot, brightness_temperature_model_points)
        ax.set_title(f"Model for Concave Pixel {pixel}")
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("Brightness Temperature (K)")
        ax.grid(True)

        plt.show()

    def plot_convex_pixel_fit(self, pixel):
        # TO DD: Implement this function
        pass


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.environ.get("DATA")

    pixel = 316
    path_to_brightness_temperature_file = os.path.join(
        DATA, "brightness_temp_per_pixel.csv"
    )
    path_to_concave_fit_parameters_file = os.path.join(DATA, "concave_pixel_fits.csv")
    path_to_convex_fit_parameters_file = os.path.join(DATA, "convex_pixel_fits.csv")
    plotter = PlotBrightnessTemperatureModel(
        path_to_brightness_temperature_file,
        path_to_concave_fit_parameters_file,
        path_to_convex_fit_parameters_file,
    )

    plotter.plot_concave_pixel_fit(pixel)
