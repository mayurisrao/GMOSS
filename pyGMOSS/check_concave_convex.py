import pandas as pd
import numpy as np


def check_concave_convex(path_to_brightness_temperature_file: str) -> pd.DataFrame:
    """
    Used the output from program "read_data_at_pixel.py" and checks the cocavity
    or convexity of frequency/brightness temperature values at each pixel.

    Parameters
    -------------
    path_to_brightness_temperature_file: str. File Name of data file containing brightness temperature per pixel.

    Returns
    -------------
    df: pandas.DataFrame of pixel number, alpha 1, alpha 2 and convexity/concavity flag.
    """

    df = pd.read_csv(path_to_brightness_temperature_file)
    n_pixels = len(df)

    convexity_list = []
    for pixel in range(n_pixels):
        convexity = {}
        convexity["PIXEL"] = pixel + 1
        print(f"Pixel Number: {pixel}")

        freq_45_in_GHz = 45 * 1e-3
        b_temp_45_list = df.loc[:, "45MHz"]
        b_temp_45 = b_temp_45_list[pixel]

        freq_150_in_GHz = 150 * 1e-3
        b_temp_150_list = df.loc[:, "150MHz"]
        b_temp_150 = b_temp_150_list[pixel]

        alpha_1 = (np.log10(b_temp_45) - np.log10(b_temp_150)) / (
            np.log10(freq_150_in_GHz) - np.log10(freq_45_in_GHz)
        )
        if alpha_1 < 2.0:
            alpha_1 = 2.0
        if alpha_1 > 3.0:
            alpha_1 = 3.0
        convexity["ALPHA_1"] = alpha_1

        freq_408_in_GHz = 408 * 1e-3
        b_temp_408_list = df.loc[:, "408MHz"]
        b_temp_408 = b_temp_408_list[pixel]

        freq_1420_in_GHz = 1420 * 1e-3
        b_temp_1420_list = df.loc[:, "1420MHz"]
        b_temp_1420 = b_temp_1420_list[pixel]

        alpha_2 = (np.log10(b_temp_408) - np.log10(b_temp_1420)) / (
            np.log10(freq_1420_in_GHz) - np.log10(freq_408_in_GHz)
        )
        if alpha_2 < 2.0:
            alpha_2 = 2.0
        if alpha_2 > 3.0:
            alpha_2 = 3.0
        convexity["ALPHA_2"] = alpha_2

        if alpha_1 < alpha_2:
            convexity["Concave/Convex"] = "Convex"
        else:
            convexity["Concave/Convex"] = "Concave"

        convexity_list.append(convexity)

    df = pd.DataFrame(convexity_list)

    return df


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()
    DATA = os.getenv("DATA")
    path_to_brightness_temperature_file = "brightness_temp_at_pixel.csv"

    df = check_concave_convex(path_to_brightness_temperature_file)
    df.to_csv(f"{DATA}convexity.csv", index=False)
