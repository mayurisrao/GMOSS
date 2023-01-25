import numpy as np
import pandas as pd
import glob
from dotenv import load_dotenv
import os

def extract_brightness_temp_at_pixel(path_to_data: str, n_pixels: int, correction_offset: float, correction_scale: float, correction_at_freq: int, tcmb: float) -> pd.DataFrame:
    """
    Takes all sky maps as input and returns a pandas DataFrame with brightness temperature at each pixel for 6 different frequencies

    Parameters:
    ---------------
    path_to_data: Path to folder where the all sky maps are stored.
    n_pixels: Number of pixels in each map
    correction_offset: Offset to be applied to brightness temperature for frequency correction_at_freq.
    correction_scale: Scaling to be applied to brightness temperature for frequency correction_at_freq.
    correction_at_freq: Frequency at which corrections to be applied to brightness temperature.

    Returns:
    ---------------
    df: pandas.DataFrame with brightness temperature values.
    """

    text_files = glob.glob(path_to_data + "*.txt")
    map_files = []
    for file in text_files:
        file_without_prepend = file.split("/")[-1]
        if file_without_prepend.startswith("map"):
            map_files.append(file)

    no_of_map_files = len(map_files)

    frequencies = np.array([])
    for file in map_files:
        file = file.split("/")[-1]
        frequencies = np.append(frequencies, float(file.split("_")[1]))

    brightness_temperature_list = []
    for pixel in range(n_pixels):
        print(f"Pixel Number: {pixel+1}")
        brightness_temperature = {}
        brightness_temperature["PIXEL"] = pixel + 1

        for i in range(no_of_map_files):
            data = np.genfromtxt(map_files[i])
            if frequencies[i] ==  22:
                if correction_at_freq == 22:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                    data[pixel] = data[pixel] - tcmb

                brightness_temperature["22MHz"] = data[pixel]

            if frequencies[i] == 45:
                if correction_at_freq == 45:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                    data[pixel] = data[pixel] - tcmb

                brightness_temperature["45MHz"] = data[pixel]

            if frequencies[i] == 150:
                if correction_at_freq == 150:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                    data[pixel] = data[pixel] - tcmb
                    
                brightness_temperature["150MHz"] = data[pixel]

            if frequencies[i] == 408:
                if correction_at_freq == 408:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                
                data[pixel] = data[pixel] - tcmb

                brightness_temperature["408MHz"] = data[pixel]

            if frequencies[i] == 1420:
                if correction_at_freq == 1420:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                
                data[pixel] = data[pixel] - tcmb

                brightness_temperature["1420MHz"] = data[pixel]

            if frequencies[i] == 23000:
                if correction_at_freq == 23000:
                    data[pixel] = (data[pixel] - correction_offset) * correction_scale
                    data[pixel] = data[pixel] - tcmb

                brightness_temperature["23000MHz"] = data[pixel]
        
        brightness_temperature_list.append(brightness_temperature)

    df = pd.DataFrame(brightness_temperature_list)

    return df


if __name__ == "__main__":
    PIXELS = 3072
    CORRECTION_150_OFFSET = 21.4
    CORRECTION_150_SCALING = 1.05
    TCMB = 2.72548

    load_dotenv()
    DATA = os.environ.get("DATA")

    df = extract_brightness_temp_at_pixel(DATA, PIXELS, CORRECTION_150_OFFSET, CORRECTION_150_SCALING, 150, TCMB)
    df.to_csv(f"{DATA}brightness_temp_per_pixel.csv", index = False)