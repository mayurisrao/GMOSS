import numpy as np
import pandas as pd
import glob
import os


def extract_brightness_temp_at_pixel(
    path_to_data: str,
    n_pixels: int,
    correction_offset: float,
    correction_scale: float,
    correction_at_freq: int,
    tcmb: float,
) -> pd.DataFrame:
    """
    Descrption:
    ---------------
    The function extract_brightness_temp_at_pixel takes the path path to folder containing the data files(all sky maps), number of pixels, offset and scaling corrections,
    frequency correction, and a value for tcmb. The function reads data from text files, applies corrections, and extracts brightness temperature data for 6
    different frequencies at each pixel. The extracted data is returned in the form of a pandas DataFrame.

    Parameters:
    ---------------
    - path_to_data : A string representing the path to the data directory containing all sky maps text files.
    - n_pixels : An integer representing the number of pixels.
    - correction_offset : A float representing the offset to be applied to brightness temperature for frequency correction_at_freq.
    - correction_scale : Scaling to be applied to brightness temperature for frequency correction_at_freq.
    - correction_at_freq : An integer representing the frequency at which corrections to be applied to brightness temperature.
    - tcmb : A float representing the value of tcmb (cosmic microwave background temparature). Recommended value to plugin 2.72548.

    Returns:
    ---------------
    - df : pandas.DataFrame with brightness temperature values.
      i.e The function returns a pandas DataFrame object containing brightness temperature data for different frequencies at each pixel.
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
            if frequencies[i] == 22:
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
    df = df.reindex(columns=["PIXEL", "22MHz", "45MHz", "150MHz", "408MHz", "1420MHz", "23000MHz"])

    return df


if __name__ == "__main__":
    from dotenv import load_dotenv
    PIXELS = 3072
    CORRECTION_150_OFFSET = 21.4
    CORRECTION_150_SCALING = 1.05
    TCMB = 2.72548

    load_dotenv()
    DATA = os.environ.get("DATA")

    df = extract_brightness_temp_at_pixel(
        DATA, PIXELS, CORRECTION_150_OFFSET, CORRECTION_150_SCALING, 150, TCMB
    )
    df.to_csv(f"{DATA}brightness_temp_per_pixel.csv", index=False)
