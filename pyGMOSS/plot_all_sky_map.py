import healpy as hp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_mollview_all_sky_map(DATA, frequency):
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
    concave_pixel_fits_df = pd.read_csv(f"{DATA}concave_pixel_fits.csv")

    PIXELS = concave_pixel_fits_df.loc[:, "PIXEL"].values
    concave_pixel_fits_df.set_index("PIXEL", inplace=True, drop=True)

    b_temp = np.array([])
    for i in range(1, 3073):
        if np.any(PIXELS == i):
            b_temp = np.append(
                b_temp,
                concave_func(
                    frequency,
                    concave_pixel_fits_df.loc[i, "FNORM1"],
                    concave_pixel_fits_df.loc[i, "FNORM2"],
                    concave_pixel_fits_df.loc[i, "ALPHA_1"],
                    concave_pixel_fits_df.loc[i, "ALPHA_2"],
                    concave_pixel_fits_df.loc[i, "T_X"],
                    concave_pixel_fits_df.loc[i, "NU_T"],
                    concave_pixel_fits_df.loc[i, "T_E"],
                ),
            )
        else:
            b_temp = np.append(b_temp, 0)

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


def concave_func(x, C_1, C_2, alpha_1, alpha_2, I_x, nu_t, T_e):
    one = np.power(x, -alpha_1)
    two = (C_2 / C_1) * np.power(x, -alpha_2)
    three = I_x * np.power(x, -2.1)
    expo = np.exp(-1 * np.power((nu_t / x), 2.1))
    eqn_one = C_1 * (one + two + three) * expo
    eqn_two = T_e * (1 - expo)
    return eqn_one + eqn_two


if __name__ == "__main__":
    from dotenv import load_dotenv
    import os

    load_dotenv()

    DATA = os.environ.get("DATA")

    frequency = 1200 * 1e-3  # GHz

    plot_mollview_all_sky_map(DATA, frequency)
