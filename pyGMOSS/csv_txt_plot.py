import numpy as np
import pandas as pd
from dotenv import load_dotenv
import os

load_dotenv()
DATA = os.environ.get("DATA")

data = np.genfromtxt(DATA + "total_model_params_11apr16_allpix.txt")

df = pd.DataFrame(
    data,
    columns=[
        "PIXEL",
        "COVEXITY",
        "COL3",
        "COL4",
        "COL5",
        "COL6",
        "COL7",
        "COL8",
        "COL9",
        "chisq",
    ],
)

convex_df = df[df.CONVEXITY == 0]

convex_df.loc[:, "COF1"] = convex_df.loc[:, "COF1"].apply(lambda x: np.power(10, x))
convex_df.rename(columns={"COF1": "FNORM"}, inplace=True)
convex_df["COF2"] = df["COF2"].apply(lambda x: np.power(10, x))
convex_df.rename(columns={"COF2": "ALPHA1"}, inplace=True)
convex_df["COF3"] = df["COF3"].apply(lambda x: np.power(10, x)) + df["COF2"].apply(
    lambda x: np.power(10, x)
)
convex_df.rename(columns={"COF3": "ALPHA2"}, inplace=True)
convex_df["COF4"] = df["COF4"].apply(lambda x: np.power(10, x) / 1e6)
convex_df.rename(columns={"COF4": "NU_BREAK_MHz"}, inplace=True)
convex_df["COF5"] = df["COF5"].apply(lambda x: np.power(10, x))
convex_df.rename(columns={"COF5": "TX"}, inplace=True)
convex_df["COF6"] = df["COF6"].apply(lambda x: np.power(10, x))
convex_df.rename(columns={"COF6": "TE"}, inplace=True)
convex_df["COF7"] = df["COF7"].apply(lambda x: np.power(10, x))
convex_df.rename(columns={"COF7": "NU_T"}, inplace=True)
