import pdb
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.special as special
import scipy.integrate as integrate
import math
import scipy as sp
import pandas as pd


class ConvexModel:
    def __init__(
        self,
        path_to_brightness_temparature_file: str,
        path_to_convexity_file: str,
        path_to_convex_model_initial_parameter_file: str,
        speed_of_light: float = 2.99792458e08,
        mass_of_electron: float = 9.1e-31,
        charge_of_electron: float = 1.6e-19,
        sine_alpha: float = 1.0,
        magnetic_field: float = 1e-9,
        GSPAN: float = 100.0,
        verbose: bool = False,
    ) -> pd.DataFrame:

        self.path_to_brightness_temparature_file = path_to_brightness_temparature_file
        self.path_to_convex_model_initial_parameter_file = (
            path_to_convex_model_initial_parameter_file
        )

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

    def func(self) -> float:  # pp is an array
        df = pd.read_csv(self.path_to_brightness_temparature_file)
        # convexity_df = pd.read_csv(self.path_to_convexity_file)
        convex_model_initial_param_df = pd.read_csv(
            self.path_to_convex_model_initial_parameter_file
        )
        pixels = convex_model_initial_param_df.loc[:, "PIXEL"].values
        fnorm_list = convex_model_initial_param_df.loc[:, "FNORM"].values
        alpha1_list = convex_model_initial_param_df.loc[:, "ALPHA1"].values
        alpha2_list = convex_model_initial_param_df.loc[:, "ALPHA2"].values
        nu_break_list = convex_model_initial_param_df.loc[:, "NU_BREAK"].values
        Tx_list = convex_model_initial_param_df.loc[:, "TX"].values
        Te_list = convex_model_initial_param_df.loc[:, "TE"].values
        nu_t_list = convex_model_initial_param_df.loc[:, "NU_T"].values

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
                m, m,

        convexity_df = convexity_df.loc[
            convexity_df.loc[:, "Concave/Convex"] == "Convex", :
        ]
        df = pd.merge(df, convexity_df, on="PIXEL", how="right")
        pixels = df.loc[:, "PIXEL"].values
        alpha_1_list = df.loc[:, "ALPHA_1"].values
        alpha_2_list = df.loc[:, "ALPHA_2"].values

        # unpacking pp
        fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
        # blowing up the
        if alpha1 < 2.0 or alpha1 > 3.0:
            alpha1 = 1000000000
        if alpha2 < 2.0 or alpha2 > 3.0:
            alpha2 = 1000000000
        if Te < 0.0 or Te > 10000:
            Te = 1000000000

    # computing chi square
    chisq = 0.0
    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu * 1e9 / GSPAN
        nu_max = nu * 1e9 * GSPAN
        gama_min = np.sqrt(nu_min / scale_gam_nu)
        gama_max = np.sqrt(nu_max / scale_gam_nu)
        gama_break = np.sqrt(nu_break / scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)

        elif xu < xb:
            C1 = alpha1
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)
        else:
            xu = xb
            C1 = alpha1
            I1, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I1 *= np.power(gama_break, 2 * C1 - 3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            I2, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha2, nu))
            I2 *= np.power(gama_break, 2 * C1 - 3)
            I = I1 + I2
        extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
        FFIT = fnorm * (
            (np.power(nu, -2) * I) + Tx * np.power(nu, -2.1)
        ) * extn + Te * (1.0 - extn)
        DDIF = (b_temp[i] - FFIT) / (b_temp[i])
        if i <= 5:
            chisq += DDIF * DDIF
    chisq /= 6.0
    print(chisq)
    return chisq


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
    func,
    x0=x_ini,
    method="Nelder-Mead",
    bounds=bounds,
    options={"verbose": 1, "maxiter": 10000},
)
print(x_ini)
print(type(result.x))
print(result)
xs = np.linspace(22e-3, 24, 100)


frequency = np.array([22.0, 45.0, 150.0, 408.0, 1420.0, 23000.0])  # x_values
frequency = np.array([np.float32(f * 10**-3) for f in frequency])
b_temp = np.array(
    [
        9.18573758e04,
        1.77507604e04,
        7.10610657e02,
        6.49989393e01,
        2.11183872e00,
        9.89014738e-04,
    ]
)  # y_values
GSPAN = 100
# adopting Dr. Rao's method and writing a function called func that has to be minimized
# therefore essentially a code that computes chi_square
NHPIX = 3072

cvel = 2.99792458e08  # m s^-1
m_e = 9.1e-31
q_e = 1.6e-19
sin_alph = 1.0
Bmag = 1e-9  # Tesla == 10 micro-Gauss
scale_gam_nu = (3.0 * q_e * Bmag * sin_alph) / (4.0 * np.pi * m_e * cvel)

df = pd.read_csv("brightness_temp_at_pixel.csv")
frequency = np.array([22.0, 45.0, 150.0, 408.0, 1420.0, 23000.0])  # x_values
frequency_string = np.array(["22", "45", "150", "408", "1420", "23000"])
b_temp = np.zeros(len(frequency_string))


df = pd.read_csv("brightness_temp_at_pixel.csv")
initial_arguments_df = pd.read_csv("convex_model_initial_params.csv")
df = pd.merge(df, initial_arguments_df, on="PIXEL", how="right")

for i, f in enumerate(frequency_string):
    b_temp[i] = df.loc[df.loc[:, "PIXEL"] == pixel, f"{f}MHz"]

# x_ini = np.array([ 9.13659851e-07,  9.66942286e+00,  2.00531364e+00,  8.86653501e+01,
#        -3.66142917e-10,  9.99999996e+03,  1.19717340e-02])
# x_ini = np.array([ -1.976754550046829e-07, 2.6728667075093107, 2.7477254162083455, 247386337.5370596/1e9, 1e-10, 8000.0, 0.001])
# x_ini = np.array([-5.761185054330454e-08, 2.6728667075093107, 2.7477254162083455, 247386337.5370596/1e9, 1e-10, 8000.0, 0.001]) #working pix 36
x_ini = np.array(
    [
        -5.766426064650115e-08,
        2.6728667075093107,
        2.7477254162083455,
        247386337.5370596 / 1e9,
        1e-10,
        8000.0,
        0.001,
    ]
)  # broken leg pix 36
# x_ini = np.array([ -6.928880, 2.6728667075093107, 2.7477254162083455, 247.3863375370596, 1e-10, 8000.0, 0.001])


def F(x):
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
        quad_term = 2429625 + 2 * x * (-1922325 + (5418382 * x) + 83221732 * (x**2))
        error_function_term = (
            1196306216
            * np.exp(x)
            * np.pi
            * (x ** (7 / 2))
            * sp.special.erfc(np.sqrt(x))
        )
        return exponential_term * ((const * quad_term) - error_function_term)


def integrand_for_convex(gama, alpha, nu):
    nu_c = (scale_gam_nu * (gama**2)) / 1e9
    x = nu / nu_c
    integrand_ = F(x) * x * np.power(gama, -1 * (2 * alpha - 3))
    return integrand_


def func(pp: np.ndarray) -> float:  # pp is an array
    # unpacking pp
    fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
    # blowing up the
    if alpha1 < 2.0 or alpha1 > 3.0:
        alpha1 = 1000000000
    if alpha2 < 2.0 or alpha2 > 3.0:
        alpha2 = 1000000000
    if Te < 0.0 or Te > 10000:
        Te = 1000000000

    # computing chi square
    chisq = 0.0
    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu * 1e9 / GSPAN
        nu_max = nu * 1e9 * GSPAN
        gama_min = np.sqrt(nu_min / scale_gam_nu)
        gama_max = np.sqrt(nu_max / scale_gam_nu)
        gama_break = np.sqrt(nu_break / scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)

        elif xu < xb:
            C1 = alpha1
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)
        else:
            xu = xb
            C1 = alpha1
            I1, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I1 *= np.power(gama_break, 2 * C1 - 3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            I2, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha2, nu))
            I2 *= np.power(gama_break, 2 * C1 - 3)
            I = I1 + I2
        extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
        FFIT = fnorm * (
            (np.power(nu, -2) * I) + Tx * np.power(nu, -2.1)
        ) * extn + Te * (1.0 - extn)
        DDIF = (b_temp[i] - FFIT) / (b_temp[i])
        if i <= 5:
            chisq += DDIF * DDIF
    chisq /= 6.0
    print(chisq)
    return chisq


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
    func,
    x0=x_ini,
    method="Nelder-Mead",
    bounds=bounds,
    options={"verbose": 1, "maxiter": 10000},
)
print(x_ini)
print(type(result.x))
print(result)
xs = np.linspace(22e-3, 24, 100)


def convex_func(frequency, pp):
    # unpacking pp
    fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
    # blowing up the
    # if alpha1 < 2.0 or alpha1 > 2.0: alpha1 = 1000000000000
    # if alpha2 < 2.0 or alpha2 > 3.0: alpha2  = 1000000000000
    # if Te < 0.0 or Te > 10000: Te = 1000000000000
    b_temps = []

    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu * 1e9 / GSPAN
        nu_max = nu * 1e9 * GSPAN
        gama_min = np.sqrt(nu_min / scale_gam_nu)
        gama_max = np.sqrt(nu_max / scale_gam_nu)
        gama_break = np.sqrt(nu_break / scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)

        elif xu < xb:
            C1 = alpha1
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I *= np.power(gama_break, 2 * C1 - 3)
        else:
            xu = xb
            C1 = alpha1
            I1, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha1, nu))
            I1 *= np.power(gama_break, 2 * C1 - 3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            I2, _ = integrate.quad(integrand_for_convex, xl, xu, args=(alpha2, nu))
            I2 *= np.power(gama_break, 2 * C1 - 3)
            I = I1 + I2
        extn = np.exp(-1.0 * np.power((nu_t / nu), 2.1))
        FFIT = fnorm * (
            (np.power(nu, -2) * I) + Tx * np.power(nu, -2.1)
        ) * extn + Te * (1.0 - extn)
        b_temps.append(FFIT)
    return b_temps


xs = np.linspace(22e-3, 24, 100)
plt.plot(xs, convex_func(xs, np.array(x_ini)), label="this is the initial guess")
plt.plot(frequency, b_temp, "r*")
plt.plot(
    frequency,
    np.array(convex_func(frequency, result.x)),
    label="the curve after fitting",
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("log Frequency[MHz]")
plt.ylabel("log Temp[K]")
plt.legend()
plt.grid()
plt.show()
# plt.plot(result.x)
