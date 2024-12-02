# pygmoss: python implementation of GMOSS (Global MOdel for the radio Sky Spectrum)

pygmoss is the python implementation of the original GMOSS from the [paper](https://export.arxiv.org/abs/1607.07453) by Dr Mayuri Rao et. al. 
It can be applied to any task that requires the generation of all sky maps in the radio frequency range. It processes input all-sky maps mainly the 6 frequencies 22MHz,45MHz, 150MHz, 406MHz, 1440MHz and 23000MHz. Additional all-sky imput maps can be added as when they are obtained.

All sky radio maps used in the above paper are given in the pyGMOSS/data directory.

## Installation

pygmoss uses version python version 3.9 and above.
To install pygmoss, use `pip`:

```bash
pip3 install pygmoss
```
(optional)

```bash
pip3 install -r requirements.txt
```
## Usage
### Step 1: Extract Brightness Temperature Values from the input all sky maps

To extract brightness temperature values from input all-sky maps, use the `extract_brightness_temperature_at_pixel` in `pygmoss.read_data_at_pixel`. It requires the path to the folder containing the input all-sky maps(requested to go through the function documentation). This function returns a DataFrame containing brightness temperature values. Save this DataFrame into a csv(using pandas.to_csv function) file(brightness tempearature file) for the next step.
```
from pygmoss.read_data_at_pixel import extract_brightness_temp_at_pixel

PIXELS = 3072
CORRECTION_150_OFFSET = 21.4
CORRECTION_150_SCALING = 1.05
TCMB = 2.72548

df = extract_brightness_temp_at_pixel(
    "/mnt/c/Users/mohit/test_gmoss",
    PIXELS,
    CORRECTION_150_OFFSET,
    CORRECTION_150_SCALING,
    150,
    TCMB,
)
df.to_csv("brightness_temp_per_pixel.csv", index=False)
```

### Step 2: Check which model is applicable to that pixel 

Use the `check_concave_convex` function in `pygmoss.check_concave_convex` to determine if pixels are concave or convex. This function takes the brightness temperature file and convexity file as the inputs and returns a DataFrame containing convexity values(concave or convex). Save this DataFrame into a csv file(convexity file) for the next step.
```
from pygmoss.check_concave_convex import check_concave_convex

load_dotenv()
DATA = os.getenv("DATA")
path_to_brightness_temperature_file = DATA + "brightness_temp_per_pixel.csv"

df = check_concave_convex(path_to_brightness_temperature_file)
df.to_csv(f"{DATA}convexity.csv", index=False)
```

### Step 3: For the pixels that use concave model

Use the `ConcaveModel`(check the documentation) class from `pygmoss.concave_model` to work with concave pixels. Provide the previously saved brightness temperature file and convexity files as inputs to initialize the class and run the `ConcaveModel.fit()` function. This function returns a DataFrame optimized concave model values. Save this DataFrame into csv file for next step.

```
from pygmoss.concave_model import ConcaveModel  


model = ConcaveModel("brightness_temp_per_pixel.csv", "convexity.csv")
df = model.fit()
df.to_csv("concave_pixel_fits.csv", index=False)
```

### Step 4: To obtain initial parameters for convex model

Use the `ConvexModelInitialParameters` class from `pygmoss.convex_model_initial_parameters` file (check the documentation) class to obtain the initial parameters file for convex model function. Provide the previously saved brightness temperature file and convexity files as inputs to initialize the class and run the `ConvexModelInitialParameters.convex_model_initial_parameter_generator()` function. This function returns a DataFrame containing the initial parameters for convex model pixels. Save this DataFrame into a file (convex model initial parameter file) for the next step.

```
from pygmoss.convex_model_initial_parameters import ConvexModelInitialParameters
convex_model_initial_parameters = ConvexModelInitialParameters(
    "brightness_temp_per_pixel.csv", "convexity.csv", verbose=True
)
df = convex_model_initial_parameters.convex_model_initial_parameter_generator()
df.to_csv("convex_model_initial_parameters.csv")
```

### Step 5: For pixels that use the convex model

Use the `ConvexModel`(check the documentation) class from `pygmoss.convex_model` to work with concave model pixels. Provide the previously saved brightness temperature file and convexity file and convex model initial parameter file as inputs to initialize the class and run the `ConvexModel.fit()` function. This function returns a DataFrame for convex pixels with optimized values. Save this DataFrame into a file for further use.

```
from pygmoss.convex_model import ConvexModel

path_to_brightness_temperature_file = "brightness_temp_per_pixel.csv"
path_to_convexity_file = "convexity.csv"
path_to_initial_parameters_file = "convex_model_initial_parameters.csv"
path_to_convex_fits_file = "convex_pixel_fits.log"
convex_model = ConvexModel(
    path_to_brightness_temperature_file,
    path_to_convexity_file,
    path_to_initial_parameters_file,
    path_to_convex_fits_file,
    verbose=True,
    start_pix=0,
    stop_pix=0,
)
convex_model.fit()
```

### Step 6: Use the obtained files to generate all sky maps 

In `utils` (`pygmoss.tools`) run the `all_sky_map()` function to generate the all sky map. Run the `btemp_vs_frequency()` to plot brightness temperature vs frequency for a specific pixel. In `utils` has help function `combine_df()` to merge dataframes by 'PIXEL' and other functions mainly concave and convex models called `concave_func()` and `convex_func()` (if working with the individual pixels is desired).

```
from pygmoss.utils import all_sky_map
from pygmoss.utils import btemp_vs_frequency
from pygmoss.utils import concave_func
from pygmoss.utils import convex_func


frequency = 1200 * 1e-3  # GHz
series = all_sky_map(
    "concave_pixel_fits.csv",
    "convex_pixel_fits.log",
    frequency,
    number_of_pixels=3072,
    plot_all_sky=True,
)
series.to_csv("all_sky_map_{frequency/1e-3}.csv", index=False, header=False)

pixel = 2060
btemp_vs_frequency(
    "concave_pixel_fits.csv",
    "convex_pixel_fits.log",
    "brightness_temp_per_pixel.csv",
    pixel_val=pixel,
)

model_concave_value = concave_func(
            1.7,
            -1.7673603339907888,
            6.209947892799811,
            2.8215667328790097,
            2.6599927540080444,
            0.16663361137887173,
            9999.790968653713,
            0.00020113422862317362
        )
model_convex_value = convex_func(
            1.7,
            8.317511113295733e-07,
            2.51968636124662,
            2.707635701950907,
            247386337.5370596,
            7338.466282609452,
            8000.0,
            0.001
        )
```

## Contributers

Dr. Mayuri Rao, Dr. Ravi Subrahmanyan, Dr Saurabh Singh, Mohith P. A.