# pygmoss: python implementation of GMOSS (Global MOdel for the radio Sky Spectrum)

pygmoss is the python implementation of the original GMOSS from the paper [paper](https://export.arxiv.org/abs/1607.07453) by Dr Mayuri Rao et. al. 
It can be applied to any task that requires the generation of all sky maps in the radio frequency range. It processes input all-sky maps mainly the 6 frequencies 22MHz,45MHz, 150MHz, 406MHz, 1440MHz and 23000MHz. Additional all-sky imput maps can be added as when they are obtained.

All sky radio maps used in the above paper are given in the pyGMOSS/data directory.

## Installation

pygmoss uses version pyhthon version 3.9 and above.
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

### Step 2: Check which model is applicable to that pixel 

Use the `check_concave_convex` function in `pygmoss.check_concave_convex` to determine if pixels are concave or convex. This function takes the brightness temperature file and convexity file as the inputs and returns a DataFrame containing convexity values(concave or convex). Save this DataFrame into a csv file(convexity file) for the next step.

### Step 3: For the pixels that use concave model

Use the `ConcaveModel`(check the documentation) class from `pygmoss.concave_model` to work with concave pixels. Provide the previously saved brightness temperature file and convexity files as inputs to initialize the class and run the `ConcaveModel.fit()` function. This function returns a DataFrame optimized concave model values. Save this DataFrame into csv file for next step.

### Step 4: To obtain initial parameters for convex model

Use the `ConvexModelInitialParameters` class from `pygmoss.convex_model_initial_parameters` file (check the documentation) class to obtain the initial parameters file for convex model function. Provide the previously saved brightness temperature file and convexity files as inputs to initialize the class and run the `ConvexModelInitialParameters.convex_model_initial_parameter_generator()` function. This function returns a DataFrame containing the initial parameters for convex model pixels. Save this DataFrame into a file (convex model initial parameter file) for the next step.

### Step 5: For pixels that use the convex model

Use the `ConvexModel`(check the documentation) class from `pygmoss.convex_model` to work with concave model pixels. Provide the previously saved brightness temperature file and convexity file and convex model initial parameter file as inputs to initialize the class and run the `ConvexModel.fit()` function. This function returns a DataFrame for convex pixels with optimized values. Save this DataFrame into a file for further use.

### Step 6: Use the obtained files to generate all sky maps 

In `tools` (`pygmoss.tools`) run the `all_sky_map()` function to generate the all sky map. `tools` has help function to merge dataframes by 'PIXEL' and other functions mainly concave and convex models (if working with the individual pixels is desired).

## Contributers

Dr. Mauri Rao, Dr. Ravi Subrahmanyan, Dr Saurabh Singh, Mohith P. A.