import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import seaborn as sns
from scipy.optimize import minimize
import time 
import glob
import os

pixels = 3072

all_files = glob.glob('*txt')
required_files = []

length_ = len(all_files)
for file in all_files:
    if file.startswith('map'):
        required_files.append(file)



root = required_files[0]
root = root.split('_')
print(root[1])