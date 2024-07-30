import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append('/users/yaguang/Onedrive/github/')
import h5py
# import asteroseismology as se
from astropy.io import ascii
from astropy.table import Table
import corner
import matplotlib.colors
import scipy
from multiprocessing import Pool
import emcee
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from scipy.optimize import curve_fit
# from grid import grid #, match_modes

# rootpath = os.getenv('WORK_DIR')+'numax-amlt/'
rootpath = '/Users/yaguang/My Drive/numax-amlt/'
sys.path.append(rootpath)
work_dir = rootpath+'numax/'