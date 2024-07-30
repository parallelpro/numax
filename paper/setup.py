from matplotlib import rcParams
rcParams["figure.dpi"] = 100
rcParams["savefig.dpi"] = 100

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors 
from matplotlib.colors import ListedColormap
import os 
import sys
# sys.path.append('/Users/yaguang/Onedrive/github/')
from astropy.io import ascii
# import asteroseismology as ast 
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from multiprocessing import get_context, Pool
import h5py
import scipy
import emcee
import corner

os.environ['TEXINPUTS'] = ".:$HOME/texmf/custom-styles//:"
os.environ['PATH'] = "/Library/Tex/texbin:" + os.getenv('PATH')

# color stuff
# red = sns.xkcd_rgb["pale red"]
# blue = sns.xkcd_rgb["denim blue"]
# green = sns.xkcd_rgb["faded green"]
# orange = sns.xkcd_rgb["amber"]
# grey = sns.xkcd_rgb["greyish"]
# darkgrey = sns.xkcd_rgb["dark grey"]
# purple = sns.xkcd_rgb["dusty purple"]
# black = sns.xkcd_rgb["black"]

red = 'indianred' #sns.xkcd_rgb["pale red"]
blue = 'steelblue' #sns.xkcd_rgb["denim blue"]
lightblue = 'lightskyblue' #sns.xkcd_rgb["light blue"]
green = 'limegreen' #sns.xkcd_rgb["faded green"]
orange = 'darkorange' #sns.xkcd_rgb["amber"]
cyan = 'cyan' #sns.xkcd_rgb["cyan"]
grey = 'lightgray' #sns.xkcd_rgb["greyish"]
darkgrey = 'darkgray' #sns.xkcd_rgb["dark grey"]
purple = 'darkmagenta' #sns.xkcd_rgb["dusty purple"]
black = 'k' #sns.xkcd_rgb["black"]

colors = ['#E69F00','#56B4E9','#009E73','#D55E00','#CC79A7','#0072B2','#F0E442',]
linestyles = ["-", "--", "-.", ":", (0, (3, 5, 1, 5, 1, 5)), "--", "-."]
markers = ['o', '^', 's', 'v', '>', '<']

# def cmap_diverging(n=10):
#     return ListedColormap(sns.diverging_palette(220, 20, n=n))

# def cmap_grey(n=10):
#     return ListedColormap(sns.color_palette("Greys", n))

# def blues(n=10):
#     return sns.color_palette("Blues", n)

rootpath = '/Users/yaguang/My Drive/numax-amlt/'
overleaf_path = '/Users/yaguang/Dropbox/Apps/Overleaf/Yaguang_numax'
work_path = rootpath
sys.path.append(work_path)

def to_overleaf(figure_file_name, subdir='figures'):
    #  return "rclone copy {:s} remote:Apps/Overleaf/Yaguang_numax/{:s}/".format(figure_file_name, subdir)
    return "cp {:s} {:s}/{:s}/".format(figure_file_name, overleaf_path, subdir)

    
# from lib.toolkit import *

fontsize = 9 # minimum fontsize is 5
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.size"] = fontsize #7.5
matplotlib.rcParams["legend.fontsize"] = fontsize#7.5
matplotlib.rcParams['text.usetex'] = True #False
matplotlib.rcParams['axes.labelsize'] = fontsize#7
matplotlib.rcParams['xtick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.direction']='out'
matplotlib.rcParams['ytick.major.size']=3.0
matplotlib.rcParams['ytick.minor.size']=2.0
matplotlib.rcParams['xtick.direction']='out'
matplotlib.rcParams['xtick.major.size']=3.0
matplotlib.rcParams['xtick.minor.size']=2.0

# matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}\usepackage{sfmath}'
# matplotlib.rcParams['font.family'] = 'sans-serif' #'sans-serif'
# matplotlib.rcParams['font.sans-serif'] = ['Helvetica'] #'Helvetica' 

# matplotlib.rcParams['text.latex.preamble'] = ''
# matplotlib.rcParams['font.family'] = 'serif' 

matplotlib.rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}\usepackage{charter}'
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ["Charter"]

# mnras size in pt
columnwidth = 240
textwidth = 504

# aas size in pt
textwidth =  513.11743
columnwidth = 242.26653

def figsize(column="one", square=False, ratio=None):
    # Thanks Dan!
    # Parameters:
    # column: "one" or "double"
    # square: True or False
    # ratio: height/width

    inches_per_pt = 1.0/72.00              # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0     # Most aesthetic ratio
    if (ratio == None): ratio = golden_mean
    if (column == "one"):
        fig_width_pt = columnwidth
    elif (column == "double"):
        fig_width_pt = textwidth
    elif (column == "triple"):
        fig_width_pt = columnwidth*3
    elif (column == "quadruple"):
        fig_width_pt = textwidth*2
    else:
        raise ValueError("column should be one of ``one'' or ``double''. ")
    fig_width = fig_width_pt*inches_per_pt # Figure width in inches
    if square:
        fig_height = fig_width
    else:
        fig_height = fig_width*ratio
    return [fig_width,fig_height]

errstyle = {'capsize':2, 'ecolor':darkgrey, 'elinewidth':1, 'capthick':1, 'linestyle':'None'}


def quantile(x, q, weights=None):
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)

    if weights is None:
        return np.percentile(x, 100.0 * q)
    else:
        weights = np.atleast_1d(weights)
        idx = np.argsort(x, axis=0)

        res = []
        for i in range(x.shape[1]):
            sw = weights[idx[:,i]]
            cdf = np.cumsum(sw)[:-1]
            cdf /= cdf[-1]
            cdf = np.append(0, cdf)
            res.append(np.interp(q, cdf, x[idx[:,i],i]))
        return np.array(res).T


def get_binned_median(xdata, ydata, bins):
    '''
    return median vals and the standand deviation of the median for given data and a bin.

    Input:
        xdata: array-like[N,]
        ydata: array-like[N,]
        bins: array-like[Nbin,]

    Output:
        xcs: array-like[Nbin,]
            center of each bin
        medians: array-like[Nbin,]
            median values of ydata in each bin
        stds: array-like[Nbin,]
            stds of the median of ydata in each bin

    '''
    
    Nbin = len(bins)-1
    medians, stds = [np.zeros(Nbin) for i in range(2)]
    for ibin in range(Nbin):
        idx = (xdata>bins[ibin]) & (xdata<=bins[ibin+1]) & (np.isfinite(xdata)) & (np.isfinite(ydata))
        if np.sum(idx)==0:
            medians[ibin], stds[ibin] = np.nan, np.nan
        else:
            medians[ibin] = np.median(ydata[idx])
            stds[ibin] = 1.253*np.std(ydata[idx])/np.sqrt(np.sum(idx))
    xcs = (bins[1:]+bins[:-1])/2.
    return xcs, medians, stds
