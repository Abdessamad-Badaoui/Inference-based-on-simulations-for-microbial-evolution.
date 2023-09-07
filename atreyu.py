import subprocess
import os
import numpy as np
import seaborn as sns
import scipy.stats
import random
import time
import multiprocessing 
import queue
import sys
import noisyopt
from matplotlib import pyplot as plt

plt.rcParams['figure.figsize'] = [16, 6]
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams["figure.autolayout"] = True
tolbright = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']
tollight = ['#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD', '#000000']
plt.rcParams["axes.prop_cycle"] = plt.cycler('color', tolbright)


def compare_dists(data1, data2, l):
    """
    Compare two distributions of the number of mutants
    (used by tests of the forward simulator)
    Plot the probability densities and cumulative densities
    And perform a Kolmogorov-Smirnov test
    """
    ks, p = scipy.stats.ks_2samp(data1, data2, mode='asymp')
    if p < 0.05:
        print('Different distributions, ks={} and p={}'.format(ks,p))
    else:
        print('Similar distributions, ks={} and p={}'.format(ks,p))

    p95 = int(max(np.percentile(data1,95), np.percentile(data2,95)))  # Cutoff x axis at 95%
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # CDF
    sns.ecdfplot(data1, ax=ax1)
    sns.ecdfplot(data2, ax=ax1)
    ax1.set_xlabel('Number of mutants')
    ax1.set_ylabel('Cumulative probability')
    ax1.set_xlim([0, p95])
    ax1.legend(l,loc='upper left')

    # PDF (raw)
    d1 = np.bincount(data1)
    d2 = np.bincount(data2)
    dw1 = d1 / sum(d1) # density from count
    dw2 = d2 / sum(d2)
    ax2.plot(dw1)
    ax2.plot(dw2)
    ax2.set_xlabel('Number of mutants')
    ax2.set_ylabel('Probability')
    ax2.set_xlim([0, p95])
    ax2.set_ylim([0, max(dw1.max(),dw2.max())*1.2])
    ax2.legend(l,loc='upper right')

    fig.tight_layout()

    return fig, p


