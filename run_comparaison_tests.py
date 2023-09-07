#!/usr/bin/env python3
"""
For a few parameter sets given in compare.test, check that the distribution obtained with our improved gillespie simulation is similar to the one obtained with the standard gillespie algorithm
"""

import subprocess
import pandas as pd
import numpy as np
import scipy.stats
import atreyu

def dist_from_output(r):
    return np.array([int(y[0]) for x in r.stdout.decode("utf-8").split("\n") if len(y:=x.split())>0])

def compare(data1, data2):
    print(np.abs(np.median(data1)-np.median(data2)) / max(np.median(data1), np.median(data2)))
    print(scipy.stats.wasserstein_distance(data1, data2))

with open('compare.test', 'r') as f:
    for i, l in enumerate(f):
        if l.startswith("#"):
            continue
        arguments = l.split('#')[0].strip()
        print(arguments, ' -> ', end='')
        dist_gillespie = dist_from_output(subprocess.run(f"./gillespie {arguments}", shell=True, capture_output=True))
        dist_atreyu_forward_simulator = dist_from_output(subprocess.run(f"./atreyu_forward_simulator {arguments}", shell=True, capture_output=True))
        fig, p = atreyu.compare_dists(dist_gillespie, dist_atreyu_forward_simulator, ['standard Gillespie','atreyu forward simulator'])
        fig.savefig(f"compare.test.{i}.pdf")
        if p<0.01:
            print(f"\033[31m Warning, with parameters {arguments}, we found different outputs for Atreyu and Gillespie \033[m")


