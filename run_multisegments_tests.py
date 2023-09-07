#!/usr/bin/env python3
"""
Check that the version with -DCONSTANT_MUTATION_RATE functions as expected
"""

import sys
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

n = 0
def check_same_output(programs, arguments):
    global n
    n += 1
    assert(len(programs)==2)
    dist1 = dist_from_output(subprocess.run(f"./{programs[0]} {arguments[0]}", shell=True, capture_output=True))
    dist2 = dist_from_output(subprocess.run(f"./{programs[1]} {arguments[1]}", shell=True, capture_output=True))
    if (dist1 != dist2).any():
        print(f"\033[31m Warning, we found different outputs for {programs[0]} with parameters {arguments[0]} than for {programs[1]} with parameters {arguments[1]} \033[m")
    else:
        print(f"  test {n} passed")

print("Testing the version with CONSTANT_MUTATION_RATE... ")
check_same_output(['atreyu_forward_simulator','atreyu_forward_simulator_ct'],['100,1E6 0.5 2e-6 0.7 0.2 1 42', '100,1E6 0.5 2e-6 0.7 0.2 1 42'])
check_same_output(['atreyu_forward_simulator','atreyu_forward_simulator_ct'],['100,1E6,1E7 0.5,0.3 2e-6,5e-7 0.7 0.2 1 42', '100,1E6,1E7 0.5,0.3 2e-6,5e-7 0.7 0.2 1 42'])
    

