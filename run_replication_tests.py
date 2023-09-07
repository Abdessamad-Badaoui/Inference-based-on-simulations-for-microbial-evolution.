#!/usr/bin/env python3
"""
Very simple tests to check that the behaviour of the simulations does not change: for a given set of parameters we check that the output is perfectly equal to the predefined one (with the specified random seed)
"""

import subprocess

for prog in ['gillespie','atreyu_forward_simulator']:
    print(f'testing program {prog}...')
    with open(f'replication.{prog}.test', 'r') as f:
        for i, l in enumerate(f):
            t = l.split(':')
            if len(t) != 2:
                continue
            arguments = t[0].strip(' "')
            expected = t[1].strip(' "\n')
            r = subprocess.run(f"./{prog} {arguments}", shell=True, capture_output=True)
            obtained = r.stdout.decode("utf-8").split("\n")[-2]
            if obtained == expected:
                print(f"  line {i+1} passed")
            else:
                print(f'\033[31m line {i+1} failed, obtained "{obtained}", expected "{expected}"\033[m')
                exit(-1)
    print('...pass')
