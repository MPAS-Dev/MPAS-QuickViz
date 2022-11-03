#!/usr/bin/env python

# note: the calls are not strictly correct because they have redirection in them!!!

import subprocess as sp
import numpy as np

# parameters
ts = np.array([25, 20, 5, 2])
te = np.array([30, 25, 10,7])

def run_160km():
    for ats, ate in zip(ts,te):
        cmd = 'time python diffusivity_binning.py -f  all_rlzn_particle_data_rads.npz -r  160000 --ts ' + str(ats) + ' --te ' + str(ate) + ' --dt '+ '2. &> log160km'+str(ats)+'-' + str(ate) + '.log'
        print cmd
        sp.call(cmd, shell=True)

def run_320km():
    for ats, ate in zip(ts,te):
        cmd = 'time python diffusivity_binning.py -f  all_rlzn_particle_data_rads.npz -r  320000 --ts ' + str(ats) + ' --te ' + str(ate) + ' --dt '+ '2. &> log320km'+str(ats)+'-' + str(ate) + '.log'
        print cmd
        sp.call(cmd, shell=True)

def run_16km():
    for ats, ate in zip(ts,te):
        cmd = 'time python diffusivity_binning.py -f  all_rlzn_particle_data_rads.npz -r  16000 --ts ' + str(ats) + ' --te ' + str(ate) + ' --dt '+ '2. &> log16km'+str(ats)+'-' + str(ate) + '.log'
        print cmd
        sp.call(cmd, shell=True)

def run_64km():
    for ats, ate in zip(ts,te):
        cmd = 'time python diffusivity_binning.py -f  all_rlzn_particle_data_rads.npz -r  64000 --ts ' + str(ats) + ' --te ' + str(ate) + ' --dt '+ '2. &> log64km'+str(ats)+'-' + str(ate) + '.log'
        print cmd
        sp.call(cmd, shell=True)

if __name__=="__main__":
    #run_160km()
    run_320km()
    run_16km()
    run_64km()
