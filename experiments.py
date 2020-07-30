#!/usr/bin/env python3
import os
import sys
import subprocess
import errno

original_dir = os.curdir

try:
    os.mkdir('temp')
except FileExistsError:
    print("WARNING: 'temp/' already exists, its contents may be overwritten.", file=sys.stdout)
finally:
    os.chdir('temp/')

network = '../src/parameters/network.mat'
experiments = '../src/parameters/experiments/'

for config in os.listdir(experiments):
    print("\nRunning experiment configuration", config, "...")
    sim = ['octave', '../src/run_experiment.m', network, experiments + config, '../src']
    proc = subprocess.run(sim, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(proc.stdout.decode())
    print(proc.stderr.decode(), file=sys.stderr)

figures = '../doc/img/'

try:
    os.mkdir(figures)
except FileExistsError:
    if input("\nWARNING: 'doc/img/' already exists, continue? [y/n]: " ) != 'y':
        exit(errno.EEXIST)

for fig in os.listdir():
    crop = ['pdfcrop', fig, figures + fig]
    proc = subprocess.run(crop, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(proc.stdout.decode())
    print(proc.stderr.decode(), file=sys.stderr)
    if proc.returncode == 0:
        os.remove(fig)
