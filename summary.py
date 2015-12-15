#!/usr/bin/env python

#color setting
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

#import public module
import sys
import math
import re
import csv
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *
import glob

name = raw_input(GREEN+'Input dir name : '+ENDC)
if os.path.exists(name):
    os.chdir(os.getcwd()+'/'+name)
    os.system('open *.png')
else:print GREEN+name+RED+' NOT EXITSTS!!'+ENDC
if os.path.exists('substances.csv'):
    f = open('substances.csv','r')
    f = csv.reader(f)
    for line in f:
        if line[1] == 0:print RED+'NG '+ENDC+': '+RED+"{:<5}".format(line[0])+ENDC+"{:>33}".format(' = ')+BLUE+line[1]+ENDC
        else:print GREEN+'OK '+ENDC+': '+GREEN+"{:<5}".format(line[0])+ENDC+"{:>33}".format(' = ')+BLUE+line[1]+ENDC
