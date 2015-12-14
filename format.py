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

#import G-Cell module
from functions import *
from proteins import *

#setup
time, location, SubList, events = Reactions().setup()
target = Reactions().Target('sequence.txt')

#Generate Monomer
#Reactions().Monomer('monomer',100,SubList,target)
Reactions().Complex('monomer',100,SubList)#Temporary placement
Reactions().Complex('complex',0,SubList)

#Events setting
Reactions().Events(Compose('complex',['monomer'],[3],1.0e-2),events)
Reactions().Events(Decompose('complex',['monomer'],[3],1),events)

#simulation
seq, mod = Reactions().Readseq(target,SubList)
logt, logd, t, tend = Showdata().logger(time, SubList, 0, 0.01)
Simulation().Run(t, tend, SubList, events, logt, logd, mod, location)

#showdata, make .png
Showdata().png(['monomer','complex'],logt, logd, SubList,'default')

#Finalize
Simulation().Save(mod,SubList)
Simulation().Makedata('default')
