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

#import whole-ecoli replication module
from functions import Reactions
from functions import Showdata
from functions import Simulation

#import functions module
from functions import Compose
from functions import Decompose
from functions import dnaA
from functions import dnaB
from functions import dnaC
from functions import dnaG
from functions import RNaseH
from functions import SSB
from functions import Topo1
from functions import SDC
from functions import DNApol1
from functions import DNApol3
from functions import DNApol3holoenzyme

print RED+"Whole-genome based simulation platform of prokaryotic DNA replication"+ENDC

#Sequence data
seq,mod,time = Reactions().Readseq('sequence.txt')

#Substance setting
SubList = []                                                  #ID
SubList.append(Reactions().Monomer('dnaA',60000))             #0
SubList.append(Reactions().Monomer('dnaB',60000))             #1
SubList.append(Reactions().Monomer('dnaC',60000))             #2
SubList.append(Reactions().Monomer('dnaG',20000))             #3
SubList.append(Reactions().Monomer('RNaseH',500))             #4
SubList.append(Reactions().Monomer('SSB',30000))              #5
SubList.append(Reactions().Monomer('Topo1',0))                #6
SubList.append(Reactions().Monomer('SDC',200000))             #7
SubList.append(Reactions().Monomer('DNApol1',0))              #8
SubList.append(Reactions().Monomer('DNApol3',100))            #9
SubList.append(Reactions().Monomer('DNApol3holoenzyme',100))  #10
SubList.append(Reactions().Monomer('OriC9',5))                #11
SubList.append(Reactions().Monomer('OriC13',3))               #12
SubList.append(Reactions().Monomer('ATP',10000))              #13
SubList.append(Reactions().Complex('dnaB','dnaC',0))          #14
SubList.append(Reactions().Complex('dnaA','dnaB/dnaC',0))     #15
SubList.append(Reactions().Complex('dnaA','ATP',0))           #16
SubList.append(Reactions().Complex('OriC9','dnaA',0))         #17
SubList.append(Reactions().Complex('OriC13','dnaA',0))        #18

#Gillespie Test
logt, logd, t, tend = Showdata().logger(time, SubList, 0, 1)
events = [Compose('dnaB','dnaC',0.01),Decompose('dnaB/dnaC',5),Compose('dnaA','dnaB/dnaC',0.01),Decompose('dnaA/dnaB/dnaC',4)]
Simulation().run(t, tend, SubList, events, logt, logd)
Showdata().figure("dnaA", logt, logd, SubList)
Showdata().figure("dnaB", logt, logd, SubList)
Showdata().figure("dnaC", logt, logd, SubList)
Showdata().figure("dnaB/dnaC", logt, logd, SubList)
Showdata().figure("dnaA/dnaB/dnaC", logt, logd, SubList)
Showdata().save("process.png")
"""
#Polymerization Process Test
r = [0]
i = 0
time = []
er = []
process = len(seq)
for location in range(process):
    dnaB().execute(location,mod,SubList)
for location in range(process):
    DNApol3().execute(location,mod,SubList,r)
    time.append(i)
    er.append(r[0])
    i += 1

#Save result
Simulation().Wcplot(time,er)
Simulation().Wcwrite(mod)
"""
Simulation().Makedata()
