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
import random

#import functions
from functions import *

class DnaA:
    def __init__(self, mod, k):
        self.mod = mod
        self.k = k
        self.OriC = [[0,1],[20,21],[30,31],[40,41],[50,51],[60,61],[70,71],[80,81]]
        self.location = random.choice(self.OriC)[0]

    def propensity(self, state, location):
        self.p = self.k * state['DnaA']*self.mod[self.location][len(state)]
        return self.p,self.location

    def execute(self, sublist, location):
        Reactions().ChangeStates('DnaA',self.location,self.mod,sublist,1)
        for l in range(self.location-5,self.location+5):
            if location-5 < 0: pass
            else: Reactions().ChangeStates('ds',l,self.mod,sublist,1)
        sublist['DnaA'] -= 1
        print "DnaA FIRE"
        return sublist

class primosome:
    def __init__(self,mod,k):
        self.mod = mod
        self.k = k

    def propensity(self, state):
        self.p = self.k * state['primosome'][1]
        return self.p

    def execute(self, state):
        if self.location >= 1 : mseq[self.location-1][Reactions().Getindex('primosome',state)] = 0
        self.mseq[self.location][Reactions().Getindex('primosome',state)] = 1
        self.mseq[self.location][Reactions().Getindex('ds',state)] = 1
        return self.location, self.mseq, state

class DNA_polymerase_III_holoenzyme:
    def __init__(self,location,mseq,k):
        pass

    def propensity(self, state):
        return

    def execute(self, state):
        return
