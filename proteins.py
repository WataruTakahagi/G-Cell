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

class YciV:#RNA/ssDNA exonuclease 5' -> 3'specific
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class YdaV:#Rac prophage; predicted DNA replication protein
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class zinc_binding_phosphatase:#zinc-binding phosphatase
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class CrfC:#Clamp-binding sister replication fork colocalization protein
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class RarA:#recombination factor
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class Hda:#regulator of DnaA that prevents premature reinitiation of DNA replication
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class DnaA_initiator_associating_factor:#DnaA initiator-associating factor for replication initiation
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class DNA_polymerase_III_holoenzyme:
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class Tus:#DNA-binding protein; inhibition of replication at Ter sites
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class DnaC:#chromosome replication; initiation and chain elongation
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class Dam:#DNA adenine methyltransferase
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class primosome:
    def __init__(self,location,mseq,k):
        self.location = location
        self.mseq = mseq
        self.k = k

    def propensity(self, state):
        self.p = self.k
        self.p = self.p * state[Reactions().Getindex('primosome',state)][1]
        return self.p

    def execute(self, state):
        if self.location >= 1 : mseq[self.location-1][Reactions().Getindex('primosome',state)] = 0
        self.mseq[self.location][Reactions().Getindex('primosome',state)] = 1
        self.mseq[self.location][Reactions().Getindex('ds',state)] = 1
        return self.location, self.mseq, state

class DnaA:
    def __init__(self, mod, k):
        self.mod = mod
        self.k = k
        self.OriC = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]

    def propensity(self, state):
        self.k = 1

    def execute(self, sublist, location):
        location = random.choice(self.OriC)[0]
        Reactions().ChangeStates('DnaA',location,self.mod,sublist,1)
        for l in range(location-5,location+5):
            if location-5 < 0: pass
            else: Reactions().ChangeStates('ds',l,self.mod,sublist,1)
        print self.mod[location,Reactions().Getindex('DnaA',sublist)]
        print self.mod[location,:]
        return sublist,location

class DNA_polymerase_III_holoenzyme:
    def __init__(self,location,mseq,k):
        pass

    def propensity(self, state):
        return

    def execute(self, state):
        return
