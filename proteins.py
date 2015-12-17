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

class General:
    def __init__(self):
        pass

    def Walk(self,name,location,mod,sublist,l):
        Reactions().ChangeStates(name,location,mod,sublist,0)
        Reactions().ChangeStates(name,location+l,mod,sublist,1)
        return location+l

    def RingOpen(self,mod,sublist,location,lg):
        for l in range(location,location+lg):
            if location < 0:Reactions().ChangeStates('ds',len(mod)+l,mod,sublist,1)
            elif len(mod) < location:Reactions().ChangeStates('ds',l-len(mod),mod,sublist,1)
            else:Reactions().ChangeStates('ds',l,mod,sublist,1)

    def Polymeraization(self,location,mod,sublist,rl):
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.23: mod[location][Reactions().Getindex('r`',sublist)]=0.76
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.26: mod[location][Reactions().Getindex('r`',sublist)]=0.73
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.24: mod[location][Reactions().Getindex('r`',sublist)]=0.75
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.25: mod[location][Reactions().Getindex('r`',sublist)]=0.74
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.76: mod[location][Reactions().Getindex('l`',sublist)]=0.23
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.73: mod[location][Reactions().Getindex('l`',sublist)]=0.26
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.75: mod[location][Reactions().Getindex('l`',sublist)]=0.24
        if mod[location][Reactions().Getindex(rl,sublist)] == 0.74: mod[location][Reactions().Getindex('l`',sublist)]=0.25

class DnaA:
    OriC = [[100,120],[120,121],[130,131],[150,151],[160,161],[170,171],[180,181],[410,411]]
    def __init__(self, mod, k):
        self.mod = mod
        self.k = k

    def propensity(self, state, location):
        self.p = self.k * state['DnaA']
        return self.p,location

    def execute(self, sublist, location):
        location = random.choice(self.OriC)[0]
        Reactions().ChangeStates('DnaA',location,self.mod,sublist,1)
        General().RingOpen(self.mod,sublist,location,5)
        sublist['DnaA'] -= 1
        return sublist

class primosome:
    def __init__(self,mod,sublist,k):
        self.mod = mod
        self.k = k
        self.life = 100
        if len(np.nonzero(self.mod.T[Reactions().Getindex('DnaA',sublist)])[0]) == 0:self.location = 0
        else: self.location = random.choice(np.nonzero(self.mod.T[Reactions().Getindex('DnaA',sublist)])[0])

    def propensity(self, sublist, location):
        self.p = self.k * sublist['primosome']
        return self.p, self.location

    def execute(self, sublist, location):
        location = self.location
        Reactions().ChangeStates('primosome',location,self.mod,sublist,1)
        General().RingOpen(self.mod,sublist,location,1)
        General().Walk('primosome',self.location,self.mod,sublist,1)
        self.life -= 1
        if self.life == 0 :
            sublist['primosome'] -= 1
            self.life = 100
        self.location += 1
        return sublist

class DNA_polymerase_III_holoenzyme:
    def __init__(self,mod,sublist,k):
        self.mod = mod
        self.k = k
        self.life = 90
        if len(np.nonzero(self.mod.T[Reactions().Getindex('ds',sublist)])[0]) == 0:self.location = 0
        else: self.location = random.choice(np.nonzero(self.mod.T[Reactions().Getindex('ds',sublist)])[0])

    def propensity(self, sublist, location):
        self.p = self.k * sublist['DNA_polymerase_III_holoenzyme'] * self.mod[self.location][Reactions().Getindex('ds',sublist)]
        return self.p, self.location

    def execute(self, sublist, location):
        location = self.location
        Reactions().ChangeStates('DNA_polymerase_III_holoenzyme',location,self.mod,sublist,1)
        General().Polymeraization(location,self.mod,sublist,'r')
        Reactions().ChangeStates('ds',location,self.mod,sublist,0)
        General().Walk('DNA_polymerase_III_holoenzyme',self.location,self.mod,sublist,1)
        self.life -= 1
        if self.life == 0 :
            sublist['DNA_polymerase_III_holoenzyme'] -= 1
            self.life = 90
        self.location += 1
        return sublist
