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

#import functions
from functions import Reactions
from functions import Showdata
from functions import Simulation
from functions import Compose
from functions import Decompose

class argP:#ArgP DNA-binding transcriptional activator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10490)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class cspD:#DNA replication inhibitor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11111)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class diaA:#DnaA initiator-associating factor for replication initiation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG12780)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dinB:#DNA polymerase IV (Y-family DNA polymerase; translesion DNA synthesis)(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=G6115)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaA:#chromosomal replication initiator protein DnaA; DNA-binding transcriptional dual regulator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10235)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaB:#replicative DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10236)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        mseq[location][0] = 1
        state[Reactions().Getindex('DnaB',state)][1] -= 1
        return location, mseq, state

class dnaC:#chromosome replication; initiation and chain elongation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10237)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaE:#DNA polymerase III, alpha subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10238)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaG:#DNA primase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10239)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaN:#DNA polymerase III, beta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10242)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaQ:#DNA polymerase III, epsilon subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10243)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaT:#primosomal protein DnaT(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10244)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class dnaX:#DNA polymerase III, gamma subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10245)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class gspB:#calcium-binding protein required for initiation of chromosome replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11263)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class gyrA:#DNA gyrase, subunit A(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10423)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class gyrB:#DNA gyrase, subunit B(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10424)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class hda:#regulator of DnaA that prevents premature reinitiation of DNA replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=G7313)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class holA:#DNA polymerase III, delta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11412)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class holB:#DNA polymerase III, delta prime subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11500)
    def __init__(self):
        pass

    def propensity(self, mseq, state, k):
        return

    def execute(self, location, mseq, state, r):
        if mseq[location][0] == 1:
            mseq[location][9] = 1
            state[9][1] -= 1
            error = 0.001#ERROR RATE(%)
            rate = rand()*100
            if mseq[location][10]=='a':
                if rate >= error:mseq[location][11]='t'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='t':
                if rate >= error:mseq[location][11]='a'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='g':
                if rate >= error:mseq[location][11]='c'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            if mseq[location][10]=='c':
                if rate >= error:mseq[location][11]='g'
                else:
                    mseq[location][11]='miss'
                    r[0] += 1
            mseq[location][0] = 0
            location += 1
        return location, mseq, state, r

class holC:#DNA polymerase III, chi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11413)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class holD:#DNA polymerase III, psi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11414)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class holE:#DNA polymerase III, theta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11505)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class mioC:#http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11199(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11199)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class parC:#dimer of topoisomerase IV subunit A(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10686)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class parE:#topoisomerase IV subunit B(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10687)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class polA:#DNA polymerase I, 5'-->3' polymerase, 5'-->3' and 3'-->5' exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10746)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class polB:#DNA polymerase II(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10747)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class priA:#primosome factor N'(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10763)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class priB:#primosomal replication protein N(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10764)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class priC:#primosomal replication protein N''(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10765)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class rarA:#recombination factor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG12690)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class rdgC:#nucleoid-associated protein RdgC(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG12158)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class recF:#ssDNA and dsDNA binding, ATP binding(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10828)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class recG:#RecG DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10829)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class recQ:#ATP-dependent DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10833)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class rep:#Rep helicase, a single-stranded DNA dependent ATPase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10837)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class rnhA:#RNase HI, degrades RNA of DNA-RNA hybrids, participates in DNA replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG10860)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class rob:#Rob DNA-binding transcriptional activator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11366)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class sbmC:#DNA gyrase inhibitor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11892)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class seqA:#SeqA, negative modulator of initiation of replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG12197)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class topA:#DNA topoisomerase I(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11013)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class topB:#DNA topoisomerase III(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11014)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class tus:#DNA-binding protein; inhibition of replication at Ter sites(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11038)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class umuC:#SOS mutagenesis and repair(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11056)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class umuD:#SOS mutagenesis; error-prone repair; processed to UmuD'; forms complex with UmuC(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11057)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return

class uvrD:#ssDNA translocase and dsDNA helicase - DNA helicase II(http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=EG11064)
    def __init__(self):
        pass

    def propensity(self, state, k):
        return

    def execute(self, location, mseq, state):
        return
