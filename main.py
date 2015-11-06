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
from functions import Compose
from functions import Decompose

#import genes
from genes import argP
from genes import scpD
from genes import diaA
from genes import dinB
from genes import dnaA
from genes import dnaB
from genes import dnaC
from genes import dnaE
from genes import dnaG
from genes import dnaN
from genes import dnaQ
from genes import dnaT
from genes import dnaX
from genes import gspB
from genes import gyrA
from genes import gyrB
from genes import hda
from genes import holA
from genes import holB
from genes import holC
from genes import holD
from genes import holE
from genes import polA
from genes import polB
from genes import priA
from genes import priB
from genes import priC
from genes import rarA
from genes import rdgC
from genes import recF
from genes import recG
from genes import recQ
from genes import rep
from genes import rnhA
from genes import rob
from genes import sbmC
from genes import seqA
from genes import topA
from genes import topB
from genes import tus
from genes import umuC
from genes import umuD
from genes import uvrD

#Sequence data
seq, mod, time, SubList = Reactions().Readseq('sequence.txt')

#Substance generate
"""
SubList.append(Reactions().Monomer('DnaA',60000))
SubList.append(Reactions().Monomer('DnaB',60000))
SubList.append(Reactions().Monomer('DnaC',60000))
SubList.append(Reactions().Monomer('DnaG',20000))
SubList.append(Reactions().Monomer('RNaseH',500))
SubList.append(Reactions().Monomer('SSB',30000))
SubList.append(Reactions().Monomer('Topo1',0))
SubList.append(Reactions().Monomer('SDC',200000))
SubList.append(Reactions().Monomer('DNApol1',0))
SubList.append(Reactions().Monomer('DNApol3',100))
SubList.append(Reactions().Monomer('DNApol3holoenzyme',100))
SubList.append(Reactions().Monomer('OriC9',5))
SubList.append(Reactions().Monomer('OriC13',3))
SubList.append(Reactions().Monomer('ATP',10000))
SubList.append(Reactions().Complex('DnaB','DnaC',0))
SubList.append(Reactions().Complex('DnaA','DnaB/DnaC',0))
"""

SubList.append(Reactions().Monomer('YciV',0))#RNA/ssDNA exonuclease 5'->3'specific(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6634-MONOMER)
SubList.append(Reactions().Monomer('YdaV',0))#Rac prophage; predicted DNA replication protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6684-MONOMER)
SubList.append(Reactions().Monomer('YcdX',0))#zinc-binding phosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6540-MONOMER)
SubList.append(Reactions().Monomer('CrfC',0))#Clamp-binding sister replication fork colocalization protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11210-MONOMER)
SubList.append(Reactions().Monomer('RarA',0))#recombination factor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12690-MONOMER)
SubList.append(Reactions().Monomer('Hda',0))#regulator of DnaA that prevents premature reinitiation of DNA replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G7313-MONOMER)
SubList.append(Reactions().Monomer('DiaA',0))#DnaA initiator-associating factor for replication initiation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=YRAO-MONOMER)
SubList.append(Reactions().Monomer('HolB',0))#DNA polymerase III, delta prime subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11500-MONOMER)
SubList.append(Reactions().Monomer('HolA',0))#DNA polymerase III, delta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11412-MONOMER)
SubList.append(Reactions().Monomer('Tus',0))#DNA-binding protein; inhibition of replication at Ter sites(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11038-MONOMER)
SubList.append(Reactions().Monomer('DnaC',60000))#chromosome replication; initiation and chain elongation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10237-MONOMER)
SubList.append(Reactions().Monomer('Dam',0))#DNA adenine methyltransferase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10204-MONOMER)
SubList.append(Reactions().Monomer('DnaG',0))#DNA primase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10239-MONOMER)
SubList.append(Reactions().Monomer('DnaN',0))#DNA polymerase III, beta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10242-MONOMER)
SubList.append(Reactions().Monomer('DnaT',0))#primosomal protein DnaT(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10244-MONOMER)
SubList.append(Reactions().Monomer('Rep',0))#Rep helicase, a single-stranded DNA dependent ATPase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10837-MONOMER)
SubList.append(Reactions().Monomer('PriB',0))#primosomal replication protein N(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10764-MONOMER)
SubList.append(Reactions().Monomer('DnaQ',0))#DNA polymerase III, epsilon subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10243-MONOMER)
SubList.append(Reactions().Monomer('DnaA',60000))#chromosomal replication initiator protein DnaA; DNA-binding transcriptional dual regulator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD03831)
SubList.append(Reactions().Monomer('MukE',0))#protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11252-MONOMER)
SubList.append(Reactions().Monomer('MukB',0))#cell division protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10618-MONOMER)
SubList.append(Reactions().Monomer('MukF',0))#MukF dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12165-MONOMER)
SubList.append(Reactions().Monomer('NrdD',0))#ribonucleoside-triphosphate reductase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=RIBONUCLEOSIDE-TRIP-REDUCT-MONOMER)
SubList.append(Reactions().Monomer('RecQ',0))#ATP-dependent DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10833-MONOMER)
SubList.append(Reactions().Monomer('DinB',0))#DNA polymerase IV (Y-family DNA polymerase; translesion DNA synthesis)(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6115-MONOMER)
SubList.append(Reactions().Monomer('NrdB',0))#ribonucleoside diphosphate reductase 1, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDB-MONOMER)
SubList.append(Reactions().Monomer('NrdE',0))#ribonucleoside-diphosphate reductase 2, alpha subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDE-MONOMER)
SubList.append(Reactions().Monomer('NrdF',0))#ribonucleoside-diphosphate reductase 2, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDF-MONOMER)
SubList.append(Reactions().Monomer('HolC',0))#DNA polymerase III, chi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11413-MONOMER)
SubList.append(Reactions().Monomer('HolD',0))#DNA polymerase III, psi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11414-MONOMER)
SubList.append(Reactions().Monomer('LigB',0))#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11334-MONOMER)
SubList.append(Reactions().Monomer('PriA',0))#primosome factor N'(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10763-MONOMER)
SubList.append(Reactions().Monomer('LigA',0))#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10534-MONOMER)
SubList.append(Reactions().Monomer('SbcC',0))#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
SubList.append(Reactions().Monomer('DnaE',0))#DNA polymerase III, alpha subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
SubList.append(Reactions().Monomer('Ssb',0))#ssDNA-binding protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10976-MONOMER)
SubList.append(Reactions().Monomer('SbcD',0))#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11094-MONOMER)
SubList.append(Reactions().Monomer('DnaB',60000))#replicative DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10236-MONOMER)
SubList.append(Reactions().Monomer('HolE',0))#DNA polymerase III, theta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11505-MONOMER)
SubList.append(Reactions().Monomer('RecF',0))#ssDNA and dsDNA binding, ATP binding(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10828-MONOMER)
SubList.append(Reactions().Monomer('LexA',0))#LexA DNA-binding transcriptional repressor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD00205)
SubList.append(Reactions().Monomer('DnaK',0))#chaperone protein DnaK(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10241-MONOMER)
SubList.append(Reactions().Monomer('DnaJ',0))#chaperone protein DnaJ(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10240-MONOMER)
SubList.append(Reactions().Monomer('MutT',0))#8-oxo-dGTP diphosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10626-MONOMER)
SubList.append(Reactions().Monomer('DnaX',0))#DNA polymerase III, tau subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10245-MONOMER)
SubList.append(Reactions().Monomer('GspB',0))#calcium-binding protein required for initiation of chromosome replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11263-MONOMER)
SubList.append(Reactions().Monomer('UvrD',0))#ssDNA translocase and dsDNA helicase - DNA helicase II(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11064-MONOMER)
SubList.append(Reactions().Monomer('PolA',0))#DNA polymerase I, 5'-->3'polymerase, 5'-->3'and3'-->5'exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10746-MONOMER)
SubList.append(Reactions().Monomer('NrdA',0))#ribonucleoside diphosphate reductase 1, alphaatom fun subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDA-MONOMER)
SubList.append(Reactions().Complex('DnaB','DnaC',0))
SubList.append(Reactions().Complex('DnaA','DnaB/DnaC',0))

"""
#Gillespie Test
logt, logd, t, tend = Showdata().logger(time, SubList, 0, 1)
events = [Compose('DnaB','DnaC',0.1),Decompose('DnaB/DnaC',2),Compose('DnaA','DnaB/DnaC',0.1),Decompose('DnaA/DnaB/DnaC',2)]
Simulation().run(t, tend, SubList, events, logt, logd)
Showdata().figure("DnaA", logt, logd, SubList)
Showdata().figure("DnaB", logt, logd, SubList)
Showdata().figure("DnaC", logt, logd, SubList)
Showdata().figure("DnaB/DnaC", logt, logd, SubList)
Showdata().figure("DnaA/DnaB/DnaC", logt, logd, SubList)
Showdata().save("process.png")

#Polymerization Process Test
r = [0]
i = 0
time = []
er = []
process = len(seq)
for location in range(process):
    dnaB().execute(location,mod,SubList)
for location in range(process):
    holB().execute(location,mod,SubList,r)
    time.append(i)
    er.append(r[0])
    i += 1

#Reading Chain5'-3'

#Ragging Chain5'-3'

#Save result
Simulation().Wcplot(time,er)
Simulation().Wcwrite(mod)
Simulation().Makedata()
"""
