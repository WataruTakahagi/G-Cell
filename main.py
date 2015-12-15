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
import cProfile

#import G-Cell module
from functions import *
from proteins import *

#setup
time, location, SubList, events = Reactions().setup()
target = Reactions().Target('mutationseq.txt')

#Generate Monomer
Reactions().Monomer('YciV',600,SubList,target)#RNA/ssDNA exonuclease 5'->3'specific(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6634-MONOMER)
Reactions().Monomer('YdaV',600,SubList,target)#Rac prophage; predicted DNA replication protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6684-MONOMER)
Reactions().Monomer('YcdX',600,SubList,target)#zinc-binding phosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6540-MONOMER)
Reactions().Monomer('CrfC',600,SubList,target)#Clamp-binding sister replication fork colocalization protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11210-MONOMER)
Reactions().Monomer('RarA',600,SubList,target)#recombination factor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12690-MONOMER)
Reactions().Monomer('Hda',600,SubList,target)#regulator of DnaA that prevents premature reinitiation of DNA replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G7313-MONOMER)
Reactions().Monomer('DiaA',600,SubList,target)#DnaA initiator-associating factor for replication initiation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=YRAO-MONOMER)
Reactions().Monomer('HolB',600,SubList,target)#DNA polymerase III, delta prime subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11500-MONOMER)
Reactions().Monomer('HolA',600,SubList,target)#DNA polymerase III, delta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11412-MONOMER)
Reactions().Monomer('Tus',600,SubList,target)#DNA-binding protein; inhibition of replication at Ter sites(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11038-MONOMER)
Reactions().Monomer('DnaC',600,SubList,target)#chromosome replication; initiation and chain elongation(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10237-MONOMER)
Reactions().Monomer('Dam',600,SubList,target)#DNA adenine methyltransferase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10204-MONOMER)
Reactions().Monomer('DnaG',600,SubList,target)#DNA primase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10239-MONOMER)
Reactions().Monomer('DnaN',600,SubList,target)#DNA polymerase III, beta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10242-MONOMER)
Reactions().Monomer('DnaT',600,SubList,target)#primosomal protein DnaT(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10244-MONOMER)
Reactions().Monomer('Rep',600,SubList,target)#Rep helicase, a single-stranded DNA dependent ATPase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10837-MONOMER)
Reactions().Monomer('PriB',600,SubList,target)#primosomal replication protein N(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10764-MONOMER)
Reactions().Monomer('DnaQ',600,SubList,target)#DNA polymerase III, epsilon subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10243-MONOMER)
Reactions().Monomer('DnaA',600,SubList,target)#chromosomal replication initiator protein DnaA; DNA-binding transcriptional dual regulator(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD03831)
Reactions().Monomer('MukE',600,SubList,target)#protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11252-MONOMER)
Reactions().Monomer('MukB',600,SubList,target)#cell division protein involved in chromosome partitioning(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10618-MONOMER)
Reactions().Monomer('MukF',600,SubList,target)#MukF dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG12165-MONOMER)
Reactions().Monomer('NrdD',600,SubList,target)#ribonucleoside-triphosphate reductase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=RIBONUCLEOSIDE-TRIP-REDUCT-MONOMER)
Reactions().Monomer('RecQ',600,SubList,target)#ATP-dependent DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10833-MONOMER)
Reactions().Monomer('DinB',600,SubList,target)#DNA polymerase IV (Y-family DNA polymerase; translesion DNA synthesis)(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=G6115-MONOMER)
Reactions().Monomer('NrdB',600,SubList,target)#ribonucleoside diphosphate reductase 1, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDB-MONOMER)
Reactions().Monomer('NrdE',600,SubList,target)#ribonucleoside-diphosphate reductase 2, alpha subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDE-MONOMER)
Reactions().Monomer('NrdF',600,SubList,target)#ribonucleoside-diphosphate reductase 2, beta subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDF-MONOMER)
Reactions().Monomer('HolC',600,SubList,target)#DNA polymerase III, chi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11413-MONOMER)
Reactions().Monomer('HolD',600,SubList,target)#DNA polymerase III, psi subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11414-MONOMER)
Reactions().Monomer('LigB',600,SubList,target)#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11334-MONOMER)
Reactions().Monomer('PriA',600,SubList,target)#primosome factor N'(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10763-MONOMER)
Reactions().Monomer('LigA',600,SubList,target)#DNA ligase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10534-MONOMER)
Reactions().Monomer('SbcC',600,SubList,target)#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
Reactions().Monomer('DnaE',600,SubList,target)#DNA polymerase III, alpha subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10927-MONOMER)
Reactions().Monomer('Ssb',600,SubList,target)#ssDNA-binding protein(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10976-MONOMER)
Reactions().Monomer('SbcD',600,SubList,target)#ATP-dependent dsDNA exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11094-MONOMER)
Reactions().Monomer('DnaB',600,SubList,target)#replicative DNA helicase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10236-MONOMER)
Reactions().Monomer('HolE',600,SubList,target)#DNA polymerase III, theta subunit(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11505-MONOMER)
Reactions().Monomer('RecF',600,SubList,target)#ssDNA and dsDNA binding, ATP binding(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10828-MONOMER)
Reactions().Monomer('RecR',600,SubList,target)#recombination and repair(http://ecocyc.org/ECOLI/NEW-IMAGE?type=POLYPEPTIDE&object=EG10834-MONOMER)
Reactions().Monomer('RecO',600,SubList,target)#protein interacts with RecR and possibly RecF proteins(http://ecocyc.org/gene?orgid=ECOLI&id=EG10832-MONOMER)
Reactions().Monomer('LexA',600,SubList,target)#LexA DNA-binding transcriptional repressor(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=PD00205)
Reactions().Monomer('DnaK',600,SubList,target)#chaperone protein DnaK(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10241-MONOMER)
Reactions().Monomer('DnaJ',600,SubList,target)#chaperone protein DnaJ(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10240-MONOMER)
Reactions().Monomer('MutT',600,SubList,target)#8-oxo-dGTP diphosphatase(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10626-MONOMER)
Reactions().Monomer('DnaX',600,SubList,target)#DNA polymerase III, tau subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10245-MONOMER)
Reactions().Monomer('GspB',600,SubList,target)#calcium-binding protein required for initiation of chromosome replication(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11263-MONOMER)
Reactions().Monomer('UvrD',600,SubList,target)#ssDNA translocase and dsDNA helicase - DNA helicase II(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG11064-MONOMER)
Reactions().Monomer('PolA',600,SubList,target)#DNA polymerase I, 5'-->3'polymerase, 5'-->3'and3'-->5'exonuclease(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10746-MONOMER)
Reactions().Monomer('NrdA',600,SubList,target)#ribonucleoside diphosphate reductase 1, alpha atom fun subunit dimer(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=NRDA-MONOMER)
Reactions().Monomer('PriC',600,SubList,target)#primosomal replication protein N''(http://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=EG10765-MONOMER)
Showdata().csv(SubList)
Reactions().Complex('zink_binding_phosphatase',0,SubList)
Reactions().Complex('DnaA_initiator_associating_factor',0,SubList)
Reactions().Complex('DNA_polymerase_III_holoenzyme',0,SubList)
Reactions().Complex('DNA_polymerase_III_core_enzyme',0,SubList)
Reactions().Complex('DNA_polymerase_III_preinitiation_complex',0,SubList)
Reactions().Complex('DNA_polymerase_III_beta_subunit',0,SubList)
Reactions().Complex('DNA_polymerase_III_tau_subunit_dimer',0,SubList)
Reactions().Complex('DNA_polymerase_III_psi_chi_subunit',0,SubList)
Reactions().Complex('replicative_DNA_helicase',0,SubList)
Reactions().Complex('DNA_primase',0,SubList)
Reactions().Complex('primosomal_protein_DnaT',0,SubList)
Reactions().Complex('primosomal_replication_protein_N',0,SubList)
Reactions().Complex('primosome',0,SubList)
Reactions().Complex('Rep_helicase',0,SubList)
Reactions().Complex('cell_division_protein',0,SubList)
Reactions().Complex('MukEF_complex',0,SubList)
Reactions().Complex('MukF_dimer',0,SubList)
Reactions().Complex('bacterial_condensin_MukBEF',0,SubList)
Reactions().Complex('ribonucleoside_triphosphate_reductase',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_1',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer',0,SubList)
Reactions().Complex('ribonucleoside_diphosphate_reductase_2',0,SubList)
Reactions().Complex('SbcCD_ATP_dependent_dsDNA_exonuclease',0,SubList)
Reactions().Complex('ssDNA_binding_protein',0,SubList)
Reactions().Complex('RecFOR_complex',0,SubList)
Reactions().Complex('LexA_DNA_binding_transcriptional_repressor',0,SubList)
Reactions().Complex('chaperone_protein_DnaJ',0,SubList)
Reactions().Complex('ssDNA_translocase_and_dsDNA_helicase',0,SubList)

#Sequence data, Make logger
seq, mod = Reactions().Readseq(target,SubList)
logt, logd, t, tend = Showdata().logger(time, SubList, 0, 0.01)

#Events setting
Reactions().Events(Compose('zink_binding_phosphatase',['YcdX'],[3],1.0e-2),events)#zink-binding phosphatase
Reactions().Events(Decompose('zink_binding_phosphatase',['YcdX'],[3],1),events)#zink-binding phosphatase
Reactions().Events(Compose('DnaA_initiator_associating_factor',['DiaA'],[4],1.0e-4),events)#DnaA initiator-associating factor
Reactions().Events(Decompose('DnaA_initiator_associating_factor',['DiaA'],[4],1),events)#DnaA initiator-associating factor
Reactions().Events(Compose('DNA_polymerase_III_core_enzyme',['DnaE','DnaQ','HolE'],[1,1,1],1.0e-2),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('DNA_polymerase_III_preinitiation_complex',['DnaX','HolB','HolA'],[3,1,1],1.0e-9),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('DNA_polymerase_III_beta_subunit',['DnaN'],[2],5.0e-1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('DNA_polymerase_III_tau_subunit_dimer',['DnaX'],[2],5.0e-1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('DNA_polymerase_III_psi_chi_subunit',['HolC','HolD'],[1,1],5.0e-1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('DNA_polymerase_III_holoenzyme',['DNA_polymerase_III_core_enzyme','DNA_polymerase_III_preinitiation_complex','DNA_polymerase_III_beta_subunit','DNA_polymerase_III_tau_subunit_dimer','DNA_polymerase_III_psi_chi_subunit'],[3,1,2,1,4],1.0e-15),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_core_enzyme',['DnaE','DnaQ','HolE'],[1,1,1],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_preinitiation_complex',['DnaX','HolB','HolA'],[3,1,1],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_beta_subunit',['DnaN'],[2],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_tau_subunit_dimer',['DnaX'],[2],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_psi_chi_subunit',['HolC','HolD'],[1,1],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Decompose('DNA_polymerase_III_holoenzyme',['DNA_polymerase_III_core_enzyme','DNA_polymerase_III_preinitiation_complex','DNA_polymerase_III_beta_subunit','DNA_polymerase_III_tau_subunit_dimer','DNA_polymerase_III_psi_chi_subunit'],[3,1,2,1,4],1),events)#DNA_polymerase_III_holoenzyme
Reactions().Events(Compose('replicative_DNA_helicase',['DnaB'],[6],1.0e-9),events)#primosome
Reactions().Events(Compose('primosomal_protein_DnaT',['DnaT'],[3],1.0e-2),events)#primosome
Reactions().Events(Compose('primosomal_replication_protein_N',['PriB'],[2],5.0e-1),events)#primosome
Reactions().Events(Compose('primosome',['replicative_DNA_helicase','primosomal_protein_DnaT','primosomal_replication_protein_N','PriA','PriC','DnaG'],[1,1,1,1,1,1],1.0e-8),events)#primosome
Reactions().Events(Decompose('replicative_DNA_helicase',['DnaB'],[6],1),events)#primosome
Reactions().Events(Decompose('primosomal_protein_DnaT',['DnaT'],[3],1),events)#primosome
Reactions().Events(Decompose('primosomal_replication_protein_N',['PriB'],[2],1),events)#primosome
Reactions().Events(Decompose('primosome',['replicative_DNA_helicase','primosomal_protein_DnaT','primosomal_replication_protein_N','PriA','PriC','DnaG'],[1,1,1,1,1,1],1),events)#primosome
Reactions().Events(Compose('Rep_helicase',['Rep'],[2],5.0e-1),events)#Rep_helicase
Reactions().Events(Decompose('Rep_helicase',['Rep'],[2],1),events)#Rep_helicase
Reactions().Events(Compose('cell_division_protein',['MukB'],[2],5.0e-1),events)#bacterial_condensin_MukBEF
Reactions().Events(Compose('MukEF_complex',['MukE','MukF_dimer'],[4,1],1.0e-9),events)#bacterial_condensin_MukBEF
Reactions().Events(Compose('MukF_dimer',['MukF'],[2],5.0e-1),events)#bacterial_condensin_MukBEF
Reactions().Events(Compose('bacterial_condensin_MukBEF',['cell_division_protein','MukEF_complex'],[1,1],5.0e-1),events)#bacterial_condensin_MukBEF
Reactions().Events(Decompose('cell_division_protein',['MukB'],[2],1),events)#bacterial_condensin_MukBEF
Reactions().Events(Decompose('MukEF_complex',['MukE','MukF_dimer'],[4,1],1),events)#bacterial_condensin_MukBEF
Reactions().Events(Decompose('MukF_dimer',['MukF'],[2],1),events)#bacterial_condensin_MukBEF
Reactions().Events(Decompose('bacterial_condensin_MukBEF',['cell_division_protein','MukEF_complex'],[1,1],1),events)#bacterial_condensin_MukBEF
Reactions().Events(Compose('ribonucleoside_triphosphate_reductase',['NrdD'],[2],5.0e-1),events)#ribonucleoside-triphosphate reductase
Reactions().Events(Decompose('ribonucleoside_triphosphate_reductase',['NrdD'],[2],1),events)#ribonucleoside-triphosphate reductase
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer',['NrdA'],[2],5.0e-1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer',['NrdB'],[2],5.0e-1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_1',['ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer'],[1,1],5.0e-1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer',['NrdA'],[2],1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer',['NrdB'],[2],1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_1',['ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer'],[1,1],1),events)#ribonucleoside diphosphate reductase 1
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer',['NrdE'],[2],5.0e-1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer',['NrdF'],[2],5.0e-1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Compose('ribonucleoside_diphosphate_reductase_2',['ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer'],[1,1],5.0e-1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer',['NrdE'],[2],1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer',['NrdF'],[2],1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Decompose('ribonucleoside_diphosphate_reductase_2',['ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer'],[1,1],1),events)#ribonucleoside-diphosphate reductase 2
Reactions().Events(Compose('SbcCD_ATP_dependent_dsDNA_exonuclease',['SbcC','SbcD'],[1,1],5.0e-1),events)#SbcCD ATP-dependent dsDNA exonuclease
Reactions().Events(Decompose('SbcCD_ATP_dependent_dsDNA_exonuclease',['SbcC','SbcD'],[1,1],1),events)#SbcCD ATP-dependent dsDNA exonuclease
Reactions().Events(Compose('ssDNA_binding_protein',['Ssb'],[4],1.0e-4),events)#ssDNA-binding protein
Reactions().Events(Decompose('ssDNA_binding_protein',['Ssb'],[4],1),events)#ssDNA-binding protein
Reactions().Events(Compose('RecFOR_complex',['RecF','RecR','RecO'],[1,1,1],1.0e-2),events)#RecFOR complex
Reactions().Events(Decompose('RecFOR_complex',['RecF','RecR','RecO'],[1,1,1],1),events)#RecFOR complex
Reactions().Events(Compose('LexA_DNA_binding_transcriptional_repressor',['LexA'],[2],5.0e-1),events)#LexA DNA-binding transcriptional repressor
Reactions().Events(Decompose('LexA_DNA_binding_transcriptional_repressor',['LexA'],[2],1),events)#LexA DNA-binding transcriptional repressor
Reactions().Events(Compose('chaperone_protein_DnaJ',['DnaJ'],[2],5.0e-1),events)#chaperone protein DnaJ
Reactions().Events(Decompose('chaperone_protein_DnaJ',['DnaJ'],[2],1),events)#chaperone protein DnaJ
Reactions().Events(Compose('ssDNA_translocase_and_dsDNA_helicase',['UvrD'],[2],5.0e-1),events)#ssDNA translocase and dsDNA helicase
Reactions().Events(Decompose('ssDNA_translocase_and_dsDNA_helicase',['UvrD'],[2],1),events)#ssDNA translocase and dsDNA helicase

#for i in range(20):
#    DnaA(mod,10).execute(SubList,location)

#simulation
#Simulation().Run(t, tend, SubList, events, logt, logd, mod, location)
cProfile.run('Simulation().Run(t, tend, SubList, events, logt, logd, mod, location)', 'profile')

#showdata, make .png
Showdata().png(['YcdX','zink_binding_phosphatase'],logt, logd, SubList,'default')
Showdata().png(['DiaA','DnaA_initiator_associating_factor'],logt, logd, SubList,'default')
Showdata().png(['DnaE','DnaQ','HolE','DNA_polymerase_III_core_enzyme'], logt, logd, SubList,'default')
Showdata().png(['DnaX','HolB','HolA','DNA_polymerase_III_preinitiation_complex'], logt, logd, SubList,'default')
Showdata().png(['DnaN','DNA_polymerase_III_beta_subunit'], logt, logd, SubList,'default')
Showdata().png(['DnaX','DNA_polymerase_III_tau_subunit_dimer'], logt, logd, SubList,'default')
Showdata().png(['HolC','HolD','DNA_polymerase_III_psi_chi_subunit'], logt, logd, SubList,'default')
Showdata().png(['DNA_polymerase_III_core_enzyme','DNA_polymerase_III_preinitiation_complex','DNA_polymerase_III_beta_subunit','DNA_polymerase_III_tau_subunit_dimer','DNA_polymerase_III_psi_chi_subunit','DNA_polymerase_III_holoenzyme'], logt, logd, SubList,'default')
Showdata().png(['DnaB','replicative_DNA_helicase'], logt, logd, SubList,'default')
Showdata().png(['DnaT','primosomal_protein_DnaT'], logt, logd, SubList,'default')
Showdata().png(['PriB','primosomal_replication_protein_N'], logt, logd, SubList,'default')
Showdata().png(['replicative_DNA_helicase','primosomal_protein_DnaT','primosomal_replication_protein_N','PriA','PriC','DnaG','primosome'], logt, logd, SubList,'default')
Showdata().png(['Rep','Rep_helicase'], logt, logd, SubList,'default')
Showdata().png(['MukB','cell_division_protein'], logt, logd, SubList,'default')
Showdata().png(['MukF','MukF_dimer'], logt, logd, SubList,'default')
Showdata().png(['MukE','MukF_dimer','MukEF_complex'], logt, logd, SubList,'default')
Showdata().png(['cell_division_protein','MukEF_complex','bacterial_condensin_MukBEF'], logt, logd, SubList,'default')
Showdata().png(['NrdD','ribonucleoside_triphosphate_reductase'], logt, logd, SubList,'default')
Showdata().png(['NrdA','ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer'], logt, logd, SubList,'default')
Showdata().png(['NrdB','ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer'], logt, logd, SubList,'default')
Showdata().png(['ribonucleoside_diphosphate_reductase_1_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_1_beta_subunit_dimer','ribonucleoside_diphosphate_reductase_1'], logt, logd, SubList,'default')
Showdata().png(['NrdE','ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer'], logt, logd, SubList,'default')
Showdata().png(['NrdF','ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer'], logt, logd, SubList,'default')
Showdata().png(['ribonucleoside_diphosphate_reductase_2_alpha_subunit_dimer','ribonucleoside_diphosphate_reductase_2_beta_subunit_dimer','ribonucleoside_diphosphate_reductase_2'], logt, logd, SubList,'default')
Showdata().png(['SbcC','SbcD','SbcCD_ATP_dependent_dsDNA_exonuclease'], logt, logd, SubList,'default')
Showdata().png(['Ssb','ssDNA_binding_protein'], logt, logd, SubList,'default')
Showdata().png(['RecF','RecR','RecO','RecFOR_complex'], logt, logd, SubList,'default')
Showdata().png(['LexA','LexA_DNA_binding_transcriptional_repressor'], logt, logd, SubList,'default')
Showdata().png(['DnaJ','chaperone_protein_DnaJ'], logt, logd, SubList,'default')
Showdata().png(['UvrD','ssDNA_translocase_and_dsDNA_helicase'], logt, logd, SubList,'default')
#Showdata().csv(SubList,'label.csv')
#Showdata().csv(logt,'time.csv')
#Showdata().csv(logd,'data.csv')

#Finalize
Showdata().csv(SubList,'summary.csv')
Simulation().Save(mod,SubList)
Simulation().Makedata('default')

#Reactions().Increase(0,'DnaB',mod,SubList)
#Reactions().Increase(0,'DnaB',mod,SubList)
#Reactions().Decrease(0,'DnaA',mod,SubList)
