#!/usr/bin/env python
# encoding: utf-8

"""
Object-oriented abstractions for GAUSSIAN input files
"""

from __future__ import print_function, division
# Python stdlib
from string import ascii_letters
import sys
import os
from collections import namedtuple, defaultdict
from datetime import datetime

if sys.version_info.major == 3:
    basestring = str

# QM_BASIS_SETS_EXT = ['d,f-Diffuse (+)', 'p,d,f-Diffuse (++)', 'd,f-Polarization (*)', 'p,d,f-Polarization (**)']
MM_FORCEFIELDS = {
    'General': ['Amber', 'GAFF', 'Dreiding', 'UFF', 'MM3'],
    'Water': ['TIP3P']
}
MEM_UNITS = ('KB', 'MB', 'GB', 'TB', 'KW', 'MW', 'GW', 'TW')
JOB_TYPES = ('SP', 'Opt', 'IRC', 'IRCMax', 'Scan', 'Freq', 'Polar', 'ADMP',
             'Force', 'Stable', 'Volume')
JOB_OPTIONS = {
    'SP': (),
    'Opt': ('Min', 'TS'),
    'IRC': (),
    'IRCMax': (),
    'Scan': (),
    'Freq': ('Raman', 'NRaman', 'NNRaman', 'NoRaman', 'VCD', 'ROA'),
    'Polar': (),
    'ADMP': (),
    'Force': (),
    'Stable': (),
    'Volume': ()
}
QM_METHODS = ('AM1', 'PM3', 'PM3MM', 'PM6', 'PDDG', 'HF', 'DFT', 'CASSCF', 'MP2',
              'MP3', 'MP4(SDQ)', 'MP4(SDTQ)', 'MP5', 'QCISD', 'CCD', 'CCSD',
              'QCISD(T)', 'QCISD(TQ)', 'BD', 'EPT', 'CBS', 'W1', 'CIS', 'TD',
              'EOM', 'ZINDO', 'DFTB', 'CI', 'GVB', 'G1', 'G2', 'G2MP2', 'G3',
              'G3MP2', 'G3B3', 'G3MP2B3', 'G4', 'G4MP2')
QM_BASIS_SETS = ('STO-3G', '3-21G', '6-21G', '4-31G', '6-31G', "6-31G(d')",
                 "6-31G(d',p')", '6-311G', 'D95V', 'D95', 'SHC', 'CEP-4G',
                 'CEP-31G', 'CEP-121G', 'LanL2MB', 'LanL2DZ', 'SDD', 'SDDAll',
                 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z', 'cc-pV6Z')
QM_BASIS_SETS_EXT = ('', '+', '++', '*', '**')
QM_BASIS_SETS_WITH_ARGS = ('SDD', 'SHC')
QM_FUNCTIONALS = {
    'Pure':      ('VSXC', 'HCTH', 'HCTH93', 'HCTH147', 'HCTH407', 'tHCTH',
                  'M06L', 'B97D', 'B97D3', 'SOGGA11', 'M11L', 'N12', 'MN12L'),
    'Hybrid':    ('B3LYP', 'B3P86', 'B3PW91', 'B1B95', 'mPW1PW91', 'mPW1LYP',
                  'mPW1PBE', 'mPW3PBE', 'B98', 'B971', 'B972', 'PBE1PBE', 'B1LYP',
                  'O3LYP', 'BHandH', 'BHandHLYP', 'BMK', 'M06', 'M06HF', 'M062X',
                  'tHCTHhyb', 'APFD', 'APF', 'SOGGA11X', 'PBEh1PBE', 'TPSSh', 'X3LYP'),
    'RS hybrid': ('HSEH1PBE', 'OHSE2PBE', 'OHSE1PBE', 'wB97XD', 'wB97', 'wB97X',
                  'LC-wPBE', 'CAM-B3LYP', 'HISSbPBE', 'M11', 'N12SX', 'MN12SX')
}
QM_FUNCTIONALS_ALL = set(f for v in QM_FUNCTIONALS.values() for f in v)

MM_ATTRIBS = ('gaffType', 'mol2type')

MM_TYPES = {
    'GAFF': {
        'C': {'element': 'C', 'MM3': ['3', '2']},
        'C1': {'element': 'C', 'MM3': ['4']},
        'C2': {'element': 'C', 'MM3': ['2', '38', '57']},
        'C3': {'element': 'C', 'MM3': ['1', '22', '56']},
        'CA': {'element': 'C', 'MM3': ['2']},
        'CP': {'element': 'C', 'MM3': ['2']},
        'CQ': {'element': 'C', 'MM3': ['2']},
        'CC': {'element': 'C', 'MM3': ['2']},
        'CD': {'element': 'C', 'MM3': ['2']},
        'CE': {'element': 'C', 'MM3': ['2']},
        'CF': {'element': 'C', 'MM3': ['2']},
        'CG': {'element': 'C', 'MM3': ['4']},
        'CH': {'element': 'C', 'MM3': ['4']},
        'CX': {'element': 'C', 'MM3': ['22']},
        'CY': {'element': 'C', 'MM3': ['56']},
        'CU': {'element': 'C', 'MM3': ['38']},
        'CV': {'element': 'C', 'MM3': ['57']},
        'CZ': {'element': 'C', 'MM3': ['2']},
        'H1': {'element': 'H', 'MM3': ['5']},
        'H2': {'element': 'H', 'MM3': ['5']},
        'H3': {'element': 'H', 'MM3': ['5']},
        'H4': {'element': 'H', 'MM3': ['5']},
        'H5': {'element': 'H', 'MM3': ['5']},
        'HA': {'element': 'H', 'MM3': ['5']},
        'HC': {'element': 'H', 'MM3': ['5', '124']},
        'HN': {'element': 'H', 'MM3': ['23', '28', '48']},
        'HO': {'element': 'H', 'MM3': ['21', '73']},
        'HP': {'element': 'H', 'MM3': ['5']},
        'HS': {'element': 'H', 'MM3': ['44']},
        'HW': {'element': 'H', 'MM3': ['21']},
        'HX': {'element': 'H', 'MM3': ['5']},
        'F': {'element': 'F', 'MM3': ['11']},
        'CL': {'element': 'Cl', 'MM3': ['12']},
        'BR': {'element': 'Br', 'MM3': ['13']},
        'I': {'element': 'I', 'MM3': ['14']},
        'N': {'element': 'N', 'MM3': ['9']},
        'N1': {'element': 'N', 'MM3': ['10']},
        'N2': {'element': 'N', 'MM3': ['37']},
        'N3': {'element': 'N', 'MM3': ['8']},
        'N4': {'element': 'N', 'MM3': ['39']},
        'NA': {'element': 'N', 'MM3': ['8', '40']},
        'NB': {'element': 'N', 'MM3': ['8']},
        'NC': {'element': 'N', 'MM3': ['37']},
        'ND': {'element': 'N', 'MM3': ['37']},
        'NE': {'element': 'N', 'MM3': ['37']},
        'NF': {'element': 'N', 'MM3': ['37']},
        'NH': {'element': 'N', 'MM3': ['8']},
        'NO': {'element': 'N', 'MM3': ['46']},
        'O': {'element': 'O', 'MM3': ['7']},
        'OH': {'element': 'O', 'MM3': ['6']},
        'OS': {'element': 'O', 'MM3': ['6', '41']},
        'OW': {'element': 'O', 'MM3': ['6']},
        'P2': {'element': 'P', 'MM3': ['60']},
        'P3': {'element': 'P', 'MM3': ['25']},
        'P4': {'element': 'P', 'MM3': ['25']},
        'P5': {'element': 'P', 'MM3': ['60']},
        'PB': {'element': 'P', 'MM3': ['60']},
        'PC': {'element': 'P', 'MM3': ['60']},
        'PD': {'element': 'P', 'MM3': ['60']},
        'PE': {'element': 'P', 'MM3': ['60']},
        'PF': {'element': 'P', 'MM3': ['60']},
        'PX': {'element': 'P', 'MM3': ['60']},
        'PY': {'element': 'P', 'MM3': ['60']},
        'S': {'element': 'S', 'MM3': ['42']},
        'S2': {'element': 'S', 'MM3': ['15']},
        'S4': {'element': 'S', 'MM3': ['15', '17']},
        'S6': {'element': 'S', 'MM3': ['18']},
        'SH': {'element': 'S', 'MM3': ['15']},
        'SS': {'element': 'S', 'MM3': ['15', '42', '104']},
        'SX': {'element': 'S', 'MM3': ['17']},
        'SY': {'element': 'S', 'MM3': ['18']}
    },
    'Amber': {
        'H': {'element': 'H'},
        'HO': {'element': 'H'},
        'HS': {'element': 'H'},
        'HC': {'element': 'H'},
        'H1': {'element': 'H'},
        'H2': {'element': 'H'},
        'H3': {'element': 'H'},
        'HP': {'element': 'H'},
        'HA': {'element': 'H'},
        'H4': {'element': 'H'},
        'H5': {'element': 'H'},
        'HW': {'element': 'H'},
        'O': {'element': 'O'},
        'O2': {'element': 'O'},
        'OW': {'element': 'O'},
        'OH': {'element': 'O'},
        'OS': {'element': 'O'},
        'CT': {'element': 'C'},
        'CA': {'element': 'C'},
        'CM': {'element': 'C'},
        'C': {'element': 'C'},
        'C*': {'element': 'C'},
        'CB': {'element': 'C'},
        'CC': {'element': 'C'},
        'CN': {'element': 'C'},
        'CK': {'element': 'C'},
        'CQ': {'element': 'C'},
        'CW': {'element': 'C'},
        'CV': {'element': 'C'},
        'CR': {'element': 'C'},
        'CY': {'element': 'C'},
        'CD': {'element': 'C'},
        'N': {'element': 'N'},
        'NA': {'element': 'N'},
        'N2': {'element': 'N'},
        'N*': {'element': 'N'},
        'NC': {'element': 'N'},
        'NB': {'element': 'N'},
        'N3': {'element': 'N'},
        'NP': {'element': 'N'},
        'NO': {'element': 'N'},
        'S': {'element': 'S'},
        'SH': {'element': 'S'},
        'P': {'element': 'P'},
        'IM': {'element': 'Cl'},
        'LI': {'element': 'Li'},
        'IP': {'element': 'Na'},
        'K': {'element': 'K'},
        'RB': {'element': 'Rb'},
        'CS': {'element': 'Cs'},
        'I': {'element': 'I'},
        'F': {'element': 'F'},
        'IB': {'element': 'Na'},
    },
    'Dreiding': {
        'H_': {'element': 'H'},
        'H_B': {'element': 'H'},
        'H_HB': {'element': 'H'},
        'B_2': {'element': 'B'},
        'B_3': {'element': 'B'},
        'C_3': {'element': 'C'},
        'C_R': {'element': 'C'},
        'C_2': {'element': 'C'},
        'C_1': {'element': 'C'},
        'N_3': {'element': 'N'},
        'N_R': {'element': 'N'},
        'N_2': {'element': 'N'},
        'N_1': {'element': 'N'},
        'O_3': {'element': 'O'},
        'O_R': {'element': 'O'},
        'O_2': {'element': 'O'},
        'O_1': {'element': 'O'},
        'F_': {'element': 'F'},
        'AL3': {'element': 'Al'},
        'SI3': {'element': 'Si'},
        'P_3': {'element': 'P'},
        'S_3': {'element': 'S'},
        'CL': {'element': 'Cl'},
        'GA3': {'element': 'Ga'},
        'GE3': {'element': 'Ge'},
        'AS3': {'element': 'As'},
        'SE3': {'element': 'Se'},
        'BR': {'element': 'Br'},
        'IN3': {'element': 'In'},
        'SN3': {'element': 'Sn'},
        'SB3': {'element': 'Sb'},
        'TE3': {'element': 'Te'},
        'I_': {'element': 'I'},
        'NA': {'element': 'Na'},
        'CA': {'element': 'Ca'},
        'FE': {'element': 'Fe'},
        'FN': {'element': 'Fe'}
    },
    'UFF': {
        'C_3': {'element': 'C', 'MM3': '1'},
        'C_22': {'element': 'C', 'MM3': '22'},
        'C_56': {'element': 'C', 'MM3': '56'},
        'C_2': {'element': 'C', 'MM3': '2'},
        'C_3C': {'element': 'C', 'MM3': '3'},
        'C_29': {'element': 'C', 'MM3': '29'},
        'C_30': {'element': 'C', 'MM3': '30'},
        'C_38': {'element': 'C', 'MM3': '38'},
        'C_57': {'element': 'C', 'MM3': '57'},
        'C_113': {'element': 'C', 'MM3': '113'},
        'C_114': {'element': 'C', 'MM3': '114'},
        'C_R': {'element': 'C', 'MM3': '2'},
        'C_50': {'element': 'C', 'MM3': '50'},
        'C_58': {'element': 'C', 'MM3': '58'},
        'C_67': {'element': 'C', 'MM3': '67'},
        'C_71': {'element': 'C', 'MM3': '71'},
        'C_1': {'element': 'C', 'MM3': '4'},
        'C_68': {'element': 'C', 'MM3': '68'},
        'C_106': {'element': 'C', 'MM3': '106'},
        'H_': {'element': 'H', 'MM3': '5'},
        'H_21': {'element': 'H', 'MM3': '21'},
        'H_23': {'element': 'H', 'MM3': '23'},
        'H_24': {'element': 'H', 'MM3': '24'},
        'H_28': {'element': 'H', 'MM3': '28'},
        'H_36': {'element': 'H', 'MM3': '36'},
        'H_44': {'element': 'H', 'MM3': '44'},
        'H_48': {'element': 'H', 'MM3': '48'},
        'H_73': {'element': 'H', 'MM3': '73'},
        'H_124': {'element': 'H', 'MM3': '124'},
        'H_B': {'element': 'H', 'MM3': '5'},
        'O_3': {'element': 'O', 'MM3': '6'},
        'O_49': {'element': 'O', 'MM3': '49'},
        'O_69': {'element': 'O', 'MM3': '69'},
        'O_70': {'element': 'O', 'MM3': '70'},
        'O_75': {'element': 'O', 'MM3': '75'},
        'O_145': {'element': 'O', 'MM3': '145'},
        'O_3_Z': {'element': 'O', 'MM3': '6'},
        'O_R': {'element': 'O', 'MM3': '6'},
        'O_41': {'element': 'O', 'MM3': '41'},
        'O_47': {'element': 'O', 'MM3': '47'},
        'O_77': {'element': 'O', 'MM3': '77'},
        'O_78': {'element': 'O', 'MM3': '78'},
        'O_79': {'element': 'O', 'MM3': '79'},
        'O_80': {'element': 'O', 'MM3': '80'},
        'O_81': {'element': 'O', 'MM3': '81'},
        'O_82': {'element': 'O', 'MM3': '82'},
        'O_148': {'element': 'O', 'MM3': '148'},
        'O_149': {'element': 'O', 'MM3': '149'},
        'O_2': {'element': 'O', 'MM3': '7'},
        'O_159': {'element': 'O', 'MM3': '159'},
        'O_1': {'element': 'O', 'MM3': '7'},
        'N_3': {'element': 'N', 'MM3': '8'},
        'N_39': {'element': 'N', 'MM3': '39'},
        'N_146': {'element': 'N', 'MM3': '146'},
        'N_150': {'element': 'N', 'MM3': '150'},
        'N_155': {'element': 'N', 'MM3': '155'},
        'N_164': {'element': 'N', 'MM3': '164'},
        'N_R': {'element': 'N', 'MM3': '37'},
        'N_40': {'element': 'N', 'MM3': '40'},
        'N_143': {'element': 'N', 'MM3': '143'},
        'N_144': {'element': 'N', 'MM3': '144'},
        'N_151': {'element': 'N', 'MM3': '151'},
        'N_2': {'element': 'N', 'MM3': '9'},
        'N_43': {'element': 'N', 'MM3': '43'},
        'N_45': {'element': 'N', 'MM3': '45'},
        'N_46': {'element': 'N', 'MM3': '46'},
        'N_72': {'element': 'N', 'MM3': '72'},
        'N_107': {'element': 'N', 'MM3': '107'},
        'N_108': {'element': 'N', 'MM3': '108'},
        'N_109': {'element': 'N', 'MM3': '109'},
        'N_110': {'element': 'N', 'MM3': '110'},
        'N_111': {'element': 'N', 'MM3': '111'},
        'N_1': {'element': 'N', 'MM3': '10'},
        'F_': {'element': 'F', 'MM3': '11'},
        'CL': {'element': 'Cl', 'MM3': '12'},
        'BR': {'element': 'Br', 'MM3': '13'},
        'I_': {'element': 'I', 'MM3': '14'},
        'S_3+2': {'element': 'S', 'MM3': '15'},
        'S_16': {'element': 'S', 'MM3': '16'},
        'S_17': {'element': 'S', 'MM3': '17'},
        'S_18': {'element': 'S', 'MM3': '18'},
        'S_104': {'element': 'S', 'MM3': '104'},
        'S_105': {'element': 'S', 'MM3': '105'},
        'S_3+4': {'element': 'S', 'MM3': '154'},
        'S_3+6': {'element': 'S', 'MM3': '154'},
        'S_R': {'element': 'S', 'MM3': '42'},
        'S_2': {'element': 'S', 'MM3': '74'},
        'P_3+3': {'element': 'P', 'MM3': '25'},
        'P_3+5': {'element': 'P', 'MM3': '25'},
        'P_60': {'element': 'P', 'MM3': '60'},
        'P_153': {'element': 'P', 'MM3': '153'},
        'P_3+Q': {'element': 'P', 'MM3': '25'},
        'B_2': {'element': 'B', 'MM3': '26'},
        'B_3': {'element': 'B', 'MM3': '27'},
        'SI3': {'element': 'Si', 'MM3': '19'},
        'SN3': {'element': 'Se', 'MM3': '32'},
        'HE4+4': {'element': 'He', 'MM3': '51'},
        'NE4+4': {'element': 'Ne', 'MM3': '52'},
        'AR4+4': {'element': 'Ar', 'MM3': '53'},
        'KR4+4': {'element': 'Kr', 'MM3': '54'},
        'XE4+4': {'element': 'Xe', 'MM3': '55'},
        'LI': {'element': 'Li', 'MM3': '163'},
        'MG3+2': {'element': 'Mg', 'MM3': '59'},
        'GE3': {'element': 'Ge', 'MM3': '31'},
        'PB3': {'element': 'Pb', 'MM3': '33'},
        'SE3+2': {'element': 'Se', 'MM3': '34'},
        'TE3+2': {'element': 'Te', 'MM3': '35'},
        'CA6+2': {'element': 'Ca', 'MM3': '125'},
        'SR6+2': {'element': 'Sr', 'MM3': '126'},
        'BA6+2': {'element': 'Ba', 'MM3': '127'},
        'FE3+2': {'element': 'Fe', 'MM3': '62'},
        'FE_61': {'element': 'Fe', 'MM3': '61'},
        'FE6+2': {'element': 'Fe', 'MM3': '62'},
        'NI4+2': {'element': 'Ni', 'MM3': '63'},
        'NI_64': {'element': 'Ni', 'MM3': '64'},
        'CO6+3': {'element': 'Co', 'MM3': '66'},
        'CO_65': {'element': 'Co', 'MM3': '65'},
        'LA3+3': {'element': 'La', 'MM3': '128'},
        'CE6+3': {'element': 'Ce', 'MM3': '129'},
        'PR6+3': {'element': 'Pr', 'MM3': '130'},
        'ND6+3': {'element': 'Nd', 'MM3': '131'},
        'PM6+3': {'element': 'Pm', 'MM3': '132'},
        'SM6+3': {'element': 'Sm', 'MM3': '133'},
        'EU6+3': {'element': 'Eu', 'MM3': '134'},
        'GD6+3': {'element': 'Gd', 'MM3': '135'},
        'TB6+3': {'element': 'Tb', 'MM3': '136'},
        'DY6+3': {'element': 'Dy', 'MM3': '137'},
        'HO6+3': {'element': 'Ho', 'MM3': '138'},
        'ER6+3': {'element': 'Er', 'MM3': '139'},
        'TM6+3': {'element': 'Tm', 'MM3': '140'},
        'YB6+3': {'element': 'Yb', 'MM3': '141'},
        'LU6+3': {'element': 'Lu', 'MM3': '142'},
        'BE3+2': {'element': 'Be', 'MM3': '165'},
        'NA': {'element': 'Na', 'MM3': '166'},
        'AL3': {'element': 'Al', 'MM3': '167'},
        'K_': {'element': 'K', 'MM3': '168'},
        'SC3+3': {'element': 'Sc', 'MM3': '169'},
        'TI3+4': {'element': 'Ti', 'MM3': '170'},
        'TI6+4': {'element': 'Ti', 'MM3': '170'},
        'V_3+5': {'element': 'V', 'MM3': '171'},
        'CR6+3': {'element': 'Cr', 'MM3': '172'},
        'MN6+2': {'element': 'Mn', 'MM3': '173'},
        'CU3+1': {'element': 'Cu', 'MM3': '174'},
        'ZN3+2': {'element': 'Zn', 'MM3': '175'},
        'GA3+3': {'element': 'Ga', 'MM3': '176'},
        'AS3+3': {'element': 'As', 'MM3': '177'},
        'RB': {'element': 'Rb', 'MM3': '178'},
        'Y_3+3': {'element': 'Y', 'MM3': '179'},
        'ZR3+4': {'element': 'Zr', 'MM3': '180'},
        'NB3+5': {'element': 'Nb', 'MM3': '181'},
        'MO6+6': {'element': 'Mo', 'MM3': '182'},
        'MO3+6': {'element': 'Mo', 'MM3': '182'},
        'TC6+5': {'element': 'Tc', 'MM3': '183'},
        'RU6+2': {'element': 'Ru', 'MM3': '184'},
        'RH6+3': {'element': 'Rh', 'MM3': '185'},
        'PD4+2': {'element': 'Pd', 'MM3': '186'},
        'AG1+1': {'element': 'Ag', 'MM3': '187'},
        'CD3+2': {'element': 'Cd', 'MM3': '188'},
        'IN3+3': {'element': 'In', 'MM3': '189'},
        'SB3+3': {'element': 'Sb', 'MM3': '190'},
        'CS': {'element': 'Cs', 'MM3': '191'},
        'HF3+4': {'element': 'Hf', 'MM3': '192'},
        'TA3+5': {'element': 'Ta', 'MM3': '193'},
        'W_6+6': {'element': 'W', 'MM3': '194'},
        'W_3+4': {'element': 'W', 'MM3': '194'},
        'W_3+6': {'element': 'W', 'MM3': '194'},
        'RE6+5': {'element': 'Re', 'MM3': '195'},
        'RE3+7': {'element': 'Re', 'MM3': '195'},
        'OS6+6': {'element': 'Os', 'MM3': '196'},
        'IR6+3': {'element': 'Ir', 'MM3': '197'},
        'PT4+4': {'element': 'Pt', 'MM3': '198'},
        'AU4+3': {'element': 'Au', 'MM3': '199'},
        'HG1+2': {'element': 'Hg', 'MM3': '200'},
        'TL3+3': {'element': 'Tl', 'MM3': '201'},
        'BI3+3': {'element': 'Bi', 'MM3': '202'},
        'PO3+2': {'element': 'Po', 'MM3': '203'},
        'AT': {'element': 'At', 'MM3': '204'},
        'RN4+4': {'element': 'Rn', 'MM3': '205'},
        'FR': {'element': 'Fr', 'MM3': '206'},
        'RA6+2': {'element': 'Ra', 'MM3': '207'},
        'AC6+3': {'element': 'Ac', 'MM3': '208'},
        'TH6+4': {'element': 'Th', 'MM3': '209'},
        'PA6+4': {'element': 'Pa', 'MM3': '210'},
        'U_6+4': {'element': 'U', 'MM3': '211'},
        'NP6+4': {'element': 'Np', 'MM3': '212'},
        'PU6+4': {'element': 'Pu', 'MM3': '213'},
        'AM6+4': {'element': 'Am', 'MM3': '214'},
        'CM6+3': {'element': 'Cm', 'MM3': '215'},
        'BK6+3': {'element': 'Bk', 'MM3': '216'},
        'CF6+3': {'element': 'Cf', 'MM3': '217'},
        'ES6+3': {'element': 'Es', 'MM3': '218'},
        'FM6+3': {'element': 'Fm', 'MM3': '219'},
        'MD6+3': {'element': 'Md', 'MM3': '220'},
        'NO6+3': {'element': 'No', 'MM3': '221'},
        'LW6+3': {'element': 'Lw', 'MM3': '222'}
    },
    'MM3': {
        '1': {'element': 'C'},
        '2': {'element': 'C'},
        '3': {'element': 'C'},
        '4': {'element': 'C'},
        '5': {'element': 'H'},
        '6': {'element': 'O'},
        '7': {'element': 'O'},
        '8': {'element': 'N'},
        '9': {'element': 'N'},
        '10': {'element': 'N'},
        '11': {'element': 'F'},
        '12': {'element': 'Cl'},
        '13': {'element': 'Br'},
        '14': {'element': 'I'},
        '15': {'element': 'S'},
        '16': {'element': 'S'},
        '17': {'element': 'S'},
        '18': {'element': 'S'},
        '19': {'element': 'Si'},
        '20': {'element': ''},
        '21': {'element': 'H'},
        '22': {'element': 'C'},
        '23': {'element': 'H'},
        '24': {'element': 'H'},
        '25': {'element': 'P'},
        '26': {'element': 'B'},
        '27': {'element': 'B'},
        '28': {'element': 'H'},
        '29': {'element': 'C'},
        '30': {'element': 'C'},
        '31': {'element': 'Ge'},
        '32': {'element': 'Sn'},
        '33': {'element': 'Pb'},
        '34': {'element': 'Se'},
        '35': {'element': 'Te'},
        '36': {'element': 'C'},
        '37': {'element': 'N'},
        '38': {'element': 'C'},
        '39': {'element': 'N'},
        '40': {'element': 'N'},
        '41': {'element': 'O'},
        '42': {'element': 'S'},
        '43': {'element': 'N'},
        '44': {'element': 'H'},
        '45': {'element': 'N'},
        '46': {'element': 'N'},
        '47': {'element': 'O'},
        '48': {'element': 'H'},
        '49': {'element': 'O'},
        '50': {'element': 'C'},
        '51': {'element': 'He'},
        '52': {'element': 'Ne'},
        '53': {'element': 'Ar'},
        '54': {'element': 'Kr'},
        '55': {'element': 'Xe'},
        '56': {'element': 'C'},
        '57': {'element': 'C'},
        '58': {'element': 'C'},
        '59': {'element': 'Mg'},
        '60': {'element': 'P'},
        '61': {'element': 'Fe'},
        '62': {'element': 'Fe'},
        '63': {'element': 'Ni'},
        '64': {'element': 'Ni'},
        '65': {'element': 'Co'},
        '66': {'element': 'Co'},
        '67': {'element': 'C'},
        '68': {'element': 'C'},
        '69': {'element': 'O'},
        '70': {'element': 'O'},
        '71': {'element': 'C'},
        '72': {'element': 'N'},
        '73': {'element': 'H'},
        '74': {'element': 'S'},
        '75': {'element': 'O'},
        '76': {'element': 'O'},
        '77': {'element': 'O'},
        '78': {'element': 'O'},
        '79': {'element': 'O'},
        '80': {'element': 'O'},
        '81': {'element': 'O'},
        '82': {'element': 'O'},
        '83': {'element': 'O'},
        '84': {'element': 'O'},
        '85': {'element': 'O'},
        '86': {'element': 'O'},
        '87': {'element': 'O'},
        '88': {'element': 'O'},
        '89': {'element': 'O'},
        '90': {'element': 'O'},
        '91': {'element': 'O'},
        '92': {'element': 'O'},
        '93': {'element': 'O'},
        '94': {'element': 'O'},
        '95': {'element': 'O'},
        '96': {'element': 'O'},
        '97': {'element': 'O'},
        '98': {'element': 'O'},
        '99': {'element': 'O'},
        '100': {'element': 'O'},
        '101': {'element': 'O'},
        '102': {'element': 'O'},
        '103': {'element': 'O'},
        '104': {'element': 'S'},
        '105': {'element': 'S'},
        '106': {'element': 'C'},
        '107': {'element': 'N'},
        '108': {'element': 'N'},
        '109': {'element': 'N'},
        '110': {'element': 'N'},
        '111': {'element': 'N'},
        '112': {'element': ''},
        '113': {'element': 'C'},
        '114': {'element': 'C'},
        '115': {'element': 'O'},
        '116': {'element': 'O'},
        '117': {'element': 'O'},
        '118': {'element': 'O'},
        '119': {'element': 'O'},
        '120': {'element': 'O'},
        '121': {'element': 'O'},
        '122': {'element': ''},
        '123': {'element': ''},
        '124': {'element': 'H'},
        '125': {'element': 'Ca'},
        '126': {'element': 'Sr'},
        '127': {'element': 'Ba'},
        '128': {'element': 'La'},
        '129': {'element': 'Ce'},
        '130': {'element': 'Pr'},
        '131': {'element': 'Nd'},
        '132': {'element': 'Pm'},
        '133': {'element': 'Sm'},
        '134': {'element': 'Eu'},
        '135': {'element': 'Gd'},
        '136': {'element': 'Tb'},
        '137': {'element': 'Dy'},
        '138': {'element': 'Ho'},
        '139': {'element': 'Er'},
        '140': {'element': 'Tm'},
        '141': {'element': 'Yb'},
        '142': {'element': 'Lu'},
        '143': {'element': 'N'},
        '144': {'element': 'N'},
        '145': {'element': 'O'},
        '146': {'element': 'N'},
        '147': {'element': ''},
        '148': {'element': 'O'},
        '149': {'element': 'O'},
        '150': {'element': 'N'},
        '151': {'element': 'N'},
        '152': {'element': ''},
        '153': {'element': 'P'},
        '154': {'element': 'S'},
        '155': {'element': 'N'},
        '156': {'element': ''},
        '157': {'element': ''},
        '158': {'element': ''},
        '159': {'element': 'O'},
        '160': {'element': 'C'},
        '161': {'element': 'C'},
        '162': {'element': 'C'},
        '163': {'element': 'Li'},
        '164': {'element': 'N'}
    }
}

"""
MM_SET_TYPES = {
    'MM3': ['Convert from gaffType attrib (GAFF)',
            'Convert from mol2type attrib (UFF)',
            'Convert from mol2type attrib (GAFF)'],
    'GAFF': ['Catch from gaffType attribute',
            'Catch from mol2type attribute']
}
"""

MM3_FROM_UFF = {
    # MM3 no metal elements (C,H,O,N,halogens,S,P,B,Si,Se,noble gases)          
    'C_3': '1', # Csp3     
    'C_22': '22',   # cyclopropane C    
    'C_56': '56',   # cyclobutane C    
    'C_2': '2', # Csp2     
    'C_3C': '3',    # C carbonyl    
    'C_29': '29',   # radical C    
    'C_30': '30',   # carbonium ion    
    'C_38': '38',   # cyclopropene C    
    'C_57': '57',   # cyclobutene C    
    'C_113': '113', # Cferrocene-H     
    'C_114': '114', # Cferrocene-C     
    'C_R': '2', # Carom     
    'C_50': '50',   # benzene (localized) C   
    'C_58': '58',   # cyclobutanone C    
    'C_67': '67',   # cyclopropanone C    
    'C_71': '71',   # ketonium C    
    'C_1': '4', # C sp    
    'C_68': '68',   # allene C    
    'C_106': '106', # ketene C    
    'H_': '5',  # normal terminal H   
    'H_21': '21',   # O-H alcohol H   
    'H_23': '23',   # NH amine/imine H   
    'H_24': '24',   # COOH carboxyl H   
    'H_28': '28',   # N-H amide H   
    'H_36': '36',   # deuterium D    
    'H_44': '44',   # SH thiol H   
    'H_48': '48',   # ammonium H    
    'H_73': '73',   # OH enol/phenol H   
    'H_124': '124', # CH acetylene H   
    'H_B': '5', # H bridge    
    'O_3': '6', # O sp3    
    'O_49': '49',   # epoxy O    
    'O_69': '69',   # amine oxide O   
    'O_70': '70',   # ketonium O    
    'O_75': '75',   # carboxyl O-H    
    'O_145': '145', # hydroxyamine O    
    'O_3_Z': '6',   # O sp3    
    'O_R': '6', # O sp3    
    'O_41': '41',   # furane O    
    'O_47': '47',   # carboxylate anion O   
    'O_77': '77',   # acid C=O    
    'O_78': '78',   # ester C=O    
    'O_79': '79',   # amide C=O    
    'O_80': '80',   # halide C=O    
    'O_81': '81',   # enone C=O    
    'O_82': '82',   # anhydride C=O    
    'O_148': '148', # anhydride (localized) O   
    'O_149': '149', # anhydride (delocalized) O   
    'O_2': '7', # O ketone or aldehyde  
    'O_159': '159', # phosphate O    
    'O_1': '7', # O sp, M-CO   
    'N_3': '8', # N sp3    
    'N_39': '39',   # ammonium N    
    'N_146': '146', # hydroxyamine N    
    'N_150': '150', # sp3 hydrazine N   
    'N_155': '155', # sp3 sulfonamide N   
    'N_164': '164', # sp3 Li-amide N   
    'N_R': '37',    # N pyr, -N=C- (deloc)  
    'N_40': '40',   # pyrrole N    
    'N_143': '143', # axoxy N (delocalized)   
    'N_144': '144', # azoxy N (delocalized)   
    'N_151': '151', # sp2 amide (delocalized)   
    'N_2': '9', # N sp2)    
    'N_43': '43',   # azoxy =N-O    
    'N_45': '45',   # azide center N   
    'N_46': '46',   # nitro N    
    'N_72': '72',   # imine N (localized)   
    'N_107': '107', # azo N (localized)   
    'N_108': '108', # oxime N    
    'N_109': '109', # azoxy N (localized)   
    'N_110': '110', # cationic N of imminium  
    'N_111': '111', # cationic N of pyridinium  
    'N_1': '10',    # N sp    
    'F_': '11', # F general    
    'CL': '12', # Cl general    
    'BR': '13', # Br general    
    'I_': '14', # I general    
    'S_3+2': '15',  # S sulfide    
    'S_16': '16',   # sulfonium S    
    'S_17': '17',   # sulfoxide S    
    'S_18': '18',   # sulfone S    
    'S_104': '104', # disulfide bridge S   
    'S_105': '105', # polysulfide S    
    'S_3+4': '154', # SX3-R, selected as sulfonamide in MM3
    'S_3+6': '154', # SO4, selected as sulfonamide in MM3
    'S_R': '42',    # sp2 thiophene S   
    'S_2': '74',    # C=S thiocarbonyl    
    'P_3+3': '25',  # general phosphine    
    'P_3+5': '25',  # general phosphine    
    'P_60': '60',   # phosphorus P(V)    
    'P_153': '153', # phosphate P    
    'P_3+Q': '25',  # general phosphine    
    'B_2': '26',    # sp2 B    
    'B_3': '27',    # sp3 B    
    'SI3': '19',    # Si general    
    'SN3': '32',    # Se general    
    'HE4+4': '51',  # He general    
    'NE4+4': '52',  # Ne general    
    'AR4+4': '53',  # Ar general    
    'KR4+4': '54',  # Kr general    
    'XE4+4': '55',  # Xe general    
    
    # MM3 representative metal and semimetal elements        
    'LI': '163',    # Li general    
    'MG3+2': '59',  # Mg general    
    'GE3': '31',    # Ge general    
    'PB3': '33',    # Pb (IV)    
    'SE3+2': '34',  # Se general    
    'TE3+2': '35',  # Te general    
    'CA6+2': '125', # Ca general    
    'SR6+2': '126', # Sr general    
    'BA6+2': '127', # Ba general    
    
    # MM3 no representative metal elements         
    'FE3+2': '62',  # Fe (III)    
    'FE_61': '61',  # Fe(II)     
    'FE6+2': '62',  # Fe (III)    
    'NI4+2': '63',  # Ni (II)    
    'NI_64': '64',  # Ni(III)     
    'CO6+3': '66',  # Co (III)    
    'CO_65': '65',  # Co(II)     
    'LA3+3': '128', # La general    
    'CE6+3': '129', # Ce general    
    'PR6+3': '130', # Pr general    
    'ND6+3': '131', # Nd general    
    'PM6+3': '132', # Pm general    
    'SM6+3': '133', # Sm general    
    'EU6+3': '134', # Eu general    
    'GD6+3': '135', # Gd general    
    'TB6+3': '136', # Tb general    
    'DY6+3': '137', # Dy general    
    'HO6+3': '138', # Ho general    
    'ER6+3': '139', # Er general    
    'TM6+3': '140', # Tm general    
    'YB6+3': '141', # Yb general    
    'LU6+3': '142', # Lu general    
   
    # Transformation of UFF atoms which don't exist in MM3      
    'BE3+2': '165', # Be     
    'NA': '166',    # Na     
    'AL3': '167',   # Al     
    'K_': '168',    # K     
    'SC3+3': '169', # Sc     
    'TI3+4': '170', # Ti     
    'TI6+4': '170', # Ti     
    'V_3+5': '171', # V     
    'CR6+3': '172', # Cr     
    'MN6+2': '173', # Mn     
    'CU3+1': '174', # Cu     
    'ZN3+2': '175', # Zn     
    'GA3+3': '176', # Ga     
    'AS3+3': '177', # As     
    'RB': '178',    # Rb     
    'Y_3+3': '179', # Y     
    'ZR3+4': '180', # Zr     
    'NB3+5': '181', # Nb     
    'MO6+6': '182', # Mo     
    'MO3+6': '182', # Mo     
    'TC6+5': '183', # Tc     
    'RU6+2': '184', # Ru     
    'RH6+3': '185', # Rh     
    'PD4+2': '186', # Pd     
    'AG1+1': '187', # Ag     
    'CD3+2': '188', # Cd     
    'IN3+3': '189', # In     
    'SB3+3': '190', # Sb     
    'CS': '191',    # Cs     
    'HF3+4': '192', # Hf     
    'TA3+5': '193', # Ta     
    'W_6+6': '194', # W     
    'W_3+4': '194', # W     
    'W_3+6': '194', # W     
    'RE6+5': '195', # Re     
    'RE3+7': '195', # Re     
    'OS6+6': '196', # Os     
    'IR6+3': '197', # Ir     
    'PT4+4': '198', # Pt     
    'AU4+3': '199', # Au     
    'HG1+2': '200', # Hg     
    'TL3+3': '201', # Tl     
    'BI3+3': '202', # Bi     
    'PO3+2': '203', # Po     
    'AT': '204',    # At     
    'RN4+4': '205', # Rn     
    'FR': '206',    # Fr     
    'RA6+2': '207', # Ra     
    'AC6+3': '208', # Ac     
    'TH6+4': '209', # Th     
    'PA6+4': '210', # Pa     
    'U_6+4': '211', # U     
    'NP6+4': '212', # Np     
    'PU6+4': '213', # Pu     
    'AM6+4': '214', # Am     
    'CM6+3': '215', # Cm     
    'BK6+3': '216', # Bk     
    'CF6+3': '217', # Cf     
    'ES6+3': '218', # Es     
    'FM6+3': '219', # Fm     
    'MD6+3': '220', # Md     
    'NO6+3': '221', # No     
    'LW6+3': '222' # Lw 
}

GAFF_DESC = {
    'c': 'Sp2 C carbonyl group',
    'c1': 'Sp C',
    'c2': 'Sp2 C', 
    'c3': 'Sp3 C',
    'ca': 'Sp2 C in pure aromatic systems',
    'cp': 'Head Sp2 C that connect two rings in biphenyl sys.',
    'cq': 'Head Sp2 C that connect two rings in biphenyl sys. identical to cp',
    'cc': 'Sp2 carbons in non-pure aromatic systems',
    'cd': 'Sp2 carbons in non-pure aromatic systems, identical to cc',
    'ce': 'Inner Sp2 carbons in conjugated systems',
    'cf': 'Inner Sp2 carbons in conjugated systems, identical to ce',
    'cg': 'Inner Sp carbons in conjugated systems',
    'ch': 'Inner Sp carbons in conjugated systems, identical to cg',
    'cx': 'Sp3 carbons in triangle systems',
    'cy': 'Sp3 carbons in square systems',
    'cu': 'Sp2 carbons in triangle systems',
    'cv': 'Sp2 carbons in square systems',
    'cz': 'Sp2 carbon in guanidine group',
    'h1': 'H bonded to aliphatic carbon with 1 electrwd. group',
    'h2': 'H bonded to aliphatic carbon with 2 electrwd. group',
    'h3': 'H bonded to aliphatic carbon with 3 electrwd. group',
    'h4': 'H bonded to non-sp3 carbon with 1 electrwd. group',
    'h5': 'H bonded to non-sp3 carbon with 2 electrwd. group',
    'ha': 'H bonded to aromatic carbon',
    'hc': 'H bonded to aliphatic carbon without electrwd. group',
    'hn': 'H bonded to nitrogen atoms',
    'ho': 'Hydroxyl group',
    'hp': 'H bonded to phosphate',
    'hs': 'Hydrogen bonded to sulphur',
    'hw': 'Hydrogen in water',
    'hx': 'H bonded to C next to positively charged group',
    'f':  'Fluorine',
    'cl': 'Chlorine',
    'br': 'Bromine',
    'i':  'Iodine',
    'n':  'Sp2 nitrogen in amide groups',
    'n1': 'Sp N',
    'n2': 'Aliphatic Sp2 N with two connected atoms',
    'n3': 'Sp3 N with three connected atoms',
    'n4': 'Sp3 N with four connected atoms',
    'na': 'Sp2 N with three connected atoms',
    'nb': 'Sp2 N in pure aromatic systems',
    'nc': 'Sp2 N in non-pure aromatic systems',
    'nd': 'Sp2 N in non-pure aromatic systems',
    'ne': 'Inner Sp2 N in conjugated systems',
    'nf': 'Inner Sp2 N in conjugated systems',
    'nh': 'Amine N connected one or more aromatic rings',
    'no': 'Nitro N',
    'o': 'Oxygen with one connected atom',
    'oh': 'Oxygen in hydroxyl group',
    'os': 'Ether and ester oxygen',
    'ow': 'Oxygen in water',
    'p2': 'Phosphate with two connected atoms',
    'p3': 'Phosphate with three connected atoms, such as PH3',
    'p4': 'Phosphate with three connected atoms, such as O=P(CH3)2',
    'p5': 'Phosphate with four connected atoms, such as O=P(OH)3',
    'pb': 'Sp2 P in pure aromatic systems',
    'pc': 'Sp2 P in non-pure aromatic systems',
    'pd': 'Sp2 P in non-pure aromatic systems, identical to pc',
    'pe': 'Inner Sp2 P in conjugated systems',
    'pf': 'Inner Sp2 P in conjugated systems, identical to pe',
    'px': 'Special p4 in conjugated systems',
    'py': 'Special p5 in conjugated systems',
    's': 'S with one connected atom',
    's2': 'S with two connected atom, involved at least one double bond',
    's4': 'S with three connected atoms', 
    's6': 'S with four connected atoms',
    'sh': 'Sp3 S connected with hydrogen',
    'ss': 'Sp3 S in thio-ester and thio-ether',
    'sx': 'Special s4 in conjugated systems',
    'sy': 'Special s6 in conjugated systems'
}

MM3_FROM_ELEMENT = {
	'Ac': ['208'],
	'Ag': ['187'],
	'Al': ['167'],
	'Am': ['214'],
	'Ar': ['53'],
	'As': ['177'],
	'At': ['204'],
	'Au': ['199'],
	'B': ['26', '27'],
	'Ba': ['127'],
	'Be': ['165'],
	'Bi': ['202'],
	'Bk': ['216'],
	'Br': ['13'],
	'C': ['1', '2', '3', '4', '22', '29', '30', '38', '50', '56', '57', '58', '67', '68', '71', '106', '113', '114', '160', '161', '162'],
	'Ca': ['125'],
	'Cd': ['188'],
	'Ce': ['129'],
	'Cf': ['217'],
	'Cl': ['12'],
	'Cm': ['215'],
	'Co': ['65', '66'],
	'Cr': ['172'],
	'Cs': ['191'],
	'Cu': ['174'],
	'D': ['36'],
	'Dy': ['137'],
	'Er': ['139'],
	'Es': ['218'],
	'Eu': ['134'],
	'F': ['11'],
	'Fe': ['61', '62'],
	'Fm': ['219'],
	'Fr': ['206'],
	'Ga': ['176'],
	'Gd': ['135'],
	'Ge': ['31'],
	'H': ['5', '21', '23', '24', '28', '44', '48', '73', '124'],
	'He': ['51'],
	'Hf': ['192'],
	'Hg': ['200'],
	'Ho': ['138'],
	'I': ['14'],
	'In': ['189'],
	'Ir': ['197'],
	'K': ['168'],
	'Kr': ['54'],
	'La': ['128'],
	'Li': ['163'],
	'Lu': ['142'],
	'Lw': ['222'],
	'Md': ['220'],
	'Mg': ['59'],
	'Mn': ['173'],
	'Mo': ['182'],
	'N': ['8', '9', '10', '37', '39', '40', '43', '45', '46', '72', '107', '108', '109', '110', '111', '143', '144', '146', '150', '151', '155', '164'],
	'Na': ['166'],
	'Nb': ['181'],
	'Nd': ['131'],
	'Ne': ['52'],
	'Ni': ['63', '64'],
	'No': ['221'],
	'Np': ['212'],
	'O': ['6', '7', '41', '47', '49', '69', '70', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '115', '116', '117', '118', '119', '120', '121', '145', '148', '149', '159'],
	'Os': ['196'],
	'P': ['25', '60', '153'],
	'Pa': ['210'],
	'Pb': ['33'],
	'Pd': ['186'],
	'Pm': ['132'],
	'Po': ['203'],
	'Pr': ['130'],
	'Pt': ['198'],
	'Pu': ['213'],
	'Ra': ['207'],
	'Rb': ['178'],
	'Re': ['195'],
	'Rh': ['185'],
	'Rn': ['205'],
	'Ru': ['184'],
	'S': ['15', '16', '17', '18', '42', '74', '104', '105', '154'],
	'Sb': ['190'],
	'Sc': ['169'],
	'Se': ['34'],
	'Si': ['19'],
	'Sm': ['133'],
	'Sn': ['32'],
	'Sr': ['126'],
	'Ta': ['193'],
	'Tb': ['136'],
	'Tc': ['183'],
	'Te': ['35'],
	'Th': ['209'],
	'Ti': ['170'],
	'Tl': ['201'],
	'Tm': ['140'],
	'U': ['211'],
	'V': ['171'],
	'W': ['194'],
	'Xe': ['55'],
	'Y': ['179'],
	'Yb': ['141'],
	'Zn': ['175'],
	'Zr': ['180']
}

MM3_FROM_GAFF = {
	'c': ['3', '2'],
	'c1': ['4'],
	'c2': ['2', '38', '57'],
	'c3': ['1', '22', '56'],
	'ca': ['2'],
	'cp': ['2'],
	'cq': ['2'],
	'cc': ['2'],
	'cd': ['2'],
	'ce': ['2'],
	'cf': ['2'],
	'cg': ['4'],
	'ch': ['4'],
	'cx': ['22'],
	'cy': ['56'],
	'cu': ['38'],
	'cv': ['57'],
	'cz': ['2'],
	'h1': ['5'],
	'h2': ['5'],
	'h3': ['5'],
	'h4': ['5'],
	'h5': ['5'],
	'ha': ['5'],
	'hc': ['5', '124'],
	'hn': ['23', '28', '48'],
	'ho': ['21', '73'],
	'hp': ['5'],
	'hs': ['44'],
	'hw': ['21'],
	'hx': ['5'],
	'f': ['11'],
	'cl': ['12'],
	'br': ['13'],
	'I': ['14'],
	'n': ['9'],
	'n1': ['10'],
	'n2': ['37'],
	'n3': ['8'],
	'n4': ['39'],
	'na': ['8', '40'],
	'nb': ['8'],
	'nc': ['37'],
	'nd': ['37'],
	'ne': ['37'],
	'nf': ['37'],
	'nh': ['8'],
	'no': ['46'],
	'o': ['7'],
	'oh': ['6'],
	'os': ['6', '41'],
	'ow': ['6'],
	'p2': ['60'],
	'p3': ['25'],
	'p4': ['25'],
	'p5': ['60'],
	'pb': ['60'],
	'pc': ['60'],
	'pd': ['60'],
	'pe': ['60'],
	'pf': ['60'],
	'px': ['60'],
	'py': ['60'],
	's': ['42'],
	's2': ['15'],
	's4': ['15', '17'],
	's6': ['18'],
	'sh': ['15'],
	'ss': ['15', '42', '104'],
	'sx': ['17'],
	'sy': ['18']
}

"""
(1)		Parameters for p4 are not enough
(2)     For bond lengty, bond angle and torsional angle parameters, hc and ha 
        are the generic names of aliphatic and aromatic hydrogen connected 
        to carbon, respectively.
(3) 	Defination priority of nitrogen: n>n3, n>n4, n>nh, n2>n, no>n, na>n
        n4>nh
(5)     Although the four atom types are same, maybe there are several different
        parameters correspond to different bond types. You may need to use 
        additional atom types to apply them correctly (with antechamber package,
        you can easily do this). 
(4) 	Polarizabilities: mg2+ 0.120, f- 0.9743
"""

class GaussianInputFile(object):

    """
    Object-oriented abstraction of a Gaussian input file.

    Represents an input file, section by section. Almost everything
    can be set directly from the class initialization, but some
    parameters will need manual assignment with special methods.

    Implemented sections are:

    - Link 0 (% commands)
    - Route (# lines)
    - Title
    - Molecule specification (see `pygaussian.GaussianAtom`)
    - Restraints (ModRedundant)
    - Extra basis sets

    It also offers some support for QM/MM jobs.

    More info about Gaussian input formats can be consulted
    in http://gaussian.com/input/.
    """

    def __init__(self, title='Untitled job', connectivity=False, *args, **kwargs):
        self.title = title
        self.connectivity = connectivity
        self._route = defaultdict(list)
        self._link = {}
        self._job = None
        self._charge = None
        self._mm_charge = None
        self._multiplicity = None
        self._mm_multiplicity = None
        self._qm_method = None
        self._qm_functional = None
        self._qm_basis_set = None
        self._qm_basis_sets_extra = []
        self._mm_forcefield = None
        self._mm_water_forcefield = None
        self._mm_forcefield_extra = []
        self._atoms = []
        self._restraints = []

        # Set and verify
        for k, v in kwargs.items():
            if hasattr(self, k) or isinstance(getattr(type(self), k, None), property):
                try:
                    setattr(self, k, v)
                except (TypeError, ValueError) as e:
                    print('! Could not set {} with value'
                          ' {} because {}'.format(k, v, e))
            else:
                print('! Keyword {} not recognized'.format(k))

    def __str__(self):
        return self.build(timestamp=False)

    def build(self, timestamp=True):
        sections = []
        if timestamp:
            sections.append(self.timestamp)
        # Link 0
        link = self.link
        if link:
            sections.append(link)
        # Route
        sections.append(self.route)
        sections.append('')
        # Title
        sections.append(self.title)
        sections.append('')
        # Charge, multiplicity and atoms
        sections.append(self.system)
        sections.append('')
        # Variables and configuration
        restraints = self.restraints
        if restraints:
            sections.append(restraints)
            sections.append('')
        if self.connectivity:
            sections.append(self.compute_connectivity())
            sections.append('')
        mm_forcefield_extra = self.mm_forcefield_extra
        if mm_forcefield_extra:
            sections.append('')
            sections.append(mm_forcefield_extra)
            sections.append('')
        qm_basis_set_extra = self.qm_basis_set_extra
        if qm_basis_set_extra:
            sections.append(qm_basis_set_extra)
            sections.append('')
        sections.append('')
        sections.append('')
        return '\n'.join(sections)

    @property
    def timestamp(self):
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        return '! Generated by pygaussian at ' + now

    # Link 0 section
    @property
    def link(self):
        s = []
        for keyword, options in self._link.items():
            if not options:
                s.append('%' + keyword)
            elif isinstance(options, (tuple, list)):
                s.append('%{}={}'.format(keyword, ','.join(options)))
            else:
                s.append('%{}={}'.format(keyword, options))
        return '\n'.join(s)

    def add_link_option(self, keyword, *options):
        self._link[keyword] = options

    @property
    def processors(self):
        return self._link.get('nprocshared')

    @processors.setter
    def processors(self, value):
        if value < 1:
            raise ValueError('processors must be greater than 0')
        self._link['nprocshared'] = int(value)

    @property
    def memory(self):
        return self._memory

    @memory.setter
    def memory(self, value):
        mem_units = 'GB'
        if isinstance(value, (list, tuple)):
            if len(value) == 1:
                value = value
            elif len(value) == 2:
                value, mem_units = value
            else:
                ValueError('Set memory with: <value>[, "units"]')
        if value <= 0:
            raise ValueError('Memory value must be greater than 0')
        if mem_units not in MEM_UNITS:
            raise ValueError('Memory unit not recognized.')
        self._memory = int(value), mem_units
        self._link['mem'] = '{}{}'.format(*self._memory)

    @property
    def checkpoint(self):
        return self._link.get('chk')

    @checkpoint.setter
    def checkpoint(self, value):
        if value:
            self._link['chk'] = value

    ########################################################
    @property
    def route(self):
        s = ['#p', self.modeling]
        for keyword, options in self._route.items():
            valid = [o for o in options if o]
            if valid:
                s.append('{}=({})'.format(keyword, ','.join(valid)))
            else:
                s.append(keyword)
        return ' '.join(s)

    def add_route_option(self, keyword, *options):
        self._route[keyword] = options

    @property
    def modeling(self):
        method = self.qm_method
        if method == 'DFT':
            method = self.qm_functional
        if not method:
            raise ValueError('Method must be set with `qm_method`. If '
                             'using DFT, `qm_functional` is also needed.')
        basis = 'gen' if self._qm_basis_sets_extra else self.qm_basis_set
        if not basis:
            raise ValueError('Basis set must be set with `qm_basis_set` '
                             'or `qm_basis_sets_extra`')
        forcefield = self.mm_forcefield
        if forcefield:
            if basis == 'gen':  # QMMM jobs usually specify ECP too
                basis = 'genecp'
            if self._mm_forcefield_extra:
                return 'oniom=({}/{}:{}/hardfirst)'.format(method, basis, forcefield)
            return'oniom=({}/{}:{})'.format(method, basis, forcefield)
        return '{}/{}'.format(method, basis)

    # Job type
    @property
    def job(self):
        job = self._job
        if job is None:
            raise ValueError('Job is not set!')
        return job

    @job.setter
    def job(self, value):
        if value not in JOB_TYPES:
            raise ValueError('Job must be either of {}'.format(JOB_TYPES))
        self._job = value
        self._route[value] = []

    @property
    def job_options(self):
        return self._route.get(self._job)

    @job_options.setter
    def job_options(self, value):
        if not isinstance(value, (tuple, list)):
            value = [value]
        self._route[self._job].extend(value)

    @property
    def freq(self):
        if 'freq' in self._route:
            freq_ = self._route.get('freq')
            if freq_:
                return freq_
            return True
        return False

    @freq.setter
    def freq(self, value):
        if self._job not in ('Opt', 'Polar'):
            raise ValueError('Freq can only be set if job is Opt or Polar.')
        if isinstance(value, basestring):
            self._route['freq'] = value
        elif value is True:
            self._route['freq'] = []
        else:
            raise ValueError('freq must be str or bool')

    # QM model
    @property
    def qm_method(self):
        return self._qm_method

    @qm_method.setter
    def qm_method(self, value):
        if value not in QM_METHODS:
            raise ValueError('Method must be either of: '
                             '{}'.format(', '.join(QM_METHODS)))
        self._qm_method = value

    @property
    def qm_functional(self):
        return self._qm_functional

    @qm_functional.setter
    def qm_functional(self, value):
        if self.qm_method != 'DFT':
            raise ValueError('Funtionals can only be set if method == DFT')
        if value not in QM_FUNCTIONALS_ALL:
            raise ValueError('functional {} not recognized'.format(value))
        self._qm_functional = value

    @property
    def qm_basis_set(self):
        if self._qm_basis_sets_extra:
            return self._qm_basis_set, self._qm_basis_sets_extra
        return self._qm_basis_set

    @qm_basis_set.setter
    def qm_basis_set(self, value):
        if value.rstrip('+*') not in QM_BASIS_SETS:
            raise ValueError('Basis set {} not recognized. '
                             'Try with one of: {}'.format(value, ', '.join(QM_BASIS_SETS)))
        self._qm_basis_set = value

    @property
    def qm_basis_set_extra(self):
        return self._qm_basis_sets_extra

    def add_extra_basis_set(self, basis_set, elements, extra_args=None, position=None):
        basis_set = CustomBasisSet(basis_set, elements, extra_args=extra_args,
                                   position=position)
        self._qm_basis_sets_extra.append(basis_set)

    # MM model
    @property
    def mm_forcefield(self):
        return self._mm_forcefield

    @mm_forcefield.setter
    def mm_forcefield(self, value):
        if value not in MM_FORCEFIELDS['General']:
            raise ValueError('Forcefield {} not recognized'.format(value))
        self._mm_forcefield = value

    @property
    def mm_water_forcefield(self):
        return self._mm_water_forcefield

    @mm_water_forcefield.setter
    def mm_water_forcefield(self, value):
        if value not in MM_FORCEFIELDS['Water']:
            raise ValueError('Forcefield {} not recognized'.format(value))
        self._mm_water_forcefield = value

    @property
    def mm_forcefield_extra(self):
        if self._mm_forcefield_extra:
            return '\n'.join(map(str, self._mm_forcefield_extra))
        else:
            return ''

    def add_mm_forcefield(self, value):
       for file in value:
            if file.endswith('.frcmod') and os.path.isfile(file):
                self._mm_forcefield_extra.extend(import_from_frcmod(file))
            else:
                raise ValueError('Supply a .frcmod file to load new parameters')

    # System
    @property
    def system(self):
        first_line = '{} {}'.format(self.charge, self.multiplicity)
        if self.mm_forcefield and None not in (self.mm_charge, self.mm_multiplicity):
            first_line += ' {} {}'.format(self.mm_charge, self.mm_multiplicity)
        return '\n'.join([first_line] + map(str, self.atoms))

    @property
    def atoms(self):
        if not self._atoms:
            raise ValueError('System does not contain any atoms! Use .add_atom()!')
        if len(self._atoms) > 250000:
            raise ValueError('Max number of atoms is 250 000.')
        if self.mm_forcefield:  # using ONIOM!
            self._compute_oniom_links(self._atoms)
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        for a in value:
            self.add_atom(atom=a)

    def add_atom(self, atom=None, *args, **kwargs):
        if atom is None:
            atom = GaussianAtom(*args, **kwargs)
        if isinstance(atom, GaussianAtom):
            self._atoms.append(atom)
        else:
            raise TypeError('Provide either `atom` or options to construct '
                            'a GaussianAtom instance.')

    # ModRedundant restraints
    @property
    def restraints(self):
        restraints = sorted(self._restraints, key=lambda r: len(r.atoms), reverse=True)
        return '\n'.join(map(str, restraints))

    def add_restraint(self, restraint):
        if self.job != 'Opt':
            raise ValueError('Restraints (ModRedundant) can only be set with job=Opt')
        self._restraints.append(restraint)

    @property
    def charge(self):
        if self._charge is None:
            raise ValueError('Please set charge or use .compute_charge()')
        charge = self.compute_charge()
        if self._charge is not None and self._charge != charge:
            print('! Registered charge ({}) does not match '
                  'computed charge from atoms ({})'.format(self._charge, charge))
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = value

    @property
    def multiplicity(self):
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, value):
        if int(value) < 1:
            raise ValueError('Multiplicity must be greater than 0')
        self._multiplicity = int(value)

    @property
    def mm_charge(self):
        return self._mm_charge

    @mm_charge.setter
    def mm_charge(self, value):
        self._mm_charge = value

    @property
    def mm_multiplicity(self):
        return self._mm_multiplicity

    @mm_multiplicity.setter
    def mm_multiplicity(self, value):
        if int(value) < 1:
            raise ValueError('Multiplicity must be greater than 0')
        self._mm_multiplicity = int(value)

    ###########################
    def compute_charge(self):
        if self._atoms:
            return sum(a.charge for a in self.atoms if a.charge is not None)

    def compute_connectivity(self):
        lines = []
        seen = set()
        for atom in self.atoms:
            if not atom.neighbors:
                continue
            seen.add(atom)
            line = []
            for neighbor, bondorder in atom.neighbors:
                if neighbor not in seen and bondorder is not None:
                    line.append('{} {:.1f}'.format(neighbor.n, bondorder))
                    seen.add(neighbor)
            if line:
                lines.append('{} {}'.format(atom.n, ' '.join(line)))

        return '\n'.join(lines)

    @staticmethod
    def _compute_oniom_links(atoms):
        by_layers = defaultdict(list)
        linked = set()
        for atom in atoms:
            atom_layer = atom.oniom_layer
            can_be_link = []
            for neighbor, order in atom.neighbors:
                if neighbor.oniom_layer != atom_layer:
                    can_be_link.append(neighbor)
            n_links = len(can_be_link)
            if n_links == 1:
                atom.oniom_link = can_be_link[0]
                can_be_link[0].oniom_link = atom
            elif n_links > 1:
                raise ValueError('Atom {} can only have one layer link. Found: {}'.format(
                    atom.atom_spec, ', '.join([a.atom_spec for a in can_be_link])))


class GaussianAtom(object):

    """
    Object-oriented abstraction of a cartesian atom specification in
    Gaussian input files.

    It offers support for QM and MM atoms. More info is available at
    http://gaussian.com/molspec/.
    """

    def __init__(self, element, coordinates, n, atom_type=None, charge=None, freeze_code=None,
                 residue_number=None, residue_name=None, pdb_name=None, fragment=None,
                 iso=None, spin=None, zeff=None, qmom=None, nmagm=None, znuc=None,
                 oniom_layer=None, oniom_link=None, oniom_bonded=None, is_link=False,
                 oniom_scale_factors=None, geometry=None):
        self.n = n
        self._element = None
        self._coordinates = None
        self._atom_type = None
        self._charge = None
        self._freeze_code = None
        self._residue_number = None
        self._residue_name = None
        self._pdb_name = None
        self._fragment = None
        self._iso = None
        self._spin = None
        self._zeff = None
        self._qmom = None
        self._nmagm = None
        self._znuc = None
        self._oniom_layer = None
        self._oniom_link = None
        self._oniom_bonded = None
        self._oniom_scale_factors = None
        self._geometry = None
        self.is_link = bool(is_link)
        self._neighbors = []

        # Set and verify
        self.element = element
        self.coordinates = coordinates
        self.atom_type = atom_type
        self.charge = charge
        self.freeze_code = freeze_code
        self.residue_number = residue_number
        self.residue_name = residue_name
        self.pdb_name = pdb_name
        self.fragment = fragment
        self.iso = iso
        self.spin = spin
        self.zeff = zeff
        self.qmom = qmom
        self.nmagm = nmagm
        self.znuc = znuc
        self.oniom_layer = oniom_layer
        self.oniom_link = oniom_link
        self.oniom_bonded = oniom_bonded
        self.oniom_scale_factors = oniom_scale_factors
        self.geometry = geometry

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, value):
        if value is None:
            self._element = None
            return
        if not value or value[0] not in ascii_letters:
            raise ValueError('Element cannot be empty and must start with a letter. '
                             'Value provided: {}'.format(value))
        self._element = value

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        try:
            if len(value) == 3:
                self._coordinates = tuple(float(v) for v in value)
                return
        except (ValueError, TypeError):
            pass
        raise ValueError('Coordinates must be 3-tuple of float (x, y, z). '
                            'Value provided: {}'.format(value))

    @property
    def atom_type(self):
        return self._atom_type

    @atom_type.setter
    def atom_type(self, value):
        if value is None:
            self._atom_type = None
            return
        if not value:
            raise ValueError('Atom_type cannot be empty.')
        self._atom_type = value

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, value):
        if value is None:
            self._charge = None
            return
        if self.atom_type is None:
            raise ValueError('Setting charges requires setting atom_type beforehand')
        try:
            self._charge = float(value)
        except TypeError:
            raise ValueError('Charge values must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def freeze_code(self):
        return self._freeze_code

    @freeze_code.setter
    def freeze_code(self, value):
        if value is None:
            self._freeze_code = None
            return
        if isinstance(value, int) and value <= 0:
            self._freeze_code = value
            return
        raise ValueError('Freeze_code must be int and <= 0. '
                         'Value provided: {}'.format(value))

    @property
    def residue_number(self):
        return self._residue_number

    @residue_number.setter
    def residue_number(self, value):
        if value is None:
            self._residue_number = None
            return
        if not isinstance(value, int):
            raise ValueError('residue_number must be int. '
                             'Value provided: {}'.format(value))
        self._residue_number = value

    @property
    def residue_name(self):
        return self._residue_name

    @residue_name.setter
    def residue_name(self, value):
        if not value:
            self._residue_name = None
            return
        if isinstance(value, basestring):
            self._residue_name = value
            return
        raise ValueError('residue_name cannot be empty and must start with a letter. '
                         'Value provided: {}'.format(value))

    @property
    def pdb_name(self):
        return self._pdb_name

    @pdb_name.setter
    def pdb_name(self, value):
        if value is None:
            self._pdb_name = None
            return
        if isinstance(value, basestring) and value and value[0] in ascii_letters:
            self._pdb_name = value
            return
        raise ValueError('pdb_name cannot be empty and must start with a letter. '
                         'Value provided: {}'.format(value))

    @property
    def fragment(self):
        return self._fragment

    @fragment.setter
    def fragment(self, value):
        if value is None:
            self._fragment = None
            return
        if not isinstance(value, int):
            raise ValueError('spin values must be int. '
                             'Value provided: {}'.format(value))
        self._fragment = value

    @property
    def iso(self):
        return self._iso

    @iso.setter
    def iso(self, value):
        if value is None:
            self._iso = None
            return
        try:
            self._iso = float(value)
        except ValueError:
            raise ValueError('iso must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def spin(self):
        return self._spin

    @spin.setter
    def spin(self, value):
        if value is None:
            self._spin = None
            return
        try:
            self._spin = float(value)
        except ValueError:
            raise ValueError('spin must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def zeff(self):
        return self._zeff

    @zeff.setter
    def zeff(self, value):
        if value is None:
            self._zeff = None
            return
        try:
            self._zeff = float(value)
        except ValueError:
            raise ValueError('zeff must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def qmom(self):
        return self._qmom

    @qmom.setter
    def qmom(self, value):
        if value is None:
            self._qmom = None
            return
        try:
            self._qmom = float(value)
        except ValueError:
            raise ValueError('qmom must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def nmagm(self):
        return self._nmagm

    @nmagm.setter
    def nmagm(self, value):
        if value is None:
            self._nmagm = None
            return
        try:
            self._nmagm = float(value)
        except ValueError:
            raise ValueError('nmagm must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def znuc(self):
        return self._znuc

    @znuc.setter
    def znuc(self, value):
        if value is None:
            self._znuc = None
            return
        try:
            self._znuc = float(value)
        except ValueError:
            raise ValueError('znuc must be float or float-like. '
                             'Value provided: {}'.format(value))

    @property
    def oniom_layer(self):
        return self._oniom_layer

    @oniom_layer.setter
    def oniom_layer(self, value):
        if value is None:
            self._oniom_layer = None
            return
        if value.upper() in ('H', 'M', 'L'):
            self._oniom_layer = value.upper()
            return
        raise ValueError('oniom_layer must be H, M, or L. '
                         'Value provided: {}'.format(value))

    @property
    def oniom_link(self):
        return self._oniom_link

    @oniom_link.setter
    def oniom_link(self, value):
        if isinstance(value, GaussianAtom) or value is None:
            self._oniom_link = value
            return
        raise ValueError('oniom_link must be a GaussianAtom instance. '
                            'Value provided: {}'.format(value))

    @property
    def oniom_bonded(self):
        return self._oniom_bonded

    @oniom_bonded.setter
    def oniom_bonded(self, value):
        if value is None:
            self._oniom_bonded = None
            return
        try:
            self._oniom_bonded = int(value)
        except ValueError:
            raise ValueError('oniom_bonded must be int or int-like. '
                             'Value provided: {}'.format(value))

    @property
    def oniom_scale_factors(self):
        return self._oniom_scale_factors

    @oniom_scale_factors.setter
    def oniom_scale_factors(self, value):
        if value is None:
            self._oniom_scale_factors = value
            return
        try:
            if len(value) <= 3:
                self._oniom_scale_factors = tuple(float(v) for v in value)
                return
        except (ValueError, TypeError):
            pass
        raise ValueError('oniom_scale_factors must be tuple of float, three values max. '
                         'Value provided: {}'.format(value))

    @property
    def geometry(self):
        return self._geometry

    @geometry.setter
    def geometry(self, value):
        if value is None:
            self._geometry = value
            return
        try:
            self._geometry = int(value)
        except (ValueError, TypeError):
            raise ValueError('geometry must be int or int-like. '
                             'Value provided: {}'.format(value))

    @property
    def neighbors(self):
        return self._neighbors

    def add_neighbor(self, neighbor, bondorder=1.0):
        self._neighbors.append((neighbor, bondorder))


    #--------------------------------------
    # Helper methods
    #--------------------------------------
    @property
    def atom_spec(self):
        """
        Summarize atom information in a single string
        """
        line = [self.element]
        if self.atom_type:
            line.append('-{}'.format(self.atom_type))
        if self.charge is not None:
            line.append('-{}'.format(self.charge))

        return ''.join(line)

    @property
    def keywords(self):
        """
        Get a dict with all the keywords
        """
        keywords = ['residue_number', 'residue_name', 'pdb_name', 'fragment',
                    'iso', 'spin', 'zeff', 'qmom', 'nmagm', 'znuc']
        return {k: getattr(self, k) for k in keywords}

    @property
    def keywords_spec(self):
        """
        Summarize all keywords in a single string
        """
        keywords = []
        if self.residue_number is not None:
            keywords.append("{}={}".format("RESNum", self.residue_number))
        if self.residue_name is not None:
            keywords.append("{}={}".format("RESName", self.residue_name))
        if self.pdb_name is not None:
            keywords.append("{}={}".format("PDBName", self.pdb_name))
        if self.fragment is not None:
            keywords.append("{}={}".format("Fragment", self.fragment))
        if self.iso is not None:
            keywords.append("{}={}".format("Iso", self.iso))
        if self.spin is not None:
            keywords.append("{}={}".format("Spin", self.spin))
        if self.zeff is not None:
            keywords.append("{}={}".format("ZEff", self.zeff))
        if self.qmom is not None:
            keywords.append("{}={}".format("QMom", self.qmom))
        if self.nmagm is not None:
            keywords.append("{}={}".format("NMagM", self.nmagm))
        if self.znuc is not None:
            keywords.append("{}={}".format("ZNuc", self.znuc))

        if keywords:
            return '({})'.format(','.join(keywords))

    @property
    def coordinates_spec(self):
        if self.coordinates is not None:
            return '{:>12.6f} {:>12.6f} {:>12.6f}'.format(*self.coordinates)

    def __str__(self):
        # Atom element, name, charge
        line = [self.atom_spec]

        # Atom keywords
        keywords = self.keywords_spec
        if keywords:
            line.append(keywords)

        # Freeze code (rigidity)
        if self.freeze_code is not None:
            line.append(' {}'.format(self.freeze_code))

        # Coordinates
        coords = self.coordinates_spec
        if coords:
            line.append(' {}'.format(coords))

        # ONIOM config
        if self.oniom_layer:
            line.append(' {}'.format(self.oniom_layer))
        link = self.oniom_link
        if link:
            line.append(' {}'.format(link.n))
            if link.oniom_bonded:
                line.append(' {}'.format(link.oniom_bonded))
            if link.oniom_scale_factors:
                line.append(' {}'.format(' '.join(link.oniom_scale_factors)))

        return ''.join(map(str, line))


class CustomBasisSet(object):

    """
    Object-oriented abstraction of an extra basis set for Gaussian
    input files.

    It includes support for ebsel local databases as well as
    online databases from BSE and Cosmo.
    """

    def __init__(self, basis_set, elements, position=0, extra_args="",):
        if basis_set.rstrip('+*') not in QM_BASIS_SETS:
            raise ValueError('Basis set {} not recognized'.format(basis_set))
        if basis_set in QM_BASIS_SETS_WITH_ARGS and not extra_args:
            raise ValueError('{} requires extra_args'.format(basis_set))
        if position and len(elements) > 1:
            raise ValueError('Position can only be set if one element is suppled')
        self.basis_set = basis_set
        self.position = position
        self.extra_args = extra_args
        self.elements = elements

    def __str__(self):
        lines = [' '.join(self.elements), self.position]
        lines.append(self.basis_set)
        lines.append(self.extra_args)
        lines.append(' ****')
        return '\n'.join(lines)

    @classmethod
    def from_database(cls, basis_set, element, database='cosmologic-services'):
        """
        Import basis_set from online database.

        TODO:
            Support https://bse.pnl.gov/bse/portal/
            Support http://cosmologic-services.de/basis-sets/basissets.php
        """
        try:
            import requests
            from bs4 import BeautifulSoup
        except ImportError:
            raise ImportError('You need to install requests and beautifulsoup4 '
                              'to import basis sets from online databases.')
        api = {
            'cosmologic-services': {
                'header': '! Downloaded from http://cosmologic-services.de/basis-sets/basissets.php',
                'url': 'http://cosmologic-services.de/basis-sets/getbasis.php',
                'bs_parser': 'html.parser',
                'parser': cls._cosmo_parser,
                'data': {'basis': basis_set, element: element,
                         'kind': 'Basis', 'format': 'Gaussian'}
            },
            'comp.chem.umn.edu': {
                'header': '! Downloaded from http://comp.chem.umn.edu/basissets/basis.cgi',
                'url': 'http://comp.chem.umn.edu/basissets/basis.cgi',
                'bs_parser': 'html.parser',
                'parser': cls._umn_parser,
                'data': {'basis_list': basis_set, 'element_list': element,
                         'format_list': 'Gaussian'}
            }
        }
        try:
            db = api[database]
        except KeyError:
            raise ValueError('Database {} not supported'.format(database))
        try:
            r = requests.post(db['url'], db['data'], timeout=10)
        except requests.exceptions.RequestException:
            raise
        if not r.ok:
            raise ValueError('Could not retrieve data from {}'.format(database))
        html = BeautifulSoup(r.content, db['bs_parser'])
        return '\n'.join([db['header'], db['parser'](html)])

    @staticmethod
    def _cosmo_parser(html):
        return html.find('pre').text

    @staticmethod
    def _umn_parser(html):
        text = str(html.find('body').contents[4])
        return text.replace('<br>', '\n').replace('</br>', '').replace(u'\xa0', u' ')

    @classmethod
    def from_bse(cls, basis_set, *elements):
        try:
            from ebsel.EMSL_local import EMSL_local
        except ImportError:
            raise ImportError('Access to Basis Set Exchange db requires ebsel package')

        db = EMSL_local(fmt="g94")
        try:
            basis = db.get_basis(basis_set, elements=elements)
        except UnboundLocalError:
            if basis_set not in [b for (b, d) in db.get_available_basis_sets()]:
                raise ValueError('Basis set {} not recognized'.format(basis_set))
            supported_elements = db.get_available_elements(basis_set)
            for element in elements:
                if element not in supported_elements:
                    raise ValueError('Basis set {} does not support element {}'.format(basis_set, element))
        else:
            # prepend a '-' to each element to prevent Gaussian errors if atom not present in system
            return '\n'.join([b.replace('****\n', '****\n-') for b in basis])



class ModRedundantRestraint(object):

    TYPES = {1: 'X', 2: 'B', 3: 'A', 4: 'D'}
    OPERATIONS = sorted('A F B K R D H S'.split())

    def __init__(self, atoms, operation, diag_elem=None, nsteps=None, stepsize=None,
                 rtype=None, min_=None, max_=None):
        if not (1 <= len(atoms) <= 4):
            raise ValueError('Provide between 1 and 4 atoms or wildcards (*).')
        if not all(a.isdigit() or a.strip() in ('', '*') for a in atoms):
            raise ValueError('atoms must be int or *')
        if operation not in self.OPERATIONS:
            raise ValueError('operation must be one of <A F B K R D H S>')
        if operation == 'S' and None in (nsteps, stepsize):
            raise ValueError('if operation=S, nsteps and stepsize must be set')
        if operation == 'H' and diag_elem is None:
            raise ValueError('if operation=H, diag_elemt must be set')

        self.atoms = [a.strip() for a in atoms if a.strip()]
        self.operation = operation
        self.diag_elem = diag_elem
        self.nsteps = nsteps
        self.stepsize = stepsize
        self.rtype = rtype
        self.min_ = min_
        self.max_ = max_

    def __str__(self):
        kwargs = dict(rtype=self.restraint_type,
                      atoms=' '.join(map(str, self.atoms)),
                      operation=self.operation,
                      args=' '.join(map(str, self._args)),
                      minmax=' '.join(map(str, self._args)))
        return '{rtype} {atoms} {operation} {args}'.format(**kwargs)

    @property
    def restraint_type(self):
        if self.rtype is not None:
            return self.rtype
        return self.TYPES[len(self.atoms)]

    @property
    def minmax(self):
        s = []
        if self.max_ is not None:
            if self.min_ is not None:
                return self.min_, self.max_
            return self.max_,

    @property
    def _args(self):
        if self.operation == 'S':
            return self.nsteps, self.stepsize
        elif self.operation == 'H':
            return self.diag_elem,
        else:
            return ()

    @property
    def atom1(self):
        try:
            return self.atoms[0]
        except IndexError:
            return ''

    @property
    def atom2(self):
        try:
            return self.atoms[1]
        except IndexError:
            return ''

    @property
    def atom3(self):
        try:
            return self.atoms[2]
        except IndexError:
            return ''

    @property
    def atom4(self):
        try:
            return self.atoms[3]
        except IndexError:
            return ''

class MmExtraDefinition(object):

    TYPES = ('VDW', 'HRMSTR1', 'HRMBND1', 'DREITRS')

    def __init__(self, dtype, mmtype, mmtype2=None, mmtype3=None, mmtype4=None,
                bond_length=None, well_depth=None, force_constant=None, 
                angle=None, idivf=None, pk=None, phase=None, pn=None):
        if dtype.upper() not in self.TYPES:
            raise ValueError('not valid definition type')

        self.dtype = dtype
        self.mmtype = mmtype
        self.mmtype2 = mmtype2
        self.mmtype3 = mmtype3
        self.mmtype4 = mmtype4
        self.bond_length = bond_length
        self.well_depth = well_depth
        self.force_constant = force_constant
        self.angle = angle
        self.idivf = idivf
        self.pk = pk
        self.phase = phase
        self.pn = pn

    def __str__(self):
        kwargs = dict(dtype=self.dtype,
                      mmtype=self.mmtype,
                      args=' '.join(map(str, self._args)))
        return '{dtype} {mmtype} {args}'.format(**kwargs)

    @property
    def _args(self):
        if self.dtype == 'VDW':
            return self.bond_length, self.well_depth
        elif self.dtype == 'HrmStr1':
            return self.mmtype2, self.force_constant, self.bond_length
        elif self.dtype == 'HrmBnd1':
            return self.mmtype2, self.mmtype3, self.force_constant, self.angle
        elif self.dtype == 'DreiTrs':
            return self.mmtype2, self.mmtype3, self.mmtype4, float(self.pk)*2, self.phase, self.pn, float(self.idivf)/2
        else:
            return ()

def import_from_frcmod(path):
    SECTIONS = ('NONB', 'BOND', 'ANGL', 'DIHE')
    section = None
    mm_definitions = []
    with open(path) as inf:
        for line in inf:
            if line.strip('\n').upper() in SECTIONS:
                section = line.strip('\n').upper()
            elif not line.strip():
                section = None
            elif section:
                    mm_definitions.append(_create_mm_definition(line, section))
    return mm_definitions

def _create_mm_definition(line, section):
    definition = None
    if section == 'NONB':
        args = line.split()
        #It is supposed that NONB section follows the 10B card type (ambermd.org/formats.html)
        definition = MmExtraDefinition('VDW', args[0], bond_length=args[1], well_depth=args[2])
    elif section == 'BOND':
        mm1, mm2 = line[:5].split('-')
        args = line[5:].split()
        #It is supposed that BOND section follows the 4 card type (ambermd.org/formats.html)
        definition = MmExtraDefinition('HrmStr1', mm1, mmtype2=mm2, force_constant=args[0], bond_length=args[1])
    elif section == 'ANGL':
        mm1, mm2, mm3 = line[:8].split('-')
        args = line[8:].split()
        #It is supposed that ANGL section follows the 5 card type (ambermd.org/formats.html)
        definition = MmExtraDefinition('HrmBnd1', mm1, mmtype2=mm2, mmtype3=mm3, force_constant=args[0], angle=args[1])
    elif section == 'DIHE':
        mm1, mm2, mm3, mm4 = line[:11].split('-')
        if mm1 == 'X ':
            mm1 = '* '
        if mm4 == 'X ':
            mm4 = '* '
        args = line[11:].split()
        definition = MmExtraDefinition('DreiTrs', mm1, mmtype2=mm2, mmtype3=mm3, mmtype4=mm4, idivf=args[0], pk=args[1], phase=args[2], pn=args[3])
    return definition

if __name__ == '__main__':
    atom = GaussianAtom(element='C', coordinates=(10, 10, 10), n=1, atom_type='CT',
                        charge=1.0, residue_number=1, oniom_layer='H')
    print(atom)
