#!/bin/env python

import os, sys
if sys.version.startswith('2'):
    print("Should use python 3")
    exit(-1)

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import Figures
import Colors
import shlex, argparse

matplotlib.use('Agg')

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lower", 
        dest="lb",
        type=float,
        default=0.0,
        help="The lower bound")
parser.add_argument("-u", "--upper", 
        dest="ub",
        type=float,
        default=0.05,
        help="The upper bound")
parser.add_argument('infile', 
    help="The input file produced by shapemapper_mutation_counter")
parser.add_argument('outPDF',
    help="The output PDF file")
args = parser.parse_args()

table = pd.read_csv(args.infile, sep="\t")

deletion    = np.array(table.iloc[:, 0:4].sum(axis=1).tolist())
insertion   = np.array(table.iloc[:, 4:8].sum(axis=1).tolist())
mistmatch   = np.array(table.iloc[:, 9:21].sum(axis=1).tolist())
multidel    = np.array(table.iloc[:, 21].tolist())
multiins    = np.array(table.iloc[:, 22].tolist())
multimis    = np.array(table.iloc[:, 23].tolist())
compdel     = np.array(table.iloc[:, 24].tolist())
compins     = np.array(table.iloc[:, 25].tolist())

effdepth    = np.array(table.loc[:, 'effective_depth'].tolist())

stackedBars = []
for i in range(table.shape[0]):
    bar = [ 
        deletion[i]/(effdepth[i]+1),
        insertion[i]/(effdepth[i]+1),
        mistmatch[i]/(effdepth[i]+1),
        multidel[i]/(effdepth[i]+1),
        multiins[i]/(effdepth[i]+1),
        multimis[i]/(effdepth[i]+1),
        compdel[i]/(effdepth[i]+1),
        compins[i]/(effdepth[i]+1)
    ]
    stackedBars.append(bar)

stackedLabels = ['del','ins','mis','multidel','multiins','multimis','compdel','compins']
stackedColors = [
    Colors.RGB['green'],
    Colors.RGB['red'],
    Colors.RGB['blue'],
    Colors.RGB['yellow'],
    Colors.RGB['teal'],
    Colors.RGB['amber'],
    Colors.RGB['gray'],
    Colors.RGB['deep_orange']
]
barLabels = range(1, table.shape[0]+1)

if table.shape[0]>100:
    width = (table.shape[0]-100)//50 + 12
else:
    width = 12

print(f"Width={width}")
plt.figure(figsize=(width,6))
Figures.stackedBarPlot( stackedBars, stackedLabels, barLabels, stackedColors)
plt.ylim(args.lb, args.ub)
plt.ylabel("mutation rate")
plt.xticks(range(1, table.shape[0], 10), range(1, table.shape[0], 10))
plt.savefig(args.outPDF)
plt.close()

