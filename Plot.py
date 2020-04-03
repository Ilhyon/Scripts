#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import argparse
import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt

df = pd.read_csv('/home/anais/Documents/Data/DemandeFondMa2020/AllEvent.tsv', sep='\t')

groups = df.groupby('Event')
for name, group in groups:
    # width of the bars
    barWidth = 0.3

    # Choose the height of the blue bars
    bars1 = group[group.Significant == 1].Mean

    # Choose the height of the cyan bars
    bars2 = group[group.Significant == 0].Mean

    # Choose the height of the error bars (bars1)
    yer1 = group[group.Significant == 1].Sd/2

    # Choose the height of the error bars (bars2)
    yer2 = group[group.Significant == 0].Sd/2

    # The x position of bars
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    # Create blue bars
    plt.bar(r1, bars1, width = barWidth, color = 'darkred', edgecolor = 'black', yerr=yer1, capsize=7, label='Significant'+name)

    # Create cyan bars
    plt.bar(r2, bars2, width = barWidth, color = 'black', edgecolor = 'black', yerr=yer2, capsize=7, label='nonSignificant')

    # general layout
    plt.xticks([r + barWidth for r in range(len(bars1))], ['cond_A', 'cond_B', 'cond_C'])
    plt.ylabel('Average densities (nb pG4/kb)')
    plt.legend()

    # Show graphic
    plt.savefig("Figure_"+name+".svg")
