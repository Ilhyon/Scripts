#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os

def main(mainPath):
    for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files
			inputfile = directory + '/' + filename
            filtered = []
            filename = mainPath+filename
            with open(filename) as f:
                lines = f.read().splitlines()
                for l in lines:
                    l = l.rstrip()
                    words = l.split('\t')
                    if words[2] >= 30:
                        filtered.append(l)
            output = open(mainPath+'Filter_'+filename, "w")
            output.write("\n".join(filtered))
            output.close()

if __name__ == '__main__':
    path = '/home/vana2406/scratch/rG4seq/GSE77282_RAW/'
    main(path)
