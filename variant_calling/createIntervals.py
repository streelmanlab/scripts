#!/usr/bin/env python3
"""

  Author: George Gruenhagen
"""
# --- Imports ---
import argparse
parser = argparse.ArgumentParser()

def parseArgs():
    global parser
    parser.add_argument("-a", help = "LG Accession Number and Length, tab delimited")
    parser.add_argument("-o", help = "Name of output file")
    parser.add_argument("-s", "--size", help = "Size of Intervals", nargs = '?', const = 1000000, type = int, default = 1000000)

    args = parser.parse_args()
    return args.a, args.o, args.size

def breakup(acNo, length, size):
    intervals = []
    for i in range(len(acNo)):
        if int(length[i]) > size:
            newIntervals = []
            # previousNum = 0
            for j in range(0, int(length[i]), size):
                if j > 1:
                    newIntervals.append(str(previousNum + 1) + "\t" + str(j))
                previousNum = j
            newIntervals.append(str(previousNum + 1) + "\t" + length[i])
            intervals.append(newIntervals)

    return intervals

def read(inputFile):
    acNum = []
    length = []
    with open(inputFile, 'r') as input:
        for line in input:
            lineSplit = line.split()
            acNum.append(lineSplit[0])
            length.append(lineSplit[1])

    return acNum, length

def write(acNo, intervals, outputFile):
    with open(outputFile, 'w+') as output:
        for i in range(len(acNo)):
            for interval in intervals[i]:
                output.write(acNo[i] + "\t" + interval + "\n")

def main():
    input, outputFile, size = parseArgs()
    size = int(size)
    
    acNo, length = read(input)
    intervals = breakup(acNo, length, size)
    write(acNo, intervals, outputFile)


main()
