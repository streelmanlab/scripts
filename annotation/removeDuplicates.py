#!/usr/bin/env python3
"""

  Author: George Gruenhagen
"""
# --- Imports ---
import argparse
parser = argparse.ArgumentParser()

def parseArgs():
    global parser
    parser.add_argument("-a", help = "File to remove duplicates from")
    parser.add_argument("-o", help = "Output file to write to")

    args = parser.parse_args()
    return args.a, args.o

def readFile(inputFile):
    with open(inputFile, 'r') as inputFile:
        lines = inputFile.readlines()

    print("Number of lines: " + str(len(lines)))
    lines = list(set(lines))
    print("Number of lines after removing duplicates: " + str(len(lines)))
    return lines

def writeFile(lines, outputFile):
    with open(outputFile, 'w') as outputFile:
        for line in lines:
            outputFile.write(line)

def main():
    inputFile, outputFile = parseArgs()
    lines = readFile(inputFile)
    writeFile(lines, outputFile)

main()
