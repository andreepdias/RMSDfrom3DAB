#!/usr/bin/python3

import re

def sequenceTo3DAB(sequence):
    new = ''

    for c in sequence:
        if c == 'I' or c == 'V' or c == 'L' or c == 'P' or c == 'C' or c == 'M' or c == 'A' or c == 'G':
            new += 'A'
        elif c == 'D' or c == 'E' or c == 'F' or c == 'H' or c == 'K' or c == 'N' or c == 'Q' or c == 'R' or c == 'S' or c == 'T' or c == 'W' or c == 'Y':
            new += 'B'
        else:
            new += c
    return new

def readProteins():
    proteins = []
    filePath = 'proteinsNames.txt'

    with open(filePath) as f:
        line = f.readline()

        while line:
            lineSplit = [x for x in re.split(r'\s{1,}', line) if x]
            proteins.append(lineSplit)

            line = f.readline()
    return proteins

def readChains():
    chains = {}
    filePath = 'proteinsChains.txt'

    with open(filePath) as f:
        line = f.readline()

        while line:
            lineSplit = [x for x in re.split(r'\s{1,}', line) if x]
            chains[lineSplit[0]] = [lineSplit[1], lineSplit[2]]

            line = f.readline()
    return chains

def getProteinInfo(protein):
    filePath = 'proteinsCIF/' + protein[0] + '.cif'

    proteinInfo = {}

    with open(filePath) as f:
        line = f.readline()

        while line:
            lineSplit = [x for x in re.split(r'\s{1,}', line) if x]

            header = lineSplit[0]
            
            if header == 'ATOM' or header == 'HETATM':
                model = lineSplit[20]

                if int(model) == 1:
                    atom = lineSplit[3]
                    chain = lineSplit[6]

                    if (header == 'ATOM' and atom == 'CA') or (header == 'HETATM' and atom == 'C') or (header == 'HETATM' and atom == 'N'):
                        if not chain in proteinInfo:
                            proteinInfo[chain] = 0
                        proteinInfo[chain] += 1

            line = f.readline()

    return [ [k, v] for k, v in proteinInfo.items() ]

def main():

    proteins = readProteins()
    chains = readChains()

    proteinsInfos = []

    for p in proteins:
        proteinInfo = getProteinInfo(p)
        proteinsInfos.append([p, proteinInfo])

    for protein in proteinsInfos:
        name = protein[0][0]
        print(name + ':')

        for c in protein[1]:
            chainLetter = c[0]
            count = c[1]
            count3DAB = int(protein[0][1])

            print(chainLetter, count, count3DAB, end='')

            if int(protein[0][1]) == int(count):
                print(' ok')
            else:
                print(' not ok')
        
        chainPDB = chains[name][1]
        chainAB = chains[name][0]

        print(chainAB)
        print(sequenceTo3DAB(chainPDB))
        print(chainPDB)

        print()

if __name__ == '__main__':
    main()