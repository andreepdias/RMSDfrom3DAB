#!/usr/bin/python3

import re

def readProteinsNames():
    proteinsNames = []
    filePath = 'proteinsNames.txt'

    with open(filePath) as f:
        line = f.readline()

        while line:
            proteinsNames.append(line.rstrip('\n'))
            line = f.readline()
            
    return proteinsNames

def getProteinInfo(proteinName):
    filePath = 'proteinsCIF/' + proteinName + '.cif'

    proteinInfo = {}

    with open(filePath) as f:
        line = f.readline()

        while line:
            lineSplit = [x for x in re.split(r'\s{1,}', line) if x]

            header = lineSplit[0]
            
            if header == 'ATOM':
                model = lineSplit[20]

                if int(model) == 1:
                    atom = lineSplit[3]
                    chain = lineSplit[6]

                    if atom == 'CA':
                        if not chain in proteinInfo:
                            proteinInfo[chain] = 0
                        proteinInfo[chain] += 1

            line = f.readline()

    return [ [k, v] for k, v in proteinInfo.items() ]

def main():
    proteinsNames = readProteinsNames()
    proteinsInfos = []


    for p in proteinsNames:
        proteinInfo = getProteinInfo(p)

        proteinsInfos.append([p, proteinInfo])

    for protein in proteinsInfos:
        name = protein[0]
        print(name + ':')

        for c in protein[1]:
            chain = c[0]
            count = c[1]
            print(chain, count)
        print()

if __name__ == '__main__':
    main()