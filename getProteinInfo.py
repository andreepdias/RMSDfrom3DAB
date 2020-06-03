#!/usr/bin/python3

import re

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
    proteinsInfos = []

    for p in proteins:
        proteinInfo = getProteinInfo(p)
        proteinsInfos.append([p, proteinInfo])

    for protein in proteinsInfos:
        name = protein[0][0]
        print(name + ':')

        for c in protein[1]:
            chain = c[0]
            count = c[1]
            print(chain, count, end='')

            if int(protein[0][1]) == int(count):
                print(' ok')
            else:
                print(' not ok')
        print()

if __name__ == '__main__':
    main()