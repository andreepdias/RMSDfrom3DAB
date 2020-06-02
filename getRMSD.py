#!/usr/bin/python3
# align 1edp, 1edp_new, cycles = 0, transform = 0

import __main__
__main__.pymol_argv = [ 'pymol', '-Qc'] # Quiet and no GUI

import sys
import re
import math
import pymol

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def calculate3DAB(angles, chain):
    proteinLength = int((len(angles) + 5) / 2)

    aminoacidPosition = [None] * proteinLength * 3

    aminoacidPosition[0] = 0.0
    aminoacidPosition[0 + proteinLength] = 0.0
    aminoacidPosition[0 + proteinLength * 2] = 0.0

    aminoacidPosition[1] = 0.0
    aminoacidPosition[1 + proteinLength] = 1.0
    aminoacidPosition[1 + proteinLength * 2] = 0.0

    aminoacidPosition[2] = math.cos(angles[0])
    aminoacidPosition[2 + proteinLength] = math.sin(angles[0]) + 1.0
    aminoacidPosition[2 + proteinLength * 2] = 0.0

    for i in range(3, proteinLength):
        aminoacidPosition[i] = aminoacidPosition[i - 1] + math.cos(angles[i - 2]) * math.cos(angles[i + proteinLength - 5])
        aminoacidPosition[i + proteinLength] = aminoacidPosition[i - 1 + proteinLength] + math.sin(angles[i - 2]) * math.cos(angles[i + proteinLength - 5])
        aminoacidPosition[i + proteinLength * 2] = aminoacidPosition[i - 1 + proteinLength * 2] + math.sin(angles[i + proteinLength - 5])
    
    obj = 0.0
    for i in range(proteinLength - 2):
        obj += (1.0 - math.cos(angles[i])) / 4.0

    for i in range(proteinLength - 2):
        for j in range(i + 2, proteinLength):
            c = -0.5
            if chain[i] == 'A' and chain[j] == 'A':
                c = 1.0
            elif chain[i] == 'B' and chain[j] == 'B':
                c = 0.5
            
            dx = aminoacidPosition[i] - aminoacidPosition[j]
            dy = aminoacidPosition[i + proteinLength] - aminoacidPosition[j + proteinLength]
            dz = aminoacidPosition[i + proteinLength * 2] - aminoacidPosition[j + proteinLength * 2]
            d = math.sqrt( (dx * dx) + (dy * dy) + (dz * dz) )

            obj += 4.0 * (1.0 / math.pow(d, 12.0) - c / math.pow(d, 6.0))
    
    return obj

def calculateNewCoordinates(angles3DAB, lengthsBetweenCAPDB):
    proteinLength = int((len(angles3DAB) + 5) / 2)
    coordinates = [ [[None] for i in range(3)] for j in range(proteinLength) ]

    coordinates[0][0] = 0.0
    coordinates[0][1] = 0.0
    coordinates[0][2] = 0.0

    coordinates[1][0] = 0.0
    coordinates[1][1] = 1.0 * lengthsBetweenCAPDB[0]
    coordinates[1][2] = 0.0

    coordinates[2][0] = math.cos(angles3DAB[0]) * lengthsBetweenCAPDB[1]
    coordinates[2][1] = coordinates[1][1] + math.sin(angles3DAB[0]) * lengthsBetweenCAPDB[1]
    coordinates[2][2] = 0.0

    for i in range(3, proteinLength):

        coordinates[i][0] = coordinates[i - 1][0] + math.cos(angles3DAB[i - 2]) * math.cos(angles3DAB[i + proteinLength - 5]) * lengthsBetweenCAPDB[i - 1]
        coordinates[i][1] = coordinates[i - 1][1] + math.sin(angles3DAB[i - 2]) * math.cos(angles3DAB[i + proteinLength - 5]) * lengthsBetweenCAPDB[i - 1]
        coordinates[i][2] = coordinates[i - 1][2] + math.sin(angles3DAB[i + proteinLength - 5]) * lengthsBetweenCAPDB[i - 1]
    
    return coordinates

def lengthsFromCoordinates(coordinates):
    distances = []

    for i in range(1, len(coordinates)):
        dx = coordinates[i][0] - coordinates[i - 1][0]
        dy = coordinates[i][1] - coordinates[i - 1][1]
        dz = coordinates[i][2] - coordinates[i - 1][2]
        d = math.sqrt( (dx * dx) + (dy * dy) + (dz * dz) )
        distances.append(d)
    
    return distances

def coordinatesFromAngles(angles):
    proteinLength = int((len(angles) + 5) / 2)

    coordinates = [ [[None] for i in range(3)] for j in range(proteinLength) ]

    coordinates[0][0] = 0.0
    coordinates[0][1] = 0.0
    coordinates[0][2] = 0.0

    coordinates[1][0] = 0.0
    coordinates[1][1] = 1.0
    coordinates[1][2] = 0.0

    coordinates[2][0] = math.cos(angles[0])
    coordinates[2][1] = math.sin(angles[0]) + 1.0
    coordinates[2][2] = 0.0

    for i in range(3, proteinLength):
        coordinates[i][0] = coordinates[i - 1][0] + math.cos(angles[i - 2]) * math.cos(angles[i + proteinLength - 5])
        coordinates[i][1] = coordinates[i - 1][1] + math.sin(angles[i - 2]) * math.cos(angles[i + proteinLength - 5])
        coordinates[i][2] = coordinates[i - 1][2] + math.sin(angles[i + proteinLength - 5])
    
    return coordinates

def readAngles3DAB(proteinName):
    filePath = 'proteinsDE/' + proteinName + '.3dab'

    angles = []

    with open(filePath) as f:
        line = f.readline().split(' ')
        angles = [(float(x) * math.pi / 180.0) for x in line] 
    
    return angles

def readChain3DAB(proteinName):
    filePath = 'proteinsDE/' + proteinName + '.3dab'

    chain = []

    with open(filePath) as f:
        line = f.readline()
        line = f.readline()
        chain = [ x for x in line ] 
    
    return chain

def readPDBFile(proteinName, chainLetter):
    filePDBpath = 'proteinsPDB/' + proteinName + '.pdb'
    
    coordinates = []

    with open(filePDBpath) as f:
        line = f.readline()

        while line:
            if 'ENDMDL' in line:
                break

            line = [x for x in re.split(r'\s{1,}', line) if x]
            header = line[0]
        
            if ((header == 'ATOM') or (header == 'HETATM')):
                atom = line[2]
                chain = line[4]
                if ((header == 'ATOM' and atom == 'CA' and chain == chainLetter) or (header == 'HETATM' and atom == 'N' and chain == chainLetter)):
                    coordinates.append([float(line[6]), float(line[7]), float(line[8])])

            line = f.readline()
    return coordinates

def plot3D(coordinates):

    X = [ x[0] for x in coordinates ]
    Y = [ x[1] for x in coordinates ]
    Z = [ x[2] for x in coordinates ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(X, Y, Z, c='r', marker='o')

    for i in range(1, len(coordinates)):
        ax.plot([X[i - 1], X[i]], [Y[i - 1], Y[i]], [Z[i - 1], Z[i]], c='r', linewidth=0.75)
    

    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')

    plt.show()

def printList(list, msg):
    print(msg)
    for i in range(len(list)):
        print(i, list[i])
    print()

def print3DMatrix(matrix, msg):
    print(msg)
    for i in range(len(matrix)):
        print('{0:0.3f} {1:0.3f} {2:0.3f}'.format(matrix[i][0], matrix[i][1], matrix[i][2]))
    print()

def getRMSD(proteinName):
    proteinPDB = proteinName
    proteinNEW = proteinName + '_new'

    proteinPDBPath = 'proteinsCIF/' + proteinPDB + '.cif'
    proteinNEWPath = 'proteinsNEW/' + proteinNEW + '.cif'

    pymol.cmd.load(proteinPDBPath)
    pymol.cmd.load(proteinNEWPath)

    x = pymol.cmd.align(proteinPDB, proteinNEW, cycles = 0, transform = 0)

    return float(x[0])


def readCIFFile(proteinName, chainLetter, middlePDBFile):
    fileCIFpath = 'proteinsCIF/' + proteinName + '.cif'
    
    headerPDBFile = ''
    tailPDBFile = ''
    coordinates = []

    buildingHeader = True

    with open(fileCIFpath) as f:

        line = f.readline()
        headerPDBFile += line

        while line:

            lineSplit = [x for x in re.split(r'\s{1,}', line) if x]
            header = lineSplit[0]
        
            if header == 'ATOM':
                buildingHeader = False

                atom = lineSplit[3]
                chain = lineSplit[6]

                if atom == 'CA' and chain == chainLetter:
                    coordinates.append([float(lineSplit[10]), float(lineSplit[11]), float(lineSplit[12])])

                    middlePDBFile.append(lineSplit)
            else:
                if buildingHeader:
                    headerPDBFile += line
                else:
                    tailPDBFile += line


            line = f.readline()

    return [ coordinates, headerPDBFile, tailPDBFile ]


def buildPDBFile(coordinates, header, tail, middle, proteinName):

    ident = [' ', '\t', ' ', '\t', ' ', ' ', ' ', ' ', '\t', ' ', ' ', ' ', '\t', ' ', ' ', ' ', '\t', ' ', ' ', '\t']

    fileText = ''
    fileText += header

    for k in range(len(middle)):
        line = str(middle[k][0])

        for i in range(1, 10):
            line += str(ident[i - 1])
            line += str(middle[k][i])

        for i in range(10, 13):
            line += str(ident[i - 1])
            line += '{:.2f}'.format(float(coordinates[k][i - 10]))

        for i in range(13, len(middle[k])):
            line += str(ident[i - 1])
            line += str(middle[k][i])

        line += '\n'
        fileText += line
    
    fileText += tail

    filePath = 'proteinsNEW/'  + proteinName + '_new.cif'
    file = open(filePath,'w')
    file.write(fileText)
    file.close()

def main():

    proteinsFilePath = 'proteinsNames.txt'
    proteins = []

    with open(proteinsFilePath) as f:
        line = f.readline()

        while line:
            lineSplit = [x for x in re.split(r'\s{1,}', line) if x] 

            name = lineSplit[0]
            chain = ''

            if len(lineSplit) > 1:
                chain = lineSplit[1]
            else:
                chain = 'A'
            
            proteins.append([name, chain])

            line = f.readline()

    for protein in proteins:
        proteinName = protein[0]
        chainLetter = protein[1]

        middlePDBFile = []

        [ coordinatesPDB, headerPDBFile, tailPDBFile ] = readCIFFile(proteinName, chainLetter, middlePDBFile)

        angles3DAB = readAngles3DAB(proteinName)
        
        lengthsBetweenCAPDB = lengthsFromCoordinates(coordinatesPDB)
        newCoordinates = calculateNewCoordinates(angles3DAB, lengthsBetweenCAPDB)

        buildPDBFile(newCoordinates, headerPDBFile, tailPDBFile, middlePDBFile, proteinName)
        rmsd = getRMSD(proteinName)

        print('Protein Name: ' + proteinName)
        print("RMSD: " + '{:.3f}'.format(rmsd))
        print()
    
    # useless stuff \/

    # chain3DAB = readChain3DAB(proteinName)
    # coordinates3DAB = coordinatesFromAngles(angles3DAB)
    # obj3DAB = calculate3DAB(angles3DAB, chain3DAB)

    # lengthsBetweenCA3DAB = lengthsFromCoordinates(coordinates3DAB)
    # newLengths = lengthsFromCoordinates(newCoordinates)

    # print('Free Energy from 3DAB: ' + str(obj3DAB) + '\n')
    # printList(lengthsBetweenCA3DAB, 'Lengths between Ca from 3DAB coordinates')
    # printList(lengthsBetweenCAPDB, 'Lengths between Ca from 3DAB coordinates')
    # printList(newLengths, 'Lengths between Ca from 3DAB coordinates and PDB lengths')
    # print3DMatrix(newCoordinates, 'New Coordinates from 3DAB coordinates and PDB lengths')

    # plot3D(coordinatesPDB)
    # plot3D(newCoordinates)


if __name__ == "__main__":
    pymol.finish_launching()
    main()
    pymol.cmd.quit()
