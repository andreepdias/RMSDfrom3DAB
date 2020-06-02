import __main__
__main__.pymol_argv = [ 'pymol', '-Qc'] # Quiet and no GUI

import pymol

pymol.finish_launching()

##
# Read User Input

proteinPDB = 'c0'
proteinNEW = 'c1'

proteinPDBPath = 'lab/' + proteinPDB + '.cif'
proteinNEWPath = 'lab/' + proteinNEW + '.cif'


pymol.cmd.load(proteinPDBPath)
pymol.cmd.load(proteinNEWPath)

x = pymol.cmd.align(proteinPDB, proteinNEW, cycles = 0, transform = 0)

print(x[0])


# Get out!
pymol.cmd.quit()