'''
(c) 2021 Anton Hanke,
    Institute for Pharmacy and Molecular Biotechnology (IPMB) Heidelberg

Licensce: MIT License
'''

from pymol import *
import numpy as np
import math as m

radians = lambda angle: angle*m.pi/180

def placePseudoAtom(radspherename, coords):
    cmd.pseudoatom(radspherename,
                   vdw=0.2,
                   mode='extent',
                   pos=coords)

def coordinates(r, phi, psi, sel):
    x = float(r) * m.cos(phi) * m.sin(psi)
    y = float(r) * m.sin(phi) * m.sin(psi)
    z = float(r) * m.cos(psi)
    return(list(np.add([x,y,z], sel)))

def radsphere(radspherename, r, selection, dens=1.0):
    """
    DESCRIPTION
        Provide command to generate a sphere of pseudoatoms
        around each atom of a selection.
        radspherename: name of object to create containing the pseudoatoms
        r: radius of spheres
        selection: selection for which to generate.
        dens: density of pseudoatoms on sphere surface (0 < dens < inf),
              default =1.0
    """
    cmd.create(radspherename, selection=selection)
    step = int(float(dens)*int(np.arccos((float(r)-1)/float(r))/(m.pi/180)))
    for atom in cmd.get_coords(selection):
        for phi in range(0, 360, step):
            for psi in range(0, 360, step):
                coords = coordinates(r=r, phi=radians(phi),
                                     psi=radians(psi),
                                     sel=atom)
                placePseudoAtom(radspherename=radspherename,
                                coords=coords)

cmd.extend('genRad', radsphere)
# radsphere(radspherename=sys.argv[1],
#           r=sys.argv[2], 
#           selection=sys.argv[3])
