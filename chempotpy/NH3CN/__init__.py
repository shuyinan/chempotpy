'''
NH3CN
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for NH3CN:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  NH3CN_VBMM:           single-state surface of NH3CN, ground state,    P/G

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The N atom closer to C will be considered as the last N atom, i.e. the N in CN.


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for NH3CN:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  NH3CN_VBMM:           single-state surface of NH3CN, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: nh3cn.f
                              special emphsize on: CN + NH3 reaction 
                              ref: J. Espinosa-Garcia, C. Rangel, M. Garcia-Chamorro,
                                   and J. C. Corchado,
                                   "Quasi-Classical Trajectory Study of the CN + NH3 
                                   Reaction Based on a Global Potential Energy Surface", 
                                   Molecules 26, 994 (2021).

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The N atom closer to C will be considered as the last N atom.

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_carbon=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for NH3CN system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for NH3CN system")
            sys.exit()
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for NH3CN system")
            sys.exit()
        if n_nitrogen!=2:
            print("number of Nitrogen atoms not equal 2 for NH3CN system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    #for NH3CN_VBMM: input Cartesian should in order of HNHHCN
    if surface=='NH3CN_VBMM':
        geom_ordered=sorted(geom, key=lambda x: ("H", "N", "C").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        #now do re-ordering. 
        r1=np.linalg.norm(xyz[3]-xyz[5])
        r2=np.linalg.norm(xyz[4]-xyz[5])
        if r1<r2:
          tmp_row=xyz[4].copy()
          xyz[4]=xyz[3]
          xyz[3]=tmp_row
          tmp_row=xyz[1].copy()
          xyz[1]=xyz[3]
          xyz[3]=tmp_row
        else:
          tmp_row=xyz[1].copy()
          xyz[1]=xyz[3]
          xyz[3]=tmp_row

    return xyz

