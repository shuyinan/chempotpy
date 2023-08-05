'''
C2O2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C2O2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  COCO_PIPNN:          single-state surface of C2O2, ground state,      P/G



        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C2O2:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  COCO_PIPNN:          single-state surface of C2O2, ground state,
                             availability: potential energy
                             functional form: PIP-NN
                             corresponding surface in POTLIB: C2O2.tar
                             special emphsize on: CO + CO reaction  
                             ref: J. Chen, J. Li, J. M. Bowman, and H. Guo, 
                                  "Energy transfer between vibrationally excited carbon 
                                  monoxide based on a highly accurate six-dimensional 
                                  potential energy surface",
                                  J. Chem. Phys. 153, 054310 (2020).



        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for C2O2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for C2O2 system")
            sys.exit()
        if n_carbon!=2:
            print("number of Nitrogen atoms not equal 2 for C2O2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for COCO_PIPNN: input Cartesian should in order of OCOC
    if surface=='COCO_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("O", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("O", "C").index(x[0]))
        tmp_row=np.copy(geom_ordered[2])
        geom_ordered[2]=geom_ordered[1]
        geom_ordered[1]=tmp_row
        tmp_row=np.copy(rank_ordered[2])
        rank_ordered[2]=rank_ordered[1]
        rank_ordered[1]=tmp_row
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    return xyz, rank_ordered


def reverse_order(v, g, d, geom, rank_ordered):

    natoms=len(geom)
    nstates=len(v)
    v_ro=np.zeros((nstates))
    g_ro=np.zeros((nstates,natoms,3))
    d_ro=np.zeros((nstates,nstates,natoms,3))

    v_ro = v
    for istate in range(nstates):
        for iatom in range(natoms):
            for idir in range(3):
                g_ro[istate][int(rank_ordered[iatom][1])][idir]=g[istate][iatom][idir]

    for istate in range(nstates):
        for jstate in range(nstates):
            for iatom in range(natoms):
                for idir in range(3):
                    d_ro[istate][jstate][int(rank_ordered[iatom][1])][idir]=d[istate][jstate][iatom][idir]


    return v_ro, g_ro, d_ro
