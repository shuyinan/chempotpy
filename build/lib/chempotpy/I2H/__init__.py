'''
I2H
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for I2H:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   I2H_LEPS_KK_1981:         single-state surface of I2H, ground state,      P/G
    2.   I2H_LEPS_PT1_1971:        single-state surface of I2H, ground state,      P/G
    3.   I2H_LEPS_PT2_1971:        single-state surface of I2H, ground state,      P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for BrH2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   I2H_LEPS_KK_1981:         single-state surface of I2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: hi2kka2001.f
                                   ref: J. A. Kaye, and A. Kuppermann,
                                        "Collinear quantum mechanical probabilities for 
                                        the I + HI → IH + I reaction using hyperspherical 
                                        coordinates",
                                        Chem. Phys. Lett. 77, 573-579 (1981).
    ============================================================================================
    2.   I2H_LEPS_PT1_1971:        single-state surface of I2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ihielepsa2001.f
                                   ref: C. A. Parr, and D. G. Truhlar,
                                        "Potential Energy Surfaces for Atom Transfer Reactions 
                                        Involving Hydrogens and Halogens"
                                        J. Phys. Chem. 75, 1844-1860 (1971).
    ============================================================================================
    3.   I2H_LEPS_PT2_1971:        single-state surface of I2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ihielepsb2001.f
                                   ref: C. A. Parr, and D. G. Truhlar,
                                        "Potential Energy Surfaces for Atom Transfer Reactions 
                                        Involving Hydrogens and Halogens"
                                        J. Phys. Chem. 75, 1844-1860 (1971).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_iodine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='I'):
                 n_iodine=n_iodine+1
        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for I2H system")
            sys.exit()
        if n_iodine!=2:
            print("number of Iodine atoms not equal 2 for I2H system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for I2H_LEPS_KK_1981, I2H_LEPS_PT1_1971, I2H_LEPS_PT2_1971: IHI    
    if surface=='I2H_LEPS_KK_1981' or surface=='I2H_LEPS_PT1_1971' or surface=='I2H_LEPS_PT2_1971':
        geom_ordered=sorted(geom, key=lambda x: ("H", "I").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "I").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        tmp_row=np.copy(xyz[0])
        xyz[0]=xyz[1]
        xyz[1]=tmp_row
        tmp_row=np.copy(rank_ordered[0])
        rank_ordered[0]=rank_ordered[1]
        rank_ordered[1]=tmp_row


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
