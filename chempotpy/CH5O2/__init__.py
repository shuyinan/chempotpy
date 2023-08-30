'''
CH5O2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH5O2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.  CH5O2_VBMM_2023:   single-state surface of CH5O2, ground state,        P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH5O2:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH5O2_VBMM_2023:   single-state surface of CH5O2, ground state,
                           availability: potential energy
                           functional form: VBMM
                           corresponding surface in POTLIB: N/A
                           special emphsize on: CH3OH + OH → CH3O/CH2OH + H2O
                           ref: J. Espinosa-Garcia, and C. Rangel,
                                "Global potential energy surface and dynamics for 
                                the OH + CH3OH reaction",
                                ChemRxiv, doi: 10.26434/chemrxiv-2022-zj8hx

        """)

    return


def check(system, surface, geom):

    n_carbon=0
    n_hydrogen=0
    n_oxygen=0
    natoms=len(geom)

    if natoms!=8:
        print("number of atoms not equal 8 for CH5O2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for CH5O2 system")
            sys.exit()
        if n_hydrogen!=5:
            print("number of Hydrogen atoms not equal 5 for CH5O2 system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH5O2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH5O2_VBMM_2023: input Cartesian should in order of COHHHHOH
    if surface=='CH5O2_VBMM_2023':
        geom_ordered=sorted(geom, key=lambda x: ("C", "O", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "O", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        rO=np.zeros((2))
        for i in range(2):
            rO[i]=np.linalg.norm(xyz[i+1]-xyz[0])
        max_idx=np.argmax(rO)
        rH=np.zeros((5))
        for i in range(5):
            rH[i]=np.linalg.norm(xyz[i+3]-xyz[max_idx+1])
        min_idx=np.argmin(rH)
        #now move the H attached to O to the last position 
        tmp_row=np.copy(xyz[7])
        xyz[7]=xyz[min_idx+3]
        xyz[min_idx+3]=tmp_row
        tmp_row=np.copy(rank_ordered[7])
        rank_ordered[7]=rank_ordered[min_idx+3]
        rank_ordered[min_idx+3]=tmp_row
        #now move the O that is not attached to C to position 2
        tmp_row=np.copy(xyz[2])
        xyz[2]=xyz[max_idx+1]
        xyz[max_idx+1]=tmp_row
        tmp_row=np.copy(rank_ordered[2])
        rank_ordered[2]=rank_ordered[max_idx+1]
        rank_ordered[max_idx+1]=tmp_row
        #and move this position 2 Oxygen to position 6
        tmp_row=np.copy(xyz[2])
        xyz[2]=xyz[6]
        xyz[6]=tmp_row
        tmp_row=np.copy(rank_ordered[2])
        rank_ordered[2]=rank_ordered[6]
        rank_ordered[6]=tmp_row

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
