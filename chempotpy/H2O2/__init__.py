'''
H2O2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H2O2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   H2O2_GEN_1998:      single-state surface of H2O2, ground state,      P/G


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for H2O2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H2O2_GEN_1998:     single-state surface of H2O2, ground state,
                            availability: potential energy, gradient
                            functional form: GEN
                            corresponding surface in POTLIB: h2o2.f
                            ref: J. Koput, S. Carter, N. C. Handy, 
                                 "Potential Energy Surface and Vibrational-
                                 Rotational Energy Levels of Hydrogen Peroxide",
                                 J. Phys. Chem. A 102, 6325-6330 (1998).


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_hydrogen=0
    natoms=len(geom)
    
    if natoms!=4:
        print("number of atoms not equal 4")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for H2O2 system")
            sys.exit()
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for H2O2 system") 
            
    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for H2O2_GEN_1998: input Cartesian should in order of HOOH 
    if surface=='H2O2_GEN_1998':
        geom_ordered=sorted(geom, key=lambda x: ("O", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("O", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

        r1=np.linalg.norm(xyz[2]-xyz[1])
        r2=np.linalg.norm(xyz[3]-xyz[1])
        if r1<r2:
          tmp_row=np.copy(xyz[0])
          xyz[0]=xyz[2]
          xyz[2]=tmp_row
          tmp_row=np.copy(rank_ordered[0])
          rank_ordered[0]=rank_ordered[2]
          rank_ordered[2]=tmp_row
        else:
          tmp_row=np.copy(xyz[0])
          xyz[0]=xyz[3]
          xyz[3]=tmp_row
          tmp_row=np.copy(rank_ordered[0])
          rank_ordered[0]=rank_ordered[3]
          rank_ordered[3]=tmp_row
          tmp_row=np.copy(xyz[2])
          xyz[2]=xyz[3]
          xyz[3]=tmp_row
          tmp_row=np.copy(rank_ordered[2])
          rank_ordered[2]=rank_ordered[3]
          rank_ordered[3]=tmp_row


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
