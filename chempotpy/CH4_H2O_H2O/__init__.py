'''
CH4_H2O_H2O
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4_H2O_H2O:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4_H2O_H2O_PIP:      single-state surface of CH4_H2O_H2O, ground state,      P

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The structure is: CH4-H2O-H2O. Because this surfaces only performs partially 
    permutationally invariant, therefore, we compute OH, CH distances to figure out
    the ordering.


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4_H2O_H2O:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4_H2O_H2O_PIP:      single-state surface of CH4_H2O_H2O, ground state,
                              availability: potential energy
                              functional form: PIP
                              corresponding surface in POTLIB: CH4-H2O-H2O.zip
                              ref: R. Conte, C. Qu, and J. M. Bowman, 
                                   "Permutationally Invariant Fitting of Many-Body, 
                                   Non-covalent Interactions with Application to 
                                   Three-Body Methane−Water−Water", 
                                   J. Chem. Theory Comput. 11, 1631-1638 (2015).

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The structure is: CH4-H2O-H2O. Because this surfaces only performs partially 
    permutationally invariant, therefore, we compute OH, CH distances to figure out
    the ordering.

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=11:
        print("number of atoms not equal 11 for CH4_H2O_H2O system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for CH4_H2O_H2O system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4_H2O_H2O system")
            sys.exit()
        if n_hydrogen!=8:
            print("number of Hydrogen atoms not equal 8 for CH4_H2O_H2O system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    ordered_rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4_H2O_H2O_PIP: input Cartesian should in order of HHHHHHHHOOC
    if surface=='CH4_H2O_H2O_PIP': 
        geom_ordered=sorted(geom, key=lambda x: ("H", "O", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "O", "C").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

        array_H=[0, 1, 2, 3, 4, 5, 6, 7]
        rO1H=np.zeros((8))
        for i in range(8):
            rO1H[i]=np.linalg.norm(xyz[i]-xyz[8])
        O1H_idx=rO1H.argsort()[:2]
        array_H.remove(O1H_idx[0])
        array_H.remove(O1H_idx[1])
        rO2H=np.zeros((8))
        for i in range(8):
            rO2H[i]=np.linalg.norm(xyz[i]-xyz[9])
        O2H_idx=rO2H.argsort()[:2]
        array_H.remove(O2H_idx[0])
        array_H.remove(O2H_idx[1])
        ordered_xyz=np.zeros((natoms,3))
        k=0
        for i in array_H:
            ordered_xyz[k]=np.copy(xyz[i])
            ordered_rank[k]=np.copy(rank_ordered[i])
            k=k+1
        ordered_xyz[4]=np.copy(xyz[O1H_idx[0]])
        ordered_xyz[5]=np.copy(xyz[O1H_idx[1]])
        ordered_xyz[6]=np.copy(xyz[O2H_idx[0]])
        ordered_xyz[7]=np.copy(xyz[O2H_idx[1]])
        ordered_xyz[8]=np.copy(xyz[8])
        ordered_xyz[9]=np.copy(xyz[9])
        ordered_xyz[10]=np.copy(xyz[10])
        xyz=np.copy(ordered_xyz)

        ordered_rank[4]=np.copy(rank_ordered[O1H_idx[0]])
        ordered_rank[5]=np.copy(rank_ordered[O1H_idx[1]])
        ordered_rank[6]=np.copy(rank_ordered[O2H_idx[0]])
        ordered_rank[7]=np.copy(rank_ordered[O2H_idx[1]])
        ordered_rank[8]=np.copy(rank_ordered[8])
        ordered_rank[9]=np.copy(rank_ordered[9])
        ordered_rank[10]=np.copy(rank_ordered[10])
        rank_ordered=np.copy(ordered_rank)


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
