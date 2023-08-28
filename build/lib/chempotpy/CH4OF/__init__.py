'''
CH4OF
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4OF:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.  CH4OF_PIPNN:       single-state surface of CH4OF, ground state,             P


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4OF:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4OF_PIPNN:       single-state surface of CH4OF, ground state,
                           availability: potential energy
                           functional form: PIP-NN
                           corresponding surface in POTLIB: CH4OF.zip
                           special emphsize on: F + CH3OH → HF + CH3O/CH2OH 
                           ref: M. L. Weichman, J. A. DeVine, M. C. Babin, J. Li,
                                L. Gup, J. Ma, H. Guo, D. M. Neumark, 
                                "Feshbach resonances in the exit channel of the
                                F + CH3OH → HF + CH3O reaction via transition
                                state spectroscopy",
                                Nat. Chem. 9, 950-955 (2017).
                                D. Lu, C. Xie, J. Li, and H. Guo,
                                "Rate Coefficients and Branching Ratio for Multi-
                                Channel HydrogenAbstractions from CH3OH by F",
                                Chinese J. Chem. Phys. 32, 84-88 (2019). 
                                D. Lu , J. Li, and H. Guo,
                                "Stereodynamical Control of Product Branching in
                                Multi-Channel Barrierless Hydrogen Abstraction of
                                CH3OH by F"
                                Chem. Sci. 10, 7994-8001 (2019).

        """)

    return


def check(system, surface, geom):

    n_carbon=0
    n_hydrogen=0
    n_oxygen=0
    n_fluorine=0
    natoms=len(geom)

    if natoms!=7:
        print("number of atoms not equal 7 for CH4OF system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for CH4OF system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4OF system")
            sys.exit()
        if n_fluorine!=1:
            print("number of Chlorine atoms not equal 1 for CH4OF system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4OF system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4OF_PIPNN: input Cartesian should in order of HHHHCFO
    if surface=='CH4OF_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("H", "C", "F", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C", "F", "O").index(x[0]))
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
