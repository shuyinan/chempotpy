'''
FH2O
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for FH2O:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.  FH2O_PIPNN:      single-state surface of FH2O, doublet, ground state,                           P
    2.  FH2O_SOC_PIPNN:  spin-orbit corrected single-state surface of FH2O, doublet, ground state,      P



        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for FH2O:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  FH2O_PIPNN:          single-state surface of FH2O, doublet, ground state,
                             availability: potential energy
                             functional form: PIP-NN
                             corresponding surface in POTLIB: soc-fh2o-nn.f
                             special emphsize on: F + H2O → HF + OH reaction 
                             ref: J. Li, R. Dawes, and H. Guo, 
                                  "An ab initio based full-dimensional global 
                                  potential energy surface for FH2O(X2A') and 
                                  dynamics for the F + H2O → HF + HO reaction",
                                  J. Chem. Phys. 137, 094304 (2012).

    ============================================================================================
    2.  FH2O_SOC_PIPNN:      spin-orbit corrected single-state surface of FH2O, doublet, ground state,
                             availability: potential energy
                             functional form: PIP-NN
                             corresponding surface in POTLIB: fh2o-pipnn.f
                             special emphsize on: F + H2O → HF + OH reaction
                             ref: J. Li, B. Jiang, and H. Guo, 
                                  "Spin-orbit corrected full-dimensional potential 
                                  energy surfaces for the two lowest-lying electronic 
                                  states of FH2O and dynamics for the F + H2O → HF 
                                  + OH reaction",
                                  J. Chem. Phys. 138, 074309 (2013).



        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_fluorine=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for FH2O system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for FH2O system")
            sys.exit()
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for FH2O system")
            sys.exit()
        if n_fluorine!=1:
            print("number of Fluorine atoms not equal 1 for FH2O system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for FH2O_PIPNN, FH2O_SOC_PIPNN: input Cartesian should in order of HHFO
    if surface=='FH2O_PIPNN' or surface=='FH2O_SOC_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("H", "F", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "F", "O").index(x[0]))
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
