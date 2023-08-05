'''
CH4OCl
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4OCl:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.  CH4OCl_PIPNN:      single-state surface of CH4OCl, ground state,           P
    2.  CH4OCl_VBMM_2023:  single-state surface of CH4OCl, ground state,           P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4OCl:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4OCl_PIPNN:      single-state surface of CH4OCl, ground state,
                           availability: potential energy
                           functional form: PIP-NN
                           corresponding surface in POTLIB: CH4OCl.zip
                           special emphsize on: Cl + CH3OH → HCl + CH3O/CH2OH 
                           ref: D. Lu, J. Li, and H. Guo, 
                                "Comprehensive Investigations of the Cl + CH3OH 
                                → HCl + CH3O/CH2OH Reaction: Validation of 
                                Experiment and Dynamic Insights",
                                CCS Chem. 2, 882-894 (2020).
    ============================================================================================
    2.  CH4OCl_VBMM_2023:  single-state surface of CH4OCl, ground state,
                           availability: potential energy
                           functional form: VBMM
                           corresponding surface in POTLIB: N/A
                           special emphsize on: CH3OH + Cl → CH3O/CH2OH + HCl
                           ref: C. Rangel, and J. Joaquin Espinosa-Garcia,
                                "Kinetics and dynamics study of the Cl(2P) + 
                                CH3OH reaction based on an analytical potential 
                                energy surface",
                                Phys. Chem. Chem. Phys. 25, 10678-10688 (2023).

        """)

    return


def check(system, surface, geom):

    n_carbon=0
    n_hydrogen=0
    n_oxygen=0
    n_chlorine=0
    natoms=len(geom)

    if natoms!=7:
        print("number of atoms not equal 7 for CH4OCl system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for CH4OCl system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4OCl system")
            sys.exit()
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for CH4OCl system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4OCl system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4OCl_PIPNN: input Cartesian should in order of HHHHCClO
    if surface=='CH4OCl_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("H", "C", "Cl", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C", "Cl", "O").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    #for CH4OCl_VBMM_2023: input Cartesian should in order of COHCl
    if surface=='CH4OCl_VBMM_2023':
        geom_ordered=sorted(geom, key=lambda x: ("C", "O", "H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "O", "H", "Cl").index(x[0]))
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
