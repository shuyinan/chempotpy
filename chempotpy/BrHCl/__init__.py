'''
BrHCl
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for BrHCl:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   BrHCl_LEPS1:              single-state surface of BrHCl, ground state,     P/G
    2.   BrHCl_LEPS2:              single-state surface of BrHCl, ground state,     P/G
    3.   BrHCl_LEPS3:              single-state surface of BrHCl, ground state,     P/G


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for BrHCl:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   BrHCl_LEPS1:              single-state surface of BrHCl, ground state,
                                   availability: potential energy, gradient
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brhcllepsk2001.f
                                   ref: D. J. Douglas, J. C. Polanyi, and J. J. Sloan,
                                        "Effect of reagent vibrational excitation on 
                                        the rate of a substantially endothermic reaction;
                                        HCl(ν′ = 1–4) + Br → Cl + HBr",
                                        J. Chem. Phys. 59, 6679-6680 (1973).
                                        J. A. Kaye, and A. Kuppermann,
                                        "Collinear quantum mechanical probabilities and 
                                        rate constants for the Br + HCl(ν = 2, 3, 4) reaction 
                                        using hyperspherical coordinates",
                                        Chem. Phys. Lett. 92, 574-580 (1982).
    ============================================================================================
    2.   BrHCl_LEPS2:              single-state surface of BrHCl, ground state,
                                   availability: potential energy, gradient
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brhcleleps2-2001.f
                                   ref: D. J. Douglas, J. C. Polanyi, and J. J. Sloan,
                                        "Effect of reagent vibrational excitation on 
                                        the rate of a substantially endothermic reaction;
                                        HCl(ν′ = 1–4) + Br → Cl + HBr",
                                        J. Chem. Phys. 59, 6679-6680 (1973).
                                        V. K. Babamov, V. Lopez, and R. A. Marcus,
                                        "Dynamics of hydrogen atom and proton transfer reactions. 
                                        Nearly degenerate asymmetric case"
                                        J. Chem. Phys. 78, 5621-5628 (1983).
    ============================================================================================
    3.   BrHCl_LEPS3:              single-state surface of BrHCl, ground state,
                                   availability: potential energy, gradient
                                   corresponding surface in POTLIB: brhcleleps2-2001.f
                                   ref: D. J. Douglas, J. C. Polanyi, and J. J. Sloan,
                                        "Effect of reagent vibrational excitation on 
                                        the rate of a substantially endothermic reaction;
                                        HCl(ν′ = 1–4) + Br → Cl + HBr",
                                        J. Chem. Phys. 59, 6679-6680 (1973).
                                        V. K. Babamov, V. Lopez, and R. A. Marcus,
                                        "Dynamics of hydrogen atom and proton transfer reactions. 
                                        Nearly degenerate asymmetric case"
                                        J. Chem. Phys. 78, 5621-5628 (1983).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_bromine=0
    n_chlorine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Br'):
                 n_bromine=n_bromine+1
            elif (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for BrHCl system")
            sys.exit()
        if n_bromine!=1:
            print("number of Bromine atoms not equal 1 for BrHCl system")
            sys.exit()
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for BrHCl system")
            sys.exit()

    xyz=np.zeros((natoms,3))
     
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='BrHCl_LEPS1' or surface=='BrHCl_LEPS2' or surface=='BrHCl_LEPS3':
        geom_ordered=sorted(geom, key=lambda x: ("Br", "H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Br", "H", "Cl").index(x[0]))
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
