'''
Cl2H
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for Cl2H:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   Cl2H_LEPS_KNPRY_1966:     single-state surface of Cl2H, ground state,      P
    2.   Cl2H_LEPS_PK3_1987:       single-state surface of Cl2H, ground state,      P  
    3.   Cl2H_LEPS_PK2_1987:       single-state surface of Cl2H, ground state,      P
    4.   Cl2H_LEPS_BCMR_1983:      single-state surface of Cl2H2, ground state,     P

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for Cl2H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   Cl2H_ELP_KNPRY_1966:      single-state surface of Cl2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: hcl2lep2001.f
                                   ref: P. J. Kuntz, E. M. Nemth, J. C. Polanyi, 
                                        S. D. Rosner, and C. E. Young, 
                                        "Energy Distribution Among Products of Exothermic 
                                        Reactions. II. Repulsive, Mixed, and Attractive 
                                        Energy Release", 
                                        J. Chem. Phys. 44, 1168-1184 (1966).
    ============================================================================================
    2.   Cl2H_LEPS_PK3_1987:       single-state surface of Cl2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clhclpk3-2001.f
                                   ref: A. Persky, and H. Kornweitz,
                                        "Oscillating Reactivity of the Light Atom Transfer 
                                        Reaction Cl + HCI —* CIH + Cl: Dependence on the 
                                        Nature of the Potential Energy Surface", 
                                        J. Phys. Chem. 91, 5496-5503 (1987).
    ============================================================================================
    3.   Cl2H_LEPS_PK2_1987:       single-state surface of Cl2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clhclpk3-2001.f
                                   ref: A. Persky, and H. Kornweitz,
                                        "Oscillating Reactivity of the Light Atom Transfer 
                                        Reaction Cl + HCI —* CIH + Cl: Dependence on the 
                                        Nature of the Potential Energy Surface", 
                                        J. Phys. Chem. 91, 5496-5503 (1987).
    ============================================================================================
    4.   Cl2H_LEPS_BCMR_1983:      single-state surface of Cl2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clhclpk3-2001.f
                                   ref: D. K. Bondi, J. N. L. Connor, J. Manz, and J. Romelt,
                                        "Exact quantum and vibrationally adiabatic quantum, 
                                        semiclassical and quasiclassical study of the collinear 
                                        reactions Cl + MuCl, Cl + HCl, Cl + DCl",
                                        Mol. Phys. 50, 467-488 (1983).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_chlorine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for Cl2H system")
            sys.exit()
        if n_chlorine!=2:
            print("number of Chlorine atoms not equal 2 for Cl2H system")
            sys.exit()

    xyz=np.zeros((natoms,3))
 
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for Cl2H_LEPS_KNPRY_1966: input Cartesian should in order of HClCl
    #for Cl2H_LEPS_PK3_1987, Cl2H_LEPS_PK2_1987, Cl2H_LEPS_BCMR_1983: 
    #    input Cartesian should in order of ClHCl
    if surface=='Cl2H_LEPS_KNPRY_1966':
        geom_ordered=sorted(geom, key=lambda x: ("H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "Cl").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='Cl2H_LEPS_PK3_1987' or surface=='Cl2H_LEPS_PK2_1987' or surface=='Cl2H_LEPS_BCMR_1983':
        geom_ordered=sorted(geom, key=lambda x: ("H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "Cl").index(x[0]))
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
