'''
IH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for IH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   IH2_LEPS_DT1_1971:        single-state surface of IH2, ground state,      P/G
    2.   IH2_LEPS_DT2_1971:        single-state surface of IH2, ground state,      P/G
    3.   IH2_LEPS_RSPTS_1970:      single-state surface of IH2, ground state,      P/G
    4.   IH2_LEPS_PPW1_1974:       single-state surface of IH2, ground state,      P/G
    5.   IH2_LEPS_KNPRY_1966:      single-state surface of IH2, ground state,      P/G

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
    1.   IH2_LEPS_DT1_1971:        single-state surface of IH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ih2rmc1d2001.f 
                                   ref: J. W. Duff, and D. G. Truhlar,
                                        "Effect of curvature of the reaction path on 
                                        dynamic effects in endothermic chemical reactions 
                                        and product energies in exothermic reactions",
                                        J. Chem. Phys. 62, 2477-2491 (1975).
    ============================================================================================
    2.   IH2_LEPS_DT2_1971:        single-state surface of IH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ih2rmc3d2001.f 
                                   ref: J. W. Duff, and D. G. Truhlar,
                                        "Effect of curvature of the reaction path on 
                                        dynamic effects in endothermic chemical reactions 
                                        and product energies in exothermic reactions",
                                        J. Chem. Phys. 62, 2477-2491 (1975).
    ============================================================================================
    3.   IH2_LEPS_RSPTS_1970:      single-state surface of IH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ih2raff2001.f
                                   ref: L. Raff, L. Stivers, R. N. Porter, D. L. Thompson, 
                                        and L. B. Sims, 
                                        "Semiempirical VB Calculation of the (H2I2) Interaction 
                                        Potential",
                                        J. Chem. Phys. 52, 3449-3457 (1970).
    ============================================================================================
    4.   IH2_LEPS_PPW1_1974:       single-state surface of IH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ih2eleps.f
                                   ref: D. S. Perry, J. C. Polanyi, and C. W. Wilson, Jr., 
                                        "Location of energy barriers. VI. The dynamics of endothermic 
                                        reactions, ab + c",
                                        Chem. Phys. 3, 317-331 (1974).
    ============================================================================================
    5.   IH2_LEPS_KNPRY_1966:      single-state surface of IH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: ih2elepsa.f
                                   ref: P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner, 
                                        and C. E. Young ,
                                        "Energy Distribution Among Products of Exothermic Reactions. 
                                        II. Repulsive, Mixed, and Attractive Energy Release",
                                        J. Chem. Phys. 44, 1168-1184 (1966).

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
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for IH2 system")
            sys.exit()
        if n_iodine!=1:
            print("number of Iodine atoms not equal 1 for IH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='IH2_LEPS_DT1_1971' or surface=='IH2_LEPS_DT2_1971' or surface=='IH2_LEPS_RSPTS_1970' or surface=='IH2_LEPS_PPW1_1974' or surface=='IH2_LEPS_KNPRY_1966':
        geom_ordered=sorted(geom, key=lambda x: ("I", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("I", "H").index(x[0]))
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
