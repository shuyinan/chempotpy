'''
F2H2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for F2H2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   F2H2_GEN_S2_1996:         single-state surface of F2H2, ground state,    P
    2.   F2H2_GEN_SQSBDE_1995:     single-state surface of F2H2, ground state,    P
    3.   F2H2_GEN_JBKKL_1990:      single-state surface of F2H2, ground state,    P

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for F2H2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   F2H2_GEN_S2_1996:         single-state surface of F2H2, ground state,
                                   availability: potential energy
                                   functional form: GEN
                                   corresponding surface in POTLIB: f2h2s2-2001.f
                                   ref: W. C. Necoechea, D. G. Truhlar,
                                        "An improved potential energy surface for 
                                        the degenerate rearrangement of (HF)2",
                                        Chem. Phys. Lett. 248, 182-188 (1996).
                                        M. Quack, and M. A. Suhm,
                                        "Potential energy surface and energy levels of 
                                        (HF)2 and its D isotopomers",
                                        Mol. Phys. 69, 791-801 (1990).
    ============================================================================================
    2.   F2H2_GEN_SQSBDE_1995:     single-state surface of F2H2, ground state,
                                   availability: potential energy
                                   functional form: GEN
                                   corresponding surface in POTLIB: f2h2sqsbde2001.f
                                   ref: M. Quack, and M. A. Suhm, 
                                        "Potential energy surface and energy levels of 
                                        (HF)2 and its D isotopomers",
                                        Mol. Phys. 69, 791-801 (1990).
                                        M. Kofranek, H. Lischka, and A. Karpfen,
                                        "Coupled pair functional study on the hydrogen 
                                        fluoride dimer. I. Energy surface and characterization 
                                        of stationary points",
                                        Chem. Phys. Lett. 121, 137-153 (1988).
                                        W. Rijks, P. E. S. Wormer, 
                                        "Correlated van der Waals coefficients. II. Dimers 
                                        consisting of CO, HF, H2O, and NH3",
                                        J. Chem. Phys. 90, 6507-6519 (1989).
    ============================================================================================
    3.   F2H2_GEN_JBKKL_1990:      single-state surface of F2H2, ground state,
                                   availability: potential energy
                                   functional form: GEN
                                   corresponding surface in POTLIB: bjkkl.f
                                   ref: R. Jensen, P. R. Bunker, A. Karpfen, M. Kofranek, 
                                        and H. Lischka, 
                                        "An ab initio calculation of the intramolecular 
                                        stretching spectra for the HF dimer and its 
                                        D‐substituted isotopic species", 
                                        J. Chem. Phys. 93, 6266-6280 (1990).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_fluorine=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for F2H2 system")
            sys.exit()
        if n_fluorine!=2:
            print("number of Fluorine atoms not equal 2 for F2H2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))
 
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    # input Cartesian should in order of H1F1H2F2. 
    if surface=='F2H2_GEN_S2_1996' or surface=='F2H2_GEN_SQSBDE_1995' or surface=='F2H2_GEN_JBKKL_1990':
        geom_ordered=sorted(geom, key=lambda x: ("H", "F").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "F").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        r1=np.linalg.norm(xyz[2]-xyz[1])
        r2=np.linalg.norm(xyz[3]-xyz[1])
        if r1<r2:
          tmp_row=np.copy(xyz[1])
          xyz[1]=xyz[3]
          xyz[3]=tmp_row
          tmp_row=np.copy(rank_ordered[1])
          rank_ordered[1]=rank_ordered[3]
          rank_ordered[3]=tmp_row
          tmp_row=np.copy(xyz[2])
          xyz[2]=xyz[3]
          xyz[3]=tmp_row
          tmp_row=np.copy(rank_ordered[2])
          rank_ordered[2]=rank_ordered[3]
          rank_ordered[3]=tmp_row
        else:
          tmp_row=np.copy(xyz[1])
          xyz[1]=xyz[2]
          xyz[2]=tmp_row
          tmp_row=np.copy(rank_ordered[1])
          rank_ordered[1]=rank_ordered[2]
          rank_ordered[2]=tmp_row

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
