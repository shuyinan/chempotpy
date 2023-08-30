'''
CH4F
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4F:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4F_LEPS_RNCG_2006:  single-state surface of CH4F, ground state,      P/G
    2.  CH4F_GEN_GC_1996:     single-state surface of CH4F, ground state,      P/G
    3.  CH4F_GEN_GCSO_1996:   single-state surface of CH4F, ground state,      P/G

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4F:
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4F_LEPS_RNCG_2006:  single-state surface of CH4F, ground state,
                              availability: potential energy, gradient
                              functional form: LEPS
                              corresponding surface in POTLIB: ch4f2006.f
                              ref: C. Rangel, M. Navarrete, J. C. Corchado, J. Garcia, 
                                   "Potential energy surface, kinetics, and dynamics 
                                   study of the Cl + CH4 → HCl + CH3 reaction",
                                   J. Chem. Phys. 124, 124306 (2006).
    ============================================================================================
    2.  CH4F_GEN_GC_1996:     single-state surface of CH4F, ground state,
                              availability: potential energy, gradient
                              functional form: GEN
                              corresponding surface in POTLIB: ch4f-noso.f 
                              ref: J. Garcia, J. C. Corchado, 
                                   "Recalibration of Two Earlier Potential Energy Surfaces 
                                   for the CH4 + H → CH3 + H2 Reaction. Application of Variational 
                                   Transition-State Theory and Analysis of the Kinetic Isotope Effects 
                                   Using Rectilinear and Curvilinear Coordinates",
                                   J. Phys. Chem. 100, 16561-16567 (1996).
    ============================================================================================
    3.  CH4F_GEN_GCSO_1996:   single-state surface of CH4F, ground state,
                              availability: potential energy, gradient
                              functional form: GEN
                              corresponding surface in POTLIB: ch4f-so.f 
                              ref: J. Garcia, J. C. Corchado, 
                                   "Recalibration of Two Earlier Potential Energy Surfaces 
                                   for the CH4 + H → CH3 + H2 Reaction. Application of Variational 
                                   Transition-State Theory and Analysis of the Kinetic Isotope Effects 
                                   Using Rectilinear and Curvilinear Coordinates",
                                   J. Phys. Chem. 100, 16561-16567 (1996).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_fluorine=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for CH4F system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_fluorine!=1:
            print("number of Fluorine atoms not equal 1 for CH4F system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4F system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4F system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4F_LEPS_RNCG_2006, CH4F_GEN_GC_1996, CH4F_GEN_GCSO_1996: input Cartesian should in order of HCHHHF
    if surface=='CH4F_LEPS_RNCG_2006' or surface=='CH4F_GEN_GC_1996' or surface=='CH4F_GEN_GCSO_1996':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H", "F").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H", "F").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[0]
        xyz[0]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[0]
        rank_ordered[0]=tmp_row

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
