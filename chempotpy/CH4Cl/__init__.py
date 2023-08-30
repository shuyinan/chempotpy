'''
CH4Cl
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4Cl:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4Cl_VBMM_2005:      single-state surface of CH4Cl, ground state,      P/G
    2.  CH4Cl_VBMM_2000:      single-state surface of CH4Cl, ground state,      P/G

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4Cl:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4Cl_VBMM_2005:      single-state surface of CH4Cl, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: ch4cl2005.f
                              special emphsize on: Cl + CH4 → HCl + CH3 reaction 
                              ref: C. Range, M. Navarrete, J. C. Corchado, 
                                   and J. Espinosa-Garcia, 
                                   "Potential energy surface, kinetics, and dynamics 
                                   study of the Cl + CH4 → HCl + CH3 reaction",
                                   J. Chem. Phys. 124, 124306 (2006).
    ============================================================================================
    1.  CH4Cl_VBMM_2000:      single-state surface of CH4Cl, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: ch4cl2001.f
                              special emphsize on: Cl + CH4 → HCl + CH3 reaction 
                              ref: J. C. Corchado, D. G. Truhlar, and J. Espinosa-Garcia,
                                   "Potential energy surface, thermal, and state-selected 
                                   rate coefficients, and kinetic isotope effects for Cl + 
                                   CH4 → HCl + CH3", 
                                   J. Chem. Phys. 112, 9375-9389 (2000).



        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_chlorine=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for CH4Cl system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for CH4Cl system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4Cl system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4Cl system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4Cl_VBMM: input Cartesian should in order of HCHHHCl
    if surface=='CH4Cl_VBMM_2005' or surface=='CH4Cl_VBMM_2000':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H", "Cl").index(x[0]))
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
