'''
CH5
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH5:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH5_VBMM_CBG_2009:    single-state surface of CH5, ground state,      P
    2.  CH5_LEPS_JJ_2002:     single-state surface of CH5, ground state,      P
    3.  CH5_GEN_J1_1987:      single-state surface of CH5, ground state,      P
    4.  CH5_GEN_J2_1987:      single-state surface of CH5, ground state,      P
    5.  CH5_LEPS_JG_1995:     single-state surface of CH5, ground state,      P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4:
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH5_VBMM_CBG_2009:    single-state surface of CH5, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: ch52008.f 
                              ref: J. C. Corchado, J. L. Bravo and J. Espinosa-Garcia,
                                   "The hydrogen abstraction reaction H + CH4⁠. I. New
                                   analytical potential energy surface based on fitting
                                   to ab initio calculations",
                                   J. Chem. Phys. 130, 184314 (2009).
    ============================================================================================
    2.  CH5_LEPS_JJ_2002:     single-state surface of CH5, ground state,
                              availability: potential energy
                              functional form: LEPS
                              corresponding surface in POTLIB: ch5jj2001.f 
                              ref: J. Espinosa-Garcia,
                                   "New analytical potential energy surface for the CH4 + H
                                   hydrogen abstraction reaction: Thermal rate constants and 
                                   kinetic isotope effects",
                                   J. Chem. Phys. 116, 10664-10673 (2002).
    ============================================================================================
    3.  CH5_GEN_J1_1987:      single-state surface of CH5, ground state,
                              availability: potential energy
                              functional form: GEN
                              corresponding surface in POTLIB: ch5j1.f
                              ref: T. Joseph, R. Steckler, and D. G. Truhlar, 
                                   "A new potential energy surface for the CH3 + H2 ↔ CH4 + H 
                                   reaction: Calibration and calculations of rate constants and 
                                   kinetic isotope effects by variational transition state theory
                                   and semiclassical tunneling calculations",
                                   J. Chem. Phys. 87, 7036-7049 (1987).
    ============================================================================================
    4.  CH5_GEN_J2_1987:      single-state surface of CH5, ground state,
                              availability: potential energy
                              functional form: GEN
                              corresponding surface in POTLIB: ch5j2.f 
                              ref: T. Joseph, R. Steckler, and D. G. Truhlar, 
                                   "A new potential energy surface for the CH3 + H2 ↔ CH4 + H 
                                   reaction: Calibration and calculations of rate constants and 
                                   kinetic isotope effects by variational transition state theory
                                   and semiclassical tunneling calculations",
                                   J. Chem. Phys. 87, 7036-7049 (1987).
    ============================================================================================
    5.  CH5_LEPS_JG_1995:     single-state surface of CH5, ground state,
                              availability: potential energy
                              functional form: LEPS
                              corresponding surface in POTLIB: ch5jg2001.f 
                              ref: M. J. T. Jordan, R. G. Gilbert,
                                   "Classical trajectory studies of the reaction CH4 + H → CH3 + H2",
                                   J. Chem. Phys. 102, 5669-5682 (1995).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for CH5 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH5 system")
            sys.exit()
        if n_hydrogen!=5:
            print("number of Hydrogen atoms not equal 5 for CH5 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH5_VBMM_CBG_2009, CH5_LEPS_JJ_2002, CH5_LEPS_JG_1995: input Cartesian should in order of HCHHHH
    if surface=='CH5_VBMM_CBG_2009' or surface=='CH5_LEPS_JJ_2002' or surface=='CH5_LEPS_JG_1995':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        r=np.zeros((5))
        for i in range(5):
            r[i]=np.linalg.norm(xyz[i+1]-xyz[0])
        max_idx=np.argmax(r)
        #first, move the furthest H to position 1
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[max_idx+1]
        xyz[max_idx+1]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[max_idx+1]
        rank_ordered[max_idx+1]=tmp_row
        #next, switch the C and position 1 H. 
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[0]
        xyz[0]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[0]
        rank_ordered[0]=tmp_row

    #for CH5_GEN_J1_1987 and CH5_GEN_J2_1987: input Cartesian should in order of HHHHCH
    if surface=='CH5_GEN_J1_1987' or surface=='CH5_GEN_J2_1987': 
        geom_ordered=sorted(geom, key=lambda x: ("H", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        r=np.zeros((5))
        for i in range(5):
            r[i]=np.linalg.norm(xyz[i+1]-xyz[0])
        max_idx=np.argmax(r)
        #first, move the furthest H to position 4
        tmp_row=np.copy(xyz[4])
        xyz[4]=xyz[max_idx]
        xyz[max_idx]=tmp_row
        tmp_row=np.copy(rank_ordered[4])
        rank_ordered[4]=rank_ordered[max_idx]
        rank_ordered[max_idx]=tmp_row
        tmp_row=np.copy(xyz[5])
        xyz[5]=xyz[4]
        xyz[4]=tmp_row
        tmp_row=np.copy(rank_ordered[5])
        rank_ordered[5]=rank_ordered[4]
        rank_ordered[4]=tmp_row

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
