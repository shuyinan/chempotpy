'''
CH4O
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4O:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4O_LEP_1995:       single-state surface of CH4O, ground state,     P/G
    2.  CH4O_LEP_1998:       single-state surface of CH4O, ground state,     P/G
    3.  CH4O_LEP_2000:       single-state surface of CH4O, ground state,     P/G
    4.  CH4O_LEP_2014:       single-state surface of CH4O, ground state,     P/G


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4O:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4O_LEP_1995:        single-state surface of CH4O, ground state,
                              availability: potential energy
                              functional form: LEP
                              corresponding surface in POTLIB: ch4o95-2001.f
                              special emphsize on: O + CH4 → OH + CH3 reaction 
                              ref: J. C. Corchado, J. Espinosa-Garcia, O. Roberto-Neto,
                                   Y.-Y. Chuang, and D. G. Truhlar,
                                   "Dual-Level Direct Dynamics Calculations of the 
                                   Reaction Rates for a Jahn−Teller Reaction:  Hydrogen 
                                   Abstraction from CH4 or CD4 by O(3P)",
                                   J. Phys. Chem. A 102, 4899-4910 (1998).
     ============================================================================================
    2.  CH4O_LEP_1998:        single-state surface of CH4O, ground state,
                              availability: potential energy
                              functional form: LEP
                              corresponding surface in POTLIB: ch4o98-2001.f
                              special emphsize on: O + CH4 → OH + CH3 reaction 
                              ref: J. C. Corchado, J. Espinosa-Garcia, O. Roberto-Neto,
                                   Y.-Y. Chuang, and D. G. Truhlar,
                                   "Dual-Level Direct Dynamics Calculations of the 
                                   Reaction Rates for a Jahn−Teller Reaction:  Hydrogen 
                                   Abstraction from CH4 or CD4 by O(3P)",
                                   J. Phys. Chem. A 102, 4899-4910 (1998).
     ============================================================================================
    3.  CH4O_LEP_2000:        single-state surface of CH4O, ground state,
                              availability: potential energy
                              functional form: LEP
                              corresponding surface in POTLIB: ch4o01-2001.f
                              special emphsize on: O + CH4 → OH + CH3 reaction 
                              ref: J. Espinosa-Garcia, J. C. Garcia-Bernandez, 
                                   "Analytical potential energy surface for the CH4 + O(3P) 
                                   → CH3 + OH reaction. Thermal rate constants and kinetic 
                                   isotope e†ects",
                                   Phys. Chem. Chem. Phys. 2, 2345-2351 (2000).
     ============================================================================================
    4.  CH4O_LEP_2014:       single-state surface of CH4O, ground state,
                              availability: potential energy
                              functional form: LEP
                              corresponding surface in POTLIB: ch4o-2014.f
                              special emphsize on: O + CH4 → OH + CH3 reaction 
                              ref: E. Gonzalez-Lavado, J.C. Corchado, Y. V. Suleimanov, 
                                   W. H. Green, and J. Espinosa-Garcia,
                                   "Theoretical Kinetics Study of the O(3P) + CH4/CD4 
                                   Hydrogen Abstraction Reaction: The Role of Anharmonicity, 
                                   Recrossing Effects, and Quantum Mechanical Tunneling",
                                   J. Phys. Chem. A 118, 3243-3252 (2014).
                                   


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for CH4O system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for CH4O system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4O system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4O system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4O_LEP_1995, CH4O_LEP_1998, CH4O_LEP_2000, CH4O_LEP_2014: input Cartesian should in order of HCHHHO
    if surface=='CH4O_LEP_1995' or surface=='CH4O_LEP_1998' or surface=='CH4O_LEP_2000' or surface=='CH4O_LEP_2014':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H", "O").index(x[0]))
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
