'''
OH3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for OH3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces: 
    1.   OH3_LEPS_SE_1986:        single-state surface of OH3, ground state,   P/G

    Multi-State Surfaces:
    1.   OH3_PIP_FFW1_2019:       multi-state surface of OH3, 1-3 state,       P/G/D
    2.   OH3_PIP_FFW2_2022:       multi-state surface of OH3, 1-3 state,       P/G/D  
    3.   OH3_PIP_FFW1_2019_DPEM:  multi-state surface of OH3, 1-3 state,       U/UG
    4.   OH3_PIP_FFW2_2022_DPEM:  multi-state surface of OH3, 1-3 state,       U/UG
 
        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for O2H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   OH3_LEPS_SE_1986:     single-state surface of OH3, ground state,
                               availability: potential energy, gradient
                               functional form: LEPS
                               corresponding surface in POTLIB: oh3.f
                               ref: G. C. Schatz, and H. Elgersma,
                                    "A quasi-classical trajectory study of product vibrational 
                                    distributions in the OH + H2 → H2O + H reaction", 
                                    Chem. Phys. Lett. 73, 21-25 (1980).

    Multi-State Surfaces:
    1.   OH3_PIP_FFW1_2019:    multi-state surface of OH3, 1-3 state,
                               availability: potential energy, gradient, nonadiabatic coupling vector
                               functional form: PIP
                               corresponding surface in POTLIB: oh3pes.zip
                               ref: Yinan Shu, Joanna Kryven, Antonio Gustavo Sampaio de Oliveira Filho,
                                    Steven L. Mielke, Linyao Zhang, Guo-Liang Song, Shaohong L. Li,
                                    Rubén Meana-Pañeda, Bina Fu, Joel M. Bowman, and Donald G. Truhlar,
                                    "Direct diabatization and analytic representation of coupled potential 
                                    energy surfaces and couplings for the reactive quenching of the excited 
                                    2Σ+ state of OH by molecular hydrogen",
                                    J. Chem. Phys. 151, 104311 (2019).
    ============================================================================================
    2.   OH3_PIP_FFW2_2022:    multi-state surface of OH3, 1-3 state,
                               availability: potential energy, gradient, nonadiabatic coupling vector
                               functional form: PIP
                               corresponding surface in POTLIB: oh3pes2022.zip
                               ref: Shanyu Han, Antonio Gustavo Sampaio de Oliveira Filho, Yinan Shu, 
                                    Donald G. Truhlar, and Hua Guo,
                                    "Semiclassical Trajectory Studies of Reactive and Nonreactive Scattering
                                    of OH(A2Sigma+) by H2 Based on an Improved Full-Dimensional Ab Initio
                                    Diabatic Potential Energy Matrix",
                                    Chem. Phys. Chem. 23, e202200039 (2022).


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_hydrogen=0
    natoms=len(geom)

    if natoms!=4:
        check=0
        print("number of atoms not equal 4")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for OH3 system")
            sys.exit()
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for OH3 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='OH3_LEPS_SE_1986':
        geom_ordered=sorted(geom, key=lambda x: ("O", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("O", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='OH3_PIP_FFW1_2019' or surface=='OH3_PIP_FFW2_2022' or surface=='OH3_PIP_FFW1_2019_DPEM' or surface=='OH3_PIP_FFW2_2022_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("H", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "O").index(x[0]))
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


def reverse_order_dpem(u, ug, geom, rank_ordered):

    natoms=len(geom)
    nstates=len(u)
    u_ro=np.zeros((nstates,nstates))
    ug_ro=np.zeros((nstates,nstates,natoms,3))

    u_ro = u

    for istate in range(nstates):
        for jstate in range(nstates):
            for iatom in range(natoms):
                for idir in range(3):
                    ug_ro[istate][jstate][int(rank_ordered[iatom][1])][idir]=ug[istate][jstate][iatom][idir]


    return u_ro, ug_ro
