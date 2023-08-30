'''
N3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for N3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    1.  N3_4App_ZV:       single-state surface of N3, quartet, A'' symmetry ground state,      P/G


    Single-State Surfaces:
    1.  N3_MBP_2017:      single-state surface of N3, ground state,                            P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for N3:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    1.  N3_4App_ZV:       single-state surface of N3, singlet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: N3_4App_MB-PIP-MEG.f90
                          ref: Z. Varga, and D. G. Truhlar,
                               "Potential energy surface for high-energy N + N2 collisions",
                               Phys. Chem. Chem. Phys. 23, 26273-26284 (2021).



    Single-State Surfaces:
    1.  N3_MBP_2017:      single-state surface of N3, ground electronic state; but spin and
                          spatial symmetry information is unknow. 
                          availability: potential energy
                          functional form: polynomial
                          corresponding surface in POTLIB: n3pes.f90
                          ref: T. K. Mankodi, U. V. Bhandarkar, B. P. Puranik,
                               "Dissociation cross sections for N2+N → 3N and O2+O → 3O using
                               the QCT method",
                               J. Chem. Phys. 146, 204307 (2017).


        """)


    return

def check(system, surface, geom):

    n_nitrogen=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3 for N3 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='N'):
                n_nitrogen=n_nitrogen+1
        if n_nitrogen!=3:
            print("number of Nitrogen atoms not equal 3 for N3 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    for iatom in range(natoms):
        for idir in range(3):
            xyz[iatom][idir]=geom[iatom][idir+1]

    return xyz, rank

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
