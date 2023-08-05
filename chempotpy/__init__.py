'''
**************************************************************************
ChemPotPy: 
A Python-based library for analytical representation of potential energy surfaces
 
ChemPotPy is based on PotLib: https://comp.chem.umn.edu/potlib/index.html 

Authors: Yinan Shu, Zoltan Varga, and Donald G. Truhlar 
         University of Minnesota
**************************************************************************
'''

import importlib
import inspect 
import os
import numpy as np
import array

print("Welcome using ChemPotPy, a library for analytical PESs")
print("To see quick introduction, use chempotpy.intro()")
print("To see the complete manual, use chempotpy.manual()")

def intro():

    print("""

    **************************************************************************
                            ChemPotPy Introduction                                 
    **************************************************************************

    ChemPotPy: A Python-based library for analytic representation of potential 
              energy surfaces for molecules and materials
  
    Authors: Yinan Shu, Zoltan Varga, and Donald G. Truhlar
             University of Minnesota, Minneapolis, MN 55455-0431, USA

    The input to ChemPotPy is ALWAYS Cartesian coordinates in angstroms.
    The output units are ALWAYS potential energies in eV, gradients in 
    eV/angstrom, and nonadiabatic coupling vectors in 1/angstrom. 

    Terminology: 
    “System” refers to the chemical formula of the system, and there may be
    more than one potential energy surface or surface set for a given system.
    Each potential energy surface for a system has its own “name”.

    ====Alphabetical List of Systems====
    To show the list of Systems, use:
      chempotpy.system()

    ====List of Surface Names for Each System===
    To show the SHORT list of surfaces for each system, use:
      chempotpy.list(system) 
    For example:
      system='O3'
      chempotpy.list(system) 

    To show the DETAILED list of surfaces for each system, use:
      chempotpy.detail(system) 

    In the detailed list, all information including the fitting 
    functional form and corresponding references are described. 

    ====Basic Usage==== 
    First, specify your system, surface, and geometry.
    Then call chempotpy working functions as follows:

    For potential energy only:
      v=chempotpy.p(system, name, geom)

    For both potential energy and ANALYTIC gradient:
      v,g=chempotpy.pg(system, name, geom)

    For both potential energy and NUMERICAL gradient:
      v,g=chempotpy.pn(system, name, geom)

    For potential energy, gradient, and nonadiabatic coupling:
      v,g,d=chempotpy.pgd(system, name, geom)

    For diabatic potential energy matrix (DPEM) element:
      u=chempotpy.u(system, name, geom)

    For DPEM element and ANALYTICAL gradient:
      u,ug=chempotpy.ug(system, name, geom)

    For DPEM element and NUMERICAL gradient:
      u,un=chempotpy.un(system, name, geom)

    For the output variables, we used v, g, and d in these examples, but the 
    user can give any names for them.


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ====!!!!!!  MUST READ  !!!!!!====
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    These surfaces require an Intel compiler:
    C3H7NO_PIP, CH4_H2O_H2O_PIP, OH3_PIP_FFW2_2022,
    C7H8S_APRP, C7H8S_III_APRP, C7H8S_S22_APRP

    The input geometry for all surfaces should be specified as Cartesians in
    angstroms. Example:
    geom=[['Y',0.8,2.0,1.0],
          ['H',1.5,1.0,0.5],
          ['R',2.6,2.0,2.0]]

    The output gradients and nonadiabatic coupling correspond to the original
    atom order provided by the user (except if the rearrangement does not
    change the atom order).

    Output energies are all in 
    energy: eV
    gradient: eV/angstrom
    nonadiabatic coupling: 1/angstrom

    ====Example==== 
    Input:
    import chempotpy    
    geom=[['Y',0.8,2.0,1.0],
          ['H',1.5,1.0,0.5],
          ['R',2.6,2.0,2.0]]
    system='YRH'
    name='YRH_LEPS_ModelSurface'
    v=chempotpy.p(system, name, geom)
    v,n=chempotpy.pn(system, name, geom)
    v,g=chempotpy.pg(system, name, geom)
    v,g,d=chempotpy.pgd(system, name, geom)
    name='YRH_LEPS_ModelSurface_DPEM'
    u=chempotpy.u(system, name, geom)
    u,ug=chempotpy.ug(system, name, geom)
    u,un=chempotpy.un(system, name, geom)

    Results:
    v[0]=0.52403335; v[1]=4.35878037
    g[0]:
    array([[-1.83816485,  2.60747969,  1.29655703],
           [ 1.88113351, -2.55666357, -1.22751566],
           [-0.04296867, -0.05081612, -0.06904137]])
    n[0]:
    array([[-1.8381639 ,  2.60747547,  1.29655704],
           [ 1.88113227, -2.55665947, -1.22751572],
           [-0.04296868, -0.0508161 , -0.06904132]])
    d[0][1]:
    array([[ 0.00047888, -0.00068238, -0.00034052],
           [-0.00091055,  0.00028886, -0.0002491 ],
           [ 0.00043167,  0.00039353,  0.00058962]])
    u:
    array([[5.24033396e-01, 4.39128619e-04],
          [4.39128619e-04, 4.35878032e+00]])
    ug[0][0]:
    array([[-1.83816523,  2.60748023,  1.2965573 ],
           [ 1.88113426, -2.55666378, -1.22751544],
           [-0.04296903, -0.05081645, -0.06904187]])
    un[0][0]:
    array([[-1.83816428,  2.60747602,  1.29655731],
           [ 1.88113302, -2.55665969, -1.22751549],
           [-0.04296904, -0.05081643, -0.06904182]])

    """)

    return


def manual():

    print("""
    **************************************************************************
                              ChemPotPy Manual                                 
    **************************************************************************

    ChemPotPy: A Python-based library for analytic representation of potential 
              energy surfaces for molecules and materials
  
    Authors: Yinan Shu, Zoltan Varga, and Donald G. Truhlar
             University of Minnesota, Minneapolis, MN 55455-0431, USA

    ChemPotPy contains subroutines for global and semiglobal representation of 
    molecular potential energies, gradients, and surface couplings as
    functions of nuclear coordinates.

    The input to ChemPotPy is ALWAYS Cartesian coordinates in angstroms.
    The output units are ALWAYS potential energies in eV, gradients in 
    eV/angstrom, and nonadiabatic coupling vectors in 1/angstrom. 
    The zero of energy (the reference geometry for which the potential is set
    to 0) is treated differently for different surfaces.

    History and references:
    The literature reference for each potential is given with that potential.
    Many of the potentials in ChemPotPy are taken from the Potlib library;
    others are taken from the literature. 
    The potentials in PotLib are mainly from previous work in the Truhlar
    group and from contributions to the PotLib library by 
        Joaquin Espinosa-Garcia, Joel Bowman, Hua Guo, Jun Li, 
        Bhalchandra P. Puranik, and David Yarkony. 
    The reference for PotLib is 
        Y. Shu, Z. Varga, A. W. Jasper, B. C. Garrett, J. Espinosa-García, 
        J. C. Corchado, R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, 
        K. R. Yang, T. C. Allison, A. F. Wagner, and D. G. Truhlar, POTLIB,
        http://comp.chem.umn.edu/potlib. 
    
    Terminology: 
    “System” refers to the chemical formula of the system, and there may be
    more than one potential energy surface or surface set for a given system.
    Each potential energy surface for a system has its own “name”.

    ====Alphabetical List of Systems====
    To show the list of Systems, use:
      chempotpy.system()

    ALC AlmHn Aln ArNO 
    BrH2 BrHCl 
    C2H2N2O C2H2O C2H4O4 C2H6Cl C2H6F C2H6H C2H6O 
    C2H6OH C2H7 COCO C3H7NO C6H5SH C6H6O C7H8S 
    CH2O CH2O2 CH2OH CH3NH2 CH3OH CH4 CH4Br CH4Cl 
    CH4CN CH4F CH4_H2O_H2O CH4O CH4OCl CH4OF CH4OH 
    CH5 CH5p Cl2H ClH2 ClH2O ClNH3 CH5O2
    F2H F2H2 FH2 FH2O 
    GeH4OH GeH5 
    H2CO H2O2 H2OBr HOBr H3 H3Cl H3plus H3S H4O2 H7p 
    H3ClO HClOH HO3 HOCO-Ar 
    I2H IH2 K2Rb2 LiFH malonaldehyde MCH MXH N2HOCp 
    N2O N2O2 N3 N4 NaFH NaH2 NH3 NH3CN NH3OH NH4 NinR 
    NO2 O2H O2H3 O3 O4 OH2 OH3 SiH4Cl SiH5 SiOH YRH
    C2N

    ====Systems Involve Multi-State Surfaces====
    C6H5SH: C6H5SH_APRP
    C6H6O: C6H6O_APRP
    C7H8S: C7H8S_APRP, C7H8S_III_APRP, C7H8S_S22_APRP
    CH3NH2: CH3NH2_APRP
    H3: H3_GEN_DMBE_1987
    LiFH: LiFH_LEPS_JHCTP_2001, LiFH_LEPS_JHTP1_2002, LiFH_LEPS_JHTP2_2002
    NaFH: NaFH_LEPS_JHCTP_2001, NaFH_LEPS_JHTP1_2002, NaFH_LEPS_JHTP2_2002
    NaH2: NaH2_LEPS_7S_2000, NaH2_LEPS_7SD_2000, NaH2_LEPS_7L_2000, 
          NaH2_LEPS_7LD_2000, NaH2_LEPS_6_2000, NaH2_LEPS_5G_1992,
          NaH2_LEPS_5F_1992
    NH3: NH3_GEN_FFW1_2006, NH3_GEN_FFW2_2007
    O3: O3_14_3Ap_2022, O3_14_3Ap_2023, O3_6_5Ap_2023
    OH3: OH3_PIP_FFW1_2019, OH3_PIP_FFW2_2022
    MCH: MCHWL_LEPS_ModelSurface, MCHWB_LEPS_ModelSurface,    
         MCHSL_LEPS_ModelSurface, MCHSB_LEPS_ModelSurface, 
         MCHTL_LEPS_ModelSurface
    MXH: MXHWL_LEPS_ModelSurface, MXHSL_LEPS_ModelSurface,
         MXHSB_LEPS_ModelSurface
    YRH: YRH_LEPS_ModelSurface

    ====List of Surface Names for Each System===
    To show the SHORT list of surfaces for each system, use:
      chempotpy.list(system) 
    For example:
      system='O3'
      chempotpy.list(system) 

    To show the DETAILED list of surfaces for each system, use:
      chempotpy.detail(system) 

    In the detailed list, all information including the fitting 
    functional form and corresponding references are described. 

    ====Basic Usage==== 
    First, specify your system, surface, and geometry.
    Then call chempotpy working functions as follows:

    For potential energy only:
      v=chempotpy.p(system, name, geom)

    For both potential energy and ANALYTIC gradient:
      v,g=chempotpy.pg(system, name, geom)

    For both potential energy and NUMERICAL gradient:
      v,g=chempotpy.pn(system, name, geom)

    For potential energy, gradient, and nonadiabatic coupling:
      v,g,d=chempotpy.pgd(system, name, geom)

    For diabatic potential energy matrix (DPEM) element:
      u=chempotpy.u(system, name, geom)

    For DPEM element and ANALYTICAL gradient:
      u,ug=chempotpy.ug(system, name, geom)

    For DPEM element and NUMERICAL gradient:
      u,un=chempotpy.un(system, name, geom)

    Note that for single-state surfaces, one will not be able to compute 
    nonadiabatic coupling vectors. Therefore, if one calls chempotpy.pdg
    working function anyway, one will receive error message. 
    For multi-state surfaces, the resulting v is an array of potential
    Energies. Suppose one has Nstates of electronic states, and Natoms of 
    atoms for the system, the v, g, d are in the following format:
    v(Nstates)
    g(Nstates,Natoms,3)
    d(Nstates,Nstates,Natoms,3)


    For the output variables, we used v, g, and d in these examples, but the 
    user can give any names for them.


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ====!!!!!!  MUST READ  !!!!!!====
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    These surfaces require an Intel compiler:
    C3H7NO_PIP, CH4_H2O_H2O_PIP, OH3_PIP_FFW2_2022,
    C7H8S_APRP, C7H8S_III_APRP, C7H8S_S22_APRP


    The input geometry for all surfaces should be specified as Cartesians in
    angstroms. Example:
    geom=[['Y',0.8,2.0,1.0],
          ['H',1.5,1.0,0.5],
          ['R',2.6,2.0,2.0]]

    When the geometry is read, the program checks whether the number of atoms
    listed in geom corresponds correctly to the number of atoms in the
    requested system. The number of atom types is also checked.

    For most of the surfaces, the user can use any arbitrary atom order;
    the program will reorder them to the required one for the given surface.

    The output gradients and nonadiabatic coupling correspond to the original
    atom order provided by the user (except if the rearrangement does not
    change the atom order).

    Output energies are all in 
    energy: eV
    gradient: eV/angstrom
    nonadiabatic coupling: 1/angstrom

    If the original subroutine output units are different from those above, 
    we use the following conversion factors:

    1 bohr = 0.529177211 angstrom 
    1 hartree = 27.211386 eV
    1 eV = 23.0609 kcal/mol
    1 kcal/mol = 349.757 cm^-1
    therefore, 
    1 hartree = 627.5190514 kcal/mol = 219479.1808631 cm^-1
    1 eV = 8065.7112013 cm^-1
    1 hartree/bohr = 27.211386 eV/0.529177211 angstrom = 51.422067 eV/angstrom


    ====A Complete Example==== 
    Input:
    import chempotpy   
    geom=[['Y',0.8,2.0,1.0],
          ['H',1.5,1.0,0.5],
          ['R',2.6,2.0,2.0]]
    system='YRH'
    name='YRH_LEPS_ModelSurface'
    v=chempotpy.p(system, name, geom)
    v,n=chempotpy.pn(system, name, geom)
    v,g=chempotpy.pg(system, name, geom)
    v,g,d=chempotpy.pgd(system, name, geom)
    name='YRH_LEPS_ModelSurface_DPEM'
    u=chempotpy.u(system, name, geom)
    u,ug=chempotpy.ug(system, name, geom)
    u,un=chempotpy.un(system, name, geom)

    Results:
    v[0]=0.52403335; v[1]=4.35878037
    g[0]:
    array([[-1.83816485,  2.60747969,  1.29655703],
           [ 1.88113351, -2.55666357, -1.22751566],
           [-0.04296867, -0.05081612, -0.06904137]])
    g[1]:
    array([[ 1.00316838, -1.39359782, -0.68143785],
           [-1.95945323,  0.49911168, -0.6449303 ],
       [ 0.95628485,  0.89448614,  1.32636815]])
    n[0]:
    array([[-1.8381639 ,  2.60747547,  1.29655704],
           [ 1.88113227, -2.55665947, -1.22751572],
       [-0.04296868, -0.0508161 , -0.06904132]])
    n[1]:
    array([[ 1.00316792, -1.39359819, -0.68143725],
           [-1.95945252,  0.49911227, -0.64493084],
       [ 0.9562846 ,  0.89448591,  1.3263681 ]])
    d[0][1]:
    array([[ 0.00047888, -0.00068238, -0.00034052],
           [-0.00091055,  0.00028886, -0.0002491 ],
       [ 0.00043167,  0.00039353,  0.00058962]])
    d[1][0]:
    array([[-0.00047888,  0.00068238,  0.00034052],
           [ 0.00091055, -0.00028886,  0.0002491 ],
       [-0.00043167, -0.00039353, -0.00058962]])

    u:
    array([[5.24033396e-01, 4.39128619e-04],
          [4.39128619e-04, 4.35878032e+00]])

    ug[0][0]:
    array([[-1.83816523,  2.60748023,  1.2965573 ],
           [ 1.88113426, -2.55666378, -1.22751544],
           [-0.04296903, -0.05081645, -0.06904187]])

    ug[0][1]:
    array([[-0.00151102,  0.0021586 ,  0.0010793 ],
           [ 0.00305192, -0.00075778,  0.00102193],
           [-0.0015409 , -0.00140082, -0.00210123]])

    ug[1][0]:
    array([[-0.00151102,  0.0021586 ,  0.0010793 ],
           [ 0.00305192, -0.00075778,  0.00102193],
           [-0.0015409 , -0.00140082, -0.00210123]])

    ug[1][1]:
    array([[ 1.00316876, -1.39359837, -0.68143812],
           [-1.95945398,  0.49911189, -0.64493053],
           [ 0.95628521,  0.89448647,  1.32636865]])


    un[0][0]:
    array([[-1.83816428,  2.60747602,  1.29655731],
           [ 1.88113302, -2.55665969, -1.22751549],
           [-0.04296904, -0.05081643, -0.06904182]])

    un[0][1]:
    array([[-0.00151102,  0.0021586 ,  0.0010793 ],
           [ 0.00305194, -0.00075777,  0.00102193],
           [-0.0015409 , -0.00140082, -0.00210124]])

    un[1][0]:
    array([[-0.00151102,  0.0021586 ,  0.0010793 ],
           [ 0.00305194, -0.00075777,  0.00102193],
           [-0.0015409 , -0.00140082, -0.00210124]])

    un[1][1]:
    array([[ 1.00316831, -1.39359873, -0.68143753],
            [-1.95945327,  0.49911249, -0.64493107],
            [ 0.95628497,  0.89448625,  1.3263686 ]])


    Notes:
    v[0]: first electronic state potential energy
    v[1]: second electronic state potential energy
    g[0]: first electronic state analytical gradient
    g[1]: second electronic state analytical gradient
    n[0]: first electronic state numerical gradient
    n[1]: second electronic state numerical gradient
    d[0][1]: nonadiabatic coupling vector between states 1 and 2
    d[1][0]: nonadiabatic coupling vector between states 2 and 1
    u: diabatic potential energy matrix U
    ug[0][0]: analytical gradient of U_11 (first row, first column element of U)
    ug[0][1]: analytical gradient of U_12
    ug[1][0]: analytical gradient of U_21
    ug[1][1]: analytical gradient of U_22
    un[0][0]: numerical gradient of U_11
    un[0][1]: numerical gradient of U_12
    un[1][0]: numerical gradient of U_21
    un[1][1]: numerical gradient of U_22


    ====Full List of Surface Names for Each System===
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    First number denotes the global index of the library
    Surface index start with S: single-state surface
    Surface index start with M: multi-state surface
    Surface index start with Z: a set of surfaces use consistent diatomic
    potentials of O2, N2 and NO by Z. Varga. 

    ALC:
    1.   S1.   ALC_LEPS_ModelSurface:     single-state, P/G
    2.   S2.   ALC5m0_LEPS_ModelSurface:  single-state, P/G
    3.   S3.   ALC6m0_LEPS_ModelSurface:  single-state, P/G
    4.   S4.   ALC7m0_LEPS_ModelSurface:  single-state, P/G
    5.   S5.   ALC7m1_LEPS_ModelSurface:  single-state, P/G
    6.   S6.   ALC7m4_LEPS_ModelSurface:  single-state, P/G
    7.   S7.   ALC7m5_LEPS_ModelSurface:  single-state, P/G
    8.   S8.   ALC7m7_LEPS_ModelSurface:  single-state, P/G
    9.   S9.   ALC7m10_LEPS_ModelSurface: single-state, P/G
    10.  S10.  ALC7p2_LEPS_ModelSurface:  single-state, P/G
    11.  S11.  ALC7p3_LEPS_ModelSurface:  single-state, P/G
    12.  S12.  ALC7p5_LEPS_ModelSurface:  single-state, P/G

    AlmHn:
    13.  S1.   AlmHn_VBO1:   single-state P
    14.  S2.   AlmHn_VBO2:   single-state P

    Aln:
    15.  S1.   Aln_BetH:           single state, P
    16.  S2.   Aln_CleR:           single state, P
    17.  S3.   Aln_CoxJM:          single state, P
    18.  S4.   Aln_deSPH_LJAT:     single state, P
    19.  S5.   Aln_deSPH_M:        single state, P
    20.  S6.   Aln_ER2AT:          single state, P
    21.  S7.   Aln_ER2BA:          single state, P
    22.  S8.   Aln_ER2CN:          single state, P
    23.  S9.   Aln_ER2CoxJM:       single state, P
    24.  S10.  Aln_ER2EACN:        single state, P
    25.  S11.  Aln_ER2EAT:         single state, P
    26.  S12.  Aln_ER2EBA:         single state, P
    27.  S13.  Aln_ER2ECN:         single state, P
    28.  S14.  Aln_ER2ESCNa:       single state, P
    29.  S15.  Aln_ER2ESCNm:       single state, P
    30.  S16.  Aln_ER2ES:          single state, P
    31.  S17.  Aln_ER2:            single state, P
    32.  S18.  Aln_ER2GEA:         single state, P
    33.  S19.  Aln_ER2MeiD:        single state, P
    34.  S20.  Aln_ER2MisFMP:      single state, P
    35.  S21.  Aln_ER2SCNa:        single state, P
    36.  S22.  Aln_ER2SCNm:        single state, P
    37.  S23.  Aln_ER2S:           single state, P
    38.  S24.  Aln_ER:             single state, P
    39.  S25.  Aln_ErkIV:          single state, P
    40.  S26.  Aln_ErkV:           single state, P
    41.  S27.  Aln_ErkVIII:        single state, P
    42.  S28.  Aln_Gol:            single state, P
    43.  S29.  Aln_HalP:           single state, P
    44.  S30.  Aln_Jac:            single state, P
    45.  S31.  Aln_MeiD:           single state, P
    46.  S32.  Aln_MisFMP:         single state, P
    47.  S33.  Aln_NP_Ad:          single state, P/G
    48.  S34.  Aln_NP_A:           single state, P
    49.  S35.  Aln_Np_Bd:          single state, P/G
    50.  S36.  Aln_NP_B:           single state, P
    51.  S37.  Aln_NP_C:           single state, P
    52.  S38.  Aln_NP_D:           single state, P
    53.  S39.  Aln_NP_E:           single state, P
    54.  S40.  Aln_PapCEP:         single state, P
    55.  S41.  Aln_PapKEP:         single state, P
    56.  S42.  Aln_PeaTHT:         single state, P
    57.  S43.  Aln_PetW:           single state, P
    58.  S44.  Aln_reCoxJM:        single state, P
    59.  S45.  Aln_redeSPH_M:      single state, P
    60.  S46.  Aln_reErkIV:        single state, P
    61.  S47.  Aln_reErkVIII:      single state, P
    62.  S48.  Aln_reGEA:          single state, P
    63.  S49.  Aln_reGol:          single state, P
    64.  S50.  Aln_reHalP:         single state, P
    65.  S51.  Aln_reJac:          single state, P
    66.  S52.  Aln_reLJAT:         single state, P
    67.  S53.  Aln_reMeiD:         single state, P
    68.  S54.  Aln_reMisFMP:       single state, P
    69.  S55.  Aln_rePetW:         single state, P
    70.  S56.  Aln_reStrM:         single state, P
    71.  S57.  Aln_reSutC:         single state, P
    72.  S58.  Aln_StrM:           single state, P
    73.  S59.  Aln_SutC:           single state, P

    ArNO:
    74.  S1.   ArNO_Ap_PIPNN:      single state, A' symmetry, P
    75.  S2.   ArNO_App_PIPNN:     single state, A'' symmetry, P

    BrH2:
    76.  S1.   BrH2_LEPS_LTBZ1_1995:   single-state, P/G
    77.  S2.   BrH2_LEPS_W_1959:       single-state, P/G
    78.  S3.   BrH2_LEPS_C_1982:       single-state, P/G
    79.  M1.   BrH2_LEPS_LTBZSO_1995:  2-state, P
    80.  M2.   BrH2_LEPS_LTBZ2_1995:   2-state, P/G/D

    BrHCl:
    81.  S1.   BrHCl_LEPS1:    single-state, P/G
    82.  S2.   BrHCl_LEPS2:    single-state, P/G
    83.  S3.   BrHCl_LEPS3:    single-state, P/G

    C2H2N2O:
    84.  S1.  C2H2N2O_PIPNN:   single-state, P

    C2H2O:
    85.  S1.  C2H2O_PIPNN:     single-state, P

    C2H6Cl:
    86.  S1.  C2H6Cl_VBMM:     single-state, P

    C2H6F:
    87.  S1.  C2H6F_VBMM:      single-state, P

    C2H6H:
    88.  S1.  C2H6H_VBMM:      single-state, P

    C2H6O:
    89.  S1.  C2H6O_VBMM:      single-state, P

    C2H7:
    90.  S1.  C2H7_VBMM:       single-state, P

    C3H7NO:
    91.  S1.  C3H7NO_PIP:      single-state, P/G

    C6H5SH:
    92.  M1.  C6H5SH_APRP:         3-state, P/G/D
    93.  M2.  C6H5SH_APRP_DPEM:    3-state, U/UG

    C6H6O:
    94.  M1.  PHOH_APRP:           3-state, P/G/D
    95.  M2.  PHOH_APRP_DPEM:      3-state, U/UG

    C7H8S:
    96.  M1.  C7H8S_APRP:          3-state, P/G/D
    97.  M2.  C7H8S_III_APRP:      3-state, P/G/D
    98.  M3.  C7H8S_S22_APRP:      3-state, P/G/D

    99.  M4.  C7H8S_APRP_DPEM:     3-state, U/UG
    100. M5.  C7H8S_III_APRP_DPEM: 3-state, U/UG
    101. M6.  C7H8S_S22_APRP_DPEM: 3-state, U/UG

    CH2O:
    102. S1.  CH2O_PIP:            single-state, P

    CH2O2:
    103. S1.  CH2O2_PIPNN:         single-state, P

    CH3NH2:
    104. M1.  CH3NH2_APRP:         2-state, P/G/D
    105. M2.  CH3NH2_APRP_DPEM:    2-state, U/UG

    CH3OH:
    106. S1.  CH3OH_PIPNN:         single-state, P

    CH4:
    107. S1.  CH4_GEN_SP_2001:     single-state, P

    CH4Br:
    108. S1.  CH4Br_VBMM_2002:     single-state, P/G

    CH4Cl:
    109. S1.  CH4Cl_VBMM_2005:     single-state, P/G
    110. S2.  CH4Cl_VBMM_2000:     single-state, P/G

    CH4CN:
    111. S1.  CH4CN_VBMM:          single-state, P/G

    CH4F:
    112. S1.  CH4F_LEPS_RNCG_2006: single-state, P/G
    113. S2.  CH4F_GEN_GC_1996:    single-state, P/G
    114. S3.  CH4F_GEN_GCSO_1996:  single-state, P/G

    CH4_H2O_H2O:
    115. S1.  CH4_H2O_H2O_PIP:     single-state, P

    CH4O:
    116. S1.  CH4O_LEP_1995:        single-state, P/G
    117. S2.  CH4O_LEP_1998:        single-state, P/G
    118. S3.  CH4O_LEP_2000:        single-state, P/G
    119. S4.  CH4O_LEP_2014:        single-state, P/G

    CH4OCl:
    120. S1.  CH4OCl_PIPNN:         single-state, P
    121. S2.  CH4OCl_VBMM_2023:     single-state, P

    CH4OF:
    122. S1.  CH4OF_PIPNN:          single-state, P

    CH4OH:
    123. S1.  CH4OH_VBMM_2015:      single-state, P/G
    124. S2.  CH4OH_VBMM_2000:      single-state, P/G

    CH5:
    125. S1.  CH5_VBMM_CBG_2009:    single-state, P
    126. S2.  CH5_LEPS_JJ_2002:     single-state, P
    127. S3.  CH5_GEN_J1_1987:      single-state, P
    128. S4.  CH5_GEN_J2_1987:      single-state, P
    129. S5.  CH5_LEPS_JG_1995:     single-state, P

    CH5O2:
    130. S1.  CH5O2_VBMM_2023:      single-state, P

    CH5p:
    131. S1.  CH5p_PIP:             single-state, P

    Cl2H:
    132. S1.  Cl2H_LEPS_KNPRY_1966: single-state, P
    133. S2.  Cl2H_LEPS_PK3_1987:   single-state, P
    134. S3.  Cl2H_LEPS_PK2_1987:   single-state, P
    135. S4.  Cl2H_LEPS_BCMR_1983:  single-state, P

    ClH2:
    136. S1.  ClH2_LEPS_GTM1_1981:      single-state, P/G
    137. S2.  ClH2_LEPS_GTM2_1981:      single-state, P/G
    138. S3.  ClH2_LEPS_GTM3_1981:      single-state, P/G
    139. S4.  ClH2_LEPS_GTM4_1982:      single-state, P/G
    140. S5.  ClH2_LEPS_PK1_1966:       single-state, P/G
    141. S6.  ClH2_LEPS_PK2_1966:       single-state, P/G
    142. S7.  ClH2_LEPS_STSBLTG1_1989:  single-state, P/G
    143. S8.  ClH2_LEPS_STSBLTG2_1989:  single-state, P/G
    144. S9.  ClH2_LEPS_KNPRY_1966:     single-state, P/G
    145. S10. ClH2_LEPS_ALTG1_1996:     single-state, P/G
    146. S11. ClH2_LEPS_ALTG2_1996:     single-state, P
    147. S12. ClH2_LEPS_SPK_1973:       single-state, P/G
    148. S13. ClH2_LEPS_LB_1980:        single-state, P/G
    149. S14. ClH2_LEPS_TTGI1_1985:     single-state, P/G
    150. S15. ClH2_LEPS_TTGI2_1985:     single-state, P/G

    ClH2O:
    151. S1.  ClH2O_PIPNN:      single-state, P

    ClNH3:
    152. S1.  ClNH3_VBMM:       single-state, P/G

    COCO:
    153. S1.  COCO_PIPNN:       single-state, P/G

    F2H:
    154. S1.  F2H_LEPS_JOT_1972:       single-state, P/G
    155. S2.  F2H_LEPS_KNPRY_1966:     single-state, P/G

    F2H2:
    156. S1.  F2H2_GEN_S2_1996:       single-state, P
    157. S2.  F2H2_GEN_SQSBDE_1995:   single-state, P
    158. S3.  F2H2_GEN_JBKKL_1990:    single-state, P

    FH2:
    159. S1.  FH2_GEN_5SEC_1991:      single-state, P/G
    160. S2.  FH2_GEN_5SECCOG1_1993:  single-state, P/G
    161. S3.  FH2_GEN_5SECCOG2_1993:  single-state, P/G
    162. S4.  FH2_GEN_5SECCOG3_1993:  single-state, P/G
    163. S5.  FH2_GEN_5SECG1_1993:    single-state, P/G
    164. S6.  FH2_GEN_5SECG2_1993:    single-state, P/G
    165. S7.  FH2_GEN_5SECG3_1993:    single-state, P/G
    166. S8.  FH2_GEN_5SECS1_1993:    single-state, P/G
    167. S9.  FH2_GEN_5SECS2_1993:    single-state, P/G
    168. S10. FH2_GEN_5SECS3_1993:    single-state, P/G
    169. S11. FH2_GEN_5SECW_1991:     single-state, P/G
    170. S12. FH2_GEN_6SEC_1993:      single-state, P/G
    171. S13. FH2_LEPS_M5_1991:       single-state, P/G
    172. S14. FH2_LEPS_M5SK_1966:     single-state, P/G
    173. S15. FH2_LEPS_NO1_1984:      single-state, P/G
    174. S16. FH2_LEPS_NO2_1984:      single-state, P/G
    175. S17. FH2_LEPS_NO2a_1984:     single-state, P/G
    176. S18. FH2_LEPS_NO3_1984:      single-state, P/G
    177. S19. FH2_LEPS_NO4_1984:      single-state, P/G
    178. S20. FH2_LEPS_NO5_1985:      single-state, P/G
    179. S21. FH2_LEPS_NO5a_1985:     single-state, P/G
    180. S22. FH2_LEPS_NO5b_1985:     single-state, P/G
    181. S23. FH2_LEPS_NO5z_1985:     single-state, P/G
    182. S24. FH2_LEPS_SK_1980:       single-state, P/G

    FH2O:
    183. S1.  FH2O_PIPNN:      single-state, P
    184. S2.  FH2O_SOC_PIPNN:  single-state, P

    GeH4OH:
    185. S1.  GeH4OH_VBMM:       single-state, P/G

    GeH5:
    186. S1.   GeH5_LEPS:        single-state, P/G

    H2CO:
    187. S1.  H2CO_PIP:          single-state, P

    H2O2:
    188. S1.  H2O2_GEN_1998:     single-state, P/G

    H2OBr:
    189. S1.  H2OBr_PIP:         single-state, P

    H3:
    190. S1.  H3_GEN_A2_2002:       single-state, P
    191. S2.  H3_GEN_A3_2002:       single-state, P
    192. S3.  H3_GEN_A4_2002:       single-state, P
    193. S4.  H3_GEN_CCI_2002:      single-state, P
    194. S5.  H3_GEN_BH_2009:       single-state, P
    195. S6.  H3_GEN_BKMP3_1996:    single-state, P
    196. S7.  H3_GEN_BKMP_1991:     single-state, P/G
    197. S8.  H3_GEN_PK2_1964:      single-state, P/G
    198. S9.  H3_GEN_TK_1984:       single-state, P/G
    199. M1.  H3_GEN_DMBE_1987:     2-state, P/G/D

    H3Cl:
    200. S1.  H3Cl_PIPNN:           single-state, P

    H3ClO:
    201. S1.  H3ClO_PIPNN:          single-state, P

    H3p:
    202. S1.  H3p_UNO_2001       single-state, P/G
    203. S2.  H3p_KBNN_2002      single-state, P/G

    H3S:
    204. S1.  H3S_PIPNN:         single-state, P/G

    H4O2:
    205. S1.  H4O2_PIP:          single-state, P

    HClOH:
    206. S1.  HClOH_PIPNN:       single-state, P

    HO3:
    207. S1.  HO3_PIPNN:         single-state, P

    HOBr:
    208. S1.  HOBr_3App_PIP:     single-state, P/G

    HOCO_Ar:
    209. S1.  HOCO_Ar_PIP:       single states, P

    I2H:
    210. S1.  I2H_LEPS_KK_1981:      single-state, P/G
    211. S2.  I2H_LEPS_PT1_1971:     single-state, P/G
    212. S3.  I2H_LEPS_PT2_1971:     single-state, P/G

    IH2:
    213. S1.  IH2_LEPS_DT1_1971:     single-state, P/G
    214. S2.  IH2_LEPS_DT2_1971:     single-state, P/G
    215. S3.  IH2_LEPS_RSPTS_1970:   single-state, P/G
    216. S4.  IH2_LEPS_PPW1_1974:    single-state, P/G
    217. S5.  IH2_LEPS_KNPRY_1966:   single-state, P/G

    K2Rb2:
    218. S1.  K2Rb2_PIPNN:           single-state, P

    LiFH:
    219. M1.  LiFH_LEPS_JHCTP_2001:       2-state, P/G/D
    220. M2.  LiFH_LEPS_JHTP1_2002:       2-state, P
    221. M3.  LiFH_LEPS_JHTP2_2002:       2-state, P/G/D
    222. M4.  LiFH_LEPS_JHCTP_2001_DPEM:  2-state, U/UG
    223. M5.  LiFH_LEPS_JHTP1_2002_DPEM:  2-state, U
    224. M6.  LiFH_LEPS_JHTP2_2002_DPEM:  2-state, U/UG

    MCH:
    225. M1.  MCHWL_LEPS_ModelSurface:      2-state, P/G/D
    226. M2.  MCHWB_LEPS_ModelSurface:      2-state, P/G/D
    227. M3.  MCHSL_LEPS_ModelSurface:      2-state, P/G/D
    228. M4.  MCHSB_LEPS_ModelSurface:      2-state, P/G/D
    229. M5.  MCHTL_LEPS_ModelSurface:      2-state, P/G/D
    230. M6.  MCHWL_LEPS_ModelSurface_DPEM: 2-state, U/UG
    231. M7.  MCHWB_LEPS_ModelSurface_DPEM: 2-state, U/UG
    232. M8.  MCHSL_LEPS_ModelSurface_DPEM: 2-state, U/UG
    233. M9.  MCHSB_LEPS_ModelSurface_DPEM: 2-state, U/UG
    234. M10. MCHTL_LEPS_ModelSurface_DPEM: 2-state, U/UG

    MXH:
    235. M1.  MXHWL_LEPS_ModelSurface:       2-state, P/G/D
    236. M2.  MXHSL_LEPS_ModelSurface:       2-state, P/G/D
    237. M3.  MXHSB_LEPS_ModelSurface:       2-state, P/G/D
    238. M4.  MXHWL_LEPS_ModelSurface_DPEM   2-state, U/UG
    239. M5.  MXHSL_LEPS_ModelSurface_DPEM   2-state, U/UG
    240. M6.  MXHSB_LEPS_ModelSurface_DPEM   2-state, U/UG

    N2HOCp:
    241. S1.  N2HOCp_PIPNN:      single-state, P

    N2O:
    242. Z1.  N2O_3Ap_ZV:        single-state, triplet, A' symmetry, P/G
    243. Z2.  N2O_3App_ZV:       single-state, triplet, A'' symmetry, P/G
    244. S1.  N2O_3Ap_PIP:       single-state, triplet, A' symmetry, P/G
    245. S2.  N2O_3App_PIP:      single-state, triplet, A'' symmetry, P/G

    N2O2:
    246. Z1.  N2O2_3A_ZV:        single-state, triplet, P/G
    247. S1.  N2O2_3A_PIP:       single-state, triplet, P/G

    N3:
    248. Z1.  N3_4App_ZV:        single-state, quartet, A'' symmetry, P/G
    249. S1.  N3_MBP_2017:       single-state, P

    N4:
    250. Z1.  N4_1A_ZV:          single-state, singlet, P/G
    251. S1.  N4_BDTC_2014:      single-state, singlet, P/G
    252. S2.  N4_1A_PIPNN:       single-state, singlet, P/G
    253. S3.  N4_1A_PIP:         single-state, singlet, P/G
    254. S4.  N4_1A_PIP_2013:    single-state, singlet, P/G

    NaFH:
    255. M1.  NaFH_LEPS_JHCTP_2001:       2-state, P/G/D
    256. M2.  NaFH_LEPS_TTYPP1_1998:      2-state, P
    257. M3.  NaFH_LEPS_TTYPP2_1998:      2-state, P/G/D 
    258. M4.  NaFH_LEPS_JHCTP_2001_DPEM:  2-state, U/UG
    259. M5.  NaFH_LEPS_TTYPP1_1998_DPEM: 2-state, U/UG
    260. M6.  NaFH_LEPS_TTYPP2_1998_DPEM: 2-state, U/UG

    NaH2:
    261. M1.  NaH2_LEPS_7S_2000:        2-state, P
    262. M2.  NaH2_LEPS_7SD_2000:       2-state, P/G/D
    263. M3.  NaH2_LEPS_7L_2000:        2-state, P
    264. M4.  NaH2_LEPS_7LD_2000:       2-state, P/G/D
    265. M5.  NaH2_LEPS_6_2000:         2-state, P/G/D
    266. M6.  NaH2_LEPS_5G_1992:        2-state, P
    267. M7.  NaH2_LEPS_5F_1992:        2-state, P
    268. M8.  NaH2_LEPS_7S_2000_DPEM:   2-state, U
    269. M9.  NaH2_LEPS_7SD_2000_DPEM:  2-state, U/UG
    270. M10. NaH2_LEPS_7L_2000_DPEM:   2-state, U
    271. M11. NaH2_LEPS_7LD_2000_DPEM:  2-state, U/UG
    272. M12. NaH2_LEPS_6_2000_DPEM:    2-state, U/UG
    273. M13. NaH2_LEPS_5G_1992_DPEM:   2-state, U
    274. M14. NaH2_LEPS_5F_1992_DPEM:   2-state, U

    NH3:
    275. M1.  NH3_GEN_FFW1_2006:      2-state, P/G/D
    276. M2.  NH3_GEN_FFW2_2007:      2-state, P/G/D
    277. M3.  NH3_GEN_FFW1_2006_DPEM: 2-state, U/UG
    278. M4.  NH3_GEN_FFW2_2007_DPEM: 2-state, U/UG

    NH3CN:
    279. S1.  NH3CN_VBMM:            single-state, P/G

    NH3OH:
    280. S1.  NH3OH_VBMM:            single-state, P/G

    NH4:
    281. S1.  NH4_VBMM_1997:         single-state, P/G
    282. S2.  NH4_VBMM_2010:         single-state, P/G

    NO2:
    283. Z1.  NO2_2Ap_ZV:         single-state, doublet, A' symmetry, P/G
    284. Z2.  NO2_4Ap_ZV:         single-state, quartet, A' symmetry, P/G
    285. Z3.  NO2_6Ap_ZV:         single-state, sextet, A' symmetry, P/G

    O2H:
    286. S1.  O2H_GEN:            single-state, P/G

    O2H3:
    287. S1.  O2H3_PIPNN:         single-state, P

    O3:
    288. Z1.  O3_1_1Ap_ZV:      single-state, singlet, A' symmetry, P/G
    289. Z2.  O3_1_1App_ZV:     single-state, singlet, A'' symmetry, P/G
    290. Z3.  O3_1_3Ap_ZV:      single-state, triplet, A' symmetry, P/G
    291. Z4.  O3_1_3App_ZV:     single-state, triplet, A'' symmetry, P/G
    292. Z5.  O3_1_5Ap_ZV:      single-state, quintet, A' symmetry, P/G
    293. Z6.  O3_1_5App_ZV:     single-state, quintet, A'' symmetry, P/G
    294. Z7.  O3_2_1Ap_ZV:      single-state, singlet, A' symmetry, 1st excited state, P/G
    295. Z8.  O3_2_3Ap_ZV:      single-state, triplet, A' symmetry, 1st excited state, P/G
    296. Z9.  O3_2_5Ap_ZV:      single-state, quintet, A' symmetry, 1st excited state, P/G
    297. S1.  O3_1_1Ap:         single-state, singlet, A' symmetry, P/G
    298. S2.  O3_1_1App:        single-state, singlet, A'' symmetry, P/G
    299. S3.  O3_1_3Ap:         single-state, triplet, A' symmetry, P/G
    300. S4.  O3_1_3App:        single-state, triplet, A'' symmetry, P/G
    301. S5.  O3_1_5Ap:         single-state, quintet, A' symmetry, P/G
    302. S6.  O3_1_5App:        single-state, quintet, A'' symmetry, P/G
    303. S7.  O3_2_1Ap:         single-state, singlet, A' symmetry, 1st excited state, P/G
    304. S8.  O3_2_3Ap:         single-state, triplet, A' symmetry, 1st excited state, P/G
    305. S9.  O3_2_5Ap:         single-state, quintet, A' symmetry, 1st excited state, P/G
    306. S10. O3_MBP_2017:      single-state surface of O3, P
    307. MZ1. O3_14_3Ap_2022:       14-state, triplet, A' symmetry, P/G/D
    308. MZ2. O3_14_3Ap_2022_DPEM:  14-state, triplet, A' symmetry, U/UG
    309. M1.  O3_14_3Ap_2023:       14-state, triplet, A' symmetry, P/G/D
    310. M2.  O3_6_5Ap_2023:        6-state, quintet, A' symmetry, P/G/D
    311. M3.  O3_14_3Ap_2023_DPEM:  14-state, triplet, A' symmetry, U/UG
    312. M4.  O3_6_5Ap_2023_DPEM:   6-state, quintet, A' symmetry, U/UG

    O4:
    313. Z1.  O4_singlet_ZV:       single-state, singlet, P/G
    314. Z2.  O4_triplet_ZV:       single-state, triplet, P/G
    315. Z3.  O4_quintet_ZV:       single-state, quintet, P/G
    316. S1.  O4_singlet:          single-state, singlet, P/G
    317. S2.  O4_triplet_v2:       single-state, triplet, P/G
    318. S3.  O4_quintet:          single-state, quintet, P/G
    319. S4.  O4_MBP_singlet_2018: single-state, singlet, P

    OH2:
    320. S1.  OH2_Ap_GEN_M21986:    single-state, A' symmetry, P/G
    321. S2.  OH2_App_GEN_M21986:   single-state, A'' symmetry, P/G
    322. S3.  OH2_App_GEN_M3n1988:  single-state, A'' symmetry, P/G
    323. S4.  OH2_App_GEN_M3na1988: single-state, A'' symmetry, P/G
    324. S5.  OH2_GEN_DIM1981:      single-state, P/G
    325. S6.  OH2_LEPS_JW1977:      single-state, P/G
    326. S7.  OH2_LEPS_JWS1985:     single-state, P/G
    327. S8.  OH2_GEN_POL1980:      single-state, P/G
    328. S9.  OH2_GEN_SL1979:       single-state, P/G

    OH3:
    329. S1.  OH3_LEPS_SE_1986:        single-state, P/G
    330. M1.  OH3_PIP_FFW1_2019:       3-state, P/G/D
    331. M2.  OH3_PIP_FFW2_2022:       3-state, P/G/D  
    332. M3.  OH3_PIP_FFW1_2019_DPEM:  3-state, U/UG
    333. M4.  OH3_PIP_FFW2_2022_DPEM:  3-state, U/UG

    SiH4Cl:
    334. S1.   SiH4Cl_VBMM:       single-state, P/G

    SiH5:
    335. S1.   SiH5_LEP_1998v1:   single-state, P/G
    336. S2.   SiH5_LEP_1998v2:   single-state, P/G

    SiOH:
    337. S1.   SiOH_GEN_1988:     single-state, P/G

    YRH:
    338. M1.   YRH_LEPS_ModelSurface:        2-state surface of YRH, P/G/D
    339. M2.   YRH_LEPS_ModelSurface_DPEM:   2-state surface of YRH, U/UG

    C2N:
    340. S1.  C2N_PIPNN_Ap:     single-state, P
    341. S1.  C2N_PIPNN_App:    single-state, P


    """)

    return

def system():
    print("""

    ====Alphabetical List of Systems====
    ALC AlmHn Aln ArNO 
    BrH2 BrHCl 
    C2H2N2O C2H2O C2H4O4 C2H6Cl C2H6F C2H6H C2H6O 
    C2H6OH C2H7 COCO C3H7NO C6H5SH C6H6O C7H8S 
    CH2O CH2O2 CH2OH CH3NH2 CH3OH CH4 CH4Br CH4Cl 
    CH4CN CH4F CH4_H2O_H2O CH4O CH4OCl CH4OF CH4OH 
    CH5 CH5p Cl2H ClH2 ClH2O ClNH3 CH5O2
    F2H F2H2 FH2 FH2O 
    GeH4OH GeH5 
    H2CO H2O2 H2OBr HOBr H3 H3Cl H3plus H3S H4O2 H7p 
    H3ClO HClOH HO3 HOCO-Ar 
    I2H IH2 K2Rb2 LiFH malonaldehyde MCH MXH N2HOCp 
    N2O N2O2 N3 N4 NaFH NaH2 NH3 NH3CN NH3OH NH4 NinR 
    NO2 O2H O2H3 O3 O4 OH2 OH3 SiH4Cl SiH5 SiOH YRH
    C2N

        """)

    return


def list(system):
    module_name = f"chempotpy.{system}"
    module = importlib.import_module(module_name)
    module.intro()

    return

def detail(system):
    module_name = f"chempotpy.{system}"
    module = importlib.import_module(module_name)
    module.intro_detail()

    return


def get_script_dir():
    return os.path.dirname(os.path.abspath(inspect.stack()[1].filename))

requires_read_file_list=[
'NO2_2Ap_PIPNN', 'NO2_4Ap_PIPNN', 'NO2_6Ap_PIPNN',
'HO3_PIPNN', 'O2H3_PIPNN', 'FH2O_PIPNN', 'FH2O_SOC_PIPNN',
'N4_1A_PIPNN', 'CH4OCl_PIPNN', 'C2H2O_PIPNN', 'CH2O2_PIPNN',
'H3Cl_PIPNN', 'CH4OF_PIPNN', 'CH2O_PIPNN', 'C2H2N2O_PIPNN',
'N2HOCp_PIPNN', 'H3S_PIPNN', 'HClOH_PIPNN', 'ArNO_Ap_PIPNN',
'ArNO_App_PIPNN', 'COCO_PIPNN', 'H3ClO_PIPNN', 'K2Rb2_PIPNN',
'HOBr_3App_PIP', 'CH3OH_PIP', 'C3H7NO_PIP', 'CH4_H2O_H2O_PIP',
'ClH2_LEPS_TTGI2_1985','CH4Br_VBMM_2002', 'H4O2_PIP',
'F2H2_GEN_JBKKL_1990','H2CO_PIP','HOCO_Ar_PIP', 'CH4_GEN_SP_2001',
'CH4F_GEN_GC_1996', 'CH4F_GEN_GCSO_1996', 'N4_BDTC_2014',
'C2N_PIPNN_Ap', 'C2N_PIPNN_App'
]

Al_list=[
'Aln_BetH', 'Aln_CleR', 'Aln_CoxJM', 'Aln_deSPH_LJAT', 
'Aln_deSPH_M', 'Aln_ER2AT', 'Aln_ER2BA', 'Aln_ER2CN', 'Aln_ER2CoxJM',
'Aln_ER2EACN', 'Aln_ER2EAT', 'Aln_ER2EBA', 'Aln_ER2ECN', 'Aln_ER2ESCNa',
'Aln_ER2ESCNm', 'Aln_ER2ES', 'Aln_ER2', 'Aln_ER2GEA', 'Aln_ER2MeiD',
'Aln_ER2MisFMP', 'Aln_ER2SCNa', 'Aln_ER2SCNm', 'Aln_ER2S', 'Aln_ER',
'Aln_ErkIV', 'Aln_ErkV', 'Aln_ErkVIII', 'Aln_Gol', 'Aln_HalP', 'Aln_Jac',
'Aln_MeiD', 'Aln_MisFMP', 'Aln_NP_Ad', 'Aln_NP_A', 'Aln_NP_Bd', 'Aln_NP_B',
'Aln_NP_C', 'Aln_NP_D', 'Aln_NP_E', 'Aln_PapCEP', 'Aln_PapKEP', 'Aln_PeaTHT',
'Aln_PetW', 'Aln_reCoxJM', 'Aln_redeSPH_M', 'Aln_reErkIV', 'Aln_reErkVIII',
'Aln_reGEA', 'Aln_reGol', 'Aln_reHalP', 'Aln_reJac', 'Aln_reLJAT', 'Aln_reMeiD',
'Aln_reMisFMP', 'Aln_rePetW', 'Aln_reStrM', 'Aln_reSutC', 'Aln_StrM', 'Aln_SutC'
]

AlmHn_list=['AlmHn_VBO1', 'AlmHn_VBO2']
SiOH_list=['SiOH_GEN_1988']

multi_state_list=['C6H5SH_APRP_DPEM', 'PHOH_APRP_DPEM', 'C7H8S_APRP_DPEM', 'C7H8S_III_APRP_DPEM', 
'C7H8S_S22_APRP_DPEM', 'CH3NH2_APRP_DPEM', 'LiFH_LEPS_JHCTP_2001_DPEM', 
'LiFH_LEPS_JHTP1_2002_DPEM', 'LiFH_LEPS_JHTP2_2002_DPEM', 'NaFH_LEPS_JHCTP_2001_DPEM', 
'NaFH_LEPS_TTYPP1_1998_DPEM', 'NaFH_LEPS_TTYPP2_1998_DPEM', 'NaH2_LEPS_7S_2000_DPEM',
'NaH2_LEPS_7SD_2000_DPEM', 'NaH2_LEPS_7L_2000_DPEM', 'NaH2_LEPS_7LD_2000_DPEM', 
'NaH2_LEPS_6_2000_DPEM', 'NaH2_LEPS_5G_1992_DPEM', 'NaH2_LEPS_5F_1992_DPEM', 
'BrH2_LEPS_LTBZSO_1995_DPEM', 'BrH2_LEPS_LTBZ2_1995_DPEM',
'NH3_GEN_FFW1_2006_DPEM', 'NH3_GEN_FFW2_2007_DPEM', 'O3_14_3Ap_2022_DPEM', 'O3_14_3Ap_2023_DPEM',
'O3_6_5Ap_2023_DPEM', 'OH3_PIP_FFW1_2019_DPEM', 'OH3_PIP_FFW2_2022_DPEM', 'MCHWL_LEPS_ModelSurface_DPEM',
'MCHWB_LEPS_ModelSurface_DPEM', 'MCHSL_LEPS_ModelSurface_DPEM', 'MCHSB_LEPS_ModelSurface_DPEM',
'MCHTL_LEPS_ModelSurface_DPEM', 'MXHWL_LEPS_ModelSurface_DPEM', 'MXHSL_LEPS_ModelSurface_DPEM',
'MXHSB_LEPS_ModelSurface_DPEM', 'YRH_LEPS_ModelSurface_DPEM']


parent_path=get_script_dir()
def p(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in requires_read_file_list:
        v, g, d = module.pes(xyz, 0, parent_path)
    elif surface in Al_list:
        natoms=int(len(geom))
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, natoms, 0)
    elif surface in AlmHn_list:
        natoms=int(len(geom))
        n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 0)
    elif surface in SiOH_list:
        natoms=int(len(geom))
        n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 0)
    else:
        v, g, d = module.pes(xyz, 0)
    v, g, d = module0.reverse_order(v, g, d, geom, rank_ordered)
    return v

def pg(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in requires_read_file_list:
        v, g, d = module.pes(xyz, 1, parent_path)
    elif surface in Al_list:
        natoms=int(len(geom))
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, natoms, 1)
    elif surface in AlmHn_list:
        natoms=int(len(geom))
        n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 1)
    elif surface in SiOH_list:
        natoms=int(len(geom))
        n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 1)
    else:
       v, g, d = module.pes(xyz, 1)
    v, g, d = module0.reverse_order(v, g, d, geom, rank_ordered)
    return v, g

def pgd(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in requires_read_file_list:
        v, g, d = module.pes(xyz, 2, parent_path)
    elif surface in Al_list:
        natoms=int(len(geom))
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, natoms, 2)
    elif surface in AlmHn_list:
        natoms=int(len(geom))
        n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 2)
    elif surface in SiOH_list:
        natoms=int(len(geom))
        n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 2)
    else:
       v, g, d = module.pes(xyz, 2)
    v, g, d = module0.reverse_order(v, g, d, geom, rank_ordered)
    return v, g, d

#function to compute numerical gradient
def pn(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in requires_read_file_list:
        v, g, d = module.pes(xyz, 0, parent_path)
    elif surface in Al_list:
        natoms=int(len(geom))
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, natoms, 0)
    elif surface in AlmHn_list:
        natoms=int(len(geom))
        n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 0)
    elif surface in SiOH_list:
        natoms=int(len(geom))
        n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
        xyz_full=np.zeros((10000, 3))
        for i in range(natoms):
            for idir in range(3):
                xyz_full[i][idir]=xyz[i][idir]
        v, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 0)
    else:
       v, g, d = module.pes(xyz, 0)

    n_atoms=len(xyz)
    n_states=len(v)
    grad=np.zeros((n_states,n_atoms,3))
    for iatom in range(n_atoms):
        for idir in range(3):
            xyz_displaced=np.copy(xyz)
            xyz_displaced[iatom][idir]=xyz_displaced[iatom][idir]+0.001
            if surface in requires_read_file_list:
                v_plus, g, d = module.pes(xyz_displaced, 0, parent_path)
            elif surface in Al_list:
                natoms=int(len(geom))
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_plus, g, d = module.pes(xyz_full, natoms, 0)
            elif surface in AlmHn_list:
                natoms=int(len(geom))
                n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_plus, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 0)
            elif surface in SiOH_list:
                natoms=int(len(geom))
                n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_plus, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 0)
            else:
                v_plus, g, d = module.pes(xyz_displaced, 0)
            xyz_displaced=np.copy(xyz)
            xyz_displaced[iatom][idir]=xyz_displaced[iatom][idir]-0.001
            if surface in requires_read_file_list:
                v_minus, g, d = module.pes(xyz_displaced, 0, parent_path)
            elif surface in Al_list:
                natoms=int(len(geom))
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_minus, g, d = module.pes(xyz_full, natoms, 0)
            elif surface in AlmHn_list:
                natoms=int(len(geom))
                n_hydrogen, n_aluminium = module0.nAlnH(system, surface, geom)
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_minus, g, d = module.pes(xyz_full, n_hydrogen, n_aluminium, natoms, 0)
            elif surface in SiOH_list:
                natoms=int(len(geom))
                n_silicon, n_oxygen, n_hydrogen = module0.nSinOnH(system, surface, geom)
                xyz_full=np.zeros((10000, 3))
                for i in range(natoms):
                    for j in range(3):
                        xyz_full[i][j]=xyz_displaced[i][j]
                v_minus, g, d = module.pes(xyz_full, n_silicon, n_oxygen, n_hydrogen, natoms, 0)
            else:
                v_minus, g, d = module.pes(xyz_displaced, 0)
            for istate in range(n_states):
                grad[istate][iatom][idir]=(v_plus[istate]-v_minus[istate])/0.002
    v, grad, d = module0.reverse_order(v, grad, d, geom, rank_ordered)
    return v, grad


def u(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in multi_state_list:
        if surface in requires_read_file_list:
            u, ug = module.dpem(xyz, 0, parent_path)
        else:
            u, ug = module.dpem(xyz, 0)
        u, ug = module0.reverse_order_dpem(u, ug, geom, rank_ordered)
        return u
    else:
        raise Exception("single-state surface, no DPEM available")

def ug(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in multi_state_list:
        if surface in requires_read_file_list:
            u, ug = module.dpem(xyz, 1, parent_path)
        else:
            u, ug = module.dpem(xyz, 1)
        u, ug = module0.reverse_order_dpem(u, ug, geom, rank_ordered)
        return u, ug
    else:
        raise Exception("single-state surface, no DPEM available")

#function to compute numerical gradient of DPEM
def un(system, surface, geom):
    module_name0 = f"chempotpy.{system}"
    module0 = importlib.import_module(module_name0)
    xyz, rank_ordered = module0.check(system, surface, geom)
    module_name = f"chempotpy.{system}.{surface}"
    module = importlib.import_module(module_name)
    if surface in multi_state_list:
        if surface in requires_read_file_list:
            u, ug = module.dpem(xyz, 0, parent_path)
        else:
            u, ug = module.dpem(xyz, 0)
        u, ug = module0.reverse_order_dpem(u, ug, geom, rank_ordered)
    else:
        raise Exception("single-state surface, no DPEM available")

    n_atoms=len(xyz)
    n_states=len(u)
    grad=np.zeros((n_states,n_states,n_atoms,3))
    for iatom in range(n_atoms):
        for idir in range(3):
            xyz_displaced=np.copy(xyz)
            xyz_displaced[iatom][idir]=xyz_displaced[iatom][idir]+0.001
            if surface in requires_read_file_list:
                u_plus, ug = module.dpem(xyz_displaced, 0, parent_path)
            else:
                u_plus, ug = module.dpem(xyz_displaced, 0)
            xyz_displaced=np.copy(xyz)
            xyz_displaced[iatom][idir]=xyz_displaced[iatom][idir]-0.001
            if surface in requires_read_file_list:
                u_minus, ug = module.dpem(xyz_displaced, 0, parent_path)
            else:
                u_minus, ug = module.dpem(xyz_displaced, 0)
            for istate in range(n_states):
                for jstate in range(n_states):
                    grad[istate][jstate][iatom][idir]=(u_plus[istate][jstate]-u_minus[istate][jstate])/0.002
    u, grad = module0.reverse_order_dpem(u, grad, geom, rank_ordered)
    return u, grad

