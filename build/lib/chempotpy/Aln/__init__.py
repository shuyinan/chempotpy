'''
Aln
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for Aln:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   Aln_BetH:           single states surface of Aln, ground state,     P
    2.   Aln_CleR:           single states surface of Aln, ground state,     P
    3.   Aln_CoxJM:          single states surface of Aln, ground state,     P
    4.   Aln_deSPH_LJAT:     single states surface of Aln, ground state,     P
    5.   Aln_deSPH_M:        single states surface of Aln, ground state,     P
    6.   Aln_ER2AT:          single states surface of Aln, ground state,     P
    7.   Aln_ER2BA:          single states surface of Aln, ground state,     P
    8.   Aln_ER2CN:          single states surface of Aln, ground state,     P
    9.   Aln_ER2CoxJM:       single states surface of Aln, ground state,     P
    10.  Aln_ER2EACN:        single states surface of Aln, ground state,     P
    11.  Aln_ER2EAT:         single states surface of Aln, ground state,     P
    12.  Aln_ER2EBA:         single states surface of Aln, ground state,     P
    13.  Aln_ER2ECN:         single states surface of Aln, ground state,     P
    14.  Aln_ER2ESCNa:       single states surface of Aln, ground state,     P
    15.  Aln_ER2ESCNm:       single states surface of Aln, ground state,     P
    16.  Aln_ER2ES:          single states surface of Aln, ground state,     P
    17.  Aln_ER2:            single states surface of Aln, ground state,     P
    18.  Aln_ER2GEA:         single states surface of Aln, ground state,     P
    19.  Aln_ER2MeiD:        single states surface of Aln, ground state,     P
    20.  Aln_ER2MisFMP:      single states surface of Aln, ground state,     P
    21.  Aln_ER2SCNa:        single states surface of Aln, ground state,     P
    22.  Aln_ER2SCNm:        single states surface of Aln, ground state,     P
    23.  Aln_ER2S:           single states surface of Aln, ground state,     P
    24.  Aln_ER:             single states surface of Aln, ground state,     P
    25.  Aln_ErkIV:          single states surface of Aln, ground state,     P
    26.  Aln_ErkV:           single states surface of Aln, ground state,     P
    27.  Aln_ErkVIII:        single states surface of Aln, ground state,     P
    28.  Aln_Gol:            single states surface of Aln, ground state,     P
    29.  Aln_HalP:           single states surface of Aln, ground state,     P
    30.  Aln_Jac:            single states surface of Aln, ground state,     P
    31.  Aln_MeiD:           single states surface of Aln, ground state,     P
    32.  Aln_MisFMP:         single states surface of Aln, ground state,     P
    33.  Aln_NP_Ad:          single states surface of Aln, ground state,     P/G
    34.  Aln_NP_A:           single states surface of Aln, ground state,     P
    35.  Aln_Np_Bd:          single states surface of Aln, ground state,     P/G
    36.  Aln_NP_B:           single states surface of Aln, ground state,     P
    37.  Aln_NP_C:           single states surface of Aln, ground state,     P
    38.  Aln_NP_D:           single states surface of Aln, ground state,     P
    39.  Aln_NP_E:           single states surface of Aln, ground state,     P
    40.  Aln_PapCEP:         single states surface of Aln, ground state,     P
    41.  Aln_PapKEP:         single states surface of Aln, ground state,     P
    42.  Aln_PeaTHT:         single states surface of Aln, ground state,     P
    43.  Aln_PetW:           single states surface of Aln, ground state,     P
    44.  Aln_reCoxJM:        single states surface of Aln, ground state,     P
    45.  Aln_redeSPH_M:      single states surface of Aln, ground state,     P
    46.  Aln_reErkIV:        single states surface of Aln, ground state,     P
    47.  Aln_reErkVIII:      single states surface of Aln, ground state,     P
    48.  Aln_reGEA:          single states surface of Aln, ground state,     P
    49.  Aln_reGol:          single states surface of Aln, ground state,     P
    50.  Aln_reHalP:         single states surface of Aln, ground state,     P
    51.  Aln_reJac:          single states surface of Aln, ground state,     P
    52.  Aln_reLJAT:         single states surface of Aln, ground state,     P
    53.  Aln_reMeiD:         single states surface of Aln, ground state,     P
    54.  Aln_reMisFMP:       single states surface of Aln, ground state,     P
    55.  Aln_rePetW:         single states surface of Aln, ground state,     P
    56.  Aln_reStrM:         single states surface of Aln, ground state,     P
    57.  Aln_reSutC:         single states surface of Aln, ground state,     P
    58.  Aln_StrM:           single states surface of Aln, ground state,     P
    59.  Aln_SutC:           single states surface of Aln, ground state,     P

 
        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for ArNO:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

 
    Single-State Surfaces:
    1.   Aln_BetH:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: BetH.f 
                             ref: B. Betz, and W. Husinsky,
                                  "Cluster bombardment of solids: A molecular dynamics study",
                                  Nucl. Instrum. Methods B 122, 311-317 (1997).
    ============================================================================================
    2.   Aln_CleR:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: CleR.f
                             ref: F. Cleri, and V. Rosato,
                                  "Tight-binding potentials for transition metals and alloys"
                                  Phys. Rev. B 48, 22-33 (1993). 
    ============================================================================================
    3.   Aln_CoxJM:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: CoxJM.f
                             ref: H. Cox, R. L. Johnston, and J. N. Murrell, 
                                  "Modelling of surface relaxation and melting of aluminium",
                                  Surf. Sci., 373, 67-84 (1997).
    ============================================================================================
    4.   Aln_deSPH_LJAT:     single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: deSPH-LJAT.f
                             ref: P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase,
                                  "Energy Transfer Dynamics in the Collision-Induced Dissociation 
                                  of Al6 and Al13 Clusters",
                                  J. Phys. Chem., 99, 8147-8161 (1995).
    ============================================================================================
    5.   Aln_deSPH_M:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: deSPH-M.f 
                             ref: P. de Sainte Claire, G. H. Peslherbe, and W. L. Hase,
                                  "Energy Transfer Dynamics in the Collision-Induced Dissociation 
                                  of Al6 and Al13 Clusters",
                                  J. Phys. Chem., 99, 8147-8161 (1995).
    ============================================================================================
    6.   Aln_ER2AT:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2AT.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    7.   Aln_ER2BA:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2BA.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    8.   Aln_ER2CN:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2CN.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    9.   Aln_ER2CoxJM:       single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2CoxJM.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    10.  Aln_ER2EACN:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2EACN.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    11.  Aln_ER2EAT:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2EAT.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    12.  Aln_ER2EBA:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2EBA.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    13.  Aln_ER2ECN:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2ECN.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    14.  Aln_ER2ESCNa:       single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2ESCNa.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    15.  Aln_ER2ESCNm:       single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2ESCNm.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    16.  Aln_ER2ES:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2ES.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    17.  Aln_ER2:            single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    18.  Aln_ER2GEA:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2GEA.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    19.  Aln_ER2MeiD:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2MeiD.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    20.  Aln_ER2MisFMP:      single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2MisFMP.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    21.  Aln_ER2SCNa:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2SCNa.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    22.  Aln_ER2SCNm:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2SCNm.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    23.  Aln_ER2S:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER2S.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    24.  Aln_ER:             single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ER.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    25.  Aln_ErkIV:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ErkIV.f
                             ref: Y. Tahtamoni, S. Erkoc,
                                  "Structural Stability and Energetics of Xn (X  = V, Cr, Nb; 
                                  n = 3 to 7) Microclusters: Empirical Many-Body Potential Energy 
                                  Funtion Calculation",
                                  Phys. Stat. Sol. (b) 161, K5-K8 (1990).
    ============================================================================================
    26.  Aln_ErkV:           single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ErkV.f
                             ref: S. Erkoc,
                                  "Empirical Potential Energy Functions Used in the Simulations 
                                  of Materials Properties",
                                  Annual Reviews of Computational Physics IX, 
                                  D. Stauffer, ed., World Scientific, Singapore 2001, pp. 1-103.
    ============================================================================================
    27.  Aln_ErkVIII:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: ErkVIII.f
                             ref: S. Erkoc,
                                  "Empirical Potential Energy Functions Used in the Simulations 
                                  of Materials Properties",
                                  Annual Reviews of Computational Physics IX, 
                                  D. Stauffer, ed., World Scientific, Singapore 2001, pp. 1-103.
    ============================================================================================
    28.  Aln_Gol:            single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: Gol.f
                             ref: H. Gollisch,
                                  "Effective binding potentials for the description of metal-metal 
                                  interaction",
                                  Surface Science 166, 87-100 (1986).
    ============================================================================================
    29.  Aln_HalP:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: HalP.f
                             ref: T. Halicioglu, and G.M. Pound,
                                  "Calculation of Potential Energy Parameters from Crystalline State 
                                  Properties",
                                  Phys. Stat. Sol. (b) 30, 619-623 (1975).
    ============================================================================================
    30.  Aln_Jac:            single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: Jac.f
                             ref: K. W. Jacobsen,
                                  Comments Cond. Mat. Phys. 14, 129 (1988). 
    ============================================================================================
    31.  Aln_MeiD:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: MeiD.f
                             ref: J. Mei, and J. W. Davenport,
                                  "Free-energy calculations and the melting point of Al",
                                  Phys. Rev. B 46, 21-25 (1992).
    ============================================================================================
    32.  Aln_MisFMP:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: MisFMP.f
                             ref: Y. Mishin, D. Farkas, M. J. Mehl, and D. A. Papaconstantopoulos,
                                  "Interatomic Potentials for Al and Ni From Experimental Data and 
                                  AB Initio Calculations",
                                  Mat. Res. Soc. Symp. Proc. 583, 535, (1999).
    ============================================================================================
    33.  Aln_NP_Ad:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-Ad.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    34.  Aln_NP_A:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-A.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    35.  Aln_Np_Bd:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-Bd.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    36.  Aln_NP_B:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-B.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    37.  Aln_NP_C:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-C.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    38.  Aln_NP_D:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-D.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    39.  Aln_NP_E:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: NP-E.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    40.  Aln_PapCEP:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: PapCEP.f
                             ref: N. I. Papanicolaou, H. Chamati, G. A. Evangelakis, D. A. Papaconstantopoulos,
                                  "Second-moment interatomic potential for Al, Ni and Ni–Al alloys, 
                                  and molecular dynamics application",
                                  Comput. Mater. Sci. 27, 191-198 (2003).
    ============================================================================================
    41.  Aln_PapKEP:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: PapKEP.f
                             ref: N. I. Papanicolaou, G. C. Kallinteris, G. A. Evangelakis, 
                                  and D. A. Papaconstantopoulos,
                                  "Second-moment interatomic potential for aluminum derived from 
                                  total-energy calculations and molecular dynamics application",
                                  Comput. Mater. Sci. 17, 224-229 (2000).
    ============================================================================================
    42.  Aln_PeaTHT:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: PeaTHT.f
                             ref: E. Pearson, T. Takai, T. Halicioglu, and W. A. Tiller,
                                  "Computer modeling of Si and SiC surfaces and surface processes 
                                  relevant to crystal growth from the vapor", 
                                  J. Cryst. Growth 70, 33-40 (1984). 
    ============================================================================================
    43.  Aln_PetW:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: PetW.f
                             ref: D. G. Pettifor, and M. A. Ward,
                                  "An analytic pair potential for simple metals",
                                  Solid State Commun. 49, 291-294 (1984). 
    ============================================================================================
    44.  Aln_reCoxJM:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reCoxJM.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    45.  Aln_redeSPH_M:      single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: redeSPH-M.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    46.  Aln_reErkIV:        single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reErkIV.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    47.  Aln_reErkVIII:      single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reErkVIII.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    48.  Aln_reGEA:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reGEA.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    49.  Aln_reGol:          single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reGol.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    50.  Aln_reHalP:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reHalP.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    51.  Aln_reJac:          single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reJac.f 
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    52.  Aln_reLJAT:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reLJAT.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    53.  Aln_reMeiD:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reMeiD.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    54.  Aln_reMisFMP:       single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reMisFMP.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    55.  Aln_rePetW:         single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: rePetW.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    56.  Aln_reStrM:         single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reStrM.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    57.  Aln_reSutC:         single states surface of Aln, ground state, 
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: reSutC.f
                             ref: A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz, 
                                  and D. G. Truhlar,
                                  "Analytic Potential Energy Functions for Aluminum Clusters", 
                                  J. Phys. Chem. B 108, 8996-9010 (2004).
    ============================================================================================
    58.  Aln_StrM:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: StrM.f
                             ref: F. H. Streitz, and J. W. Mintmire, 
                                  "Electrostatic potentials for metal-oxide surfaces and interfaces",
                                  Phys. Rev. B, 50, 11996-12003 (1994).
    ============================================================================================
    59.  Aln_SutC:           single states surface of Aln, ground state,
                             availability: potential energy
                             functional form: GEN
                             corresponding surface in POTLIB: SutC.f
                             ref: A. P. Sutton, and J. Chen,
                                  "Long-range Finnis–Sinclair potentials",
                                  Philos. Mag. Lett. 61, 139-146 (1990).


        """)

    return

def check(system, surface, geom):

    n_aluminium=0
    natoms=len(geom)

    if natoms==1:
        print("number of atoms can not equal 1")
        sys.exit()
    if natoms>10000:
        print("number of atoms can not exceed 10000")
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
