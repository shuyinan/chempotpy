****************************************************************************
* The PIP-NN PESs for ground C2N(^2A',^2A") system at MC-PDFT/AVTZ level.
* Ref: Junxiang Zuo et al., J. Chem. Theory Comput., 2022, 18 (12), 7121-7131.
****************************************************************************

These files are for the PESs of C2N system, including
  1. Main fortran files containing subroutines:
	(1) 'c2n_2a_prime.f' for C2N(^2A');
        (2) 'c2n_2a_dprime.f' for C2N(^2A").
  2. Parameter files:
        (1) c2n_2a_prime_para_a.txt
        (2) c2n_2a_prime_para_b.txt
        (3) c2n_2a_prime_para_c.txt
        (4) c2n_2a_dprime_para_a.txt
        (5) c2n_2a_dprime_para_b.txt
        (6) c2n_2a_dprime_para_c.txt
----------------------------------------------------------------------------
  To use C2N(^2A') PES, please
        call c2n_2a_prime(rNC,rCC,theta,Vpes)
  and for C2N(^2A") PES, please
        call c2n_2a_dprime(rNC,rCC,theta,Vpes)
  where,
        rNC  : The bond length of NC in angstrom;
	rCC  : The bond length of CC in angstrom;
        theta: The angle of N-C-C in degree;
        Vpes : total energy in eV.
****************************************************************************
