      subroutine pes(xyz,natoms,igrad,p,g,d)
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, intent(in) :: natoms
      double precision, intent(in) :: xyz(10000,3)
      integer, intent(in) :: igrad
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: cx(natoms), cy(natoms), cz(natoms)
      double precision :: dvdx(natoms), dvdy(natoms), dvdz(natoms)
      double precision :: v
      integer :: istate, iatom, i

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
         cx(iatom)=xyz(iatom,1)/0.529177211
         cy(iatom)=xyz(iatom,2)/0.529177211
         cz(iatom)=xyz(iatom,3)/0.529177211
      enddo
 
      call pot(cx, cy, cz, v, natoms, 10000)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine
      subroutine pot(x,y,z,v,natom,maxatom)

C   System:                        Aln
C   Functional form:               Extended Rydberg with Screening and Coordination Number
C   Common name:                   NP-A
C   Number of derivatives:         1
C   Number of bodies:              variable
C   Number of electronic surfaces: 1
C   Interface:                     HO-MM-1
C
c  Notes:  Many-body aluminum potential energy function.  The functional
c           form is from Ref. 2.  The parameters were re-optimized in Ref. 1
c           against a data set of energies for aluminum clusters and
c           nanoparticles and bulk data.  This PEF has a cutoff at 6.5 A.
c
c  References: (1) A. W. Jasper, N. E. Schultz, and D. G. Truhlar, "Analytic
c              Potential Energy Functions for Simulating Aluminum
c              Nanoparticles," in preparation.
c              (2) A. W. Jasper, P. Staszewski, G. Staszewska, N. E. Schultz,
c              and D. G. Truhlar, "Analytic Potential Energy Functions
c              for Aluminum Clusters," J. Phys. Chem. B 108, 8996(2004).
C
C  Units:
C       energies    - hartrees
C       coordinates - bohrs
C
C  v       --- energy of the system
C  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
C                    of the system
C  dvdx,dvdy,dvdz   --- one-dimensional arrays of the gradients with respect
C                       to the Cartesian coordinates of the system (output)
C  natom   --- number of atoms in the system
C  maxatom --- maximal number of atoms
C
C  Note:  Derivatives were generated using ADIFOR, whose disclaimer follows.
C
C                           DISCLAIMER
C
C   This file was generated on 09/13/04 by the version of
C   ADIFOR compiled on June, 1998.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
        implicit double precision(a-h, o-z)
        parameter (autoev = 27.2113961d0)
        parameter (autoang = 0.529177249d0)
        dimension x(maxatom),y(maxatom),z(maxatom)
        double precision gxs(maxatom),gys(maxatom),gzs(maxatom),
     &  gxgcorr(maxatom,maxatom),gygcorr(maxatom,maxatom),
     &  gzgcorr(maxatom,maxatom),gcorr(maxatom)
        double precision dvdy(maxatom),dvdz(maxatom),dvdx(maxatom)

c Nanoparticle parameters (Ref. 2)
c       two-body interaction
        de = 1.71013678553981441 / autoev
        re = 5.08182706399609163
        a1 = 1.24074255007327805
        a2 = 0.551880801172447422
        a3 = 0.129970688812896917
        au2 = 0.143243771372740580
        bu2 = 6.50000000000000000 / autoang
c       screening
        xk1= 4.24002677622442103
        xk2= 0.117656503960960862
        xk3= 4.78063179546451522
        au23= 1.63973192904916298
c       coordination number
        g1= 0.708483373073205747
        du2= 1.13286279334603357
        gu2= 0.663930057862113232
        gzero= 8.54498572971970027
        g2= 5.39584023677170066/autoang
c       effective two-body
        deb= 1.42526928794948882/autoev
        reb= 4.87735706664722812
        a1b= 1.20666644170640880
        a2b= 0.728296669115275908
        a3b= 0.215461507389864804
        au2b= 0.138211749991007299
        bu2b= 6.50000000000000000/autoang

c Initialize arrays
        do 99998 i = 1, natom
        dvdx(i) = 0.0d0
        dvdy(i) = 0.0d0
        dvdz(i) = 0.0d0
        v = 0.d0
99998   continue

C do S term, needs xk1,xk2,xk3
        do 510 i = 1, natom
          do 515 j = i + 1, natom
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            d1_w = dx * dx + dy * dy + dz * dz
            d2_v = sqrt(d1_w)
            rij = d2_v
            if (rij .ge. bu2) goto 515
            do k = 1, natom
            gxs(k) = 0.0d0
            gys(k) = 0.0d0
            gzs(k) = 0.0d0
            enddo
            s = 0.d0
            do 512 l = 1, natom
              if (l .eq. j .or. l .eq. i) goto 512
              dx = x(i) - x(l)
              dy = y(i) - y(l)
              dz = z(i) - z(l)
              d1_w = dx * dx + dy * dy + dz * dz
              d2_v = sqrt(d1_w)
              ril = d2_v
              dx = x(j) - x(l)
              dy = y(j) - y(l)
              dz = z(j) - z(l)
              d1_w = dx * dx + dy * dy + dz * dz
              d2_v = sqrt(d1_w)
              rjl = d2_v
              if (rjl .ge. bu2) goto 512
              if (ril .ge. bu2) goto 512

              rrr = (ril + rjl) / rij
              d3_v = 1d0 - rjl / bu2
              d6_b1 = -au23 / (d3_v**2 * bu2)
              ac = au23 * (1d0 - 1.d0/d3_v)
              d4_v = 1d0 - rij / bu2
              d5_v = 1d0 / d4_v
              d8_b1 = (-((-au23) * ((-d5_v) / d4_v))) * (1.0d0 / bu2)
              ac = ac + au23 * (1d0 - d5_v)
              d4_v = 1d0 - ril / bu2
              d5_v = 1d0 / d4_v
              d8_b2 = (-((-au23) * ((-d5_v) / d4_v))) * (1.0d0 / bu2)
              ac = ac + au23 * (1d0 - d5_v)
              term = xk1 * dexp(-xk2 * rrr ** xk3 + ac)
              d5_b = -xk2 * xk3 * rrr ** ( xk3 - 1.0d0)
              d6_b = d5_b * (1.0d0 / rij) * term
              d7_b = -d5_b * (rrr / rij) * term

              gxs(i) = gxs(i) + (term*d8_b1+d7_b)*(x(i)-x(j))/rij
              gxs(i) = gxs(i) + (term*d8_b2+d6_b)*(x(i)-x(l))/ril
              gxs(j) = gxs(j) + (term*d6_b1+d6_b)*(x(j)-x(l))/rjl
              gxs(j) = gxs(j) + (term*d8_b1+d7_b)*(x(j)-x(i))/rij
              gxs(l) = gxs(l) + (term*d8_b2+d6_b)*(x(l)-x(i))/ril
              gxs(l) = gxs(l) + (term*d6_b1+d6_b)*(x(l)-x(j))/rjl

              gys(i) = gys(i) + (term*d8_b1+d7_b)*(y(i)-y(j))/rij
              gys(i) = gys(i) + (term*d8_b2+d6_b)*(y(i)-y(l))/ril
              gys(j) = gys(j) + (term*d6_b1+d6_b)*(y(j)-y(l))/rjl
              gys(j) = gys(j) + (term*d8_b1+d7_b)*(y(j)-y(i))/rij
              gys(l) = gys(l) + (term*d8_b2+d6_b)*(y(l)-y(i))/ril
              gys(l) = gys(l) + (term*d6_b1+d6_b)*(y(l)-y(j))/rjl

              gzs(i) = gzs(i) + (term*d8_b1+d7_b)*(z(i)-z(j))/rij
              gzs(i) = gzs(i) + (term*d8_b2+d6_b)*(z(i)-z(l))/ril
              gzs(j) = gzs(j) + (term*d6_b1+d6_b)*(z(j)-z(l))/rjl
              gzs(j) = gzs(j) + (term*d8_b1+d7_b)*(z(j)-z(i))/rij
              gzs(l) = gzs(l) + (term*d8_b2+d6_b)*(z(l)-z(i))/ril
              gzs(l) = gzs(l) + (term*d6_b1+d6_b)*(z(l)-z(j))/rjl

              s = s + term

512           continue

            d2_v = tanh (s)
            d1_p = 1.0d0 - ( d2_v *  d2_v)
            fs = d2_v

c effective two-body
            rr = rij
            rhob = rr-reb
            polyb=1.d0+a1b*rhob+a2b*rhob**2+a3b*rhob**3
            eeb=dexp(-a1b*rhob)
            cob=0.d0
            if (rr.le.bu2b) cob=dexp(au2b*(1.d0-1.d0/(1.d0-rr/bu2b)))
            tmpb=deb*eeb*polyb*cob
            dcobdrr=0.d0
            if (rr.le.bu2b) dcobdrr=(-cob/bu2b)*au2b/(1.d0-rr/bu2b)**2
            dpolybdrr=a1b+2.d0*a2b*rhob+3.d0*a3b*rhob**2
            deebdrr=-a1b*eeb
            dtmpbdrr=deb*(eeb*polyb*dcobdrr+eeb*dpolybdrr*cob
     &         +deebdrr*polyb*cob)

            dvdx(i) = dvdx(i) + dtmpbdrr*(x(i)-x(j))/rr*fs
            dvdx(j) = dvdx(j) + dtmpbdrr*(x(j)-x(i))/rr*fs
            dvdy(i) = dvdy(i) + dtmpbdrr*(y(i)-y(j))/rr*fs
            dvdy(j) = dvdy(j) + dtmpbdrr*(y(j)-y(i))/rr*fs
            dvdz(i) = dvdz(i) + dtmpbdrr*(z(i)-z(j))/rr*fs
            dvdz(j) = dvdz(j) + dtmpbdrr*(z(j)-z(i))/rr*fs

            do k = 1, natom
            dvdx(k) = dvdx(k) + d1_p*tmpb*gxs(k)
            dvdy(k) = dvdy(k) + d1_p*tmpb*gys(k)
            dvdz(k) = dvdz(k) + d1_p*tmpb*gzs(k)
            enddo

            v = v + tmpb * d2_v
515       continue
510     continue

C CALCULATE CN FUNCTION
        g1e = dexp(g1)
        do 20 i = 1, natom
          do j = 1, natom
            gxgcorr(j,i) = 0.0d0
            gygcorr(j,i) = 0.0d0
            gzgcorr(j,i) = 0.0d0
          enddo
          gcorr(i) = 0.d0
          do 21 j = 1, natom
            if (i .eq. j) goto 21
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            d1_w = dx * dx + dy * dy + dz * dz
            d2_v = sqrt(d1_w)
            rij = d2_v
            if (rij .ge. g2) goto 21
            d3_v = 1d0 - rij / g2
            d4_v = (-g1) / d3_v
            d4_b = (-((-d4_v) / d3_v)) * (1.0d0 / g2)
            d1_w = d4_v
            d2_v = dexp(d1_w)
            d1_p =  d2_v
            dum = g1e*d1_p*d4_b/rij
            term = d2_v
            gxgcorr(i,i) =  dum*(x(i)-x(j)) + gxgcorr(i,i)
            gxgcorr(j,i) = -dum*(x(i)-x(j)) + gxgcorr(j,i)
            gygcorr(i,i) =  dum*(y(i)-y(j)) + gygcorr(i,i)
            gygcorr(j,i) = -dum*(y(i)-y(j)) + gygcorr(j,i)
            gzgcorr(i,i) =  dum*(z(i)-z(j)) + gzgcorr(i,i)
            gzgcorr(j,i) = -dum*(z(i)-z(j)) + gzgcorr(j,i)
            gcorr(i) = gcorr(i) + term
21          continue
          gcorr(i) = gcorr(i) * g1e
20      continue

C     then calculate CN term
        do 30 i = 1, natom
          do 31 j = 1, natom
            if (i.eq.j) goto 31
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            d1_w = dx * dx + dy * dy + dz * dz
            d2_v = sqrt(d1_w)
            rij = d2_v
            if (rij.ge.bu2b) goto 31

            term = 0.d0
            gdum = 0.d0
            if (rij.ge.g2) goto 32
            d3_v = 1.d0-rij/g2
            d4_v = 1.d0/d3_v
            d6_b = -g1*d4_v/(d3_v*g2)
            d1_w = g1 * (1.d0 - d4_v)
            d2_v = dexp(d1_w)
            d1_p =  d2_v
            gdum = d1_p * d6_b /rij
            term = d2_v
32          continue

            cz1 = (gcorr(i) - term) / gzero
            if (cz1.gt.1.d-20) then
               d2_v = cz1 ** ( gu2 - 2.0d0)
               d2_v =  d2_v * cz1
               d1_p =  gu2 *  d2_v
               d2_v =  d2_v * cz1
            else
               cz1 = 0.d0
               d2_v = 0.d0
               d1_p = 0.d0
            endif
            d3_v = 1d0 + d2_v
            d4_ba = -d1_p / d3_v**2
            es1 = 1d0 / d3_v

            cz2 = (gcorr(j) - term) / gzero
            if (cz2.gt.1.d-20) then
               d2_v = cz2 ** ( gu2 - 2.0d0)
               d2_v =  d2_v * cz2
               d1_p =  gu2 *  d2_v
               d2_v =  d2_v * cz2
            else
               cz2 = 0.d0
               d2_v = 0.d0
               d1_p = 0.d0
            endif
            d3_v = 1d0 + d2_v
            d4_v = 1d0 / d3_v
            d4_bb = (-d4_v) / d3_v * d1_p
            ec1 = d4_v

            tt1 = du2 * ec1 * d4_ba
            tt2 = du2 * es1 * d4_bb

c modified two-body
              rr = rij
              rhob = rr-reb
              polyb=1.d0+a1b*rhob+a2b*rhob**2+a3b*rhob**3
              eeb=dexp(-a1b*rhob)
              cob=0.d0
              if (rr.le.bu2b) cob=dexp(au2b*(1.d0-1.d0/(1.d0-rr/bu2b)))
              tmpb=0.5d0*deb*eeb*polyb*cob
              dcobdrr=0.d0
              if (rr.le.bu2b) dcobdrr=(-cob/bu2b)*au2b/(1.d0-rr/bu2b)**2
              dpolybdrr=a1b+2.d0*a2b*rhob+3.d0*a3b*rhob**2
              deebdrr=-a1b*eeb
              dtmpbdrr=0.5d0*deb*(eeb*polyb*dcobdrr+eeb*dpolybdrr*cob
     &                +deebdrr*polyb*cob)

            do k = 1, natom
              dvdx(k) = dvdx(k) -
     &            tmpb*(tt2*gxgcorr(k,j)+tt1*gxgcorr(k,i))/gzero
              dvdy(k) = dvdy(k) -
     &            tmpb*(tt2*gygcorr(k,j)+tt1*gygcorr(k,i))/gzero
              dvdz(k) = dvdz(k) -
     &            tmpb*(tt2*gzgcorr(k,j)+tt1*gzgcorr(k,i))/gzero
            enddo
              dvdx(i)=dvdx(i)+tmpb*(tt2+tt1)*gdum/gzero*(x(i)-x(j))
              dvdx(j)=dvdx(j)-tmpb*(tt2+tt1)*gdum/gzero*(x(i)-x(j))
              dvdy(i)=dvdy(i)+tmpb*(tt2+tt1)*gdum/gzero*(y(i)-y(j))
              dvdy(j)=dvdy(j)-tmpb*(tt2+tt1)*gdum/gzero*(y(i)-y(j))
              dvdz(i)=dvdz(i)+tmpb*(tt2+tt1)*gdum/gzero*(z(i)-z(j))
              dvdz(j)=dvdz(j)-tmpb*(tt2+tt1)*gdum/gzero*(z(i)-z(j))

            fcn = du2 * (es1 * ec1 - 1d0)
            dvdx(i) = dvdx(i) - dtmpbdrr*(x(i)-x(j))/rr*fcn
            dvdx(j) = dvdx(j) - dtmpbdrr*(x(j)-x(i))/rr*fcn
            dvdy(i) = dvdy(i) - dtmpbdrr*(y(i)-y(j))/rr*fcn
            dvdy(j) = dvdy(j) - dtmpbdrr*(y(j)-y(i))/rr*fcn
            dvdz(i) = dvdz(i) - dtmpbdrr*(z(i)-z(j))/rr*fcn
            dvdz(j) = dvdz(j) - dtmpbdrr*(z(j)-z(i))/rr*fcn

            v = v - tmpb * fcn

31          continue
30        continue

C pairwise term
        do 99989 i = 1, natom
          do 99990 j = i + 1, natom
            dx = x(i) - x(j)
            dy = y(i) - y(j)
            dz = z(i) - z(j)
            d1_w = dx * dx + dy * dy + dz * dz
            d2_v = sqrt(d1_w)

            d1_p = 1.0d0 / (2.0d0 *  d2_v)
            g_rr = d1_p * 2.d0
            rr = d2_v
            if (rr .gt. bu2) goto 99990
            g_rho = g_rr
            rho = rr - re
            d4_v = rho * rho
            d2_p = 2.0d0 * rho
            d7_v = rho ** ( 3 - 2)
            d7_v =  d7_v * rho
            d1_p =  3 *  d7_v
            d7_v =  d7_v * rho
            d5_b = a3 * d1_p + a2 * d2_p + a1
            g_poly = d5_b * g_rho
            poly = 1.d0 + a1 * rho + a2 * d4_v + a3 * d7_v
            d1_w = (-a1) * rho
            d2_v = dexp(d1_w)
            d1_p =  d2_v
            g_ee = d1_p * (-a1) * g_rho
            ee = d2_v
            g_co = 0.0d0
            co = 0.d0
            if (rr .le. bu2) then
              d3_v = 1.d0 - rr / bu2
              d4_v = 1.d0 / d3_v
              d6_b = (-((-au2) * ((-d4_v) / d3_v))) * (1.0d0 / bu2)
              d1_w = au2 * (1.d0 - d4_v)
              d2_v = dexp(d1_w)
              d1_p =  d2_v
              g_co = d1_p * d6_b * g_rr
              co = d2_v
            endif
            d2_v = de * ee
            d4_v = d2_v * poly
            d5_b = co * d2_v
            d6_b = co * poly * de
            g_tmp = d4_v * g_co + d5_b * g_poly + d6_b * g_ee
            tmp = d4_v * co
            dum =  - g_tmp
            dvdx(i) = dvdx(i) + dum * dx
            dvdx(j) = dvdx(j) - dum * dx
            dvdy(i) = dvdy(i) + dum * dy
            dvdy(j) = dvdy(j) - dum * dy
            dvdz(i) = dvdz(i) + dum * dz
            dvdz(j) = dvdz(j) - dum * dz
            v = v - tmp
99990     continue
99989   continue
        return
      end
