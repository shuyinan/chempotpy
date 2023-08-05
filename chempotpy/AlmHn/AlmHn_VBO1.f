      subroutine pes(xyz,nAl,nH,natoms,igrad,p,g,d)
      implicit none
      ! number of electronic state
      integer, parameter :: nstates=1
      integer, intent(in) :: nAl, nH, natoms
      double precision, intent(in) :: xyz(10000,3)
      integer, intent(in) :: igrad
      double precision, intent(out) :: p(nstates), g(nstates,natoms,3)
      double precision, intent(out) :: d(nstates,nstates,natoms,3)

      double precision :: cx(natoms), cy(natoms), cz(natoms)
      double precision :: v
      character*2 :: symbol(natoms)
      integer :: istate, iatom

      !initialize 
      v=0.d0
      g=0.d0
      d=0.d0

      do iatom=1,natoms
         cx(iatom)=xyz(iatom,1)/0.529177211
         cy(iatom)=xyz(iatom,2)/0.529177211
         cz(iatom)=xyz(iatom,3)/0.529177211
      enddo

      do iatom=1,nAl
        symbol(iatom)='Al'
      enddo
      do iatom=1,nH
        symbol(iatom+nAl)='H'
      enddo 

      call alhpot(cx, cy, cz, v, natoms, 10000, symbol)

      if (igrad==0) then
        do istate=1,nstates
          p(istate)=v*27.211386
        enddo
      else
        write (*,*) 'Only energy is available'
      endif

      endsubroutine

      subroutine alhpot(x,y,z,v,natom,maxatom,symbol)
c   System:                        Al, H, Al-H
c   Functional form:               Valence-Bond Order
c   Common name:                   VBO.beta1
c   Number of derivatives:         0
c   Number of bodies:              variable
c   Number of electronic surfaces: 1
c   Interface:                     HO-MM-0
c
c   Notes:  Valence-Bond Order potential energy funciton. The functional
c           form is from Ref. 1. 
c
c  References: (1)"valence-Bond Order(VBO): A New Approach to Modeling
c                 Reactive Potential Energy Surfaces for Materials and Nanoparticles"
c                 Meiyu Zhao, Mark A. Iron, Przemys?aw Staszewski, Nathan E. Schultz,
c                Rosendo Valero, and Donald G. Truhlar, manuscript in preparation, Oct.18, 2008.
c
c  Units:
c       energies    - hartrees
c       coordinates - bohrs
c
c  v       --- energy of the system (output)
c  x,y,z   --- one-dimensional arrays of the Cartesian coordinates
c              of the system (input)
c  natom   --- number of atoms in the system (input)
c  maxatom --- maximal number of atoms (input)
c  symbol  --- character*2 array of dimension holds the atomic symbols.

C This version of the code corrects two typos around lines 64 and 65. January 5, 2009 by Meiyu Zhao. 

      implicit real*8(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autoang=0.529177249d0)
      character*2 symbol(maxatom)
      dimension x(maxatom),y(maxatom),z(maxatom),r(maxatom,maxatom)
      dimension id(maxatom)
      double precision b(4,maxatom,maxatom),sigma(4,maxatom,maxatom)
      double precision cn(4,maxatom,maxatom), gam(4,maxatom,maxatom)

c   gam1_h to gam4_h are gammaij in Eq. (1) of Ref. 1
c   c1_h to c4_h are ci,1 to ci,4 in Eqs. (6) and (7)
c   esn_h is the exponent of Vi in Eq. (6)

           gam1_h=       2.7292
           gam2_h=       8.8100
           gam3_h=      16.1462
           gam4_h=      16.1330
           c1_h=         2.0112/autoev
           c2_h=         0.8084/autoev
           esn_h=        0.1536
	   c3_h=         0.00001
	   c4_h=         0.00001
c           c3_h=     23378.9026/autoev**(1/esn_h)
c           c4_h=     22212.1589/autoev**(1/esn_h)

c   The same definitions apply for Aluminum

           gam1_al=      0.3599
           gam2_al=      4.2538
           gam3_al=      0.5684
           gam4_al=      2.9562
           c1_al=        0.2248/autoev
           c2_al=        0.2172/autoev 
           esn_al=       0.7931
	   c3_al=         0.0133
	   c4_al=         0.0068
c           c3_al=        0.8576/autoev**(1/esn_al)
c           c4_al=        0.4396/autoev**(1/esn_al)

c   gam1_hal to gam4_hal and gam3_alh, gam4_alh are gammaij in Eq. (1)
c   c1_hal and c2_hal are ci,1 and ci,2 in Eqs. (9) and (10)
c   esn_h is the exponent of Vi in Eq. (9)
c   sig3_hal to sig4_alh are the sigma in Eq. (10)
   

           gam1_hal=     7.4891  
           gam2_hal=     1.8993
           gam3_hal=    13.5669
           gam4_hal=    13.6902
           gam3_alh=     1.6911
           gam4_alh=     2.8674
           sig3_hal=     0.00001
           sig4_hal=     0.0191
           sig3_alh=     2.9306
           sig4_alh=     0.00008
           c1_hal=       0.3509/autoev
           c2_hal=       0.5174/autoev

c   Determine which atoms in the data file are Aluminum and which Hydrogen

      do i=1,natom
	if((symbol(i).eq.'H ').or.(symbol(i).eq.'h ').or.
     &   (symbol(i).eq.' H').or.(symbol(i).eq.' h'))then
	  id(i)=1
	elseif((symbol(i).eq.'AL').or.(symbol(i).eq.'Al').or.
     &   (symbol(i).eq.'al').or.(symbol(i).eq.'aL'))then
	  id(i)=2
	else
        write(13,*)'Atomic symbol ',symbol(i),' undefined. Stopping.'
	 stop "Undefined atomic symbol."
	endif
      enddo

c    Calculate all internuclear distances

        do 1 i=1,natom-1
          r(i,i)=0.d0
        do 1 j=i+1,natom
          dx=x(i)-x(j)
          dy=y(i)-y(j)
          dz=z(i)-z(j)
          r(i,j)=dsqrt(dx*dx+dy*dy+dz*dz)
1        r(j,i)=r(i,j)
         r(natom,natom)=0.d0

c    Calculate en1 which corresponds to the first term in parenthesis in Eq. (6) (pure Al or pure H) and Eq. (9) (AlH).

      en1=0d0
      en2=0d0
       do 2 i=1,natom
      v1=0d0
      do 3 j=1,natom
        if(i.eq.j)go to 3
	if((id(i).eq.1).and.(id(j).eq.1))then
           del1=5.29/autoang

           gam(1,i,j)=gam1_h
           gam(2,i,j)=gam2_h
           gam(3,i,j)=gam3_h
           gam(4,i,j)=gam4_h
           en=0.5
           c1= c1_h
           c2= c2_h

           r1=0.7414/autoang

c    Calculate Nij according to Eq. (2)

           den1=1d0-(r1/del1)**en
           cn(1,i,j)=dexp(gam(1,i,j)/den1)
           den2=1d0-(r1/del1)**en
           cn(2,i,j)=dexp(gam(2,i,j)/den2)
           den3=1d0-(r1/del1)**en
           cn(3,i,j)=dexp(gam(3,i,j)/den3)
           den4=1d0-(r1/del1)**en
           cn(4,i,j)=dexp(gam(4,i,j)/den4)

	   sigma(3,i,j)=1.0
	   sigma(4,i,j)=1.0

	else if((id(i).eq.2).and.(id(j).eq.2))then
           del1= 6.88/autoang

	   gam(1,i,j)=gam1_al
	   gam(2,i,j)=gam2_al
	   gam(3,i,j)=gam3_al
	   gam(4,i,j)=gam4_al
           en=0.5
           c1=c1_al
           c2=c2_al

	   r1= 2.7306/autoang

           den1=1d0-(r1/del1)**en
           cn(1,i,j)=dexp(gam(1,i,j)/den1)
           den2=1d0-(r1/del1)**en
           cn(2,i,j)=dexp(gam(2,i,j)/den2)
           den3=1d0-(r1/del1)**en
           cn(3,i,j)=dexp(gam(3,i,j)/den3)
           den4=1d0-(r1/del1)**en
           cn(4,i,j)=dexp(gam(4,i,j)/den4)

	   sigma(3,i,j)=1.0
	   sigma(4,i,j)=1.0

	else
           del1=6.88/autoang

	    gam(1,i,j)=gam1_hal
	    gam(2,i,j)=gam2_hal
	   if((id(i).eq.1).and.(id(j).eq.2))then
	    gam(3,i,j)=gam3_hal
	    gam(4,i,j)=gam4_hal
	    sigma(3,i,j)=sig3_hal
	    sigma(4,i,j)=sig4_hal
	   else if((id(i).eq.2).and.(id(j).eq.1)) then
	    gam(3,i,j)=gam3_alh
	    gam(4,i,j)=gam4_alh
	    sigma(3,i,j)=sig3_alh
	    sigma(4,i,j)=sig4_alh
	   end if

	   en=0.5
	   c1=c1_hal
	   c2=c2_hal

           r1= 1.6637/autoang

           den1=1d0-(r1/del1)**en
           cn(1,i,j)=dexp(gam(1,i,j)/den1)
           den2=1d0-(r1/del1)**en
           cn(2,i,j)=dexp(gam(2,i,j)/den2)
           den3=1d0-(r1/del1)**en
           cn(3,i,j)=dexp(gam(3,i,j)/den3)
           den4=1d0-(r1/del1)**en
           cn(4,i,j)=dexp(gam(4,i,j)/den4)

	endif

c   Calculate bij according to Eq. (1)

        rij=r(i,j)
        if(rij.lt.del1)then
	  b(1,i,j)=dexp(-gam(1,i,j)/(1d0-(rij/del1)**en))*cn(1,i,j)
	  b(2,i,j)=dexp(-gam(2,i,j)/(1d0-(rij/del1)**en))*cn(2,i,j)
	  b(3,i,j)=dexp(-gam(3,i,j)/(1d0-(rij/del1)**en))*cn(3,i,j)
	  b(4,i,j)=dexp(-gam(4,i,j)/(1d0-(rij/del1)**en))*cn(4,i,j)
          fa=c1*b(1,i,j)+c2*b(2,i,j)
          v1=v1+fa
	else
	  do k=1,4
	    b(k,i,j)=0.0d0
	  enddo
	endif
   3  continue   
      en1=en1+v1
   2  continue 

c    Calculate en2 which corresponds to the last term in Eq. (6) (pure Al or pure H) and Eq. (9) (AlH)
 
   
      do 4 i=1,natom
        if(id(i).eq.1)then
          c3  = c3_h
          c4  = c4_h
          esn = esn_h
	elseif(id(i).eq.2)then
          c3  = c3_al
          c4  = c4_al
          esn = esn_al
      endif
	Vi3=0.0d0
	Vi4=0.0d0
        do 5 j=1,natom
	   if(i.eq.j)goto 5
	     Vi3=Vi3+c3*sigma(3,i,j)*b(3,i,j)
	     Vi4=Vi4+c4*sigma(4,i,j)*b(4,i,j)
 5      continue
        add=(Vi3+Vi4)**esn
        en2=en2-add
 4    continue

      v=en1+en2
      end      


