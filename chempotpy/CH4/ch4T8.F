C   System:          CH4
C   Functional form:
C   Common name:     Schwenke-Partridge T8
C   Number of derivatives: 0
C   Number of bodies: 5
C   Number of electronic surfaces: 1
C   Interface: Section-2
C   Data file:	ch4T8datafile.tar
C
C   References:      D. W. Schwenke and H. Partridge, Spectrochim. Acta A 57, (2001) 887
C
C   Notes:    
C      The data file (ch4T8datafile.tar) contains two files:
C      ch4pes is a data file read once by the PES routine
C      getpot_ch4t8.f90 is a calling routine that accepts cartesian coordinates, converts them to Radau coordinates 
C      (the required input for the T8 PES), calls the surface, zeroes it so that the ZOE is at the global minimum, 
C      and returns the energy.  
C      There are no derivatives.  One should use nuclear rather than atomic masses for the C and H atoms.
C

      subroutine vibpot(rij,v,n,path)
      implicit real*8 (a-h,o-z)
      character(len=1024), intent(in) :: path
      character(len=1024) :: file_path1
      character*80 title
c
c     ch4 potential based on Radau coordinates
c
c     rij(i,1) - ch1
c     rij(i,2) - ch2
c     rij(i,3) - ch3
c     rij(i,4) - ch4
c     rij(i,5) - h1ch2
c     rij(i,6) - h1ch3
c     rij(i,7) - h1ch4
c     rij(i,8) - h2ch3
c     rij(i,9) - h2ch4
c     rij(i,10) - h3ch4
c                  ALL UNITS ARE IN au,i.e., bohr, radian and hartree. YU
      parameter (id=2000)
      dimension rij(n,10),v(n,1),cs(6),cb(6),                           7d10s01
     $          fmat(9,9),nx(8),ixp(id,8,8),ixw(id,8,8),ccc(id,8)       5d30s00
      dimension ixa(id,5),ixr(id,4),cands(id,id),rmat(id),              7d27s00
     $          ixam(id,5)                                              7d27s00
      save
      data ifirst /0/

      if(ifirst.eq.0)then
      tocm=1d0/4.556335d-6                                              5d23s91
      bohr=0.529177249d0
      bohri=1d0/bohr

       file_path1 = trim(path)//"/CH4/ch4pes"

       open(unit=1,file=file_path1,form='formatted',status='old')
      read(1,1492)title                                                 8d1s00
 1492 format(a80)                                                       8d1s00
      read(1,*)re0,ex0,beta0
      read(1,*)r01,r02,r03,r04
      read(1,*)a02,a03,a04,a06,a08                                      8d3s00
      re0=re0*bohri                                                     8d1s00
      r01=r01*6d0
      r02=r02*6d0
      r03=r03*6d0
      r04=r04*6d0
      a02=a02*4d0
      a03=a03*4d0
      a04=a04*4d0
c      call flush(6)
       call setupvp(nx,ixp,ixw,ccc,id,reau,ex,thrd,exa,srh,sr12,cs,cb,
     $              beta,const,ixa,nexpa,ixr,cands,nexpr,r01,r02,       7d27s00
     $              r03,r04,a02,a03,a04,a06,a08)                        10d6s00
       do 510 j=1,5                                                     7d27s00
        do 1676 i=1,nexpa                                               7d27s00
         ixam(i,j)=ixa(i,j)-1
 1676   continue                                                        7d27s00
  510  continue                                                         7d27s00
       ifirst=1
      end if

      do 100 ig=1,n
       t1=cos(rij(ig,5))+thrd                                           7d24s00
       t2=cos(rij(ig,6))+thrd                                           7d24s00
       t3=cos(rij(ig,7))+thrd                                           7d24s00
       t4=cos(rij(ig,8))+thrd                                           7d24s00
       t5=cos(rij(ig,9))+thrd                                           7d24s00
       t6=cos(rij(ig,10))+thrd                                          7d24s00
       rr1=rij(ig,1)-re0                                                7d24s00
       rr2=rij(ig,2)-re0                                                7d24s00
       rr3=rij(ig,3)-re0                                                7d24s00
       rr4=rij(ig,4)-re0                                                7d24s00
       d1=exp(-beta0*(rr1**2+rr2**2))
       d2=exp(-beta0*(rr1**2+rr3**2))
       d3=exp(-beta0*(rr1**2+rr4**2))
       d4=exp(-beta0*(rr2**2+rr3**2))
       d5=exp(-beta0*(rr2**2+rr4**2))
       d6=exp(-beta0*(rr3**2+rr4**2))
       y1=1d0-exp(-ex0*rr1)
       y2=1d0-exp(-ex0*rr2)
       y3=1d0-exp(-ex0*rr3)
       y4=1d0-exp(-ex0*rr4)
       t12=t1*t1
       t22=t2*t2
       t32=t3*t3
       t42=t4*t4
       t52=t5*t5
       t62=t6*t6
       vs=y1*(r01+y1*(r02+y1*(r03+y1*r04)))
     $            +y2*(r01+y2*(r02+y2*(r03+y2*r04)))
     $            +y3*(r01+y3*(r02+y3*(r03+y3*r04)))
     $            +y4*(r01+y4*(r02+y4*(r03+y4*r04)))
     $            +d1*t12*(a02+t1*(a03+t1*(a04+t12*(a06+t12*a08))))     8d3s00
     $            +d2*t22*(a02+t2*(a03+t2*(a04+t22*(a06+t22*a08))))     8d3s00
     $            +d3*t32*(a02+t3*(a03+t3*(a04+t32*(a06+t32*a08))))     8d3s00
     $            +d4*t42*(a02+t4*(a03+t4*(a04+t42*(a06+t42*a08))))     8d3s00
     $            +d5*t52*(a02+t5*(a03+t5*(a04+t52*(a06+t52*a08))))     8d3s00
     $            +d6*t62*(a02+t6*(a03+t6*(a04+t62*(a06+t62*a08))))     8d3s00
              d1=rij(ig,1)-reau                                                   5d8s00
              d2=rij(ig,2)-reau                                                   5d8s00
              d3=rij(ig,3)-reau                                                   5d8s00
              d4=rij(ig,4)-reau                                                   5d8s00
       if(ex.gt.0d0)then                                                10d6s00
        d1=1d0-exp(-ex*d1)                                              10d11s00
        d2=1d0-exp(-ex*d2)                                              10d11s00
        d3=1d0-exp(-ex*d3)                                              10d11s00
        d4=1d0-exp(-ex*d4)                                              10d11s00
       end if                                                           10d6s00
              damp=exp(-beta*(d1*d1+d2*d2+d3*d3+d4*d4))                 7d24s00
              s1=0.5d0*(d1+d2+d3+d4)
              s2a=sr12*(2d0*t1-t2-t3-t4-t5+2d0*t6)
              s2b=0.5d0*(t2-t3-t4+t5)
              s3x=0.5d0*(d1-d2+d3-d4)
              s3y=0.5d0*(d1-d2-d3+d4)
              s3z=0.5d0*(d1+d2-d3-d4)
              s4x=srh*(t5-t2)
              s4y=srh*(t4-t3)
              s4z=srh*(t6-t1)
              fmat(1,1)=1d0
              fmat(1,2)=1d0
              fmat(1,3)=1d0
              fmat(1,4)=1d0
              fmat(1,5)=1d0
              fmat(1,6)=1d0
              fmat(1,7)=1d0
              fmat(1,8)=1d0
              fmat(1,9)=1d0
              fmat(2,1)=s1
              fmat(2,2)=s2a
              fmat(2,3)=s2b
              fmat(2,4)=s3x
              fmat(2,5)=s3y
              fmat(2,6)=s3z
              fmat(2,7)=s4x
              fmat(2,8)=s4y
              fmat(2,9)=s4z
              do 4353 i=3,9
               fmat(i,1)=fmat(i-1,1)*fmat(2,1)
               fmat(i,2)=fmat(i-1,2)*fmat(2,2)
               fmat(i,3)=fmat(i-1,3)*fmat(2,3)
               fmat(i,4)=fmat(i-1,4)*fmat(2,4)
               fmat(i,5)=fmat(i-1,5)*fmat(2,5)
               fmat(i,6)=fmat(i-1,6)*fmat(2,6)
               fmat(i,7)=fmat(i-1,7)*fmat(2,7)
               fmat(i,8)=fmat(i-1,8)*fmat(2,8)
               fmat(i,9)=fmat(i-1,9)*fmat(2,9)
 4353         continue
              do 67 i=1,nx(1)                                            5d30s00
               vs=vs+ccc(i,1)*fmat(ixp(i,1,1),ixw(i,1,1))*damp               5d30s00
   67         continue                                                  5d30s00
              do 68 i=1,nx(2)                                            5d30s00
               vs=vs+ccc(i,2)*fmat(ixp(i,1,2),ixw(i,1,2))*damp               5d30s00
     $                       *fmat(ixp(i,2,2),ixw(i,2,2))               5d30s00
   68         continue                                                  5d30s00
              do 69 i=1,nx(3)                                            5d30s00
               vs=vs+ccc(i,3)*fmat(ixp(i,1,3),ixw(i,1,3))*damp               5d30s00
     $                       *fmat(ixp(i,2,3),ixw(i,2,3))               5d30s00
     $                       *fmat(ixp(i,3,3),ixw(i,3,3))               5d30s00
   69         continue                                                  5d30s00
              do 70 i=1,nx(4)                                            5d30s00
               vs=vs+ccc(i,4)*fmat(ixp(i,1,4),ixw(i,1,4))*damp               5d30s00
     $                       *fmat(ixp(i,2,4),ixw(i,2,4))               5d30s00
     $                       *fmat(ixp(i,3,4),ixw(i,3,4))               5d30s00
     $                       *fmat(ixp(i,4,4),ixw(i,4,4))               5d30s00
   70         continue                                                  5d30s00
              do 71 i=1,nx(5)                                            5d30s00
               vs=vs+ccc(i,5)*fmat(ixp(i,1,5),ixw(i,1,5))*damp               5d30s00
     $                       *fmat(ixp(i,2,5),ixw(i,2,5))               5d30s00
     $                       *fmat(ixp(i,3,5),ixw(i,3,5))               5d30s00
     $                       *fmat(ixp(i,4,5),ixw(i,4,5))               5d30s00
     $                       *fmat(ixp(i,5,5),ixw(i,5,5))               5d30s00
   71         continue                                                  5d30s00
              do 72 i=1,nx(6)                                            5d30s00
               vs=vs+ccc(i,6)*fmat(ixp(i,1,6),ixw(i,1,6))*damp               5d30s00
     $                       *fmat(ixp(i,2,6),ixw(i,2,6))               5d30s00
     $                       *fmat(ixp(i,3,6),ixw(i,3,6))               5d30s00
     $                       *fmat(ixp(i,4,6),ixw(i,4,6))               5d30s00
     $                       *fmat(ixp(i,5,6),ixw(i,5,6))               5d30s00
     $                       *fmat(ixp(i,6,6),ixw(i,6,6))               5d30s00
   72         continue                                                  5d30s00
              do 73 i=1,nx(7)                                            5d30s00
               vs=vs+ccc(i,7)*fmat(ixp(i,1,7),ixw(i,1,7))*damp               5d30s00
     $                       *fmat(ixp(i,2,7),ixw(i,2,7))               5d30s00
     $                       *fmat(ixp(i,3,7),ixw(i,3,7))               5d30s00
     $                       *fmat(ixp(i,4,7),ixw(i,4,7))               5d30s00
     $                       *fmat(ixp(i,5,7),ixw(i,5,7))               5d30s00
     $                       *fmat(ixp(i,6,7),ixw(i,6,7))               5d30s00
     $                       *fmat(ixp(i,7,7),ixw(i,7,7))               5d30s00
   73         continue                                                  5d30s00
              do 74 i=1,nx(8)                                            5d30s00
               vs=vs+ccc(i,8)*fmat(ixp(i,1,8),ixw(i,1,8))*damp               5d30s00
     $                       *fmat(ixp(i,2,8),ixw(i,2,8))               5d30s00
     $                       *fmat(ixp(i,3,8),ixw(i,3,8))               5d30s00
     $                       *fmat(ixp(i,4,8),ixw(i,4,8))               5d30s00
     $                       *fmat(ixp(i,5,8),ixw(i,5,8))               5d30s00
     $                       *fmat(ixp(i,6,8),ixw(i,6,8))               5d30s00
     $                       *fmat(ixp(i,7,8),ixw(i,7,8))               5d30s00
     $                       *fmat(ixp(i,8,8),ixw(i,8,8))               5d30s00
   74         continue                                                  5d30s00
              vs=vs+const*damp                                          8d3s00
              v(ig,1)=vs                                                7d25s00
  100 continue
      return
      end

      subroutine setupvp(nx,ixp,ixw,ccc,idx,reau,ex,thrd,exa,srh,sr12,
     $           cs,cb,beta,const,ixa,nexpa,ixr,cands,nexpr,r01,        7d27s00
     $           r02,r03,r04,a02,a03,a04,a06,a08)                       10d6s00
      implicit real*8 (a-h,o-z)
      parameter (id=6400,idu=20,idf=100)
      dimension nx(8),ixp(idx,8,8),ixw(idx,8,8),ccc(idx,8),             7d25s00
     $          ixa(idx,5),ixr(idx,4),ixar(id,2),cands(idx,1)           7d27s00
      dimension cs(6),cb(6),
     $          ix(id,9),cc(id),iperm(4,24),scart(3,4),
     $          iaperm(6,24),ca(6),stry(8,24),rx(4),try(24),ihit(24),
     $          ixs(id,9),ccs(id),jnu(idu),inux(idf,9),cnu(idu,idf),    6d10s00
     $          cnua(idu,idu),cnub(idu,idf),cnuc(idu),ipvt(idu),        6d10s00
     $          cnuu(idu,idf),                                          6d10s00
     $          icnuc(idu),cnuaa(idu,idu),cnud(idu,idf),icnud(idu)      6d10s00
      data iperm/1,2,3,4, 1,2,4,3, 1,3,2,4, 1,3,4,2, 1,4,2,3, 1,4,3,2,  5d8s00
     $           2,1,3,4, 2,1,4,3, 2,3,1,4, 2,3,4,1, 2,4,1,3, 2,4,3,1,  5d8s00
     $           3,1,2,4, 3,1,4,2, 3,2,1,4, 3,2,4,1, 3,4,1,2, 3,4,2,1,  5d8s00
     $           4,1,2,3, 4,1,3,2, 4,2,1,3, 4,2,3,1, 4,3,1,2, 4,3,2,1/  5d8s00
      data iaperm/1,2,3,4,5,6, 1,3,2,5,4,6, 2,1,3,4,6,5, 2,3,1,6,4,5,
     $            3,1,2,5,6,4, 3,2,1,6,5,4, 1,4,5,2,3,6, 1,5,4,3,2,6,
     $            4,1,5,2,6,3, 4,5,1,6,2,3, 5,1,4,3,6,2, 5,4,1,6,3,2,
     $            2,4,6,1,3,5, 2,6,4,3,1,5, 4,2,6,1,5,3, 4,6,2,5,1,3,
     $            6,2,4,3,5,1, 6,4,2,5,3,1, 3,5,6,1,2,4, 3,6,5,2,1,4,
     $            5,3,6,1,4,2, 5,6,3,4,1,2, 6,3,5,2,4,1, 6,5,3,4,2,1/
       read(1,*)re,ex,exa,beta,iradau,lm
       tocm=1d0/4.556335d-6                                              5d23s91
       sr12=sqrt(1d0/12d0)
       srh=sqrt(0.5d0)
       nexp=0
       ix1=0
       ix2=0
c       call flush(6)
    2  continue
        read(1,*,end=4)i1,i2,i3,i4,i5,i6,i7,i8,i9,coef
        in=i3+i4+i5+i6+i7+i8+i9
        j1=in+i2
        j2=in+i1
        if(j1.eq.0.and.lm.ne.0)then                                     6d23s00
         cs(i1)=coef*6d0
         ix1=max(ix1,i1)
         coef=0d0                                                       6d13s00
        else if(j2.eq.0.and.lm.ne.0)then                                6d23s00
         cb(i2)=coef*4d0
         ix2=max(ix2,i2)
         coef=0d0                                                       6d13s00
        end if                                                          6d13s00
         nexp=nexp+1
         if(nexp.gt.id)stop 'id'
         ix(nexp,1)=i1
         ix(nexp,2)=i2
         ix(nexp,3)=i3
         ix(nexp,4)=i4
         ix(nexp,5)=i5
         ix(nexp,6)=i6
         ix(nexp,7)=i7
         ix(nexp,8)=i8
         ix(nexp,9)=i9
         cc(nexp)=coef
c$$$         cc(nexp)=cc(nexp)*4.556335d-6                                  7d24s00
    3    format(i5,9i2,1pe22.14)
c$$$        end if
        go to 2
    4  continue
       energy=1d0/(1d-18*6.0221367d23*1d-3*3.808798d-4)
       reau=re*1.889725989d0
       do 5 i=1,nexp
c$$$        cc(i)=cc(i)*energy
        do 6 j=1,9
         ix(i,j)=ix(i,j)+1
    6   continue
    5  continue
       thrd=1d0/3d0
       close(unit=1)
       nx(1)=0                                                             5d30s00
       nx(2)=0                                                             5d30s00
       nx(3)=0                                                             5d30s00
       nx(4)=0                                                             5d30s00
       nx(5)=0
       nx(6)=0
       nx(7)=0
       nx(8)=0
       nxzero=0                                                         7d24s00
       do 600 i=1,nexp                                                  5d30s00
        nnz=0                                                           5d30s00
        do 601 j=1,9                                                    5d30s00
         if(ix(i,j).gt.1)nnz=nnz+1                                      5d30s00
  601   continue                                                        5d30s00
        if(nnz.eq.0)then                                                7d24s00
         nxzero=nxzero+1                                                7d24s00
         const=cc(i)                                                    7d24s00
        else                                                            7d24s00
         if(nnz.eq.7)then
 1492     format(9i2)
         end if
        nx(nnz)=nx(nnz)+1                                                 5d30s00
        if(nx(nnz).gt.idx)stop 'idx'                                    7d24s00
        ccc(nx(nnz),nnz)=cc(i)                                           5d30s00
        ii=1                                                            5d30s00
        do 602 j=1,9                                                    5d30s00
         if(ix(i,j).gt.1)then                                           5d30s00
          ixp(nx(nnz),ii,nnz)=ix(i,j)                                    5d30s00
          ixw(nx(nnz),ii,nnz)=j                                          5d30s00
          ii=ii+1                                                       5d30s00
         end if                                                         5d30s00
  602   continue                                                        5d30s00
        end if                                                          7d24s00
  600  continue                                                         5d30s00
       return
       end
