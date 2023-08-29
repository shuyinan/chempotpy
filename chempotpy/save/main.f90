program main
  implicit none
  character(len=1024) :: system, surface, path
  integer, parameter :: natoms=3, nstates=6
  integer :: igrad
  double precision :: geom(natoms,3)
  double precision :: p(nstates)
  double precision :: g(nstates,natoms,3)
  double precision :: d(nstates,nstates,natoms,3)
  double precision :: u(nstates,nstates)
  double precision :: ug(nstates,nstates,natoms,3)

  igrad=2
  geom(1,1)=-0.506
  geom(1,2)=0.0
  geom(1,3)=-0.094
  geom(2,1)=-1.789
  geom(2,2)=0.0
  geom(2,3)=-0.902
  geom(3,1)=-1.467
  geom(3,2)=0.0
  geom(3,3)=0.0

  p=0.d0
  g=0.d0
  d=0.d0
  u=0.d0
  ug=0.d0

  path="/home/truhlard/shuxx055/test/chempotpy"
  system="O3"
  surface="O3_6_5Ap_2023"

  call chempotpy(system,surface,path,geom,nstates,natoms,igrad,p,g,d,u,ug)

  write(*,*) "adiabats"
  write(*,*) p(1), p(2), p(3)
  write(*,*) p(4), p(5), p(6)
  write(*,*) "gradients of state 1"
  write(*,*) g(1,1,1), g(1,1,2), g(1,1,3)
  write(*,*) g(1,2,1), g(1,2,2), g(1,2,3)
  write(*,*) g(1,3,1), g(1,3,2), g(1,3,3)
  write(*,*) "NAC of state 1 and state 2"
  write(*,*) d(1,2,1,1), d(1,2,1,2), d(1,2,1,3)
  write(*,*) d(1,2,2,1), d(1,2,2,2), d(1,2,2,3)
  write(*,*) d(1,2,3,1), d(1,2,3,2), d(1,2,3,3)

contains

subroutine chempotpy(system,surface,path,geom,nstates,natoms,igrad,p,g,d,u,ug)
  use, intrinsic :: iso_c_binding
  implicit none
  character(len=1024), intent(in) :: system, surface, path
  integer, intent(in) :: nstates, natoms, igrad
  double precision, intent(in) :: geom(natoms,3)
  double precision, intent(inout) :: p(nstates)
  double precision, intent(inout) :: g(nstates,natoms,3)
  double precision, intent(inout) :: d(nstates,nstates,natoms,3)
  double precision, intent(inout) :: u(nstates,nstates)
  double precision, intent(inout) :: ug(nstates,nstates,natoms,3)

  interface
    subroutine call_pes(so_name, geom, igrad, p, g, d) bind(C, name="call_pes")
      use iso_c_binding, only: c_char, c_int, c_double
      character(kind=c_char), intent(in) :: so_name(*)
      integer(kind=c_int), intent(in) :: igrad
      real(kind=c_double), dimension(:,:), intent(in) :: geom
      real(kind=c_double), dimension(:), intent(inout) :: p
      real(kind=c_double), dimension(:,:,:), intent(inout) :: g
      real(kind=c_double), dimension(:,:,:,:), intent(inout) :: d
    end subroutine call_pes

    subroutine call_dpem(so_name, geom, igrad, u, ug) bind(C, name="call_dpem")
      use iso_c_binding, only: c_char, c_int, c_double
      character(kind=c_char), intent(in) :: so_name(*)
      integer(kind=c_int), intent(in) :: igrad
      real(kind=c_double), dimension(:,:), intent(in) :: geom
      real(kind=c_double), dimension(:,:), intent(inout) :: u
      real(kind=c_double), dimension(:,:,:,:), intent(inout) :: ug
    end subroutine call_dpem
  end interface

  character(len=1024) :: so_path
  so_path = trim(path)// '/' //trim(system) // '/' // trim(surface) // '.so'

  write(*,*) "geometry"
  write(*,*) geom(1,1), geom(1,2), geom(1,3)
  write(*,*) geom(2,1), geom(2,2), geom(2,3)
  write(*,*) geom(3,1), geom(3,2), geom(3,3)  

  call call_pes(trim(so_path), geom, igrad, p, g, d)

end subroutine

end program
