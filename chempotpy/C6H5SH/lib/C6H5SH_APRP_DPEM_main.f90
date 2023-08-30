program C6H5SH_APRP_DPEM_main
  implicit none
  integer :: nstates, natoms, igrad
  double precision, allocatable :: geom(:,:), u(:,:), ug(:,:,:,:)
  character(len=1024) :: buffer,path
  integer :: id, i, j, k, l
  
  call get_command_argument(1, buffer)
  read(buffer,*) natoms
  call get_command_argument(2, buffer)
  read(buffer,*) nstates
  call get_command_argument(3, buffer)
  read(buffer,*) igrad
  allocate(geom(natoms, 3))
  id=4
  do i=1,natoms
  do j=1,3
    call get_command_argument(id,buffer)
    read(buffer,*) geom(i,j)
    id=id+1
  enddo
  enddo
  allocate(u(nstates,nstates))
  do i=1,nstates
  do j=1,nstates
    call get_command_argument(id,buffer)
    read(buffer,*) u(i,j)
    id=id+1
  enddo
  enddo
  allocate(ug(nstates, nstates, natoms, 3))
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    call get_command_argument(id,buffer)
    read(buffer,*) ug(i,j,k,l)
    id=id+1
  enddo
  enddo
  enddo
  enddo
  call get_command_argument(id, path)
  path=trim(path)
  call dpem(geom,igrad,u,ug)
  open(unit=10, file="shared_data.bin", access="stream", status="replace")
  do i=1,nstates
  do j=1,nstates
    write(10) u(i,j)
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    write(10) ug(i,j,k,l)
  enddo
  enddo
  enddo
  enddo
  close(10)
  deallocate(geom)
  deallocate(u)
  deallocate(ug)
end program
