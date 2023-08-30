subroutine chempotpy(system,surface,path,geom,nstates,natoms,igrad,p,g,d,u,ug)
  implicit none
  character(len=1024), intent(in) :: system, surface, path
  integer, intent(in) :: nstates, natoms, igrad
  double precision, intent(in) :: geom(natoms,3)
  double precision, intent(inout) :: p(nstates)
  double precision, intent(inout) :: g(nstates,natoms,3)
  double precision, intent(inout) :: d(nstates,nstates,natoms,3)
  double precision, intent(inout) :: u(nstates,nstates)
  double precision, intent(inout) :: ug(nstates,nstates,natoms,3)
  character(len=100), dimension(1000) :: List_of_surface_need_path
  integer :: pos
  
  pos = INDEX(surface, "_DPEM")
  if (pos > 0) then
    call dpem_interface(system,surface,path,geom,nstates,natoms,igrad,u,ug)
  else
    call adiabats_interface(system,surface,path,geom,nstates,natoms,igrad,p,g,d)
  endif
end subroutine chempotpy
  
subroutine adiabats_interface(system,surface,path,geom,nstates,natoms,igrad,p,g,d)
  character(len=1024), intent(in) :: system, surface, path
  integer, intent(in) :: nstates, natoms, igrad
  double precision, intent(in) :: geom(natoms,3)
  double precision, intent(inout) :: p(nstates)
  double precision, intent(inout) :: g(nstates,natoms,3)
  double precision, intent(inout) :: d(nstates,nstates,natoms,3)
  character(len=65536) :: serialized_data, commend
  character(len=1024) :: main_name, buffer
  integer :: i, j, k, l
  
  serialized_data=" "
  write(buffer,*) natoms
  serialized_data=trim(serialized_data)//" "// trim(buffer)
  write(buffer,*) nstates
  serialized_data=trim(serialized_data)//" "//trim(buffer)
  write(buffer,*) igrad
  serialized_data=trim(serialized_data)//" "//trim(buffer)
  do i=1,natoms
  do j=1,3
    write(buffer,*) geom(i,j)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  do i=1,nstates
    write(buffer,*) p(i)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  do i=1,nstates
  do j=1,natoms
  do k=1,3
    write(buffer,*) g(i,j,k)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    write(buffer,*) d(i,j,k,l)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  enddo
  enddo
  serialized_data=trim(serialized_data)//" "//trim(path)
  serialized_data = trim(serialized_data)
  main_name=trim(surface)//"_main"
  commend = trim(trim(path)//"/"//trim(system)//"/"//"lib/"//trim(main_name)//" "//serialized_data)
  call EXECUTE_COMMAND_LINE(commend)
  open(unit=10, file="shared_data.bin", access="stream", form="unformatted", status="old")
  do i=1,nstates
    read(10) p(i)
  enddo
  do i=1,nstates
  do j=1,natoms
  do k=1,3
    read(10) g(i,j,k)
  enddo
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    read(10) d(i,j,k,l)
  enddo
  enddo
  enddo
  enddo
  close(10)
end subroutine adiabats_interface
  
subroutine dpem_interface(system,surface,path,geom,nstates,natoms,igrad,u,ug)
  character(len=1024), intent(in) :: system, surface, path
  integer, intent(in) :: nstates, natoms, igrad
  double precision, intent(in) :: geom(natoms,3)
  double precision, intent(inout) :: u(nstates,nstates)
  double precision, intent(inout) :: ug(nstates,nstates,natoms,3)
  character(len=65536) :: serialized_data, commend
  character(len=1024) :: main_name, buffer
  integer :: i, j, k, l
  
  serialized_data=" "
  write(buffer,*) natoms
  serialized_data=trim(serialized_data)//" "// trim(buffer)
  write(buffer,*) nstates
  serialized_data=trim(serialized_data)//" "//trim(buffer)
  write(buffer,*) igrad
  serialized_data=trim(serialized_data)//" "//trim(buffer)
  do i=1,natoms
  do j=1,3
    write(buffer,*) geom(i,j)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
    write(buffer,*) u(i,j)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    write(buffer,*) ug(i,j,k,l)
    serialized_data=trim(serialized_data)//" "//trim(buffer)
  enddo
  enddo
  enddo
  enddo
  serialized_data=trim(serialized_data)//" "//trim(path)
  serialized_data = trim(serialized_data)
  main_name=trim(surface)//"_main"
  commend = trim(trim(path)//"/"//trim(system)//"/"//"lib/"//trim(main_name)//" "//serialized_data)
  call EXECUTE_COMMAND_LINE(commend)
  open(unit=10, file="shared_data.bin", access="stream", form="unformatted", status="old")
  do i=1,nstates
  do j=1,nstates
    read(10) u(i,j)
  enddo
  enddo
  do i=1,nstates
  do j=1,nstates
  do k=1,natoms
  do l=1,3
    read(10) ug(i,j,k,l)
  enddo
  enddo
  enddo
  enddo
  close(10)
end subroutine dpem_interface
