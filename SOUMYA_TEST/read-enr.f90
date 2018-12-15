
!---------------------------------------------------------------------------
! This is fortran code for reading the energy values
!---------------------------------------------------------------------------

implicit none
integer :: 
real(kind=8), allocatable :: energy(:),fgr(:),sgr(:)
read(*,*)nstate,fgrid,sgrid

allocate(energy(nstate))
allocate(fgr(fgrid))
allocate(sgr(sgr))

do ii = 1,fgrid
  do jj = 1,sgrid
    do i = 1,4
      read(*,*)
    enddo
    do j = 1,nstate
      read(*,*)enr(j)
    enddo
    write(*,*)fgr(ii),sgr(jj),enr
  enddo
  write(*,*)
enddo



 







