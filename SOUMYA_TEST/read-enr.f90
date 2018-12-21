
!---------------------------------------------------------------------------
! This is fortran code for reading the energy values
!---------------------------------------------------------------------------

!use hdf5

program input

implicit none
integer :: i,j,ii,jj,nstate,ncoupling,fgrid,sgrid
real(kind=8), allocatable :: energy(:),nactddr(:),fgr(:),sgr(:)
character(len=100) :: molecule,nacttype,filename1,filename2,filename3,filename4
read(*,*)nstate,fgrid,sgrid,molecule,nacttype

filename1=trim(molecule)//'_energy.res'
filename2=trim(molecule)//'_nact.res'
filename3=trim(molecule)//'_energy.dat'
filename4=trim(molecule)//'_nact.dat'

open(11,file=filename1,status='old')
open(12,file=filename2,status='old')
open(13,file=filename3,status='unknown')
open(14,file=filename4,status='unknown')

ncoupling = nstate*(nstate-1)/2
 
!character(len=9), parameter :: filename = "energy.h5"  ! File name
!character(len=8), parameter :: dsetname = "IntArray" 

!integer(hid_t) :: file_id       ! File identifier
!integer(hid_t) :: dset_id       ! Dataset identifier
!integer(hid_t) :: dataspace     ! Dataspace identifier
!integer(hid_t) :: memspace
!
!integer(hsize_t), dimension(3) :: dimsm = (/fgrid,sgrid,3/) ! Dataset dimensions in memory
!integer(hsize_t), dimension(2) :: dimsf = (/fgrid,sgrid/)   ! Dataset dimensions 

allocate(energy(nstate))
allocate(fgr(fgrid))
allocate(sgr(sgrid))

do ii = 1,fgrid
  do jj = 1,sgrid
    do i = 1,4
      read(*,*)
    enddo
    do j = 1,nstate
      read(11,*)energy(j)
    enddo
    write(13,*)fgr(ii),sgr(jj),energy
  enddo
  write(13,*)
enddo

select case(nacttype)

  case('ddr')
   do ii = 1,fgrid
     do jj = 1,sgrid
       do i = 1,4
         read(*,*)
       enddo
       read(12,*)aa,bb,(nactddr(j),j=1,2*ncoupling
       write(14,*)fgr(ii),sgr(jj),nactddr
     enddo
     write(14,*)
   enddo

        

end select

deallocate(energy)
deallocate(fgr)
deallocate(sgr)

stop

end program input
 







