
!---------------------------------------------------------------------------
! This is fortran code for reading the energy and nact values
!---------------------------------------------------------------------------

!use hdf5

module parameters

 implicit none
 real*8, parameter :: pi = 4.0d0*datan(1.0d0)
 real*8, parameter :: hcross = 0.063508d0
 real*8, parameter :: angtobohr = 1.889725989d0

end module parameters

program input

use parameters

implicit none
integer :: i,j,l,ii,jj,kk,nstate,natom,ncoupling,fgrid,sgrid
real(kind=8), allocatable :: energy(:),nactddr(:),fgr(:),sgr(:),coeff1(:),coeff2(:),mass(:)
real(kind=8), allocatable :: taur(:,:,:),taup(:,:,:),grad(:,:,:,:)
real(kind=8) :: omega1,omega2,ph,aa,bb
character(len=100) :: molecule,nacttype,filename1,filename2,filename3,filename4
read(*,*)nstate,fgrid,sgrid,molecule,nacttype

ncoupling = nstate*(nstate-1)/2

allocate(coeff1(natom))
allocate(coeff2(natom))
allocate(mass(natom))
allocate(energy(nstate))
allocate(fgr(fgrid))
allocate(sgr(sgrid))
allocate(taur(fgrid,sgrid,natom))
allocate(taup(fgrid,sgrid,natom))
allocate(grad(fgrid,sgrid,ncoupling,3*natom))

if(nacttype.eq.'analytic') then 
   read(*,*)omega1,omega2
   read(*,*)mass
   read(*,*)coeff1
   read(*,*)coeff2 
endif

filename1=trim(molecule)//'_energy.res'
filename2=trim(molecule)//'_nact.res'
filename3=trim(molecule)//'_energy.dat'
filename4=trim(molecule)//'_nact.dat'

open(11,file=filename1,status='old')
open(12,file=filename2,status='old')
open(13,file=filename3,status='unknown')
open(14,file=filename4,status='unknown')
 
!character(len=9), parameter :: filename = "energy.h5"  ! File name
!character(len=8), parameter :: dsetname = "IntArray" 

!integer(hid_t) :: file_id       ! File identifier
!integer(hid_t) :: dset_id       ! Dataset identifier
!integer(hid_t) :: dataspace     ! Dataspace identifier
!integer(hid_t) :: memspace
!
!integer(hsize_t), dimension(3) :: dimsm = (/fgrid,sgrid,3/) ! Dataset dimensions in memory
!integer(hsize_t), dimension(2) :: dimsf = (/fgrid,sgrid/)   ! Dataset dimensions 

do ii = 1,fgrid
  do jj = 1,sgrid
    do i = 1,4
      read(*,*)
    enddo
    do j = 1,nstate
      read(11,*)fgr(ii),sgr(jj),energy(j)
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
         read(12,*)
       enddo
       read(12,*)aa,bb,(nactddr(j),j=1,2*ncoupling)
       write(14,*)fgr(ii),sgr(jj),nactddr
     enddo
     write(14,*)
   enddo

   case('analytic')
    coeff1=coeff1*dsqrt(hcross/mass/omega1)*angtobohr
    coeff2=coeff2*dsqrt(hcross/mass/omega2)*angtobohr
    taur = 0.0d0
    taup = 0.0d0
    do ii = 1,fgrid
      do jj = 1,sgrid
        ph = sgr(jj)*pi/180.0d0
        do kk = 1,ncoupling
          do i = 1,4
            read(12,*)
          enddo
          do i = 1,natom
            read(12,*)(grad(ii,jj,kk,l),l=3*i-2,3*i)
          enddo
          do i = 1,3*natom
            taur(ii,jj,kk) = taur(ii,jj,kk)+coeff1(i)*grad(ii,jj,kk,i)*dcos(ph)+coeff2(i)*grad(ii,jj,kk,i)*dsin(ph)
            taup(ii,jj,kk) = taup(ii,jj,kk)-coeff1(i)*grad(ii,jj,kk,i)*dsin(ph)+coeff2(i)*grad(ii,jj,kk,i)*dcos(ph)
          enddo 
        enddo
        write(14,*)fgr(ii),sgr(jj),taur,taup
      enddo
      write(14,*)
    enddo

end select

deallocate(energy)
deallocate(fgr)
deallocate(sgr)
deallocate(coeff1)
deallocate(coeff2)
deallocate(mass)
deallocate(taur)
deallocate(taup)
deallocate(grad)

stop

end program input
 







