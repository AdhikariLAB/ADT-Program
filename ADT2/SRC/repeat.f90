
program repeatation

  implicit none
  integer :: nstate,ngrid1,ngrid2
  real(kind=8), allocatable :: energy(:)
  open(1,file='input.dat',status='old')
  open(2,file='xx.res',status='old')

  read(1,*)nstate,ngrid1,ngrid2
  allocate energy(nstate,ngrid1,ngrid2),  

  read(2,*) 
end program repeatation 
