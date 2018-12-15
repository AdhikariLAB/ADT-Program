
program adt

  implicit none
  integer :: mm1,mm2,n,nn,path,i,j,k,ii,jj
  real(kind=8) :: r,p,dr,dp,h,toutr,toutp
  real(kind=8), allocatable :: rr(:),ph(:),taur(:,:,:),taup(:,:,:),uu(:,:,:),umat(:,:),ss(:)
  real(kind=8), allocatable :: y(:),yy(:,:),yyy(:,:,:),aa(:,:),wa(:,:),a(:,:,:),rr1(:),ph1(:),taur_1(:,:,:),taup_1(:,:,:)
  real(kind=8), dimension(:), allocatable :: y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,ak10,ak11
  real(kind=8) :: s,cx1,cx2,cx3,cy21,cy31,cy32,cy41,cy42,cy43,cy51,cy52,cy53,cy61,cy62,cy63,cy64,cy71,cy72,cy73,cy74,cy75,cy81,cy82
  real(kind=8) :: cy83,cy84,cy91,cy92,cy93,cy94,cy95,cy101,cy102,cy103,cy104,cy105,cy106,cy111,cy112,cy113,cy114,cy115,cy116
  
  read(1,*) mm1,mm2,n,path
  nn = n*(n-1)/2

  allocate(rr(mm1))
  allocate(ph(mm2))
  allocate(taur(mm1,mm2,nn))
  allocate(taup(mm1,mm2,nn))
  allocate(uu(mm1,mm2,n))
  allocate(umat(n,n))
  allocate(ss(nn))
  allocate(y(nn))
  allocate(yyy(mm1,mm2,nn))
  allocate(aa(n,n))
  allocate(wa(n,n))
  allocate(a(nn,n,n))
  allocate(rr1(mm1+2))
  allocate(ph1(mm2+2))
  allocate(taur_1(mm1+2,mm2+2,nn))
  allocate(taup_1(mm1+2,mm2+2,nn))
  allocate(y1(nn))
  allocate(y2(nn))
  allocate(y3(nn))
  allocate(y4(nn))
  allocate(y5(nn))
  allocate(y6(nn))
  allocate(y7(nn))
  allocate(y8(nn))
  allocate(y9(nn))
  allocate(y10(nn))
  allocate(y11(nn))
  allocate(ak1(nn))
  allocate(ak2(nn))
  allocate(ak3(nn))
  allocate(ak4(nn))
  allocate(ak5(nn))
  allocate(ak6(nn))
  allocate(ak7(nn))
  allocate(ak8(nn))
  allocate(ak9(nn))
  allocate(ak10(nn))
  allocate(ak11(nn))

  select case(path)
  case(1:4)
     allocate(yy(mm1,nn))
  case(5:8)
     allocate(yy(mm2,nn))
  end select
 
! include 'openfile.h'

  dr = rr(2)-rr(1)
  dp = ph(2)-ph(1)

  rr1(2:mm1+1) = rr(:)
  rr1(1) = rr(2)-dr
  rr1(mm1+2) = rr1(mm1+1)+dr

  ph1(2:mm2+1) = ph(:)
  ph1(1) = ph1(2)-dp
  ph1(mm2+2) = ph1(mm2+1)+dp

  taur_1(2:mm1+1,2:mm2+1,:) = taur(:,:,:)
  taup_1(2:mm1+1,2:mm2+1,:) = taup(:,:,:)
 
  taur_1(1,2,:) = taur_1(2,2,:)
  taur_1(mm1+2,2,:) = taur_1(mm1+1,2,:)
  taur_1(2,mm2+2,:) = taur_1(2,mm2+1,:)
  taur_1(mm1+1,mm2+2,:) = taur_1(mm1+1,mm2+1,:)

  taup_1(1,2,:) = taup_1(2,2,:)
  taup_1(mm1+2,2,:) = taup_1(mm1+1,2,:)
  taup_1(2,mm2+2,:) = taup_1(2,mm2+1,:)
  taup_1(mm1+1,mm2+2,:) = taup_1(mm1+1,mm2+1,:)

  taur_1(:,1,:) = taur_1(:,2,:)
  taup_1(:,1,:) = taup_1(:,2,:)

  taur_1(1,:,:) = taur_1(2,:,:)
  taup_1(1,:,:) = taup_1(2,:,:)

  taur_1(mm1+2,:,:) = taur_1(mm1+1,:,:)
  taup_1(mm1+2,:,:) = taup_1(mm1+1,:,:)

  taur_1(:,mm2+2,:) = taur_1(:,mm2+1,:)
  taup_1(:,mm2+2,:) = taup_1(:,mm2+1,:)

! include 'integration.h'
 
  deallocate(rr)
  deallocate(ph)
  deallocate(taur)
  deallocate(taup)
  deallocate(uu)
  deallocate(umat)
  deallocate(ss)
  deallocate(y)
  deallocate(yy)
  deallocate(yyy)
  deallocate(aa)
  deallocate(wa)
  deallocate(a)
  deallocate(rr1)
  deallocate(ph1)
  deallocate(taur_1)
  deallocate(taup_1)
  deallocate(y1)
  deallocate(y2)
  deallocate(y3)
  deallocate(y4)
  deallocate(y5)
  deallocate(y6)
  deallocate(y7)
  deallocate(y8)
  deallocate(y9)
  deallocate(y10)
  deallocate(y11)
  deallocate(ak1)
  deallocate(ak2)
  deallocate(ak3)
  deallocate(ak4)
  deallocate(ak5)
  deallocate(ak6)
  deallocate(ak7)
  deallocate(ak8)
  deallocate(ak9)
  deallocate(ak10)
  deallocate(ak11)

  stop

end program adt



subroutine amat(y,aa,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: y(n*(n-1)/2)
  real(kind=8), intent(out) :: aa(n,n)
  integer :: i,j,k,k1
  real(kind=8) :: aa1(n,n),b(n,n),a(n*(n-1)/2,n,n)

  a = 0.0d0
  do i = 1,n
    a(:,i,i) = 1.0d0
  enddo

  k1 = 0
  do j = 2,n
    do i = 1,j-1
      k1 = k1+1
      a(k1,i,i) = dcos(y(k1))
      a(k1,j,j) = dcos(y(k1))
      a(k1,i,j) = dsin(y(k1))
      a(k1,j,i) = -dsin(y(k1))
    enddo
  enddo

  aa = 0.0d0
  do i = 1,n
    aa(i,i) = 1.0d0
  enddo

  do k1 = 1,n*(n-1)/2
     b(:,:) = a(k1,:,:)
     aa1 = matmul(aa,b)
     aa = aa1
  enddo

  return

end subroutine amat



subroutine amatspl(y,aa,afor,aback,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: y(n*(n-1)/2)
  real(kind=8), intent(out) :: aa(n,n),afor(0:n*(n-1)/2,n,n),aback(n*(n-1)/2+1,n,n)
  integer :: i,j,k,k1
  real(kind=8) :: a(n*(n-1)/2,n,n),b(n,n),aa1(n,n)

  a = 0.0d0
  afor = 0.0d0
  aback = 0.0d0
  do i = 1,n
    a(:,i,i) = 1.0d0
    afor(0,i,i) = 1.0d0
    aback(n*(n-1)/2+1,i,i) = 1.0d0
  enddo

  k1 = 0
  do j = 2,n
    do i = 1,j-1
      k1 = k1+1
      a(k1,i,i) = dcos(y(k1))
      a(k1,j,j) = dcos(y(k1))
      a(k1,i,j) = dsin(y(k1))
      a(k1,j,i) = -dsin(y(k1))
    enddo
  enddo

  aa = 0.0d0
  do i = 1,n
    aa(i,i) = 1.0d0
  enddo

  do k1 = 1,n*(n-1)/2
     b(:,:) = a(k1,:,:)
     aa1 = matmul(aa,b)
     afor(k1,:,:) = aa1(:,:)
     aa = aa1
  enddo

  aa = 0.0d0
  do i = 1,n
    aa(i,i) = 1.0d0
  enddo

  do k1 = n*(n-1)/2,1,-1
     b(:,:) = a(k1,:,:)
     aa1 = matmul(aa,b)
     aback(k1,:,:) = aa1(:,:)
     aa = aa1
  enddo

  return

end subroutine amatspl



subroutine fexp(c1,c2,rr,ph,taup,y,fp,n,mm1,mm2)

  implicit none
  integer, intent(in) :: n,mm1,mm2
  real(kind=8), intent(in) :: c1,c2,rr(mm1+2),ph(mm2+2),taup(mm1+2,mm2+2,n*(n-1)/2),y(n*(n-1)/2)
  real(kind=8), intent(out) :: fp(n*(n-1)/2)
  integer :: i,j,nn
  real(kind=8) :: dsec,tau(n*(n-1)/2),sol(n*(n-1)/2),g(n*(n-1)/2,n*(n-1)/2),gi(n*(n-1)/2,n*(n-1)/2),a(n*(n-1)/2,n,n)

  nn = n*(n-1)/2
  call interpol(mm1,mm2,nn,rr,ph,taup,c1,c2,tau)
  call res(y,tau,fp,n)

  return
      
end subroutine fexp 



subroutine fexr(c1,c2,rr,ph,taur,y,fr,n,mm1,mm2)

  implicit none
  integer, intent(in) :: n,mm1,mm2
  real(kind=8), intent(in) :: c1,c2,rr(mm1+2),ph(mm2+2),taur(mm1+2,mm2+2,n*(n-1)/2),y(n*(n-1)/2)
  real(kind=8), intent(out) :: fr(n*(n-1)/2)
  integer :: i,j,nn
  real(kind=8) :: dsec,tau(n*(n-1)/2),sol(n*(n-1)/2),g(n*(n-1)/2,n*(n-1)/2),gi(n*(n-1)/2,n*(n-1)/2),a(n*(n-1)/2,n,n)

  nn = n*(n-1)/2
  call interpol(mm1,mm2,nn,rr,ph,taur,c1,c2,tau)
  call res(y,tau,fr,n)

  return
      
end subroutine fexr



subroutine gradcomat(y,afor,aback,g,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: afor(0:n*(n-1)/2,n,n),aback(n*(n-1)/2+1,n,n),y(n*(n-1)/2)
  real(kind=8), intent(out) :: g(n*(n-1)/2,n*(n-1)/2)
  integer :: i,j,k,k1,counter
  real(kind=8) :: adiff(n*(n-1)/2+1,n,n),aa1(n,n),aa2(n,n),b1(n,n),b2(n,n),b3(n,n)

  adiff = 0.0d0
  
  k1 = 0
  do j = 2,n
    do i = 1,j-1
      k1 = k1+1
      adiff(k1,i,i) = -dsin(y(k1))
      adiff(k1,j,j) = -dsin(y(k1))
      adiff(k1,i,j) = dcos(y(k1))
      adiff(k1,j,i) = -dcos(y(k1))
    enddo
  enddo

  do k1 = n*(n-1)/2,1,-1
     b1(:,:) = afor(k1-1,:,:)
     b2(:,:) = adiff(k1,:,:)
     b3(:,:) = aback(k1+1,:,:)
     aa1 = matmul(b1,b2)
     aa2 = matmul(aa1,b3)
     counter = 0
     do i = 2,n
       do j = 1,i-1
         counter = counter+1
         g(counter,k1) = aa2(i,j)
       enddo
     enddo
  enddo

  return

end subroutine gradcomat



subroutine inverse(g,n,gi)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: g(n,n)
  real(kind=8), intent(out) :: gi(n,n)
  integer :: i,j,ii,jj,irank,irow
  real(kind=8) :: zero,tmp,fc,a(n,n),en(n,n) 
  
  zero = 1.0d-20

  a = g 
  en = 0.0d0
  do i = 1,n
    en(i,i) = 1.0d0
  enddo

  gi = en
  irank = n
 
1 irow = n+1-irank

  if(dabs(a(irow,irow)) .le. zero) then
    j = 0
    do i = irow+1,n
      if(dabs(a(i,irow)) .gt. zero) then
        j = i
        exit
      endif
    enddo

    if(j .eq. 0) then
      write(*,*) "Inverse of the given matrix does not exist!"
      gi = 0.0d0
      goto 2
    else
      do i = 1,n
        tmp = a(irow,i)
        a(irow,i) = a(j,i)
        a(j,i) = tmp
        tmp = gi(irow,i)
        gi(irow,i) = gi(j,i)
        gi(j,i) = tmp
      enddo
    endif 
  endif 
   
  fc = a(irow,irow)

  a(irow,:) = a(irow,:)/fc
  gi(irow,:) = gi(irow,:)/fc

  do i = 1,irow-1
    fc = a(i,irow)
    a(i,:) = a(i,:)-a(irow,:)*fc 
    gi(i,:) = gi(i,:)-gi(irow,:)*fc
  enddo

  do i = irow+1,n
    fc = a(i,irow)
    a(i,:) = a(i,:)-a(irow,:)*fc 
    gi(i,:) = gi(i,:)-gi(irow,:)*fc
  enddo

  if(irank .eq. 1) goto 2

  irank = irank-1
  goto 1

2 return
  
end subroutine inverse



subroutine negtau(tau,taumat,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: tau(n*(n-1)/2)
  real(kind=8), intent(out) :: taumat(n,n)
  integer :: i,j,counter
 
  taumat = 0.0d0
  counter = 0
  do i = 2,n
    do j = 1,i-1
      counter = counter+1
      taumat(i,j) = tau(counter)
    enddo
  enddo

  do i = 1,n-1
    do j = i+1,n
      taumat(i,j) = -taumat(j,i)
    enddo
  enddo

  return
  
end subroutine negtau



subroutine res(ang,nact,f,n)

  implicit none
  integer, intent(in) :: n
  real(kind=8), intent(in) :: ang(n*(n-1)/2), nact(n*(n-1)/2)
  real(kind=8), intent(out) :: f(n*(n-1)/2)
  integer :: i,j,counter
  real(kind=8) :: val(n*(n-1)/2),g(n*(n-1)/2,n*(n-1)/2),gi(n*(n-1)/2,n*(n-1)/2),amat(n,n),tmat(n,n),prod(n,n)
  real(kind=8) :: a(n*(n-1)/2,n,n),afor(0:n*(n-1)/2,n,n),aback(n*(n-1)/2+1,n,n)

  call amatspl(ang,amat,afor,aback,n)
  call gradcomat(ang,afor,aback,g,n)
  call inverse(g,n*(n-1)/2,gi)
  call negtau(nact,tmat,n)
  
  prod = matmul(tmat,amat)

  counter = 0
  do i=2,n
    do j=1,i-1
      counter = counter+1
      val(counter) = prod(i,j)
    enddo
  enddo

  f = matmul(gi,val)

  return

end subroutine res



subroutine wmat(n,umat,aa,wa)

  implicit none
  integer, intent(in) :: 
end subroutine wmat 
