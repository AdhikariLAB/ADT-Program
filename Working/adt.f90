
!###########################################################################################################################
! Authors are Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
!###########################################################################################################################

! Compile this module using f2py as
! f2py -c adt.f90 -m adt_module --f90flags='-fopenmp' -lgomp only: get_angle amat


module adt
!$ use omp_lib
implicit none
! 'e' means expanded
real(8) ,allocatable, dimension(:,:,:) :: taur, taup, etaur, etaup
real(8) ,allocatable, dimension(:)     :: gridr, gridp, egridr, egridp
integer(8)                             :: ngridr, ngridp, ntau, nstate
real(8)                                :: wt(16,16), gridr_val, gridp_val, s,cx(3), cy(117)

data wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,&
        &4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,&
        &0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,&
        &5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

private :: wt, gridr_val, gridp_val, s,cx, cy
contains


! This subroutine returns the parameters of 8th order Runge-Kutta method  

subroutine init()
    s=dsqrt(21.0d0)

    cx(1)=0.5d0
    cx(2)=(7.0d0+s)/14.0d0
    cx(3)=(7.0d0-s)/14.0d0

    cy(21) =0.5d0
    cy(31) =0.25d0
    cy(32) =0.25d0
    cy(41) =1.0d0/7.0d0
    cy(42) =(-7.0d0-3.0d0*s)/98.0d0
    cy(43) =(21.0d0+5.0d0*s)/49.0d0
    cy(51) =(11.0d0+s)/84.0d0
    cy(52) =(18.0d0+4.0d0*s)/63.0d0
    cy(53) =(21.0d0-s)/252.0d0
    cy(61) =(5.0d0+s)/48.0d0
    cy(62) =(9.0d0+s)/36.0d0
    cy(63) =(231.0d0+14.0d0*s)/360.0d0
    cy(64) =(63.0d0-7.0d0*s)/80.0d0
    cy(71) =(10.0d0-s)/42.0d0
    cy(72) =(-432.0d0+92.0d0*s)/315.0d0
    cy(73) =(633.0d0-145.0d0*s)/90.0d0
    cy(74) =(-504.0d0+115.0d0*s)/70.0d0
    cy(75) =(63.0d0-13.0d0*s)/35.0d0
    cy(81) =1.0d0/14.0d0
    cy(82) =(14.0d0-3.0d0*s)/126.0d0
    cy(83) =(13.0d0-3.0d0*s)/63.0d0
    cy(84) =1.0d0/9.0d0
    cy(91) =1.0d0/32.0d0
    cy(92) =(91.0d0-21.0d0*s)/576.0d0
    cy(93) =11.0d0/72.0d0
    cy(94) =(-385.0d0-75.0d0*s)/1152.0d0
    cy(95) =(63.0d0+13.0d0*s)/128.0d0
    cy(101)=1.0d0/14.0d0
    cy(102)=1.0d0/9.0d0
    cy(103)=(-733.0d0-147.0d0*s)/2205.0d0
    cy(104)=(515.0d0+111.0d0*s)/504.0d0
    cy(105)=(-51.0d0-11.0d0*s)/56.0d0
    cy(106)=(132.0d0+28.0d0*s)/245.0d0
    cy(111)=(-42.0d0+7.0d0*s)/18.0d0
    cy(112)=(-18.0d0+28.0d0*s)/45.0d0
    cy(113)=(-273.0d0-53.0d0*s)/72.0d0
    cy(114)=(301.0d0+53.0d0*s)/72.0d0
    cy(115)=(28.0d0-28.0d0*s)/45.0d0
    cy(116)=(49.0d0-7.0d0*s)/18.0d0
end subroutine init


! This subroutine solves coupled differential equations by Runge-Kutta method   

subroutine rungeKutta8(fun,x,xy,y,dx,n)
    !takes y at x returns y at x+dx
    !fun = dy/dx
    external fun
    integer(8), intent(in):: n
    real(8), intent(in)   :: x,dx,xy
    real(8), intent(inout):: y(n)

    real(8) ,dimension(n) :: tmp,fo,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11


    call fun(x,xy, y, fo,n)
    a1=fo*dx
    tmp=y+a1*cy(21)

    call fun(x+cx(1)*dx,xy, tmp, fo,n)
    a2=fo*dx
    tmp=y+a1*cy(31)+a2*cy(32)

    call fun(x+cx(1)*dx,xy, tmp, fo,n)
    a3=fo*dx
    tmp=y+a1*cy(41)+a2*cy(42)+a3*cy(43)

    call fun(x+cx(2)*dx,xy, tmp, fo,n)
    a4=fo*dx
    tmp=y+a1*cy(51)+a3*cy(52)+a4*cy(53)

    call fun(x+cx(2)*dx,xy, tmp, fo,n)
    a5=fo*dx
    tmp=y+a1*cy(61)+a3*cy(62)+a4*cy(63)+a5*cy(64)

    call fun(x+cx(1)*dx,xy, tmp, fo,n)
    a6=fo*dx
    tmp=y+a1*cy(71)+a3*cy(72)+a4*cy(73)+a5*cy(74)+a6*cy(75)

    call fun(x+cx(3)*dx,xy, tmp, fo,n)
    a7=fo*dx
    tmp=y+a1*cy(81)+a5*cy(82)+a6*cy(83)+a7*cy(84)

    call fun(x+cx(3)*dx,xy, tmp, fo,n)
    a8=fo*dx
    tmp=y+a1*cy(91)+a5*cy(92)+a6*cy(93)+a7*cy(94)+a8*cy(95)

    call fun(x+cx(1)*dx,xy, tmp, fo,n)
    a9=fo*dx
    tmp=y+a1*cy(101)+a5*cy(102)+a6*cy(103)+a7*cy(104)+a8*cy(105)+a9*cy(106)

    call fun(x+cx(2)*dx,xy, tmp, fo,n)
    a10=fo*dx
    tmp=y+a5*cy(111)+a6*cy(112)+a7*cy(113)+a8*cy(114)+a9*cy(115)+a10*cy(116)

    call fun(x+dx,xy, tmp,fo,n)
    a11=fo*dx

    y=y+(9.0d0*a1+49.0d0*a8+64.0d0*a9+49.0d0*a10+9.0d0*a11)/180.0d0
end subroutine rungeKutta8


! This subroutine returns ADT angles over a 2D grid of geometries

subroutine get_angle(full_angle, ngridr, ngridp, ntau, path)    
    integer(8),  intent(in):: ngridr, ngridp, ntau, path
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)

    call init()

    select case (path)

        case (1) 
            call path1(full_angle, ngridr, ngridp, ntau)

        case (2)
            call path2(full_angle, ngridr, ngridp, ntau)

        case (3) 
            call path3(full_angle, ngridr, ngridp, ntau)

        case (4)
            call path4(full_angle, ngridr, ngridp, ntau)

        case (5)
            call path5(full_angle, ngridr, ngridp, ntau)

        case (6)
            call path6(full_angle, ngridr, ngridp, ntau)

        case (7)
            call path7(full_angle, ngridr, ngridp, ntau)

        case (8)
            call path8(full_angle, ngridr, ngridp, ntau)

        case default
            stop "This line should never be executed."

    end select
end subroutine get_angle



! This subroutine returns ADT angles over a 2D grid of geometries along path1

subroutine path1(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridp, ntau)
    integer(8)             :: i,j

    call init()
    h1 = gridr(2) - gridr(1)
    h2 = gridp(2) - gridp(1)

    angle=0.0d0

    gridp_val = gridp(1)-0.5d0*h2
    gridr_val = gridr(1)

    do i=1,ngridp
        call rungeKutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
        fangle(i,:) = angle
        gridp_val = gridp(i)
    enddo 

    !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
    do i=1,ngridp
        angle = fangle(i,:)
        gridp_val = gridp(i)
        gridr_val = gridr(1)-0.5d0*h1
        do j=1,ngridr
            call rungeKutta8(funcr, gridr_val, gridp_val, angle, h1, ntau)
            gridr_val = gridr(j)
            ! print *,omp_get_thread_num()
            full_angle(j,i,:) = angle
        enddo 
    enddo
    !$omp end parallel do

end subroutine path1



! This subroutine returns ADT angles over a 2D grid of geometries along path2

subroutine path2(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridp, ntau)
    integer(8)             :: i,j

    h1 = gridr(2) - gridr(1)
    h2 = -gridp(2) + gridp(1)
    angle=0.0d0


    gridp_val = gridp(ngridp)-0.5d0*h2
    gridr_val = gridr(1)

    do i=ngridp,1,-1
        call rungeKutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
        fangle(i,:) = angle
        gridp_val = gridp(i)
    enddo




    do i=ngridp,1,-1
        angle = fangle(i,:)
        gridr_val = gridr(1)-0.5d0*h1
        gridp_val = gridp(i)

        do j=1,ngridr
            call rungeKutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
            full_angle(j,i,:)=angle
            gridr_val = gridr(j)
        enddo
    enddo
end subroutine path2



! This subroutine returns ADT angles over a 2D grid of geometries along path3

subroutine path3(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridp, ntau)
    integer(8)             :: i,j

    h1 = -gridr(2) + gridr(1)
    h2 = gridp(2) - gridp(1)
    angle=0.0d0




    gridp_val = gridp(1)-0.5d0*h2
    gridr_val = gridr(ngridr)
    do i=1,ngridp
        call rungeKutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
        fangle(i,:) = angle
        gridp_val = gridp(i)
    enddo



    do i=1,ngridp
        angle = fangle(i,:)
        gridr_val = gridr(ngridr)-0.5d0*h1
        gridp_val = gridp(i)

        do j=ngridr,1,-1
            call rungeKutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
            full_angle(j,i,:)=angle
            gridr_val = gridr(j)
        enddo
    enddo
end subroutine path3



! This subroutine returns ADT angles over a 2D grid of geometries along path4

subroutine path4(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridp, ntau)
    integer(8)             :: i,j

    h1 = -gridr(2)+ gridr(1)
    h2 = -gridp(2)+ gridp(1)
    angle=0.0d0


    gridp_val = gridp(ngridp)-0.5d0*h2
    gridr_val = gridr(ngridr)

    do i=ngridp,1,-1
        call rungeKutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
        fangle(i,:) = angle
        gridp_val = gridp(i)

    enddo


    do i=ngridp,1,-1
        angle = fangle(i,:)
        gridr_val = gridr(ngridr)-0.5d0*h1
        gridp_val = gridp(i)

        do j=ngridr,1,-1
            call rungeKutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
            full_angle(j,i,:)=angle
            gridr_val = gridr(j)
        enddo
    enddo
end subroutine path4



! This subroutine returns ADT angles over a 2D grid of geometries along path5

subroutine path5(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridr, ntau)
    integer(8)             :: i,j

    h1 = gridr(2) - gridr(1)
    h2 = gridp(2) - gridp(1)
    angle=0.0d0


    gridr_val = gridr(1)-0.5d0*h1
    gridp_val = gridp(1)
    do i=1,ngridr
        call rungeKutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
        fangle(i,:) = angle
        gridr_val = gridr(i)
    enddo


    do i=1,ngridr
        angle = fangle(i,:)
        gridp_val = gridp(1)-0.5d0*h2
        gridr_val = gridr(i)

        do j=1,ngridp
            call rungeKutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
            full_angle(i,j,:)=angle
            gridp_val = gridp(j)
        enddo
    enddo
end subroutine path5



! This subroutine returns ADT angles over a 2D grid of geometries along path6

subroutine path6(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridr, ntau)
    integer(8)             :: i,j

    h1 = -gridr(2)+ gridr(1)
    h2 = gridp(2) - gridp(1)
    angle=0.0d0


    gridr_val = gridr(ngridr)-0.5d0*h1
    gridp_val = gridp(1)

    do i=ngridr,1,-1
        call rungeKutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
        fangle(i,:) = angle
        gridr_val = gridr(i)
    enddo

    do i=ngridr,1,-1
        angle = fangle(i,:)
        gridp_val = gridp(1)-0.5d0*h2
        gridr_val = gridr(i)

        do j=1,ngridp
            call rungeKutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
            full_angle(i,j,:)=angle
            gridp_val = gridp(j)
        enddo
    enddo
end subroutine path6



! This subroutine returns ADT angles over a 2D grid of geometries along path7

subroutine path7(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridr, ntau)
    integer(8)             :: i,j

    h1 = gridr(2) - gridr(1)
    h2 = -gridp(2) + gridp(1)
    angle=0.0d0

    gridr_val = gridr(1)-0.5d0*h1
    gridp_val = gridp(ngridp)

    do i=1,ngridr
        call rungeKutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
        fangle(i,:) = angle
        gridr_val = gridr(i)
    enddo

    do i=1,ngridr
        angle = fangle(i,:)
        gridp_val = gridp(ngridp)-0.5d0*h2
        gridr_val = gridr(i)
        do j=ngridp,1,-1
            call rungeKutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
            full_angle(i,j,:)=angle
            gridp_val = gridp(j)
        enddo
    enddo

end subroutine path7



! This subroutine returns ADT angles over a 2D grid of geometries along path8

subroutine path8(full_angle, ngridr, ngridp, ntau)     
    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
    real(8)                :: angle(ntau), h1, h2, fangle(ngridr, ntau)
    integer(8)             :: i,j

    h1 = -gridr(2) + gridr(1)
    h2 = gridp(2) - gridp(1)
    angle=0.0d0

    gridr_val = gridr(ngridr)-0.5d0*h1
    gridp_val = gridp(ngridp)

    do i=ngridr,1,-1
        call rungeKutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
        fangle(i,:) = angle
        gridr_val = gridr(i)
    enddo


    do i=ngridr,1,-1
        angle = fangle(i,:)
        gridp_val = gridp(ngridp)-0.5d0*h2
        gridr_val = gridr(i)

        do j=ngridp,1,-1
            call rungeKutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
            full_angle(i,j,:)=angle
            gridp_val = gridp(j)
        enddo
    enddo

end subroutine path8



!*****************************************************************************************
!This subroutine returns the magnitude of $^{N}$C$_{2}$ numbers of ADT angles at a 
!specific nuclear geometry using a pre-calculated/predefined set of ADT angles (initial 
!value). This routine is implemented during integration along first coordinate when the 
!second coordinate is fixed at a definite value 
!*****************************************************************************************

subroutine funcr(gridr_val,gridp_val, ini_angle, output,ntau)


    integer(8),  intent(in):: ntau
    real(8), intent(in)    :: gridr_val, ini_angle(ntau), gridp_val
    real(8), intent(out)   :: output(ntau)
    real(8) :: tau_val(ntau)


    call interpol(etaur, gridr_val, gridp_val, tau_val, ngridr, ngridp, ntau)
    call res(ini_angle,tau_val,output, ntau, nstate)
end subroutine funcr



!*****************************************************************************************
!This subroutine returns the magnitude of $^{N}$C$_{2}$ numbers of ADT angles at a 
!specific nuclear geometry using a pre-calculated/predefined set of ADT angles (initial 
!value). This routine is implemented during integration along second coordinate when the 
!first coordinate is fixed at a definite value 
!*****************************************************************************************

subroutine funcp(gridp_val,gridr_val, ini_angle, output,ntau)

    integer(8),  intent(in):: ntau
    real(8), intent(in) :: gridp_val, ini_angle(ntau), gridr_val
    real(8), intent(out) :: output(ntau)
    real(8) :: outval(ntau)

    call interpol(etaup, gridr_val, gridp_val, outval, ngridr, ngridp, ntau)
    call res(ini_angle,outval,output, ntau, nstate)
end subroutine funcp



!***********************************************************************************************************
!      Subroutine 'INTERPOL' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
!      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
!      Delhi-110040, 2000.')
!***********************************************************************************************************

!*****************************************************************************************************
!While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
!intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
!is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used to  
!perform this interpolation.
!*****************************************************************************************************

subroutine interpol(tau,x1,y1,tout,ngridr, ngridp, ntau)

    integer(8),  intent(in):: ngridr, ngridp, ntau
    real(8),intent(in) :: tau(ngridr+2, ngridp+2, ntau), x1,y1
    real(8),intent(out):: tout(ntau)
    real(8)            :: yf(4),yf1(4),yf2(4),yf12(4), dx,dy,x1l,x1u,x2l,x2u,ansy!,ansy1,ansy2        
    integer(8)         :: ii1,ii2,jj1,jj2,k



    dx=egridr(2)-egridr(1)
    dy=egridp(2)-egridp(1)
    ! call locate(egridr,ngridr+2,x1,ii1)
    ii1 = int((x1-egridr(1))/dx)+1

    ii2=ii1+1
    x1l=egridr(ii1)
    x1u=egridr(ii2)
    ! call locate(egridp,ngridp+2,y1,jj1)
    jj1 = int((y1-egridp(1))/dy)+1
    jj2=jj1+1
    x2l=egridp(jj1)
    x2u=egridp(jj2)


    do k=1,ntau
        !zeroth derivatives
        yf(1)=tau(ii1,jj1,k)
        yf(2)=tau(ii2,jj1,k)
        yf(3)=tau(ii2,jj2,k)
        yf(4)=tau(ii1,jj2,k)
        !first derivatives

        if(ii1.eq.1)then
            !forward differencing(f.d)
            yf1(1)=(tau(ii1+1,jj1,k)-tau(ii1,jj1,k))/dx
            yf1(4)=(tau(ii1+1,jj2,k)-tau(ii1,jj2,k))/dx
        else
            !central differencing(c.d)
            yf1(1)=0.5d0*(tau(ii1+1,jj1,k)-tau(ii1-1,jj1,k))/dx
            yf1(4)=0.5d0*(tau(ii1+1,jj2,k)-tau(ii1-1,jj2,k))/dx
        endif

        if(ii2.eq.ngridr+2)then
            !backward differencing(b.d)
            yf1(2)=(tau(ii2,jj1,k)-tau(ii2-1,jj1,k))/dx
            yf1(3)=(tau(ii2,jj2,k)-tau(ii2-1,jj2,k))/dx
          else
            !c.d
            yf1(2)=0.5d0*(tau(ii2+1,jj1,k)-tau(ii2-1,jj1,k))/dx
            yf1(3)=0.5d0*(tau(ii2+1,jj2,k)-tau(ii2-1,jj2,k))/dx
        endif

        if(jj1.eq.1)then
            !f.d
            yf2(1)=(tau(ii1,jj1+1,k)-tau(ii1,jj1,k))/dy
            yf2(2)=(tau(ii2,jj1+1,k)-tau(ii2,jj1,k))/dy
          else
            !c.d
            yf2(1)=0.5d0*(tau(ii1,jj1+1,k)-tau(ii1,jj1-1,k))/dy
            yf2(2)=0.5d0*(tau(ii2,jj1+1,k)-tau(ii2,jj1-1,k))/dy
        endif

        if(jj2.eq.ngridp+2)then
            !b.d
            yf2(3)=(tau(ii2,jj2,k)-tau(ii2,jj2-1,k))/dy
            yf2(4)=(tau(ii1,jj2,k)-tau(ii1,jj2-1,k))/dy
          else
            !c.d
            yf2(3)=0.5d0*(tau(ii2,jj2+1,k)-tau(ii2,jj2-1,k))/dy
            yf2(4)=0.5d0*(tau(ii1,jj2+1,k)-tau(ii1,jj2-1,k))/dy
        endif

        !second derivatives
        if(ii1.eq.1.and.jj1.eq.1)then
            !f.d & f.d
            yf12(1)=((tau(ii1+1,jj1+1,k)-tau(ii1,jj1+1,k))-(tau(ii1+1,jj1,k)-tau(ii1,jj1,k)))/dx/dy
          elseif(ii1.eq.1.and.jj1.ne.1)then
            !f.d & c.d
            yf12(1)=0.5d0*((tau(ii1+1,jj1+1,k)-tau(ii1,jj1+1,k))-(tau(ii1+1,jj1-1,k)-tau(ii1,jj1-1,k)))/dx/dy
          elseif(ii1.ne.1.and.jj1.eq.1)then
            !c.d & f.d
            yf12(1)=0.5d0*((tau(ii1+1,jj1+1,k)-tau(ii1-1,jj1+1,k))-(tau(ii1+1,jj1,k)-tau(ii1-1,jj1,k)))/dx/dy
          else
            !c.d & c.d
            yf12(1)=0.25d0*((tau(ii1+1,jj1+1,k)-tau(ii1-1,jj1+1,k))-(tau(ii1+1,jj1-1,k)-tau(ii1-1,jj1-1,k)))/dx/dy
        endif

        if(ii2.eq.ngridr+2.and.jj1.eq.1)then
            !b.d & f.d
            yf12(2)=((tau(ii2,jj1+1,k)-tau(ii2-1,jj1+1,k))-(tau(ii2,jj1,k)-tau(ii2-1,jj1,k)))/dx/dy
          elseif(ii2.eq.ngridr+2.and.jj1.ne.1)then
            !b.d & c.d
            yf12(2)=0.5d0*((tau(ii2,jj1+1,k)-tau(ii2-1,jj1+1,k))-(tau(ii2,jj1-1,k)-tau(ii2-1,jj1-1,k)))/dx/dy
          elseif(ii2.ne.ngridr+2.and.jj1.eq.1)then
            !c.d & f.d
            yf12(2)=0.5d0*((tau(ii2+1,jj1+1,k)-tau(ii2-1,jj1+1,k))-(tau(ii2+1,jj1,k)-tau(ii2-1,jj1,k)))/dx/dy
          else
            !c.d & c.d
            yf12(2)=0.25d0*((tau(ii2+1,jj1+1,k)-tau(ii2-1,jj1+1,k))-(tau(ii2+1,jj1-1,k)-tau(ii2-1,jj1-1,k)))/dx/dy
        endif

        if(ii2.eq.ngridr+2.and.jj2.eq.ngridp+2)then
            !b.d & b.d
            yf12(3)=((tau(ii2,jj2,k)+tau(ii2-1,jj2,k))-(tau(ii2,jj2-1,k)-tau(ii2-1,jj2-1,k)))/dx/dy
          elseif(ii2.eq.ngridr+2.and.jj2.ne.ngridp+2)then
            !b.d & c.d
            yf12(3)=0.5d0*((tau(ii2,jj2+1,k)+tau(ii2-1,jj2+1,k))-(tau(ii2,jj2-1,k)-tau(ii2-1,jj2-1,k)))/dx/dy
          elseif(ii2.ne.ngridr+2.and.jj2.eq.ngridp+2)then
            !c.d & b.d
            yf12(3)=0.5d0*((tau(ii2+1,jj2,k)+tau(ii2-1,jj2,k))-(tau(ii2+1,jj2-1,k)-tau(ii2-1,jj2-1,k)))/dx/dy
          else
            !c.d & c.d
            yf12(3)=0.25d0*((tau(ii2+1,jj2+1,k)+tau(ii2-1,jj2+1,k))-(tau(ii2+1,jj2-1,k)-tau(ii2-1,jj2-1,k)))/dx/dy
        endif

        !f.d & b.d
        if(ii1.eq.1.and.jj2.eq.ngridp+2)then
            yf12(4)=((tau(ii1+1,jj2,k)-tau(ii1,jj2,k))-(tau(ii1+1,jj2-1,k)-tau(ii1,jj2-1,k)))/dx/dy
          elseif(ii1.eq.1.and.jj2.ne.ngridp+2)then
            !f.d & c.d
            yf12(4)=0.5d0*((tau(ii1+1,jj2+1,k)-tau(ii1,jj2+1,k))-(tau(ii1+1,jj2-1,k)-tau(ii1,jj2-1,k)))/dx/dy
          elseif(ii1.ne.1.and.jj2.eq.ngridp+2)then
            !c.d & b.d
            yf12(4)=0.5d0*((tau(ii1+1,jj2,k)-tau(ii1-1,jj2,k))-(tau(ii1+1,jj2-1,k)-tau(ii1-1,jj2-1,k)))/dx/dy
          else
            !c.d & c.d
            yf12(4)=0.25d0*((tau(ii1+1,jj2+1,k)-tau(ii1-1,jj2+1,k))-(tau(ii1+1,jj2-1,k)-tau(ii1-1,jj2-1,k)))/dx/dy
        endif

        call bcuint(yf,yf1,yf2,yf12,x1l,x1u,x2l,x2u,x1,y1,ansy)
        tout(k)=ansy
    enddo
end subroutine interpol



! subroutine locate(xx,n,x,jj)

!     integer(8), intent(in) :: n
!     real(8),    intent(in) :: xx(n),x
!     integer(8), intent(out):: jj
!     integer(8)             :: jl,ju,jm

!     jl=0
!     ju=n+1
!     do while((ju-jl).gt.1)
!         jm=(ju+jl)/2
!         if((xx(n).gt.xx(1)).eqv.(x.ge.xx(jm)))then
!             jl=jm
!           else
!             ju=jm
!         endif
!     enddo
!     jj=jl
! end subroutine locate



!***********************************************************************************************************
!      Subroutine 'BCUINT' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
!      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
!      Delhi-110040, 2000.')
!***********************************************************************************************************

!*****************************************************************************************************
!While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
!intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
!is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used as a 
!secondary subroutine in 'INTERPOL'.
!*****************************************************************************************************

subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy)

    real(8), intent(in) :: x1,x1l,x1u,x2,x2l,x2u,y(4),y1(4),y12(4),y2(4)
    real(8), intent(out)::ansy!,ansy1,ansy2

    integer(8):: i
    real(8):: t,u,c(4,4)
    call bcucof (y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
    if(x1u.eq.x1l .or. x2u.eq.x2l) stop 'bad input in bcuint'
    t=(x1-x1l)/(x1u-x1l)
    u=(x2-x2l)/(x2u-x2l)
    ansy=0.0d0
    do i=4,1,-1
        ansy =t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
    enddo
end subroutine bcuint



!***********************************************************************************************************
!      Subroutine 'BCUCOF' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery, 
!      Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela, 
!      Delhi-110040, 2000.')
!***********************************************************************************************************

!*****************************************************************************************************
!While integrating the differential equations by 8th order Runge-Kutta method, NACT values for the 
!intermediate geometries between two grid points are required and therefore, bi-cubic interpolation 
!is adopted to dig out the magnitude of NACTs at those unknown points. This subroutine is used as a 
!secondary subroutine in 'INTERPOL'.
!*****************************************************************************************************

subroutine bcucof (y,y1,y2,y12,d1,d2,c)

    real(8), intent(in) :: d1,d2,y(4),y1(4),y12(4),y2(4)
    real(8), intent(out):: c(4,4)
    real(8)             :: x(16)
    x=(/y, y1*d1, y2*d2, y12*d1*d2/)
    c = transpose(reshape(matmul(wt,x),(/4,4/)))
end subroutine bcucof



!**********************************
!      Subroutine 'RES'
!**********************************

!**************************************************************************************
!Magnitudes of gradient of ADT angles at every grid point in the nuclear
!CS is returned by this subroutine
!**************************************************************************************

subroutine res(ang,nact,f, ntau, nstate)


    integer(8),  intent(in)::  ntau, nstate
    real(8), intent(in) :: ang(ntau), nact(ntau)
    real(8), intent(out) :: f(ntau)
    integer(8) :: i,j,counter
    real(8) :: val(ntau),g(ntau,ntau),gi(ntau,ntau),amat(nstate,nstate),tmat(nstate,nstate),prod(nstate,nstate)
    real(8) :: afor(0:ntau,nstate,nstate),aback(ntau+1,nstate,nstate)


    call amatspl(ang,amat,afor,aback, ntau, nstate)
    call gradcomat(ang,afor,aback,g, ntau, nstate)
    call inverse(g,gi,ntau)
    call negtau(nact,tmat, ntau, nstate)

    prod = matmul(tmat,amat)

    counter = 0
    do i=2,nstate
        do j=1,i-1
            counter = counter+1
            val(counter) = prod(i,j)
        enddo
    enddo

    f = matmul(gi,val)
end subroutine res



!**********************************
!      Subroutine 'AMATSPL'
!**********************************

!**********************************************************************************************************
!This subroutine is implemented not only to generate the complete ADT matrix (AA), but also two sets of 
!partially multiplied ADT matrices (AFOR, ABACK). Adiabatic to diabatic transformation matrix
!is generated by multiplying elementary rotation matrices in a definite order, but partial ADT matrices 
!can be constructed by collecting one or more [2,3,4,.....,(N-1); N = number of elementary matrices] 
!matrices from that set and multiplying them in the same order as the parent one.
!**********************************************************************************************************

subroutine amatspl(y,aa,afor,aback, ntau, nstate)


    integer(8),  intent(in)::  ntau, nstate
    real(8), intent(in) :: y(ntau)
    real(8), intent(out) :: aa(nstate,nstate),afor(0:ntau,nstate,nstate),aback(ntau+1,nstate,nstate)
    integer(8) :: i,j,k
    real(8) :: a(ntau,nstate,nstate),aa1(nstate,nstate)

    a = 0.0d0
    aa=0d0
    aa1=0d0
    afor = 0.0d0
    aback = 0.0d0
    do i = 1,nstate
        a(:,i,i) = 1.0d0
        aa(i,i) = 1.0d0
        aa1(i,i) = 1.0d0
        afor(0,i,i) = 1.0d0
        aback(ntau+1,i,i) = 1.0d0
    enddo

    k = 0
    do j = 2,nstate
        do i = 1,j-1
            k = k+1
            a(k,i,i) = dcos(y(k))
            a(k,j,j) = dcos(y(k))
            a(k,i,j) = dsin(y(k))
            a(k,j,i) = -dsin(y(k))
        enddo
    enddo


    do k = 1,ntau
        aa1 = matmul(aa1,a(k,:,:))
        afor(k,:,:) = aa1
    enddo

    do k = ntau,1,-1 
        aa = matmul(a(k,:,:),aa)
        aback(k,:,:) = aa
    enddo
end subroutine amatspl



!**********************************
!      Subroutine 'INVERSE'
!**********************************

!********************************************************************************
!The inverse of coefficient matrix of gradient of ADT angles is
!evaluated by Gauss-Jordon Method.
!********************************************************************************
      
subroutine inverse(g,gi, ntau)


    integer(8),  intent(in)::  ntau
    real(8), intent(in) :: g(ntau,ntau)
    real(8), intent(out) :: gi(ntau,ntau)
    integer(8) :: i,j,irank,irow
    real(8) :: zero,tmp,fc,a(ntau,ntau),en(ntau,ntau)

    zero = 1.0d-20
    a = g
    en = 0.0d0

    forall(i=1:ntau) en(i,i) = 1.0d0


    gi = en
    irank = ntau
    ! 1 irow = ntau+1-irank
    do while(irank .ne. 0)
        irow = ntau+1-irank


        if(dabs(a(irow,irow)) .le. zero) then
            j = 0
            do i = irow+1,ntau
                if(dabs(a(i,irow)) .gt. zero) then
                    j = i
                    exit
                endif
            enddo

            if(j .eq. 0) then
                write(*,*) "Inverse of the given matrix does not exist!"
                gi = 0.0d0
                return
              else
                do i = 1,ntau
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

        do i = irow+1,ntau
            fc = a(i,irow)
            a(i,:) = a(i,:)-a(irow,:)*fc
            gi(i,:) = gi(i,:)-gi(irow,:)*fc
        enddo

        irank = irank-1
    enddo
end subroutine inverse



!**********************************
!      Subroutine 'NEGTAU'
!**********************************

!*****************************************************************************************
!This subprogram accepts the NACM as input and returns the negative form
!of the matrix at each and every grid point in nuclear CS.         
!*****************************************************************************************

subroutine negtau(tau,taumat, ntau, nstate)


    integer(8),  intent(in):: ntau, nstate
    real(8), intent(in) :: tau(ntau)
    real(8), intent(out) :: taumat(nstate,nstate)
    integer(8) :: i,j,counter

    taumat = 0.0d0
    counter = 0
    do i = 2,nstate
        do j = 1,i-1
            counter = counter+1
            taumat(i,j) = tau(counter)
        enddo
    enddo

    do i = 1,nstate-1
        do j = i+1,nstate
            taumat(i,j) = -taumat(j,i)
        enddo
    enddo
end subroutine negtau



!**********************************
!      Subroutine 'GRADCOMAT'
!**********************************

!**************************************************************************************
!This subprogram is implemented for calculating the elements of coefficient matrix 
!of gradient of ADT angles
!**************************************************************************************

subroutine gradcomat(y,afor,aback,g, ntau, nstate)


    integer(8),  intent(in):: ntau, nstate
    real(8), intent(in) :: afor(0:ntau,nstate,nstate),aback(ntau+1,nstate,nstate),y(ntau)
    real(8), intent(out) :: g(ntau,ntau)
    integer(8) :: i,j,k,counter
    real(8) :: adiff(ntau+1,nstate,nstate),aa1(nstate,nstate),aa2(nstate,nstate),b1(nstate,nstate),b2(nstate,nstate),&
                    &b3(nstate,nstate)

    adiff = 0.0d0

    k = 0
    do j = 2,nstate
        do i = 1,j-1
            k = k+1
            adiff(k,i,i) = -dsin(y(k))
            adiff(k,j,j) = -dsin(y(k))
            adiff(k,i,j) = dcos(y(k))
            adiff(k,j,i) = -dcos(y(k))
        enddo
    enddo

    do k = ntau,1,-1
        b1 = afor(k-1,:,:)
        b2 = adiff(k,:,:)
        b3 = aback(k+1,:,:)
        aa1 = matmul(b1,b2)
        aa2 = matmul(aa1,b3)
        counter = 0
        do i = 2,nstate
            do j = 1,i-1
                counter = counter+1
                g(counter,k) = aa2(i,j)
            enddo
        enddo
    enddo
end subroutine gradcomat



!**********************************
!      Subroutine 'AMAT'
!**********************************

!*********************************************************************
!This subroutine returns numerically calculated ADT matrix (AA) for a
!particular set of ADT angles (Y).
!*********************************************************************

subroutine amat(y,aa, ntau, nstate)
    integer(8),intent(in):: ntau, nstate
    integer(8):: i,j,k
    real(8), intent(in)::y(ntau)
    real(8), intent(out):: aa(nstate,nstate)
    real(8)::a(ntau,nstate,nstate)
    !f2py intent(hide) :: ntau = shape(y,0)

    a=0.d0
    aa=0.d0
    do i=1,nstate
        a(:,i,i)=1.0d0
        aa(i,i) =1.0d0
    enddo

    k=0

    do j=2,nstate
        do i=1,j-1
            k=k+1 !i.e.a12 is denoted as a1, a13 is denoted as a2, a23 is denoted as a3 etc.

            a(k,i,i)= dcos(y(k))
            a(k,j,j)= a(k,i,i)          
            a(k,i,j)= dsin(y(k))
            a(k,j,i)=-a(k,i,j)
        enddo
    enddo

    do i=1,ntau
        aa = matmul(aa,a(i,:,:))
    enddo
end subroutine amat


end module adt
