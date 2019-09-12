
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                      !
!    This fortran file contains all necessary fortran subroutines to solve the adiabatic to diabatic transformation    !
!    (ADT) equations along eight different paths of integration. The stiff differential equations are solved by        !
!    employing 8th order Runge-Kutta method. While performing the initialization step, 'f2py' generates the python     !
!    module, 'adt_module.so' according to user specified compiler flags. In order to carry out this step, see the      !
!    instructions in 'init.sh'.                                                                                        !
!                                                                                                                      !
!    Written by Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Satyam Ravi, Saikat Mukherjee, Subhankar Sardar and !
!    Satrajit Adhikari                                                                                                 !
!                                                                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module adt
    !$ use omp_lib
    implicit none
    ! 'e' means expanded
    real(8) ,allocatable, dimension(:,:,:) :: taur, taup, etaur, etaup
    real(8) ,allocatable, dimension(:)     :: gridr, gridp, egridr, egridp, grid
    integer(8)                             :: ngridr, ngridp, ngrid, ntau, nstate
    real(8) ,allocatable, dimension(:,:)   :: tau
    integer(8),allocatable,dimension(:)    :: order
    ! tau and grid varaible are only here just for the 1D ADT case
    real(8)                                :: wt(16,16), gridr_val, gridp_val, s,cx(3), cy(117)

    data wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,&
            &4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,&
            &0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,&
            &5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

    private :: wt, gridr_val, gridp_val, s,cx, cy
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine init()

    !This subroutine returns the parameters of 8th order Runge-Kutta method

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_angle(full_angle, ngridr, ngridp, ntau, path)

        ! This subroutine returns ADT angles over a 2D grid of geometries along any one of the eight paths of integration

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path1(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path1

        integer(8),  intent(in):: ngridr, ngridp, ntau
        real(8),    intent(out):: full_angle(ngridr,ngridp,ntau)
        real(8)                :: angle(ntau), h1, h2, fangle(ngridp, ntau)
        integer(8)             :: i,j

        h1 = gridr(2) - gridr(1)
        h2 = gridp(2) - gridp(1)

        angle=0.0d0

        gridp_val = gridp(1)-0.5d0*h2
        gridr_val = gridr(1)

        do i=1,ngridp
            call rungekutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
            fangle(i,:) = angle
            gridp_val = gridp(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=1,ngridp
            angle = fangle(i,:)
            gridp_val = gridp(i)
            gridr_val = gridr(1)-0.5d0*h1
            do j=1,ngridr
                call rungekutta8(funcr, gridr_val, gridp_val, angle, h1, ntau)
                gridr_val = gridr(j)
                ! print *,omp_get_thread_num()
                full_angle(j,i,:) = angle
            enddo
        enddo
        !$omp end parallel do

    end subroutine path1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path2(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path2

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
            call rungekutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
            fangle(i,:) = angle
            gridp_val = gridp(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=ngridp,1,-1
            angle = fangle(i,:)
            gridr_val = gridr(1)-0.5d0*h1
            gridp_val = gridp(i)

            do j=1,ngridr
                call rungekutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
                full_angle(j,i,:)=angle
                gridr_val = gridr(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path3(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path3

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
            call rungekutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
            fangle(i,:) = angle
            gridp_val = gridp(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=1,ngridp
            angle = fangle(i,:)
            gridr_val = gridr(ngridr)-0.5d0*h1
            gridp_val = gridp(i)

            do j=ngridr,1,-1
                call rungekutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
                full_angle(j,i,:)=angle
                gridr_val = gridr(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path3

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path4(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path4

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
            call rungekutta8(funcp, gridp_val,gridr_val,angle,h2,ntau)
            fangle(i,:) = angle
            gridp_val = gridp(i)

        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=ngridp,1,-1
            angle = fangle(i,:)
            gridr_val = gridr(ngridr)-0.5d0*h1
            gridp_val = gridp(i)

            do j=ngridr,1,-1
                call rungekutta8(funcr, gridr_val,gridp_val, angle, h1, ntau)
                full_angle(j,i,:)=angle
                gridr_val = gridr(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path5(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path5

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
            call rungekutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
            fangle(i,:) = angle
            gridr_val = gridr(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=1,ngridr
            angle = fangle(i,:)
            gridp_val = gridp(1)-0.5d0*h2
            gridr_val = gridr(i)

            do j=1,ngridp
                call rungekutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
                full_angle(i,j,:)=angle
                gridp_val = gridp(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path5

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path6(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path6

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
            call rungekutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
            fangle(i,:) = angle
            gridr_val = gridr(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=ngridr,1,-1
            angle = fangle(i,:)
            gridp_val = gridp(1)-0.5d0*h2
            gridr_val = gridr(i)

            do j=1,ngridp
                call rungekutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
                full_angle(i,j,:)=angle
                gridp_val = gridp(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path6

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path7(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path7

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
            call rungekutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
            fangle(i,:) = angle
            gridr_val = gridr(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=1,ngridr
            angle = fangle(i,:)
            gridp_val = gridp(ngridp)-0.5d0*h2
            gridr_val = gridr(i)
            do j=ngridp,1,-1
                call rungekutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
                full_angle(i,j,:)=angle
                gridp_val = gridp(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path7

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine path8(full_angle, ngridr, ngridp, ntau)

        ! This subroutine returns ADT angles over a 2D grid of geometries along path8

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
            call rungekutta8(funcr, gridr_val,gridp_val,angle,h1,ntau)
            fangle(i,:) = angle
            gridr_val = gridr(i)
        enddo

        !$omp parallel do default(shared) private(angle, gridp_val, gridr_val, i,j)
        do i=ngridr,1,-1
            angle = fangle(i,:)
            gridp_val = gridp(ngridp)-0.5d0*h2
            gridr_val = gridr(i)

            do j=ngridp,1,-1
                call rungekutta8(funcp, gridp_val,gridr_val, angle, h2, ntau)
                full_angle(i,j,:)=angle
                gridp_val = gridp(j)
            enddo
        enddo
        !$omp end parallel do

    end subroutine path8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine rungekutta8(fun,x,xy,y,dx,n)

        ! This subroutine solves coupled differential equations (ADT equations) by 8th order Runge-Kutta method

        ! takes y at x returns y at x+dx
        ! fun = dy/dx

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

    end subroutine rungekutta8

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine funcr(gridr_val,gridp_val, ini_angle, output,ntau)

        ! This subroutine returns the magnitude of ADT angles at a specific nuclear geometry during integration along
        ! first coordinate

        integer(8),  intent(in):: ntau
        real(8), intent(in)    :: gridr_val, ini_angle(ntau), gridp_val
        real(8), intent(out)   :: output(ntau)
        real(8) :: tau_val(ntau)


        call interpol(etaur, gridr_val, gridp_val, tau_val, ngridr, ngridp, ntau)
        call res(ini_angle,tau_val,output, ntau, nstate)

    end subroutine funcr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine funcp(gridp_val,gridr_val, ini_angle, output,ntau)

        ! This subroutine returns the magnitude of ADT angles at a specific nuclear geometry during integration along
        ! second coordinate

        integer(8),  intent(in):: ntau
        real(8), intent(in) :: gridp_val, ini_angle(ntau), gridr_val
        real(8), intent(out) :: output(ntau)
        real(8) :: outval(ntau)

        call interpol(etaup, gridr_val, gridp_val, outval, ngridr, ngridp, ntau)
        call res(ini_angle,outval,output, ntau, nstate)

    end subroutine funcp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !***********************************************************************************************************
    !   Subroutine 'interpol' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery,
    !   Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela,
    !   Delhi-110040, 2000.')
    !***********************************************************************************************************

    subroutine locate(xx,n,x,jj)

        integer(8), intent(in) :: n
        real(8),    intent(in) :: xx(n),x
        integer(8), intent(out):: jj
        integer(8)             :: jl,ju,jm
    
        jl=0
        ju=n+1
        do while((ju-jl).gt.1)
            jm=(ju+jl)/2
            if((xx(n).gt.xx(1)).eqv.(x.ge.xx(jm)))then
                jl=jm
              else
                ju=jm
            endif
        enddo
        jj=jl
    end subroutine locate
    

    subroutine interpol(tau,x1,y1,tout,ngridr, ngridp, ntau)

        ! This subroutine is used to perform bi-cubic interpolation to evaluate the magnitude of nonadiabatic coupling terms
        ! (NACTs) at intermediate geometries between two grid points

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !***********************************************************************************************************
    !   Subroutine 'bcucof' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery,
    !   Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela,
    !   Delhi-110040, 2000.')
    !***********************************************************************************************************

    subroutine bcucof (y,y1,y2,y12,d1,d2,c)

        ! This subroutine is used as a secondary subroutine in 'interpol'

        real(8), intent(in) :: d1,d2,y(4),y1(4),y12(4),y2(4)
        real(8), intent(out):: c(4,4)
        real(8)             :: x(16)
        x=(/y, y1*d1, y2*d2, y12*d1*d2/)
        c = transpose(reshape(matmul(wt,x),(/4,4/)))
    end subroutine bcucof

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !***********************************************************************************************************
    !   Subroutine 'bcuint' (taken from 'W. H. Press, S. A. Teukolsky, W.T. Vetterling, B. P. Flannery,
    !   Numerical Recipes in FORTRAN, Cambridge University Press, A-229, DSIDC Industrial Park, Narela,
    !   Delhi-110040, 2000.')
    !***********************************************************************************************************

    subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy)

        ! This subroutine is used as a secondary subroutine in 'interpol'

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !
    subroutine get_angle1d(ngrid, ntau,full_angle)

        ! This subroutine returns ADT angles over a 1D grid of geometries 
        implicit none
        integer(8),  intent(in):: ngrid, ntau

        real(8),    intent(out):: full_angle(ngrid,ntau)
        real(8)                ::  h1, diff(ngrid,ntau), angle(ntau)
        integer(8)             :: i

        call init()
        h1 = grid(2) - grid(1)

        angle=0.0d0

        ! store the derivatives in diff variable
        do i=1,ntau
            call spline(grid, tau(:,i), ngrid, diff(:,i) )
        enddo

        full_angle(1,:) = 0
        do i=2,ngrid
            call rungekutta81( grid(i-1), angle, h1, diff,ngrid,ntau)
            full_angle(i,:) = angle
        enddo

    end subroutine get_angle1d




    subroutine rungekutta81(x,y,dx,diff, m,n)

        ! This subroutine solves coupled differential equations (ADT equations) by 8th order Runge-Kutta method


        integer(8), intent(in):: m,n
        real(8), intent(in)   :: x,dx, diff(m,n)
        real(8), intent(inout):: y(n)

        real(8) ,dimension(n) :: tmp,fo,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11


        call func( diff, y, x, m,n, fo)
        a1=fo*dx
        tmp=y+a1*cy(21)

        call func( diff, tmp, x+cx(1)*dx, m,n,fo)
        a2=fo*dx
        tmp=y+a1*cy(31)+a2*cy(32)

        call func( diff, tmp, x+cx(1)*dx,  m,n,fo)
        a3=fo*dx
        tmp=y+a1*cy(41)+a2*cy(42)+a3*cy(43)

        call func( diff, tmp, x+cx(2)*dx,  m,n,fo)
        a4=fo*dx
        tmp=y+a1*cy(51)+a3*cy(52)+a4*cy(53)

        call func( diff, tmp, x+cx(2)*dx,  m,n,fo)
        a5=fo*dx
        tmp=y+a1*cy(61)+a3*cy(62)+a4*cy(63)+a5*cy(64)

        call func( diff, tmp, x+cx(1)*dx,  m,n,fo)
        a6=fo*dx
        tmp=y+a1*cy(71)+a3*cy(72)+a4*cy(73)+a5*cy(74)+a6*cy(75)

        call func( diff, tmp, x+cx(3)*dx,  m,n,fo)
        a7=fo*dx
        tmp=y+a1*cy(81)+a5*cy(82)+a6*cy(83)+a7*cy(84)

        call func( diff, tmp, x+cx(3)*dx,  m,n,fo)
        a8=fo*dx
        tmp=y+a1*cy(91)+a5*cy(92)+a6*cy(93)+a7*cy(94)+a8*cy(95)

        call func( diff, tmp, x+cx(1)*dx,  m,n,fo)
        a9=fo*dx
        tmp=y+a1*cy(101)+a5*cy(102)+a6*cy(103)+a7*cy(104)+a8*cy(105)+a9*cy(106)

        call func( diff, tmp, x+cx(2)*dx,  m,n,fo)
        a10=fo*dx
        tmp=y+a5*cy(111)+a6*cy(112)+a7*cy(113)+a8*cy(114)+a9*cy(115)+a10*cy(116)

        call func( diff, tmp, x+dx, m,n,fo)
        a11=fo*dx

        y=y+(9.0d0*a1+49.0d0*a8+64.0d0*a9+49.0d0*a10+9.0d0*a11)/180.0d0

    end subroutine rungekutta81



    subroutine func(diff, ini_angle, x_val, m,n, output)


        ! y are the tau_val at point x
        ! This subroutine returns the magnitude of ADT angles at a specific nuclear geometry during integration along
        ! first coordinate
        ! m is the grid points and n is number of nact

        integer(8),  intent(in):: m,n
        real(8), intent(in)    :: diff(m,n), ini_angle(n), x_val
        real(8), intent(out)   :: output(n)
        real(8)                :: y_val(n)
        integer(8)             :: i

        do i=1,n
            call splint(grid, tau(:,i), diff(:,i), m, x_val, y_val(i))
        enddo
        call res(ini_angle,y_val,output, n, nstate)

    end subroutine func



    subroutine spline(x,y,n,y2)
        implicit none
        integer(8),  intent(in):: n
        real(8), intent(in)    :: x(n),y(n)

        real(8), intent(out)   :: y2(n)
        integer(8)             :: i ,k
        real(8)                :: u(n),sig, p, qn, un, yp1, ypn

        yp1 = 1.0D30
        ypn = 1.0D30

        if (yp1.gt..99e30) then
            y2(1)=0.
            u(1)=0.
        else
            y2(1)=-0.5
            u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do i=2,n-1
            sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
            p=sig*y2(i-1)+2.
            y2(i)=(sig-1.)/p
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo
        if (ypn.gt..99e30) then
            qn=0.
            un=0.
        else
            qn=0.5
            un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
        do k=n-1,1,-1
            y2(k)=y2(k)*y2(k+1)+u(k)
        enddo

    end subroutine spline

    subroutine splint(xa,ya,y2a,n,x,y)
        implicit none
        integer(8),  intent(in):: n
        real(8), intent(in)    :: xa(n),y2a(n),ya(n),x
        real(8), intent(out)   :: y
        integer(8)             :: k, khi, klo
        real(8)                :: h,a,b


        klo=1
        khi=n
        do while ((khi-klo).gt.1)
            k=(khi+klo)/2
            if(xa(k).gt.x)then
                khi=k
            else
                klo=k
            endif
        enddo
        h=xa(khi)-xa(klo)
        if (h.eq.0.) stop "bad xa input in splint"
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+ ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

    end subroutine splint




    subroutine res(ang,nact,f, ntau, nstate)

        ! This subroutine returns magnitude of ADT angles at every grid point in the nuclear configuration space (CS)

        integer(8),  intent(in)::  ntau, nstate
        real(8), intent(in) :: ang(ntau), nact(ntau)
        real(8), intent(out) :: f(ntau)
        integer(8) :: i,j,counter
        real(8) :: val(ntau),g(ntau,ntau),gi(ntau,ntau),amat(nstate,nstate),tmat(nstate,nstate),prod(nstate,nstate)
        real(8) :: afor(0:ntau,nstate,nstate),aback(ntau+1,nstate,nstate), tmp(ntau)

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

        tmp = f 
        ! f(1) = tmp(3)
        ! f(2) = tmp(2)
        ! f(3) = tmp(1)

        ! order the result in correct order
        do i=1,ntau
            do  j = 1,ntau 
                if (order(j) .eq. i) then
                    f(i) = tmp(j)
                    exit
                endif
            enddo
        enddo


    end subroutine res

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine amatspl(y,aa,afor,aback, ntau, nstate)

        ! This subroutine generates complete ADT matrix (aa) and two sets of partially multiplied ADT matrices (afor,aback)

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

        call reordermatrix(a, ntau, nstate)

        do k = 1,ntau
            aa1 = matmul(aa1,a(k,:,:))
            afor(k,:,:) = aa1
        enddo

        do k = ntau,1,-1
            aa = matmul(a(k,:,:),aa)
            aback(k,:,:) = aa
        enddo

    end subroutine amatspl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine gradcomat(y,afor,aback,g, ntau, nstate)

        ! This subroutine returns the elements of coefficient matrix of gradient of ADT angles

        integer(8),  intent(in):: ntau, nstate
        real(8), intent(in) :: afor(0:ntau,nstate,nstate),aback(ntau+1,nstate,nstate),y(ntau)
        real(8), intent(out) :: g(ntau,ntau)
        integer(8) :: i,j,k,counter
        real(8) :: adiff(ntau,nstate,nstate),aa1(nstate,nstate),aa2(nstate,nstate),b1(nstate,nstate),b2(nstate,nstate),&
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

        call reordermatrix(adiff, ntau, nstate)
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine inverse(g,gi, ntau)

        ! This subroutine returns the inverse of coefficient matrix of gradient of ADT angles by Gauss-Jordon method

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine negtau(tau,taumat, ntau, nstate)

        ! This subroutine returns the negative form of nonadiabatic coupling matrix (NACM) at each and every grid point
        ! in nuclear CS

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine amat(y,aa, ntau, nstate)

        ! This subroutine returns numerically calculated ADT matrix for a particular set of ADT angles.

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
        call reordermatrix(a, ntau, nstate)
        do i=1,ntau
            aa = matmul(aa,a(i,:,:))
        enddo

    end subroutine amat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine reordermatrix(mat, m, n)
        integer(8), intent(in):: m, n       !<<<--- m=ntau; n=nstate
        real(8),intent(inout) :: mat(m,n,n)
        real(8) ::  mattmp(m,n,n)
        integer(8):: i
        ! reorders the array of rotation matrix in a given order
        
        mattmp = mat
        do i=1,m 
            mat(i,:,:) = mattmp(order(i),:,:)
        enddo
    end subroutine reordermatrix

   end module adt
