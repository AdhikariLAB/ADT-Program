module expand
    implicit none
    contains
    subroutine expandGrid(x,y)
        real(kind=8), intent(in) :: x(:)
        real(kind=8), intent(out) :: y(:)
        integer :: n 
        n = size(x)
        if((size(x)+2)/=size(y)) stop "invalid array dimension"
        y(2:n+1) = x 

        ! extend the grid
        y(1)   = x(1)- (x(2)-x(1))
        y(n+2) = x(n)+ (x(n) - x(n-1))
    end

    subroutine expandGrid3d(x,y) ! expansion infirst two axis
        real(kind=8), intent(in) :: x(:,:,:)
        real(kind=8), intent(out) :: y(:,:,:)
        integer :: n1,n2
        n1 = size(x,1); n2 = size(x,2)

        if((n1+2/=size(y,1)).or.(n2+2/=size(y,2))) stop "invalid array dimension"
        y(2:n1+1,2:n2+1,:) = x
        y(2:n1+1,1,:) = y(2:n1+1,2,:)
        y(2:n1+1,n2+2,:) = y(2:n1+1,n2+1,:)
        y(1,:,:) = y(2,:,:)
        y(n1+2,:,:) = y(n1+1,:,:)
    end
end




program name
    use expand
    implicit none
    real(8) :: x(10), y(12)
    x = [1,2,3,4,5,6,7,8,9,10]
    call expandGrid(x,y)
    print *, all([1,2]==[1,3])
    print *, shape([1,2,3,4])
end program name