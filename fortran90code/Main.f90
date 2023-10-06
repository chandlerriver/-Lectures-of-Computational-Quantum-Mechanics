program main 
    use Comphy_Findroot
    implicit none 
    real*8 ::x
    integer :: imethod, niter
    external :: myfunc

    write(*, "(/, 'Input initial x: ', $)");read(*, *) x
    write(*, "(/, 'Input iteration number: ', $)");read(*, *) niter 
    do 
        print *
        print "(a)", "All methods"
        print *, "1: Bisection method"

        write(*, "(/, 'Select method: ', $)");read(*, "(i8)") imethod
        print *

        select case(imethod)
        case(1)
            call Bisection(myfunc, niter, x)

        case default
            exit 
        end select 
    end do 

end program 

subroutine myfunc(x, fx, dfx, gx, dgx)
    implicit none 
    real*8,intent(in) :: x
    real*8,intent(out) :: fx 
    real*8,intent(out),optional :: dfx, gx, dgx 

    fx = x**2 - 2.0d0 * x - 4.0
end subroutine

