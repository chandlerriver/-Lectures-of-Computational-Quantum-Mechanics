program main 
    use Comphy_ODE
    implicit none 
    real*8 :: h
    integer :: nstep, istep, imethod 
    real*8 :: xn, yn, fn, yreal, error, inxin, inyin
    real*8 :: preyn, prefn, xn2, yn2, preyn2, prefn2, tpy 
    real*8 :: mfn(-3:1),mfn2(-3:1)
    interface
        subroutine getODEinfo(x, y, f)
            real*8,intent(in) :: x 
            real*8,intent(inout) :: y 
            real*8,intent(out),optional :: f 
        end subroutine
    end interface

    print "(a)", "Ordinary Differential Equation dy/dx=exp(x)+y"
    write(*, "(/, 'Input step size: ', $)");read(*, *)h
    write(*, "(/, 'Input step number: ', $)");read(*, *)nstep
    inxin = 0.0d0;inyin = 0.0d0; error = 0.0d0 

    do 
        print *
        print "(a)", "All methods"
        print *, "1:Basic Euler Method"
        print *, "2:Improved Euler Method"
        print *, "3:Leap Frog Method"
        print *, "4:Simpson Method"
        print *, "5:Runge-Kutta Method(four order)"
        print *, "6:Multistep Method(four step)"
        write(*, "(/, 'Selected method:', $)");read(*, "(i8)")imethod 
        print *
        print "(a)", "    x        y        yn        y-yn        posteriori"

        xn = inxin;yn = inyin;error = 0.0d0 
    do istep=1, nstep
        xn2 = xn;yn2 = yn 

        select case(imethod)
        
        case(1)
            call ODEBasicEuler(xn, yn, h)
            call ODEBasicEuler(xn2, yn2, 0.5d0*h)
            call ODEBasicEuler(xn2, yn2, 0.5d0*h)

            error = error + 2.0d0 * (yn2 - yn)

        case default
            exit 
        end select 

        call getODEinfo(xn, yreal)

        if (mod(istep, 20)==0) print "(3f8.3, 2e11.1e2)", xn, yreal, yn, yreal-yn, error
    end do 
    end do 

end program 

subroutine getODEinfo(x, y, f)
    real*8,intent(in) :: x 
    real*8,intent(inout) :: y 
    real*8,intent(out),optional :: f 
    if (.not. present(f)) then 
        y = x * exp(p)
    else 
        f = exp(x) + y 
    end if 
end subroutine