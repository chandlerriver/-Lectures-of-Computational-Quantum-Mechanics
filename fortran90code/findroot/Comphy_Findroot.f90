module Comphy_Findroot
implicit none 
contains

subroutine Bisection(func, niter, x)
    implicit none 
    real*8, intent(inout) :: x
    integer,intent(in) :: niter 
    integer :: iter
    real*8 :: xup, xlow, fx, flow, fup
    interface
        subroutine func(x, fx, dfx, gx, dgx)
            real*8,intent(in) :: x
            real*8,intent(out) :: fx
            real*8,intent(out),optional :: dfx, gx, dgx 
        end subroutine 
    end interface 

    xlow = x - 1.5d0
    xup = x + 1.5d0 
    
    call func(xlow, flow);call func(xup, fup)
    
    if (flow*fup>0.0d0) then 
        print *, "There is no root in the interval!";return 
    end if 

    print "(a)", "iter   x   f(x)   |dx|"
    do iter = 1,niter 
        x = (xlow+xup)/2.0d0
        call func(x, fx)
        print "(i5, 3f10.3)", iter, x, fx, abs(xup-xlow)
        if (fx*flow < 0.0d0) then 
            xup = x
            fup = fx
        else
            xlow = x
            flow = fx 
        end if 
    end do 
end subroutine



end module 

