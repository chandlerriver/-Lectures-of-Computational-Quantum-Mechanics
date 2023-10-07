module Comphy_ODE
implicit none 
contains 

subroutine ODEBasicEuler(xn, yn, h)
    real*8,intent(in) :: h 
    real*8,intent(inout) :: xn, yn
    real*8 :: fn  

    call getODEinfo(xn, yn, fn)

    yn = yn + h * fn 
    xn = xn + h 
end subroutine 

end module 