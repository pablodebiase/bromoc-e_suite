subroutine cmapspline(x,y,n,y2)

implicit none
! INPUT:
! x(1)<x(2)<...<x(n) -> grid points
! y(1:n)             -> tabulated function
! OUPUT:
! y2(1:n)            -> second derivative of the interpolating function at the grid points
integer n
real x(n),y(n),y2(n)
! local variables
integer i
real p,sig,u(n-1)

! The lower boundary condition is set to be "natural"
y2(1) = 0.0
u(1) = 0.0
! Descomposition loop of the tridiagonal algortihm. y2 and u are used 
! for temporary storage of the descomposed factors
do i = 2,n-1
  sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
  p = sig*y2(i-1) + 2.0
  y2(i) = (sig-1.0)/p
  u(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) &
       - sig*u(i-1))/p 
end do
! The upper boundary condition is set to be "natural
y2(n) = 0.0
! Backsubstitution loop of the tridiagonal algorithm
do i = n-1,1,-1
  y2(i) = y2(i)*y2(i+1) + u(i)
end do

return
end subroutine
