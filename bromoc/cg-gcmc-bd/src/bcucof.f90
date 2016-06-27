subroutine bcucof(y,y1,y2,y12,d,c)

use explatmod 

implicit none
! INPUT:
! d           -> length of the grid cell 
! y,y1,y2,y12 -> function, gradients, and cross derivatives at the four grid points
!                of a rectangular grid cell (numbered counterclockwise from the lower left)
! OUPUT:
! c           -> coefficients for the bicubic interpolation
real d,c(16),y(4),y1(4),y2(4),y12(4)
! local variables
integer i,j
real d2,xx,x(16)

d2 = d*d
! Pack a temporary vector x
do i = 1,4
  x(i) = y(i)
  x(i+4) = y1(i)*d
  x(i+8) = y2(i)*d
  x(i+12) = y12(i)*d2
end do
! Matrix multiply by the stored table
do i = 1,16
  xx = 0.0
  do j = 1,16
    xx = xx + wt(i,j)*x(j)
  end do
  c(i) = xx
end do

return
end subroutine
