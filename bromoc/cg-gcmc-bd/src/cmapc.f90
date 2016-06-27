subroutine cmapc()

use explatmod

implicit none
! local variables
integer icmap,n,n2,i,j,j1,j2,j3,j4,k,l,m
real y(4),y1(4),y2(4),y12(4),d,c(16)

do icmap = 1,ncmap
  n = cmap(icmap)
  n2 = n*n
  j1 = 0
  do i = 1,n
    do j = 1,n
      ! Grid points for a rectangular grid cell (numbered counter clockwise from the lower left)
      j1 = j1 + 1
      j2 = j1 + n
      j3 = j2 + 1
      j4 = j1 + 1
      ! Symmetry considerations for grid points
      if (j.eq.n) then
        j3 = j3 - n
        j4 = j4 - n
      end if
      if (i.eq.n) then
        j2 = j2 - n2
        j3 = j3 - n2
      end if
      ! Function at the four grid points
      y(1) = fcmap(j1,icmap)
      y(2) = fcmap(j2,icmap)   
      y(3) = fcmap(j3,icmap)
      y(4) = fcmap(j4,icmap) 
      ! Gradients at the four grid points 
      y1(1) = ftcmap(j1,icmap)
      y1(2) = ftcmap(j2,icmap)
      y1(3) = ftcmap(j3,icmap)
      y1(4) = ftcmap(j4,icmap)
      y2(1) = fpcmap(j1,icmap)
      y2(2) = fpcmap(j2,icmap)
      y2(3) = fpcmap(j3,icmap)
      y2(4) = fpcmap(j4,icmap)
      ! Cross derivative at the four grid points
      y12(1) = ftpcmap(j1,icmap)
      y12(2) = ftpcmap(j2,icmap)
      y12(3) = ftpcmap(j3,icmap)
      y12(4) = ftpcmap(j4,icmap)
      ! Obtain the CMAP coefficient for the bicubic interpolation
      d = gscmap(icmap)
      call bcucof(y,y1,y2,y12,d,c)
      m = 0
      do k = 1,4
        do l = 1,4
          m = m + 1
          ccoef(m,j1,icmap) = c(m)
        end do
      end do
    end do ! next j
  end do ! next i
end do ! next icmap

return
end subroutine
