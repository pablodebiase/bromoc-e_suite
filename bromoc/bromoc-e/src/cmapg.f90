subroutine cmapg()

use explatmod

implicit none
! local variables
integer icmap,n,n2,i,i1,j,j1
real dang,idang
real, parameter :: cte1 = 1.0/3.0, cte2 = 1.0/6.0
real, allocatable :: xp(:),yp(:),y2p(:)
real, allocatable :: xt(:),yt(:),y2t(:)
real, allocatable :: csstmp1(:),csstmp2(:)

do icmap = 1,ncmap
  ! total number of grid points in each direction
  n = cmap(icmap)
  n2 = n*n
  allocate (xp(n+1),yp(n+1),y2p(n+1))
  allocate (xt(n+1),yt(n+1),y2t(n+1))
  allocate (csstmp1(n2),csstmp2(n2))
  ! grid spacing
  dang = gscmap(icmap)
  idang = 1.0/dang
  ! *** Gradient calculation
  ! Psi angle
  do i = 1,n
    i1 = n*(i-1)
    do j = 1,n
      j1 = i1 + j
      ! Grid points
      xp(j) = -180.0 + (j-1)*dang
      ! Interpolating function at the grid points 
      yp(j) = fcmap(j1,icmap)
    end do
    xp(n+1) = 180.0
    yp(n+1) = yp(1)
    ! Second derivative of the interpolating function at the grid points
    call cmapspline(xp,yp,n+1,y2p)
    ! Psi gradient
    do j = 1,n
      j1 = i1 + j
      fpcmap(j1,icmap) = (yp(j+1)-yp(j))*idang - y2p(j)*dang*cte1 - y2p(j+1)*dang*cte2
    end do
  end do ! next i
  ! Theta angle
  do j = 1,n
    do i = 1,n
      i1 = (i-1)*n + j
      ! Grid points
      xt(i) = -180.0 + (i-1)*dang
      ! Interpolating function at the grid points
      yt(i) = fcmap(i1,icmap)
    end do
    xt(n+1) = 180.0
    yt(n+1) = yt(1)
    ! Second derivative of the interpolating function at the grid points
    call cmapspline(xt,yt,n+1,y2t)
    ! Theta gradient
    do i = 1,n 
      i1 = (i-1)*n + j
      ftcmap(i1,icmap) = (yt(i+1)-yt(i))*idang - y2t(i)*dang*cte1 - y2t(i+1)*dang*cte2
    end do
  end do ! next j
  ! *** Cross derivative calculation
  ! Cross derivative (theta-psi)
  do i = 1,n
    i1 = n*(i-1)
    do j = 1,n
      j1 = i1 + j
      ! Interpolating function at the grid points 
      yp(j) = ftcmap(j1,icmap) 
    end do
    yp(n+1) = yp(1)
    ! Second derivative of the interpolating function at the grid points
    call cmapspline(xp,yp,n+1,y2p)
    ! Cross derivative 
    do j = 1,n
      j1 = i1 + j
      csstmp1(j1) = (yp(j+1)-yp(j))*idang - y2p(j)*dang*cte1 - y2p(j+1)*dang*cte2
    end do
  end do ! next i
  ! Cross derivative (psi-theta)
  do j = 1,n
    do i = 1,n
      i1 = (i-1)*n + j
      ! Interpolating function at the grid points
      yt(i) = fpcmap(i1,icmap)
    end do
    yt(n+1) = yt(1)
    ! Second derivative of the interpolating function at the grid points
    call cmapspline(xt,yt,n+1,y2t)
    ! Cross derivative
    do i = 1,n 
      i1 = (i-1)*n + j
      csstmp1(i1) = (yt(i+1)-yt(i))*idang - y2t(i)*dang*cte1 - y2t(i+1)*dang*cte2
    end do
  end do ! next j
  ! Cross derivative (average)
  j1 = 0
  do i = 1,n
    do j = 1,n
      j1 = j1 + 1
      ftpcmap(j1,icmap) = (csstmp1(j1)+csstmp2(j1))*0.5
    end do
  end do   
  deallocate (xp,yp,y2p)
  deallocate (xt,yt,y2t)
  deallocate (csstmp1,csstmp2)
end do ! next icmap

return
end subroutine
