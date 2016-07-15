!    BROMOC  -  CG-GCMC-BD
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Fortran common block for Grand Canonical Monte Carlo
! Dimension of arrays
module listmod
implicit none

! TYPES DEFINITIONS
! Cartesian Type
type :: car
    real :: x
    real :: y
    real :: z
end type car

! Particle Properties Type
type :: parpro
    integer :: ptyp                            !! Identifies the Particle Type with an Index number
    integer :: sr                              !! Starting in the position vector
    integer :: ne                              !! Number of element in the particle
end type parpro

! Particle Type
type :: partype
    integer :: ne                              !! Number of Elements in Particle
    integer,allocatable,dimension(:) :: etyp   !! Particle Element Types Vector :: Size of ne
    type(car),allocatable,dimension(:) :: r    !! Position Vector of Particle Elements :: Size of ne
end type partype

! LIST DEFINITIONS
!! Particles List
integer                                 :: npar       !! Number of Particles 
type(parpro), allocatable, dimension(:) :: parl        !! Particles List Description :: Size of Number of Particles (npar)

!! Particle-Type List
integer                                  :: ntyp      !! Number of Types
type(partype), allocatable, dimension(:) :: ptypl    !! Type of Particle list :: Size of ntyp

!! Force and position vectors
integer nele                                !! Total Number of elements
type(car), allocatable, dimension(:) :: r  !! Elements Position Vector
type(car), allocatable, dimension(:) :: f  !! Elements Force Vector

contains

subroutine resizeptypl(newsize)
implicit none
integer vdim,newsize
type(partype), allocatable, dimension(:) :: ptypltmp    !! Type of Particle list :: Size of ntyp
if (newsize .lt. 1) newsize = 1
! Measure Vector Size
vdim=size(ptypl)
! Allocate double of original size in temp vectors
allocate (ptypltmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(ptypl) .and. vdim .gt. 0) ptypltmp(1:vdim)=ptypl(1:vdim)
! Deallocate original Vector
if (allocated(ptypl)) deallocate (ptypl)
! Move allocations
call move_alloc(ptypltmp,ptypl)
end subroutine

subroutine resizeparl(newsize)
implicit none
integer vdim,newsize
type(parpro), allocatable, dimension(:) :: parltmp        !! Particles List Description :: Size of Number of Particles (npar)
if (newsize .lt. 1) newsize = 1
! Measure Vector Size
vdim=size(parl)
! Allocate double of original size in temp vectors
allocate (parltmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(parl) .and. vdim .gt. 0) parltmp(1:vdim)=parl(1:vdim)
! Deallocate original Vector
if (allocated(parl)) deallocate (parl)
! Move allocations
call move_alloc(parltmp,parl)
end subroutine

subroutine resizecvec(newsize)
implicit none
integer vdim,newsize
type(car), allocatable, dimension(:) :: rtmp  !! Elements Position Vector
type(car), allocatable, dimension(:) :: ftmp  !! Elements Force Vector
if (newsize .lt. 1) newsize = 1
! Measure Vector Size
vdim=size(r)
! Allocate double of original size in temp vectors
allocate (rtmp(newsize),ftmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(r) .and. vdim .gt. 0) rtmp(1:vdim)=r(1:vdim)
if (allocated(f) .and. vdim .gt. 0) ftmp(1:vdim)=f(1:vdim)
! Deallocate original Vector
if (allocated(r)) deallocate (r)
if (allocated(f)) deallocate (f)
! Move allocations
call move_alloc(rtmp,r)
call move_alloc(ftmp,f)
end subroutine

subroutine addpar(ptype)
implicit none
integer ptype,lastnele,i
if (ptype .gt. 0 .and. ptype .le. ntyp) then
  npar=npar+1
  lastnele=nele
  nele=nele+ptypl(ptype)%ne
  if (nele .gt. size(r)) call resizecvec(nele)
  if (npar .gt. size(parl)) call resizeparl(npar)
  parl(npar)%ptyp=ptype     ! Save Particle Type Idx in Particle List
  parl(npar)%sr=lastnele    ! Save Previous Position of last Added Element in Cartesian Vectors to Particle List
  parl(npar)%ne=ptypl(ptype)%ne  ! Save Number of Elements of this particle in the Particle List
  ! Copy Coordinates from Particle Type Template 
  do i=1,parl(npar)%ne
    r(lastnele+i)=ptypl(ptype)%r(i)
  enddo
endif
end subroutine

subroutine movepar(parn,rcent)
implicit none
type(car) :: rcent,arcent  ! Particle Centroid
integer ne,sr,parn ! Particle Number 
real ine
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Compute Actual Centroid
call getcentroid(arcent,parn)
! Compute Displacement: Substract Actual Centroid to New Centroid
call subcar2par(parn,arcent)
! Add Displacement to all elements position
call addcar2par(parn,rcent)
end subroutine

! Randomly rotates a particle
subroutine uranrot(parn)
implicit none
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510
real, parameter :: twopi=2.0*pi
integer parn
real,external :: rndm
real phi, theta, psi, cosp, sinp, ocosp
real x,y,z,rot(3,3)
type(car) :: cent
if (parl(parn)%ne .lt. 2) return
theta=acos(2.0*rndm()-1.0) ! from 0 to pi
phi=twopi*rndm()           ! from 0 to 2*pi
psi=twopi*rndm()           ! from 0 to 2*pi
! Make Random Vector
sinp=sin(theta)
x=sinp*cos(phi)
y=sinp*sin(phi)
z=cos(theta)
! Build Rotation Matrix
cosp=cos(psi)
ocosp=1.0-cosp
sinp=sin(psi)
! First Row
rot(1,1)=cosp+x*x*ocosp
rot(1,2)=x*y*ocosp-z*sinp
rot(1,3)=x*z*ocosp+y*sinp
! Second Row
rot(2,1)=x*y*ocosp+z*sinp
rot(2,2)=cosp+y*y*ocosp
rot(2,3)=y*z*ocosp-x*sinp
! Third Row
rot(3,1)=x*z*ocosp-y*sinp
rot(3,2)=y*z*ocosp+x*sinp
rot(3,3)=cosp+z*z*ocosp
! Get Centroid
call getcentroid(cent,parn)
! Remove Centroid
call subcar2par(parn,cent)
! Rotate Particle
call rotatepar(parn,rot)
! Restore Centroid
call addcar2par(parn,cent)
end subroutine

! Given a rotation matrix, rotates particle
subroutine rotatepar(parn,rot)
implicit none
integer i,ne,sr,parn ! Particle Number 
real ine,rot(3,3)
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call rotatecar(r(sr+i),rot)
enddo
end subroutine

subroutine rotatecar(rc,mat)
implicit none
type(car) :: rc  ! Particle Centroid
real mat(3,3),x,y,z
x=mat(1,1)*rc%x+mat(1,2)*rc%y+mat(1,3)*rc%z
y=mat(2,1)*rc%x+mat(2,2)*rc%y+mat(2,3)*rc%z
z=mat(3,1)*rc%x+mat(3,2)*rc%y+mat(3,3)*rc%z
rc%x=x
rc%y=y
rc%z=z
end subroutine

subroutine getcentroid(rc,parn)
implicit none
type(car) :: rc  ! Particle Centroid
integer i,ne,sr,parn ! Particle Number 
real ine
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Compute Actual Centroid
! Zero vector
call zero(rc)
! Add each position vector to rc
do i=1,ne
   call addcar(rc,r(sr+i))
enddo
! divide rc by number of elements
call mulcar(rc,ine)
end subroutine

subroutine addcar2par(parn,rc)
implicit none
type(car) :: rc  ! Cartesian coordinates
integer i,ne,sr,parn ! Particle Number 
real ine
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call addcar(r(sr+i),rc)
enddo
end subroutine 

subroutine subcar2par(parn,rc)
implicit none
type(car) :: rc  ! Cartesian coordinates
integer i,ne,sr,parn ! Particle Number 
real ine
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call subcar(r(sr+i),rc)
enddo
end subroutine

function getcar(x,y,z)
implicit none
type(car) :: getcar
real :: x,y,z
getcar%x=x
getcar%y=y
getcar%z=z
end function

subroutine setcar(r,x,y,z)
implicit none
type(car) :: r
real :: x,y,z
r%x=x
r%y=y
r%z=z
end subroutine

subroutine zero(r)
implicit none
type(car) :: r
call setcar(r,0.0,0.0,0.0)
end subroutine

! Add Second Cartesian Type to first
subroutine addcar(r1,r2)
implicit none
type(car) :: r1,r2
r1%x=r1%x+r2%x
r1%y=r1%y+r2%y
r1%z=r1%z+r2%z
end subroutine

! Subtract Cartesian Type
subroutine subcar(r1,r2)
implicit none
type(car) :: r1,r2
r1%x=r1%x-r2%x
r1%y=r1%y-r2%y
r1%z=r1%z-r2%z
end subroutine

! Dot Product of Cartesian Types
function dotcar(r1,r2)
implicit none
real dotcar
type(car) :: r1,r2
dotcar=r1%x*r2%x+r1%y*r2%y+r1%z*r2%z
end function

function dotcarvec(r,v)
implicit none
real dotcarvec,v(3)
type(car) :: r
dotcarvec=v(1)*r%x+v(2)*r%y+v(3)*r%z
end function

! Scalar Multiplication of Cartesian Types
subroutine mulcar(r,sc)
implicit none
type(car) :: r
real sc
r%x=sc*r%x
r%y=sc*r%y
r%z=sc*r%z
end subroutine

subroutine delpar(parn)
implicit none
integer ptype,i,n,ne,sr,parn
! If parn is the last particle just remove
if (parn .eq. npar) then
    if (npar .eq. 1) then
        nele=0
        npar=0
    else
        nele=nele-parl(parn)%ne
        npar=npar-1
    endif 
    return
endif
! Look for last particle of same type
ptype=parl(parn)%ptyp
ne=parl(parn)%ne
n=npar
do while (parl(n)%ptyp .ne. ptype .and. n.gt.parn)
    n=n-1
enddo
! If there is another of the same type closer to the end copy that one to parn's place
if (n .gt. parn) r(parl(parn)%sr+1:parl(parn)%sr+parl(parn)%ne)=r(parl(n)%sr+1:parl(n)%sr+parl(n)%ne)
! Move all coordinates backwards and remove that n particle from the coordinate vector
sr=parl(n)%sr
do i=sr+1,nele-ne
    r(i)=r(i+ne)
enddo
nele=nele-ne
! Move all particles backwards and remove one particle from particle list
do i=n,npar-1
   parl(i)%ptyp=parl(i+1)%ptyp
   parl(i)%sr=parl(i+1)%sr-ne
   parl(i)%ne=parl(i+1)%ne
enddo
npar=npar-1
end subroutine

end module
