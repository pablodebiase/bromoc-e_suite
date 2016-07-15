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
use typesmod
implicit none

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

subroutine increaseptypl(newsize)
implicit none
integer vdim,newsize
type(partype), allocatable, dimension(:) :: ptypltmp    !! Type of Particle list :: Size of ntyp

! Measure Vector Size
vdim=size(ptypl)

! Allocate double of original size in temp vectors
allocate (ptypltmp(newsize))
! Copy Vectors
ptypltmp(1:vdim)=ptypl
! Deallocate original Vector
deallocate (ptypl)
! Move allocations
call move_alloc(ptypltmp,ptypl)
end subroutine

subroutine increaseparl(newsize)
implicit none
integer vdim,newsize
type(parpro), allocatable, dimension(:) :: parltmp        !! Particles List Description :: Size of Number of Particles (npar)

! Measure Vector Size
vdim=size(parl)

! Allocate double of original size in temp vectors
allocate (parltmp(newsize))
! Copy Vectors
parltmp(1:vdim)=parl
! Deallocate original Vector
deallocate (parl)
! Move allocations
call move_alloc(parltmp,parl)

end subroutine

subroutine increasecvec(newsize)
implicit none
integer cvdim,newsize
type(car), allocatable, dimension(:) :: rtmp  !! Elements Position Vector
type(car), allocatable, dimension(:) :: ftmp  !! Elements Force Vector

! Measure Vector Size
cvdim=size(r)

! Allocate double of original size in temp vectors
allocate (rtmp(newsize),ftmp(newsize))
! Copy Vectors
rtmp(1:cvdim)=r
ftmp(1:cvdim)=f
! Deallocate original Vector
deallocate (r,f)
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
  if (nele .gt. size(r)) call increasecvec(nele)
  if (npar .gt. size(parl)) call increaseparl(npar)
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
integer i,parn ! Particle Number 
! Compute Actual Centroid
   call getcentroid(arcent,parn)
! Compute Displacement: Substract Actual Centroid to New Centroid
   ! to do
! Add Displacement to all elements position
   ! to do
end subroutine

subroutine getcentroid(rc,parn)
implicit none
type(car) :: rc  ! Particle Centroid
integer i,ne,sr,parn ! Particle Number 
real ine
! Compute Actual Centroid
   call zero(rc)
   ne=parl(parn)%ne
   ine=1.0/ne
   sr=parl(parn)%sr
   do i=1,ne
      call addcar(rc,r(sr+i))
   enddo
   mulcar(rc,ine)
end subroutine

subroutine zero(r)
implicit none
type(car) :: r
r%x=0.0
r%y=0.0
r%z=0.0
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
real dotcar
dotcar=r1%x*r2%x+r1%y*r2%y+r1%z*r2%z
end function

! Scalar Multiplication of Cartesian Types
subroutine mulcar(r,sc)
type(car) :: r
real sc
r%x=sc*r%x
r%y=sc*r%y
r%z=sc*r%z
end subroutine

subroutine delpar(parn)
implicit none
integer ptype
! Look for 
npar=npar-1
end subroutine

end module
