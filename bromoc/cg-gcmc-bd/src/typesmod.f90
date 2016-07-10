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
module typesmod
implicit none

! TYPES DEFINITIONS
! Position Type
type :: pos
    real :: x
    real :: y
    real :: z
end type pos

! Particle Properties
type :: parpro
    integer :: ptyp                            !! Identifies the Particle Type with an Index number
    integer :: ne                              !! Number of element in the particle
    integer :: sr                              !! Starting in the position vector
end type parpro

! Particle Types
type :: partype
    integer :: ne                              !! Number of Elements in Particle
    integer,allocatable,dimension(:) :: etyp   !! Particle Element Types Vector :: Size of ne
    type(pos),allocatable,dimension(:) :: r    !! Position Vector of Particle Elements :: Size of ne
end type partype

! Particle Type Vectors
type :: parvec
    integer,allocatable,dimension(:)   :: etyp     !! Element Types Vector :: Size of Estimated Maximum number of particles of each Type
    type(pos),allocatable,dimension(:) :: r        !! Position Vector :: Size of Estimated Maximum number of particles of each Type
    type(pos),allocatable,dimension(:) :: f        !! Force Vector :: Size of Estimated Maximum number of particles of each Type
end type parvec


!! Particles
integer                                 :: npar       !! Number of Particles 
type(parpro), allocatable, dimension(:) :: par        !! Particles List Description :: Size of Number of Particles (npar)
type(parvec), allocatable, dimension(:) :: pvec       !! Particle Vectors :: Size of Number of Particles Types (ntyp)

!! Types
integer                                  :: ntyp      !! Number of Types
type(partype), allocatable, dimension(:) :: partyp    !! Type of Particle list :: Size of ntyp


contains

subroutine addpar()
implicit none
end subroutine

subroutine delpar()
implicit none
npar=npar-1
end subroutine

subroutine addptype(nelem)
implicit none
integer nelem
type(partype) ptyp
type(partype), allocatable, dimension(:) :: partyptmp

! Create New Particle Type
ptyp%ne = nelem
allocate (ptyp%etyp(nelem),ptyp%r(nelem))

! Add new Particle Type to partyp
ntyp=ntyp+1
! Increase vector is needed
if (size(partyp) < ntyp) then
   allocate(partyptmp(ntyp))
   partyptmp(1:size(partyp)) = partyp
   deallocate(partyp)
   call move_alloc(partyptmp, partyp)
end if

partyp(ntyp) = ptyp
end subroutine

end module
