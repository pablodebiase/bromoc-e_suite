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
! Position Type
type :: pos
    real :: x
    real :: y
    real :: z
end type pos
! Particle Properties
type :: parpro
    integer :: typidx                          !! Identifies the Particle Type with an Index number
    integer :: ne                              !! Number of element components of the particle
    type(pos), dimension(:), pointer :: r      !! Pointer to Position Vector of each element
    integer, dimension(:), pointer :: eleidx   !! Pointer to Elements Type Index Vector
end type parpro

integer                   :: npar
type(parpro), allocatable :: par(:)
type(pos), allocatable, target :: rt(:)

! allocate rt(parpro%ne)
! parpro => rt(
end module
