!    BROMOC-E
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

module efpmod
implicit none

type :: coef3
    real :: a,b,c
end type coef3

type :: efpot
    real                                 :: xl           !! square of x lower limit of the potential after the head (ex dmi)
    real                                 :: xl2,xu2      !! square of x upper limit of the potential before the tail(ex dm2)
    integer                              :: n            !! number of points for the cf
    type(coef3)                          :: sc           !! Head and Tail Coefficients
    type(coef3),allocatable,dimension(:) :: ep           !! Effective Potential Square Coefficients
end type efpot

type(efpot),allocatable,dimension(:)     :: efp          !! effective potentials
real                                     :: res,ires     !! Resolution and inverse resolution
real,allocatable,dimension(:)            :: fct          !! Coulombic Factor
logical*1,allocatable :: Qcol(:)
end module

