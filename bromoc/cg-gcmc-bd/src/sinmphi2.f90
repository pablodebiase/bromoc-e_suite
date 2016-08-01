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

subroutine sinmphi2(j,m,cc0,s0,as)
!------------------------------------------------------------------------
!sin(M*phi) calculation (M > 0)
use grandmod
implicit none

integer   m,i,j
real    cc0,s0,as(nele,0:m)

as(j,0) = 0.0
as(j,1) = s0
as(j,2) = 2.0*cc0*s0
as(j,3) = 2.0*cc0*as(j,2)-s0
do i = 4, m
  as(j,i) = 2.0*cc0*as(j,i-1)-as(j,i-2)
enddo
 
return
end
