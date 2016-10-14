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

subroutine dynamics0()
use grandmod
use listmod
use constamod

implicit none
integer itype, i
real,external :: rgauss
real delx, dely, delz
!$omp parallel private(i,itype,delx,dely,delz)
!$omp do
do i = nelenuc+1, nele
  itype = et(i)
  delx = f(i)%x*fact1a(itype)
  dely = f(i)%y*fact1a(itype)
  delz = f(i)%z*fact1a(itype)
  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)
  r(i)%x = r(i)%x + delx + fact2a(itype)*rgauss()
  r(i)%y = r(i)%y + dely + fact2a(itype)*rgauss()
  r(i)%z = r(i)%z + delz + fact2a(itype)*rgauss()
enddo
!$omp end do
!$omp end parallel
end subroutine
