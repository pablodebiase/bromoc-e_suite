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

subroutine dynamics3()
use grandmod
use listmod
use constamod
use sevalmod

implicit none
integer itype,i
real fact1, fact2
real rgauss
external rgauss
real delx, dely, delz, delDx, delDy, delDz
real idiffusion, didiffusion
real ebet,idist

!$omp parallel private(i,itype,ebet,idiffusion,didiffusion,idist,fact1,fact2,delDx,delDy,delDz,delx,dely,delz)
!$omp do
do i=nelenuc+1, nele
  itype = et(i)
  if (dids(1,i).gt.0.0.and.dids(1,i).le.diffcutoff) then
! dna proximity-dependent diffusion constant
    ebet=exp(-dids(1,i)*ibeta)
    idiffusion = etypl(itype)%dif-(etypl(itype)%dif-diff0)*ebet
    didiffusion = (etypl(itype)%dif-diff0)*ibeta*ebet
    idist=1.0/dids(2,i)
    fact1 = idiffusion*kBTdt
    fact2 = sqrt(2.0*dt*idiffusion)
    delDx = didiffusion*dt*dids(3,i)*idist
    delDy = didiffusion*dt*dids(4,i)*idist
    delDz = didiffusion*dt*dids(5,i)*idist
  else
    fact1=fact1a(itype)
    fact2=fact2a(itype)
!    fact1=diff0*kBTdt
!    fact2=fact2pd
    delDx=0.0
    delDy=0.0
    delDz=0.0
  endif
  delx = f(i)%x*fact1
  dely = f(i)%y*fact1
  delz = f(i)%z*fact1
  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)
  r(i)%x=r(i)%x+delx+fact2*rgauss()+delDx
  r(i)%y=r(i)%y+dely+fact2*rgauss()+delDy
  r(i)%z=r(i)%z+delz+fact2*rgauss()+delDz
enddo
!$omp end do
!$omp end parallel
end subroutine
