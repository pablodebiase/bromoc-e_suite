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

subroutine dynamics1()
use grandmod
use listmod
use constamod

implicit none
integer itype,i
real fact1, fact2
real rgauss
external rgauss
real delx, dely, delz, delDz
real sw, dsw, idiffusion, didiffusion
real zz, pp3, pore1, pore2,ipp3

!$omp parallel private(i,itype,delx,dely,delz,zz,pp3,pore1,pore2,sw,dsw,ipp3,idiffusion,didiffusion,fact1,fact2,delDz)
!$omp do
do i = nelenuc+1, nele
  itype=et(i)
! space-dependent diffusion constant
  zz = r(i)%z
  pp3 = p3(itype)
  pore1 = -plength2 + pcenter
  pore2 =  plength2 + pcenter
  if (zz.gt.(pore1-pp3) .and. zz.lt.(pore2+pp3)) then
    if (zz.ge.pore1 .and. zz.le.pore2) then
      sw  = 0.0
      dsw = 0.0
    else
      ipp3=1.0/pp3
      if (zz.gt.pore2) then
        delz = (zz-pore2)
        sw = 2.0*(delz*ipp3)**3-3.0*(delz*ipp3)**2
        dsw = 6.0*((delz*ipp3)**2-delz*ipp3)*ipp3
      else if (zz.lt.pore1) then
        delz = (zz-pore1)
        sw = -2.0*(delz*ipp3)**3-3.0*(delz*ipp3)**2
        dsw = -6.0*((delz*ipp3)**2+delz*ipp3)*ipp3
      endif
      sw = -sw
      dsw = -dsw
    endif
  else
    sw = 1.0
    dsw = 0.0
  endif
  idiffusion = etypl(itype)%dif*(ampl3(itype)+(1.0-ampl3(itype))*sw)
  didiffusion = etypl(itype)%dif*(1.0-ampl3(itype))*dsw
  fact1 = idiffusion*kBTdt
  delx  = f(i)%x*fact1
  dely  = f(i)%y*fact1
  delz  = f(i)%z*fact1
  fact2 = sqrt(2.0*dt*idiffusion)
  delDz = didiffusion*dt

  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)

  r(i)%x = r(i)%x + delx + fact2*rgauss()
  r(i)%y = r(i)%y + dely + fact2*rgauss()
  r(i)%z = r(i)%z + delz + fact2*rgauss() + delDz

enddo
!$omp end do
!$omp end parallel
end subroutine
