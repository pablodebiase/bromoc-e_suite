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

subroutine dynamics2()
use grandmod
use listmod
use constamod
use splinemod
use sevalmod

implicit none
integer itype,i
real fact1, fact2
real rgauss
external rgauss
real delx, dely, delz, delDz
real sw,dsw, idiffusion, didiffusion
real zz

do i=nelenuc+1, nele
  itype = et(i)
! space-dependent diffusion constant
  zz = r(i)%z
  if(xs(1).lt.xs(nspline)) then 
    if(zz.lt.xs(1)) then
      sw=ys(1)
      dsw=0.0
    elseif(zz.gt.xs(nspline)) then
      sw=ys(nspline)
      dsw=0.0
    else
      sw=seval(nspline,zz,xs,ys,b,c,d)
      dsw=sevald(nspline,zz,xs,b,c,d)
    endif
  elseif(xs(1).gt.xs(nspline)) then
    if(zz.gt.xs(1)) then
      sw=ys(1)
      dsw=0.0
    elseif(zz.lt.xs(nspline)) then
      sw=ys(nspline)
      dsw=0.0
    else
      sw=seval(nspline,zz,xs,ys,b,c,d)
      dsw=sevald(nspline,zz,xs,b,c,d)
    endif
  endif
!   sw  = seval (nspline,zz,xs,ys,b,c,d)
!   dsw = sevald(nspline,zz,xs,b,c,d)
  idiffusion = etypl(itype)%dif*sw
  didiffusion = etypl(itype)%dif*dsw
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
return
end
