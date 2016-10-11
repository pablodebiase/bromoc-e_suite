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

subroutine fixcoor()
use grandmod
use listmod
implicit none
integer i,j,sr,ne
real  dx,dy,dz
real  dxmax,dymax,dzmax
real  dxmin,dymin,dzmin
real  r1,r2,lzrad
type(car) :: rr

do i=1+nparnuc,npar
  sr=parl(i)%sr
  ne=parl(i)%ne
  dxmax=0.0
  dymax=0.0
  dzmax=0.0
  dxmin=0.0
  dymin=0.0
  dzmin=0.0
  do j=sr+1,ne+sr
    ! spherical system
    if (Qsphere) then
      r2 = (r(j)%x-cx)**2+(r(j)%y-cy)**2+(r(j)%z-cz)**2
      if (r2.gt.Rsphe2) then
        lzrad=2.0*rsphe/sqrt(r2)-1.0
        r(j)%x = lzrad*(r(j)%x-cx)+cx
        r(j)%y = lzrad*(r(j)%y-cy)+cy
        r(j)%z = lzrad*(r(j)%z-cz)+cz
      endif
    ! cyllindrical system
    elseif (Qecyl) then
      r2=((r(j)%x-cx)*iecx)**2+((r(j)%y-cy)*iecy)**2
      if (r2.gt.1.0) then
         r1=0.999999999999999/sqrt(r2)
         r(j)%x = r1*(r(j)%x-cx)+cx
         r(j)%y = r1*(r(j)%y-cy)+cy
      endif
      if (r(j)%z.lt.lz2m) then
        r(j)%z = lz2m
      else if (r(j)%z.gt.lz2p) then
        r(j)%z = lz2p
      endif
    ! rectangular system
    else
      if (ne.eq.1) then
        if (r(j)%x.lt.lx2m) then
          r(j)%x = lx2m
        else if (r(j)%x.gt.lx2p) then !ge
          r(j)%x = lx2p
        endif
        if (r(j)%y.lt.ly2m) then
          r(j)%y = ly2m
        else if (r(j)%y.gt.ly2p) then !ge
          r(j)%y = ly2p
        endif
        if (r(j)%z.lt.lz2m) then
          r(j)%z = lz2m
        else if (r(j)%z.gt.lz2p) then
          r(j)%z = lz2p
        endif
      else 
        if (r(j)%x.lt.lx2m) then
          dx = lx2m - r(j)%x
          dxmax = max(dxmax,dx)
        else if (r(j)%x.gt.lx2p) then !ge
          dx = lx2p - r(j)%x
          dxmin = min(dxmin,dx)
        endif
        if (r(j)%y.lt.ly2m) then
          dy = ly2m - r(j)%y
          dymax = max(dymax,dy)
        else if (r(j)%y.gt.ly2p) then !ge
          dy = ly2p - r(j)%y
          dymin = min(dymin,dy)
        endif
        if (r(j)%z.lt.lz2m) then
          dz = lz2m - r(j)%z
          dzmax = max(dzmax,dz)
        else if (r(j)%z.gt.lz2p) then
          dz = lz2p - r(j)%z
          dzmin = min(dzmin,dz)
        endif
      endif
    endif
  enddo
  if (.not.Qsphere.and..not.Qecyl.and.ne.gt.1) then
    if (abs(dxmax).gt.abs(dxmin)) then
      rr%x=dxmax
    else
      rr%x=dxmin
    endif
    if (abs(dymax).gt.abs(dymin)) then
      rr%y=dymax
    else
      rr%y=dymin
    endif
    if (abs(dzmax).gt.abs(dzmin)) then
      rr%z=dzmax
    else
      rr%z=dzmin
    endif
    if (rr%x.ne.0.0.or.rr%y.ne.0.0.or.rr%z.ne.0.0) call addcar2par(i,rr)
  endif
enddo
end subroutine

subroutine fixcoornuc()
use nucleotmod
use grandmod
use listmod
implicit none
integer j
real  r1,r2,lzrad
if (Qdnafree) return
do j=1,nelenuc
  ! spherical system
  if (Qsphere) then
    r2 = (r(j)%x-cx)**2+(r(j)%y-cy)**2+(r(j)%z-cz)**2
    if (r2.gt.Rsphe2) then
      lzrad=2.0*rsphe/sqrt(r2)-1.0
      r(j)%x = lzrad*(r(j)%x-cx)+cx
      r(j)%y = lzrad*(r(j)%y-cy)+cy
      r(j)%z = lzrad*(r(j)%z-cz)+cz
    endif
  ! cyllindrical system
  elseif (Qecyl) then
    r2=((r(j)%x-cx)*iecx)**2+((r(j)%y-cy)*iecy)**2
    if (r2.gt.1.0) then
       r1=0.999999999999999/sqrt(r2)
       r(j)%x = r1*(r(j)%x-cx)+cx
       r(j)%y = r1*(r(j)%y-cy)+cy
    endif
    if (r(j)%z.lt.lz2m) then
      r(j)%z = lz2m
    else if (r(j)%z.gt.lz2p) then
      r(j)%z = lz2p
    endif
  ! rectangular system
  else
    if (r(j)%x.lt.lx2m) then
      r(j)%x = lx2m
    else if (r(j)%x.gt.lx2p) then !ge
      r(j)%x = lx2p
    endif
    if (r(j)%y.lt.ly2m) then
      r(j)%y = ly2m
    else if (r(j)%y.gt.ly2p) then !ge
      r(j)%y = ly2p
    endif
    if (r(j)%z.lt.lz2m) then
      r(j)%z = lz2m
    else if (r(j)%z.gt.lz2p) then
      r(j)%z = lz2p
    endif
  endif
enddo
end subroutine
