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

subroutine chrden
!-----------------------------------------------------------------------
!This subroutine computes the electrostatic force resulting from the 
!static field generated by the solvent at each atom.

!Static field forces - RXNF

use constamod
use grandmod
use nucleotmod
use gsbpmod
use charfuncmod, only: sng   !command parser

!local variables
implicit none
integer ncyz,i,ix,iy,iz,in3,itype
real  chi,xi,yi,zi,chitp,chitn,ichitn,ichitp
logical*1 ok

ncyz = ncly4*nclz4

ichitp=1.0
ichitn=1.0

if (Qchdenorm) then
  chitp=0.0
  chitn=0.0
  do i=1,nele
    if (i.le.nelenuc.or.i.gt.nelenuc) then
      ok = .false.
      if (i.le.nelenuc) then
        if (namsite(i).eq.'P ') then
          chi = cgnuc
          ok = .true.
        endif
      else if (i.gt.nelenuc) then
        itype = abs(typei(i))
        chi = cg(itype)
        ok = .true.
      endif
      if (ok) ok=x(i).le.xbcen4+tranx4.and.x(i).ge.xbcen4-tranx4.and. &
                 y(i).le.ybcen4+trany4.and.y(i).ge.ybcen4-trany4.and. &
                 z(i).le.zbcen4+tranz4.and.z(i).ge.zbcen4-tranz4
      if (ok) then
        if (chi.gt.0.0) then
           chitp=chitp+chi
        else
           chitn=chitn+chi
        endif
      endif
    endif
  enddo
  ichitp=1.0/chitp
  ichitn=1.0/abs(chitn)
endif

do i = 1, nele
  if (i.le.nelenuc .or. i.gt.nelenuc) then 
    ok = .false.
    if (i.le.nelenuc) then
      if (namsite(i).eq.'P ') then
        chi = cgnuc
        ok = .true.
      endif
    else if (i.gt.nelenuc) then
      itype = abs(typei(i))
      chi = cg(itype)
      ok = .true.
    endif
    if (ok) ok=x(i).le.xbcen4+tranx4.and.x(i).ge.xbcen4-tranx4.and. &
               y(i).le.ybcen4+trany4.and.y(i).ge.ybcen4-trany4.and. &
               z(i).le.zbcen4+tranz4.and.z(i).ge.zbcen4-tranz4
    if (ok) then
!      ion cartesian coordinates in the local grid system  
      xi = x(i) + tranx4-xbcen4
      yi = y(i) + trany4-ybcen4
      zi = z(i) + tranz4-zbcen4
      ix = int(xi*idcel4) 
      iy = int(yi*idcel4) 
      iz = int(zi*idcel4)
!     Atom charge distribution by 8 adjacent grid points
      in3 = ix*ncyz + iy*nclz4 + iz + 1
      if (chi.gt.0.0) then 
        chden(in3)=chden(in3)+sng(chi*ichitp)
      else
        chden(in3)=chden(in3)+sng(chi*ichitn)
      endif
    endif ! ok  
  endif  
enddo ! i = 1,...,nele

return
end
