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
use listmod
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
    ok = .true.
    itype=et(i)
    chi=q(i)
    if (chi.eq.0.0) ok = .false.
    if (ok) ok=r(i)%x.le.xbcen4+tranx4.and.r(i)%x.ge.xbcen4-tranx4.and. &
               r(i)%y.le.ybcen4+trany4.and.r(i)%y.ge.ybcen4-trany4.and. &
               r(i)%z.le.zbcen4+tranz4.and.r(i)%z.ge.zbcen4-tranz4
    if (ok) then
      if (chi.gt.0.0) then
         chitp=chitp+chi
      else
         chitn=chitn+chi
      endif
    endif
  enddo
  ichitp=1.0/chitp
  ichitn=1.0/abs(chitn)
endif

do i = 1, nele
  ok = .true.
  itype=et(i)
  chi=q(i)
  if (chi.eq.0.0) ok = .false.
  if (ok) ok=r(i)%x.le.xbcen4+tranx4.and.r(i)%x.ge.xbcen4-tranx4.and. &
             r(i)%y.le.ybcen4+trany4.and.r(i)%y.ge.ybcen4-trany4.and. &
             r(i)%z.le.zbcen4+tranz4.and.r(i)%z.ge.zbcen4-tranz4
  if (ok) then
!    ion cartesian coordinates in the local grid system  
    xi = r(i)%x + tranx4-xbcen4
    yi = r(i)%y + trany4-ybcen4
    zi = r(i)%z + tranz4-zbcen4
    ix = int(xi*idcel4) 
    iy = int(yi*idcel4) 
    iz = int(zi*idcel4)
!   Atom charge distribution by 8 adjacent grid points
    in3 = ix*ncyz + iy*nclz4 + iz + 1
    if (chi.gt.0.0) then 
      chden(in3)=chden(in3)+sng(chi*ichitp)
    else
      chden(in3)=chden(in3)+sng(chi*ichitn)
    endif
  endif ! ok  
enddo ! i = 1,...,nele

return
end
