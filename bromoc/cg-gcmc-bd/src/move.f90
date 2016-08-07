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

subroutine move(iat,rr)
use grandmod
use nucleotmod
use listmod

implicit none
integer iat
real radius,nremovablepar
logical*1 endok
type(car) :: rr, rs

!Pick one atom randomly and move it (if it moves outside the 
!limits, pick a different atom)

endok = .true.  
nremovablepar = float(npar - nparnuc)
do while (endok)
  iat = nparnuc + int(nremovablepar*rndm()) + 1 ! [nelenuc+1,nele]
  call getcentroid(iat, rs)
  rr%x = rs%x + mcmax*(rndm()-0.5)
  rr%y = rs%y + mcmax*(rndm()-0.5)
  rr%z = rs%z + mcmax*(rndm()-0.5)
  if (Qsphere) then
    radius = (rr%x-cx)**2+(rr%y-cy)**2+(rr%z-cz)**2
    endok = radius.gt.Rsphe2
  elseif (Qecyl) then
    endok = (((rr%x-cx)*iecx)**2+((rr%y-cy)*iecy)**2).gt.1.0.or.rr%z.lt.lz2m.or.rr%z.gt.lz2p
  else
    endok = rr%x.lt.lx2m.or.rr%x.gt.lx2p.or.rr%y.lt.ly2m.or.rr%y.gt.ly2p.or.rr%z.lt.lz2m.or.rr%z.gt.lz2p
  endif
enddo 
  
return
end subroutine
