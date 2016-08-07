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

subroutine find(ib,ip,iat)
use listmod
implicit none
!Find particle number ip from buffer number ib
!INPUTS : ib, ip
!OUTPUT : iat
integer ib,ip,iat,ipx

ipx = 0 
do iat = nparnuc+1, npar
   if (parl(iat)%ibuf.eq.ib) then
      ipx = ipx + 1
      if (ip.eq.ipx) return
   endif
enddo

return
end subroutine
