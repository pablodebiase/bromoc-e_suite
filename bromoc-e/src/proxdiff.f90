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
subroutine proxdiff(i,j,is,dist)
use grandmod
use listmod
use efpmod
integer i,j,is
real dist,didstmp
if (j.gt.nelenuc) return
if (dist.gt.0.0) then
  if (j.eq.1) then
    dids(1,i)=dist-efp(is)%xl
    dids(2,i)=dist
    dids(3,i)=r(i)%x-r(j)%x
    dids(4,i)=r(i)%y-r(j)%y
    dids(5,i)=r(i)%z-r(j)%z
  else
    didstmp=dist-efp(is)%xl
    if (didstmp.lt.dids(1,i)) then
      dids(1,i)=didstmp
      dids(2,i)=dist
      dids(3,i)=r(i)%x-r(j)%x
      dids(4,i)=r(i)%y-r(j)%y
      dids(5,i)=r(i)%z-r(j)%z
    endif
  endif
endif
end subroutine
