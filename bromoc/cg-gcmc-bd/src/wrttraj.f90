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

subroutine wrttraj
use grandmod
use stdiomod
use nucleotmod
use charfuncmod, only: sng   !command parser

implicit none
integer i

write(iuntrj) runtime                    ! simulation time in pico-second (before was nano)
write(iuntrj) nele                       ! total number of ions and sites in motion
write(iuntrj) (et(i),i=1,nele)       ! ion and nucleotides types
write(iuntrj) (sng(r(i)%x),i=1,nele)       ! x coordinates
write(iuntrj) (sng(r(i)%y),i=1,nele)       ! y coordinates
write(iuntrj) (sng(r(i)%z),i=1,nele)       ! z coordinates 
return
end
