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

subroutine sphe_svpol(nmpol,lstpl,lstpm)
!-----------------------------------------------------------------------
!store the basis-set number 
!in LSTPL and LSTPM array for Spherical Harmonics

implicit none
integer nmpol,lstpl(*),lstpm(*)
!local
integer l,m,norder

!always the same order in spherical harmonics
norder=0
do l=0,nmpol-1
   norder=norder+1 
   lstpl(norder) =l
   lstpm(norder) =0
   do m=1,l
      norder=norder+1 
      lstpl(norder) =l
      lstpm(norder) =m
      norder=norder+1 
      lstpl(norder) =l
      lstpm(norder) =-m
   enddo
enddo

return
end subroutine
