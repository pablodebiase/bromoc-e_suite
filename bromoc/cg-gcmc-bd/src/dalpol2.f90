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

subroutine dalpol2(j,lmax,mmax,xc,ap,adp)
!------------------------------------------------------------------------
!Derivatives of Associate Legendre Polynomials (From Smythe's BOOK)
!This is different from FUNCTION DALPOL because we don't consider xc > 1 here.

use grandmod
implicit none

integer lmax,mmax,j
real  xc,ap(ntot,0:lmax,0:mmax+1),adp(ntot,0:lmax,0:mmax)
!local
integer l,m
real  fact

if(xc.eq.1.0.or.xc.eq.-1.0) then
   do l = 0,lmax
      adp(j,l,0) = 0.0
      if(xc.eq. 1.0) adp(j,l,0) = l*(l+1.0)*0.5
      if(xc.eq.-1.0) adp(j,l,0) = (-1.0)**(l+1)*l*(l+1.0)*0.5
      do m = 1, mmax
         adp(j,l,m) = 0.0
      enddo
   enddo
   return
endif

do l = 0, lmax
   do m = 0, mmax
      fact = 1.0/sqrt(1.0-xc*xc)
      adp(j,l,m) = fact*(ap(j,l,m+1)-fact*m*xc*ap(j,l,m))
   enddo
enddo
 
return
end subroutine
