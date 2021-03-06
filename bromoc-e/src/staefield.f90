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

subroutine staefield(xj,yj,zj,efx,efy,efz)
!-----------------------------------------------------------------------
!This subroutine computes the electrostatic force resulting from the 
!static field generated by the solvent at each atom.

!Static field forces - RXNF

use constamod
use gsbpmod
!local variables
implicit none
integer ncyz,ix,iy,iz,n1,n2,n3,in3
real  efx,efy,efz,one,xj,yj,zj
real  xi,yi,zi,ai,bi,ci,fi
real  aisign,bisign,cisign,prefac
logical*1 ok

one=1.0
ncyz = ncly1*nclz1

!      initialization of efield 
efx = 0.0
efy = 0.0
efz = 0.0

ok=xj.le.xbcen1+tranx1.and.xj.ge.xbcen1-tranx1.and. &
   yj.le.ybcen1+trany1.and.yj.ge.ybcen1-trany1.and. &
   zj.le.zbcen1+tranz1.and.zj.ge.zbcen1-tranz1
if (ok) then
!      ion cartesian coordinates in the local grid system  
  xi = xj + tranx1-xbcen1
  yi = yj + trany1-ybcen1
  zi = zj + tranz1-zbcen1
  ix = int(xi*idcel1) 
  iy = int(yi*idcel1) 
  iz = int(zi*idcel1)
  if (ix.eq.nclx1-1) ix=nclx1-2
  if (iy.eq.ncly1-1) iy=ncly1-2
  if (iz.eq.nclz1-1) iz=nclz1-2
 !    Atom charge distribution by 8 adjacent grid points
  do n1 = ix, ix+1
    ai = xi - n1*dcel1
    aisign = sign(one,ai)
    ai = one - abs(ai)*idcel1
    do n2 = iy, iy+1
      bi = yi - n2*dcel1
      bisign = sign(one,bi)
      bi = one - abs(bi)*idcel1
      do n3 = iz, iz+1
        ci = zi - n3*dcel1
        cisign = sign(one,ci)
        ci = one - abs(ci)*idcel1
        fi = ai*bi*ci ! fraction of the charge assigned to a grid point
        in3 = n1*ncyz + n2*nclz1 + n3 + 1
    ! Electrostatic field
        prefac = phix(in3)*celec2*idcel1
        if ((ai.lt.(1.0-rsmall)).and.(ai.gt.rsmall)) efx = efx + aisign*bi*ci*prefac
        if ((bi.lt.(1.0-rsmall)).and.(bi.gt.rsmall)) efy = efy + bisign*ai*ci*prefac
        if ((ci.lt.(1.0-rsmall)).and.(ci.gt.rsmall)) efz = efz + cisign*ai*bi*prefac
      enddo ! n3
    enddo ! n2
  enddo ! n1
endif ! ok  

return
end subroutine
