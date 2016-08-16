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

subroutine staticf1
!-----------------------------------------------------------------------
!This subroutine computes the electrostatic force resulting from the 
!static field generated by the solvent at each atom.

!Static field forces - RXNF

use constamod
use grandmod
use listmod
use gsbpmod
!local variables
implicit none
integer ncyz,i,ix,iy,iz,n1,n2,n3,in3,itype
real  rxnfx,rxnfy,rxnfz,one
real  chi,xi,yi,zi,ai,bi,ci,fi
real  aisign,bisign,cisign,prefac
logical*1 ok

one=1.0
ncyz = ncly1*nclz1
estaticf = 0.0

!Main loop by atoms


do i = 1, nele
  ok = .true.
  itype=et(i)
  chi=etypl(itype)%chg
  if (chi.eq.0.0) cycle
  if (.not.(r(i)%x.le.xbcen1+tranx1.and.r(i)%x.ge.xbcen1-tranx1.and. &
            r(i)%y.le.ybcen1+trany1.and.r(i)%y.ge.ybcen1-trany1.and. &
            r(i)%z.le.zbcen1+tranz1.and.r(i)%z.ge.zbcen1-tranz1)) cycle
!      initializations forces       
  rxnfx = 0.0
  rxnfy = 0.0
  rxnfz = 0.0
!    ion cartesian coordinates in the local grid system  
  xi = r(i)%x + tranx1-xbcen1
  yi = r(i)%y + trany1-ybcen1
  zi = r(i)%z + tranz1-zbcen1
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
        fi = ai*bi*ci ! fraction of the charge assigned to a
                      ! grid point
        in3 = n1*ncyz + n2*nclz1 + n3 + 1
  ! Electrostatic forces
        if (Qforces) then 
          prefac = phix(in3)*celec*chi*idcel1
          if ((ai.lt.(1.0-rsmall)).and.(ai.gt.rsmall)) then
            rxnfx = rxnfx + aisign*bi*ci*prefac
          endif
          if ((bi.lt.(1.0-rsmall)).and.(bi.gt.rsmall)) then
             rxnfy = rxnfy + bisign*ai*ci*prefac
          endif
          if ((ci.lt.(1.0-rsmall)).and.(ci.gt.rsmall)) then
             rxnfz = rxnfz + cisign*ai*bi*prefac
          endif
        endif
  ! Electrostatic Energy 
         estaticf = estaticf + (fi*chi*phix(in3)*celec)
      enddo ! n3
    enddo ! n2
  enddo ! n1
  if (Qforces) then 
    f(i)%x = f(i)%x + rxnfx
    f(i)%y = f(i)%y + rxnfy
    f(i)%z = f(i)%z + rxnfz
  endif
enddo ! i = 1,...,nele

return
end subroutine
