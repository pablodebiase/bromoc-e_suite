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

subroutine vdwgd1spln
! -----------------------------------------------------------------------
!This subroutine computes the repulsive potental energy and forces
!based on the grid

!Repulsive forces - VDWF

use ioxmod
use constamod
use stdiomod
use grandmod
use listmod
use gsbpmod
!local variables
implicit none
integer ncyz,ncel3,i,ix,iy,iz,ifir
integer k,l,m,ipx,ipy,ipz
real  vdwfx,vdwfy,vdwfz,esvdw
real  xi,yi,zi,ai,bi,ci,fi,dai,dbi,dci,prefac
real  xc,yc,zc,m3,dm3,evdwgdloc
real  phisum,phivi,phivsv
logical*1 ok

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
evdwgd = 0.0
ifir = 0
esvdw = svdw

!Main loop by atoms
!$omp parallel private(i,ok,xi,yi,zi,ifir,esvdw,vdwfx,vdwfy,vdwfz,ix,iy,iz,phisum,k,ipx,xc,ai,dai,l,ipy,yc,bi,dbi,m,ipz,zc,ci,dci,fi,phivi,phivsv,prefac,evdwgdloc)
evdwgdloc=0.0
!$omp do
do i = 1, nele
  ok=r(i)%x.le.xbcen2+tranx2-dcel2.and.r(i)%x.ge.xbcen2-tranx2+dcel2.and. &
     r(i)%y.le.ybcen2+trany2-dcel2.and.r(i)%y.ge.ybcen2-trany2+dcel2.and. &
     r(i)%z.le.vzmax              .and.r(i)%z.ge.vzmin
  if (.not.ok) cycle
! ion cartesian coordinates in the local grid system
  xi = r(i)%x + tranx2-xbcen2
  yi = r(i)%y + trany2-ybcen2
  zi = r(i)%z + tranz2-zbcen2
  ifir=0
  if (Qnmcden) then
    ifir = (et(i)-1)*ncel3
  else
    if (Qsvdw) esvdw = svdw * scal(et(i))
  endif
  !initializations forces
  vdwfx = 0.0
  vdwfy = 0.0
  vdwfz = 0.0
  !integer counter for ion cartesian coordinates
  ix = nint(xi*idcel2)
  iy = nint(yi*idcel2)
  iz = nint(zi*idcel2)
  if(ix.eq.0) ix=1
  if(iy.eq.0) iy=1
  if(iz.eq.0) iz=1
  if(ix.eq.nclx2-1)ix=nclx2-2
  if(iy.eq.ncly2-1)iy=ncly2-2
  if(iz.eq.nclz2-1)iz=nclz2-2
  !Atom charge distribution by 27 adjacent grid points
  phisum = 0.0
  do k = ix-1, ix+1
    ipx = k*ncyz
    xc = k*dcel2
    ai = 1.5 -(xi-xc)*idcel2
    dai = dm3(ai)
    ai = m3(ai)
    if (ai.ne.0.0.or.dai.ne.0.0) then
      do l = iy-1, iy+1
        ipy= l*nclz2 + ipx
        yc = l*dcel2
        bi = 1.5 -(yi-yc)*idcel2
        dbi = dm3(bi)
        bi = m3(bi)
        if (bi.ne.0.0.or.dbi.ne.0.0) then
          do m = iz-1, iz+1
            ipz = m + ipy + 1
            zc = m*dcel2
            ci = 1.5 -(zi-zc)*idcel2
            dci = dm3(ci)
            ci = m3(ci)
            if (ci.ne.0.0.or.dci.ne.0.0) then
              fi = ai*bi*ci
              phivi=phiv(ipz+ifir)
              phivsv=phivi*esvdw
              phisum = phisum + phivi
              !Repulsive Forces
              if (Qforces) then
                prefac = phivsv*idcel2
                vdwfx = vdwfx + dai*bi*ci*prefac
                vdwfy = vdwfy + ai*dbi*ci*prefac
                vdwfz = vdwfz + ai*bi*dci*prefac
              endif
              !Repulsive Energy 
              evdwgdloc = evdwgdloc + fi*phivsv
            endif
          enddo
        endif
      enddo
    endif
  enddo
  if (Qforces) then 
    f(i)%x = f(i)%x + vdwfx
    f(i)%y = f(i)%y + vdwfy
    f(i)%z = f(i)%z + vdwfz
  endif
enddo
!$omp end do
!$omp critical
evdwgd = evdwgd + evdwgdloc
!$omp end critical
!$omp end parallel
end subroutine
