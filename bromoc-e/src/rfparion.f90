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


!==============================================================================
subroutine rfparionj(parn,energy)
!==============================================================================
!
!     2006
!
!     This subroutine calculates the self reaction field energy parameter
!     srfej and the effective radius parameter reffj for ion j from the 
!     grids gsrfen and greff.
!
!     Input:   xj, yj, zj (coordinates of ion j)
!              jtype (ion type of ion j)
!              nclx, ncly, nclz (number of grid points in x, y, and z 
!                 direction)
!              dcel (grid spacing)
!              tranx, trany, tranz (half of the length covered by the
!                 grid in x, y, and z direction, e.g. 
!                 tranx=0.5*(nclx-1)*dcel)
!              xbcen, ybcen, zbcen (position of grid center)
!              gsrfen (grids with self reaction field energy parameters)
!              greff (grids with effective radius parameters)
!              reffac (factor for the effective radius parameters)
!     Outpout: srfej (self reaction field energy parameter of ion j)
!              reffj (effective radius parameter of ion j)
!     Calls:   rfparion              
!
use constamod
use listmod
use grandmod
use gsbpmod
implicit none
real srfe(nele)
real reff(nele)
integer ncyz,ncel3,i,j,ix,iy,iz,n1,n2,n3,in3,ifir,parn
real gaux1,gaux2,xj,yj,zj
real tau,dist2,aux1,aux2,reffij
real energy
real xi,yi,zi,ai,bi,ci,fi
real aisign,bisign,cisign
logical*1 ok
ncyz=ncly3*nclz3
ncel3=nclx3*ncyz
energy=0.0
j=parl(parn)%sr+1
if (q(j).eq.0.0) return
xj=r(j)%x
yj=r(j)%y
zj=r(j)%z
!     Main loop by atoms
if (.not.(xj.le.xbcen3+tranx3.and.xj.ge.xbcen3-tranx3.and. &
          yj.le.ybcen3+trany3.and.yj.ge.ybcen3-trany3.and. &
          zj.le.zbcen3+tranz3.and.zj.ge.zbcen3-tranz3)) return
srfe=0.0
reff=0.0
do i=1,nele
  if (q(i).eq.0.0) cycle
  xi=r(i)%x+tranx3-xbcen3
  yi=r(i)%y+trany3-ybcen3
  zi=r(i)%z+tranz3-zbcen3
  ok=r(i)%x.le.xbcen3+tranx3.and.r(i)%x.ge.xbcen3-tranx3.and. &
     r(i)%y.le.ybcen3+trany3.and.r(i)%y.ge.ybcen3-trany3.and. &
     r(i)%z.le.zbcen3+tranz3.and.r(i)%z.ge.zbcen3-tranz3
  if (.not.ok) cycle
  if (Qrfpsin) then
    ifir=0
  else
    ifir=(et(i)-netnuc-1)*ncel3
  endif
  aux1=0.0
  aux2=0.0
  ix=int(xi*idcel3)
  iy=int(yi*idcel3)
  iz=int(zi*idcel3)
  if (ix.eq.nclx3-1) ix=nclx3-2
  if (iy.eq.ncly3-1) iy=ncly3-2
  if (iz.eq.nclz3-1) iz=nclz3-2
  ! Calculate GB radius from 8 next neighbor grid point values    
  do n1=ix,ix+1
    ai=xi-n1*dcel3
    aisign=sign(1.0,ai)
    ai=1.0-abs(ai)*idcel3
    do n2=iy,iy+1
      bi=yi-n2*dcel3
      bisign=sign(1.0,bi)
      bi=1.0-abs(bi)*idcel3
      do n3=iz,iz+1
        ci=zi-n3*dcel3
        cisign=sign(1.0,ci)
        ci=1.0-abs(ci)*idcel3
        fi=ai*bi*ci
        in3=n1*ncyz+n2*nclz3+n3+1
        gaux1=gsrfen(in3+ifir)
        gaux2=sqrfac*greff(in3+ifir)
        ! Local reaction field parameters     
        aux1=aux1+fi*gaux1
        aux2=aux2+fi*gaux2
      enddo
    enddo
  enddo
  srfe(i)=aux1
  reff(i)=aux2
enddo

if (srfe(j).eq.0.0) return
tau = celec*q(j)
! self reaction field energy minus Born energy
energy = energy+0.5*tau*q(j)*srfe(j)**2
if (reff(j).eq.0.0) return
do i=1,nele
  if (i.eq.j) cycle
  if (q(i).eq.0.0) cycle
  if (srfe(i).eq.0.0) cycle
  if (reff(i).eq.0.0) cycle
  dist2 = dist2car(r(i),r(j))
  reffij = reff(j)*reff(i)
  ! reaction field energy 
  energy=energy+tau*q(i)*reffij*srfe(j)*srfe(i)/sqrt(reffij**2+dist2)
enddo
end subroutine
      
!==============================================================================
subroutine rfparion
!(nele,x,y,z,type,nclx,ncly,nclz,dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,gsrfen,greff,reffac,srfedx,srfedy,srfedz,srfe,reffdx,reffdy,reffdz,reff,qforces)
!==============================================================================
!
!     2006
!
!     This subroutine calculates the self reaction field energy parameters
!     srfe, the effective radius parameters reff, and the corresponding
!     derivatives for all ions from the grids gsrfen and greff.
!
!     Input:   nele (number of ions)
!              x, y, z (coordinates of the ions)
!              type (ion types)
!              nclx, ncly, nclz (number of grid points in x, y, and z 
!                 direction)
!              dcel (grid spacing)
!              tranx, trany, tranz (half of the length covered by the
!                 grid in x, y, and z direction, e.g. 
!                 tranx=0.5*(nclx-1)*dcel)
!              xbcen, ybcen, zbcen (position of grid center)
!              gsrfen (grids with self reaction field energy parameters)
!              greff (grids with effective radius parameters)
!              reffac (factor for the effective radius parameters)
!              Qforces (calculate forces if Qforces is true)
!     Outpout: srfedx,srfedy,srfedz (derivation of the self reaction 
!                 field energy parameters of the ions)
!              srfe (self reaction field energy parameters of the ions)
!              reffdx,reffdy,reffdz (derivation of the effective radius
!                 parameters of the ions)
!              reff (effective radius parameters of the ions)
!     Calls:                 
!
use constamod
use grandmod
use listmod
use gsbpmod
implicit none

real srfedx(nele),srfedy(nele),srfedz(nele),srfe(nele)
real reffdx(nele),reffdy(nele),reffdz(nele),reff(nele)
integer ncyz,ncel3,i,j,ix,iy,iz,n1,n2,n3,in3,ifir
real aux1dx,aux1dy,aux1dz,gaux1
real aux2dx,aux2dy,aux2dz,gaux2
real aux,tau,de,dist2,rfdn,rfcf,aux0,aux1,aux2,aux3,srfeij,reffij
real prefa1,prefa2
real xi,yi,zi,ai,bi,ci,fi
real aisign,bisign,cisign

ncyz=ncly3*nclz3
ncel3=nclx3*ncyz
erfpar=0.0
srfe=0.0
reff=0.0
srfedx=0.0
srfedy=0.0
srfedz=0.0
reffdx=0.0
reffdy=0.0
reffdz=0.0

!     Main loop by atoms
do i=1,nele
  if (q(i).eq.0.0) cycle
  if (.not.(r(i)%x.le.xbcen3+tranx3.and.r(i)%x.ge.xbcen3-tranx3.and. &
      r(i)%y.le.ybcen3+trany3.and.r(i)%y.ge.ybcen3-trany3.and. &
      r(i)%z.le.zbcen3+tranz3.and.r(i)%z.ge.zbcen3-tranz3)) cycle
  if (Qrfpsin) then
    ifir=0
  else
    ifir=(et(i)-netnuc-1)*ncel3
  endif
  aux1=0.0
  aux1dx=0.0
  aux1dy=0.0
  aux1dz=0.0
  aux2=0.0
  aux2dx=0.0
  aux2dy=0.0
  aux2dz=0.0
  xi=r(i)%x+tranx3-xbcen3
  yi=r(i)%y+trany3-ybcen3
  zi=r(i)%z+tranz3-zbcen3
  ix=int(xi*idcel3)
  iy=int(yi*idcel3)
  iz=int(zi*idcel3)
  if (ix.eq.nclx3-1) ix=nclx3-2
  if (iy.eq.ncly3-1) iy=ncly3-2
  if (iz.eq.nclz3-1) iz=nclz3-2
  !     Calculate GB radius from 8 next neighbor grid point values    
  do n1=ix,ix+1
    ai=xi-n1*dcel3
    aisign=sign(1.0,ai)
    ai=1.0-abs(ai)*idcel3
    do n2=iy,iy+1
      bi=yi-n2*dcel3
      bisign=sign(1.0,bi)
      bi=1.0-abs(bi)*idcel3
      do n3=iz,iz+1
        ci=zi-n3*dcel3
        cisign=sign(1.0,ci)
        ci=1.0-abs(ci)*idcel3
        fi=ai*bi*ci
        in3=n1*ncyz+n2*nclz3+n3+1
        gaux1=gsrfen(in3+ifir) 
        gaux2=sqrfac*greff(in3+ifir)
        prefa1=gaux1*idcel3
        prefa2=gaux2*idcel3
!     Local reaction field parameters     
        aux1=aux1+fi*gaux1
        aux2=aux2+fi*gaux2
        if (Qforces) then  
!     Local reaction field parameter derivatives     
          if ((ai.lt.(1.0-rsmall)).and.(ai.gt.rsmall)) then
            aux=aisign*bi*ci
            aux1dx=aux1dx-aux*prefa1
            aux2dx=aux2dx-aux*prefa2
          end if
          if ((bi.lt.(1.0-rsmall)).and.(bi.gt.rsmall)) then
            aux=bisign*ai*ci
            aux1dy=aux1dy-aux*prefa1
            aux2dy=aux2dy-aux*prefa2
          end if
          if ((ci.lt.(1.0-rsmall)).and.(ci.gt.rsmall)) then
            aux=cisign*ai*bi
            aux1dz=aux1dz-aux*prefa1
            aux2dz=aux2dz-aux*prefa2
          end if 
        end if   
      enddo
    enddo
  enddo
  srfe(i)=aux1
  reff(i)=aux2
  if (Qforces) then  
    srfedx(i)=aux1dx
    srfedy(i)=aux1dy
    srfedz(i)=aux1dz
    reffdx(i)=aux2dx
    reffdy(i)=aux2dy
    reffdz(i)=aux2dz            
  endif
  tau = celec*q(i)
! self reaction field energy minus Born energy
  aux = tau*q(i)*srfe(i)
  ! reaction field energy 
  erfpar = erfpar + 0.5*aux*srfe(i)
! forces related to the reaction field energy                  
  if(Qforces)then
    ! forces caused by variation of srfe(j)
    de = -aux
    f(i)%x = f(i)%x + de*srfedx(i)
    f(i)%y = f(i)%y + de*srfedy(i)
    f(i)%z = f(i)%z + de*srfedz(i)
  endif
  do j=1,i-1
    if (q(j).eq.0.0) cycle 
    dist2 = (r(j)%x-r(i)%x)**2+(r(j)%y-r(i)%y)**2+(r(j)%z-r(i)%z)**2
    srfeij = srfe(j)*srfe(i)
    reffij = reff(j)*reff(i)
    rfdn = 1.0/sqrt(reffij*reffij+dist2)
    rfcf = reffij*rfdn
    aux0 = tau*q(j)
    aux1 = aux0*rfcf
    ! reaction field energy 
    erfpar = erfpar + aux1*srfeij
    ! forces related to the reaction field energy                  
    if(Qforces)then
      aux2 = aux0*srfeij*rfdn**3
      ! forces due to variation of srfe(i) and srfe(j)
      de = -aux1*srfe(i)
      f(j)%x = f(j)%x + de*srfedx(j)
      f(j)%y = f(j)%y + de*srfedy(j)
      f(j)%z = f(j)%z + de*srfedz(j)
      de = -aux1*srfe(j)
      f(i)%x = f(i)%x + de*srfedx(i)
      f(i)%y = f(i)%y + de*srfedy(i)
      f(i)%z = f(i)%z + de*srfedz(i)
      ! forces due to variation of reff(i) and reff(j)
      aux3=-aux2*(dist2+reff(j)*reff(i)*0.5)
      de = aux3*reff(i)
      f(j)%x = f(j)%x + de*reffdx(j)
      f(j)%y = f(j)%y + de*reffdy(j)
      f(j)%z = f(j)%z + de*reffdz(j)
      de = aux3*reff(j)
      f(i)%x = f(i)%x + de*reffdx(i)
      f(i)%y = f(i)%y + de*reffdy(i)
      f(i)%z = f(i)%z + de*reffdz(i)
      ! forces due to variation of ion positions 
      de = aux2*reffij
      f(j)%x = f(j)%x + de*(r(j)%x-r(i)%x)
      f(j)%y = f(j)%y + de*(r(j)%y-r(i)%y)
      f(j)%z = f(j)%z + de*(r(j)%z-r(i)%z)
      f(i)%x = f(i)%x - de*(r(j)%x-r(i)%x)
      f(i)%y = f(i)%y - de*(r(j)%y-r(i)%y)
      f(i)%z = f(i)%z - de*(r(j)%z-r(i)%z)
    endif
  enddo
enddo

end subroutine

