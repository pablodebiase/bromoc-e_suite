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
subroutine readrfpar(unitfi,unn,outu,adjust)
!==============================================================================
!
!     2006
!     
!     This subroutine reads reaction field parameter files.
!
!     Input:   unitfi (file unit)
!              outu (output unit)
!              ntype (number of ion types)
!              lx, lx, lz (system size)
!              atnam (names of the ions)
!     Outpout: nclx, ncly, nclz (number of grid points in x, y, and z 
!                 direction)
!              dcel (grid spacing)
!              xbcen, ybcen, zbcen (position of grid center)
!              tranx, trany, tranz (half of the length covered by the
!                 grid in x, y, and z direction, e.g. 
!                 tranx=0.5*(nclx-1)*dcel)
!              radion (Born radii of the ions)
!              gsrfen (grids with self reaction field energy parameters)
!              greff (grids with effective radius parameters)
!     Calls:                 
!
use gsbpmod
use grandmod
use errormod
use listmod

implicit none
integer unn
integer unitfi(unn),outu
integer*4 nclx,ncly,nclz
integer*4 nclxa,nclya,nclza
real*8 dcel,xbcen,ybcen,zbcen,aux
real*8 dcela,xbcena,ybcena,zbcena
real*4,allocatable :: gsrfent(:),grefft(:)      
real*8,allocatable :: radion(:)
real*8 xmmn,xmmm,ymmn,ymmm,zmmn,zmmm
integer i,j,k,itype,ncel3,ifir,ilas,iunit,ncel3p
integer x1, x2, y1, y2, z1, z2, niont
integer*4 nxo, nyo, nzo, pos, posn
logical*1 adjust,first


first=.true.
if (Qrfpsin) then
  if (unn.ne.1) call error('rfpar','Number maps must be one if rfpsingle is used',faterr)
  niont=1
else
  niont=netyp-netnuc
  if (niont.ne.unn) call error('rfpar','Number of ion types differ from number of RFPAR types',faterr)
endif
if (.not.allocated(radion)) allocate (radion(niont))

do itype=1,niont
  iunit=unitfi(itype)
  read(iunit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
  read(iunit) radion(itype),aux,aux,aux,aux,aux
  
  if (first) then 
    tranx3 = 0.5*(nclx-1)*dcel
    trany3 = 0.5*(ncly-1)*dcel
    tranz3 = 0.5*(nclz-1)*dcel
    idcel3=1.0/dcel
    ncel3=nclx*ncly*nclz
    allocate (gsrfent(ncel3*niont),grefft(ncel3*niont))
    first=.false.
  endif

  write(outu,*)
  if (Qrfpsin) then
    write(outu,'(6x,a)') 'Single Map will be considered for all types'
    write(outu,'(6x,a,i0,a)') 'Data from unit ',iunit,' is assigned to all particle types'
  else
    write(outu,103) 'Data from unit ',iunit,' is assigned to particle type ',etypl(itype+netnuc)%nam
  endif
  write(outu,102) 'Born radius                        = ',radion(itype)
  write(outu,101) 'Number of grid point in X   (NCLX) = ',nclx 
  write(outu,101) 'Number of grid point in Y   (NCLY) = ',ncly 
  write(outu,101) 'Number of grid point in Z   (NCLZ) = ',nclz
  write(outu,102) 'Grid spacing                (DCEL) = ',dcel
  write(outu,102) 'Center of grid box in X     (XBCEN)= ',xbcen
  write(outu,102) 'Center of grid box in Y     (YBCEN)= ',ybcen
  write(outu,102) 'Center of grid box in Z     (ZBCEN)= ',zbcen
  write(outu,102) 'Grid box in X from ',xbcen-tranx3,' to ',xbcen+tranx3
  write(outu,102) 'Grid box in Y from ',ybcen-trany3,' to ',ybcen+trany3
  write(outu,102) 'Grid box in Z from ',zbcen-tranz3,' to ',zbcen+tranz3

  if (itype.gt.1) then
    if (nclxa.ne.nclx) call error('readrfpar','Different NCLX values in reaction field parameter files !',faterr)
    if (nclya.ne.ncly) call error('readrfpar','Different NCLY values in reaction field parameter files !',faterr)
    if (nclza.ne.nclz) call error('readrfpar','Different NCLZ values in reaction field parameter files !',faterr)
    if (dcela.ne.dcel) call error('readrfpar','Different DCEL values in reaction field parameter files !',faterr)
    if (xbcena.ne.xbcen) call error('readrfpar','Different XBCEN values in reaction field parameter files !',faterr)
    if (ybcena.ne.ybcen) call error('readrfpar','Different YBCEN values in reaction field parameter files !',faterr) 
    if (zbcena.ne.zbcen) call error('readrfpar','Different ZBCEN values in reaction field parameter files !',faterr) 
  end if

  nclxa=nclx
  nclya=ncly
  nclza=nclz
  dcela=dcel
  xbcena=xbcen
  ybcena=ybcen
  zbcena=zbcen

  ifir=(itype-1)*ncel3+1
  ilas=itype*ncel3
  read(iunit) (gsrfent(i),i=ifir,ilas)
  read(iunit) (grefft(i),i=ifir,ilas)
enddo

nclx3=nclx
ncly3=ncly
nclz3=nclz
dcel3=dcel
xbcen3=xbcen
ybcen3=ybcen
zbcen3=zbcen

write(outu,'(6x,/A)') '         SYSTEM-SIZE                        MAP-DIMENSION'
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'x ',lx2m,' - ',lx2p,xbcen3-tranx3,' - ',xbcen3+tranx3
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'y ',ly2m,' - ',ly2p,ybcen3-trany3,' - ',ybcen3+trany3
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'z ',lz2m,' - ',lz2p,zbcen3-tranz3,' - ',zbcen3+tranz3
if (lx2m.lt.xbcen3-tranx3.or.lx2p.gt.xbcen3+tranx3.or.ly2m.lt.ybcen3-trany3.or.ly2p.gt.ybcen3+trany3.or.lz2m.lt.zbcen3-tranz3.or.lz2p.gt.zbcen3+tranz3) then
  write(outu,'(6x,a)') 'SYSTEM-DIMENSION does not fit into MAP-DIMENSION'
else 
  write(outu,'(6x,a)') 'SYSTEM-DIMENSION fits into MAP-DIMENSION'
endif
if (lx2m.gt.xbcen3-tranx3.or.lx2p.lt.xbcen3+tranx3.or.ly2m.gt.ybcen3-trany3.or.ly2p.lt.ybcen3+trany3.or.lz2m.gt.zbcen3-tranz3.or.lz2p.lt.zbcen3+tranz3) then
  write(outu,'(6x,a)') 'MAP-DIMENSION does not fit into SYSTEM-DIMENSION'
else
  write(outu,'(6x,a)') 'MAP-DIMENSION fits into SYSTEM-DIMENSION'
endif

if (adjust) then
  write(outu,'(/6x,a)') 'Resizing map to the system limits to optimize memory usage'
  write(outu,'(6x,a)') 'New boundaries:'
  x1=int((lx2m+tranx3-xbcen3)/dcel3)
  y1=int((ly2m+trany3-ybcen3)/dcel3)
  z1=int((lz2m+tranz3-zbcen3)/dcel3)
  x2=int((lx2p+tranx3-xbcen3)/dcel3)+2
  y2=int((ly2p+trany3-ybcen3)/dcel3)+2
  z2=int((lz2p+tranz3-zbcen3)/dcel3)+2
  if (x1.lt.1) x1=1
  if (y1.lt.1) y1=1
  if (z1.lt.1) z1=1
  if (x2.gt.nclx3) x2=nclx3
  if (y2.gt.ncly3) y2=ncly3
  if (z2.gt.nclz3) z2=nclz3
  xmmn=(x1-1)*dcel3-tranx3+xbcen3
  ymmn=(y1-1)*dcel3-trany3+ybcen3
  zmmn=(z1-1)*dcel3-tranz3+zbcen3
  xmmm=(x2-1)*dcel3-tranx3+xbcen3
  ymmm=(y2-1)*dcel3-trany3+ybcen3
  zmmm=(z2-1)*dcel3-tranz3+zbcen3
  nxo=nclx3
  nyo=ncly3
  nzo=nclz3
  nclx3=x2-x1+1
  ncly3=y2-y1+1
  nclz3=z2-z1+1
  ncel3p=ncel3
  ncel3 = nclx3*ncly3*nclz3
  tranx3 = 0.5*(nclx3-1)*dcel3
  trany3 = 0.5*(ncly3-1)*dcel3
  tranz3 = 0.5*(nclz3-1)*dcel3
  xbcen3 = 0.5*(xmmm+xmmn)
  ybcen3 = 0.5*(ymmm+ymmn)
  zbcen3 = 0.5*(zmmm+zmmn)
  write(outu,101) 'Number of grid point in X   (nclx) = ',nclx3
  write(outu,101) 'Number of grid point in Y   (ncly) = ',ncly3
  write(outu,101) 'Number of grid point in Z   (nclz) = ',nclz3
  write(outu,*)
  write(outu,102) 'Center of box in X          (xbcen)= ',xbcen3
  write(outu,102) 'Center of box in Y          (ybcen)= ',ybcen3
  write(outu,102) 'Center of box in Z          (zbcen)= ',zbcen3
  write(outu,*)
  write(outu,102) 'Box in X from ',xbcen3-tranx3,' to ',xbcen3+tranx3
  write(outu,102) 'Box in Y from ',ybcen3-trany3,' to ',ybcen3+trany3
  write(outu,102) 'Box in Z from ',zbcen3-tranz3,' to ',zbcen3+tranz3
  allocate (gsrfen(ncel3*niont),greff(ncel3*niont))
  posn=0
  do itype=1,niont
    ifir=(itype-1)*ncel3p
    do i=x1,x2
      do j=y1,y2
        do k=z1,z2
          pos=(i-1)*nyo*nzo + (j-1)*nzo + k + ifir
          posn=posn+1
          gsrfen(posn)=gsrfent(pos)
          greff(posn)=grefft(pos)
        enddo
      enddo
    enddo
  enddo
  deallocate (gsrfent,grefft)
else
  allocate (gsrfen(ncel3*niont),greff(ncel3*niont))
  gsrfen=gsrfent
  greff=grefft
  deallocate (gsrfent,grefft)
endif
            
101  format(6X,A,I6)
102  format(6X,A,F8.3,A,F8.3)
103  format(6X,A,I3,A,A)
      
end subroutine


 
