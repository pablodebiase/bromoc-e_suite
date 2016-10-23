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

subroutine par_membrane(parn, energy)
use listmod
use grandmod
use constamod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
integer ne,sr
real r2, charge
integer iat, itype
real sw1, dsw1, sw2, dsw2
energy = 0.0
ne = parl(parn)%ne
sr = parl(parn)%sr
do iat = 1+sr, ne+sr
  itype=et(iat)
  charge=q(iat)
  ! Nernst transmembrane potential
  if (voltage.ne.0.0.and.charge.ne.0.0) then
    if (r(iat)%z.lt.zmemb1) then ! REGION 1: z=z(iat)-zmemb1 < 0, lim{z->-inf} pot(1) = 0
      energy = energy + charge*afact*exp(ikappa*(r(iat)%z-zmemb1)) ! Eq. (31) paper
    elseif ((r(iat)%z.ge.zmemb1) .and. (r(iat)%z.le.zmemb2)) then ! REGION 2
      energy = energy + charge*afact*(ceps*ikappa*(r(iat)%z-zmemb1) + 1.0) ! Eq. (31) paper
    elseif (r(iat)%z.gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
      energy = energy + charge*(voltage-afact*exp(-ikappa*(r(iat)%z-zmemb2))) ! Eq. (31) paper
    endif
  endif ! voltage
  ! membrane boundary
  if (ampl1(itype).gt.0.0) then
    call switch1(sw1,dsw1,r(iat)%z,p1(1,itype),p1(2,itype),zmemb1,zmemb2)
    if(sw1.ne.0.0) then ! if particle is not in the bulk
      if (Qpore) then ! cylindrical channel
        r2 = r(iat)%x**2+r(iat)%y**2
        call switch2(sw2,dsw2,r2,rcylinder(itype),p2(1,itype),p2(2,itype))
        if (sw2.ne.0.0) then ! if particle is not inside the pore
          if (sw1.eq.1.0) then ! if particle is inside membrane or in pore wall
            if (sw2.eq.1.0) then ! particle in membrane
              energy=energy+ampl2(itype)
            else        ! particle in pore wall
              energy=energy+ampl2(itype)*sw2
            endif
          else ! if particle is in membrane wall or in membrane+pore wall
            if (sw2.eq.1.0) then ! particle in membrane wall
              energy=energy+ampl1(itype)*sw1
            else ! particle in membrane and pore wall
              energy = energy + 0.5*(ampl2(itype)*sw2+ampl1(itype)*sw1)
            endif
          endif
        endif
      else ! no cylindrical channel
        energy = energy + ampl1(itype)*sw1
      endif
    endif
  endif
enddo
end subroutine

subroutine par_rfparion(parn, energy)
use constamod
use grandmod
use listmod
use gsbpmod
implicit none
integer parn
real energy,energyloc
integer ne,sr
real srfe(nele)
real reff(nele)
integer ncyz,ncel3,i,j,k,ix,iy,iz,n1,n2,n3,in3,ifir
real gaux1,gaux2
real dist2,rfdn,aux1,aux2,reffij
real xi,yi,zi,ai,bi,ci,fi
real aisign,bisign,cisign
logical ok(nele)
integer li(nele),pnele,ll,l
type(pair) llu(nele*(nele-1)/2)

ncyz=ncly3*nclz3
ncel3=nclx3*ncyz
l=nele*(nele-1)/2
energy=0.0
srfe=0.0
reff=0.0
ne = parl(parn)%ne
sr = parl(parn)%sr
ok=.false.

pnele=0
do i=1,nele
  if (q(i).eq.0.0) cycle
  if (.not.(r(i)%x.le.xbcen3+tranx3.and.r(i)%x.ge.xbcen3-tranx3.and. &
            r(i)%y.le.ybcen3+trany3.and.r(i)%y.ge.ybcen3-trany3.and. &
            r(i)%z.le.zbcen3+tranz3.and.r(i)%z.ge.zbcen3-tranz3)) cycle
  ok(i)=.true.
  pnele=pnele+1
  li(pnele)=i
enddo

!     Main loop by atoms
!$omp parallel private(i,j,k,ifir,aux1,aux2,xi,yi,zi,ix,iy,iz,n1,ai,aisign,n2,bi,bisign,n3,ci,cisign,fi,in3,gaux1,gaux2,dist2,reffij,rfdn,energyloc)
energyloc=0.0
!$omp do
do k=1,pnele
  i=li(k)
  if (Qrfpsin) then
    ifir=0
  else
    ifir=(et(i)-netnuc-1)*ncel3
  endif
  aux1=0.0
  aux2=0.0
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
      ! Local reaction field parameters     
        aux1=aux1+fi*gaux1
        aux2=aux2+fi*gaux2
      enddo
    enddo
  enddo
  srfe(i)=aux1
  reff(i)=aux2
  ! self reaction field energy minus Born energy
  if (i.gt.sr.and.i.le.ne+sr) energyloc=energyloc+0.5*celec*(q(i)*srfe(i))**2
enddo
!$omp end do
!$omp single
ll=0
do k=1,l
  i=lu(k)%a
  j=lu(k)%b
  if (.not.ok(i)) cycle
  if (.not.ok(j)) cycle
  if (srfe(i).eq.0.0) cycle
  if (srfe(j).eq.0.0) cycle
  if (reff(i).eq.0.0) cycle
  if (reff(j).eq.0.0) cycle
  if (.not.((i.gt.sr.and.i.le.ne+sr).or.(j.gt.sr.and.j.le.ne+sr))) cycle
  ll=ll+1
  llu(ll)=lu(k)
enddo
!$omp end single

! reaction field energy 
!$omp do
do k=1,ll
  i=llu(k)%a
  j=llu(k)%b
  dist2 = dist2car(r(i),r(j))
  reffij = reff(i)*reff(j)
  rfdn = 1.0/sqrt(reffij*reffij+dist2)
  ! reaction field energy 
  energyloc = energyloc + celec*q(i)*q(j)*reffij*rfdn*srfe(i)*srfe(j)
enddo
!$omp end do 
!$omp critical
energy=energy+energyloc
!$omp end critical
!$omp end parallel
end subroutine

subroutine par_staticf(parn, energy)
!-----------------------------------------------------------------------
!This subroutine computes only the external static field energy
!for one particle used in subroutine INTERACT in simul.f   
use constamod
use grandmod
use listmod
use gsbpmod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
integer ne,sr
!input variables
integer i 
real  xj, yj, zj
!local variables
integer ncyz,ix,iy,iz,n1,n2,n3,in3
real  chi,xi,yi,zi,ai,bi,ci,fi
logical*1 ok

ncyz = ncly1*nclz1
energy = 0.0
ne = parl(parn)%ne
sr = parl(parn)%sr

do i=sr+1,sr+ne
  chi = q(i)
  if (chi.eq.0.0) cycle
  xj=r(i)%x
  yj=r(i)%y
  zj=r(i)%z
  ok=xj.le.xbcen1+tranx1.and.xj.ge.xbcen1-tranx1
  if (.not.ok) cycle
  ok=yj.le.ybcen1+trany1.and.yj.ge.ybcen1-trany1
  if (.not.ok) cycle
  ok=zj.le.zbcen1+tranz1.and.zj.ge.zbcen1-tranz1
  if (.not.ok) cycle
  !  ion cartesian coordinates in the local grid system              
  xi = xj + tranx1-xbcen1
  yi = yj + trany1-ybcen1
  zi = zj + tranz1-zbcen1
  !  integer*4 counter for ion cartesian coordinates        
  ix = int(xi*idcel1)
  iy = int(yi*idcel1)
  iz = int(zi*idcel1)
  if (ix.eq.nclx1-1) ix=nclx1-2
  if (iy.eq.ncly1-1) iy=ncly1-2
  if (iz.eq.nclz1-1) iz=nclz1-2
  !Atom charge distribution by 8 adjacent grid points
  do n1 = ix, ix+1
    ai = xi - n1*dcel1
    ai = 1.0 - abs(ai)*idcel1
    do n2 = iy, iy+1
      bi = yi - n2*dcel1
      bi = 1.0 - abs(bi)*idcel1
      do n3 = iz, iz+1
        ci = zi - n3*dcel1
        ci = 1.0 - abs(ci)*idcel1
        fi = ai*bi*ci
        in3 = n1*ncyz + n2*nclz1 + n3 + 1
  !Electrostatic Energy 
        energy = energy + (fi*chi*phix(in3)*celec)
      enddo ! n3
    enddo ! n2
  enddo ! n1
enddo
end subroutine

!-----------------------------------------------------------------------
!This subroutine computes only the repulsive potential energy
!for one particle used in subroutine INTERACT in simul.f
!Repulsive forces - VDWF
subroutine par_vdwgd_trln(parn, energy)
use constamod
use grandmod
use listmod
use gsbpmod     
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
!local
integer ne,sr
integer i, itype
real  xj, yj, zj
logical*1 ok
integer ncyz,ncel3,ix,iy,iz,n1,n2,n3,in3,ifir
real  xi,yi,zi,ai,bi,ci,fi,esvdw
real  phisum,phis
ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz

energy = 0.0
ifir = 0
ne = parl(parn)%ne
sr = parl(parn)%sr

do i=sr+1,sr+ne
  itype=et(i)
  xj=r(i)%x
  yj=r(i)%y
  zj=r(i)%z
  ok=xj.le.xbcen2+tranx2.and.xj.ge.xbcen2-tranx2.and. &
     yj.le.ybcen2+trany2.and.yj.ge.ybcen2-trany2.and. &
     zj.le.vzmax.and.zj.ge.vzmin
  if (ok) then
    !ion cartesian coordinates in the local grid system              
    xi = xj + tranx2-xbcen2
    yi = yj + trany2-ybcen2
    zi = zj + tranz2-zbcen2
    esvdw = svdw
    if (Qnmcden) then
      ifir = (itype-1)*ncel3
    else
      if (Qsvdw) esvdw = svdw * scal(itype)
    endif
    !integer counter for ion cartesian coordinates        
    ix = int(xi*idcel2)
    iy = int(yi*idcel2)
    iz = int(zi*idcel2)
    if (ix.eq.nclx2-1) ix=nclx2-2
    if (iy.eq.ncly2-1) iy=ncly2-2
    if (iz.eq.nclz2-1) iz=nclz2-2
    !Atom charge distribution by 8 adjacent grid points
    phisum = 0.0
    do n1 = ix, ix+1
      ai = xi - n1*dcel2
      ai = 1.0 - abs(ai)*idcel2
      do n2 = iy, iy+1
        bi = yi - n2*dcel2
        bi = 1.0 - abs(bi)*idcel2
        do n3=iz,iz+1
          ci = zi - n3*dcel2
          ci = 1.0 - abs(ci)*idcel2
          fi = ai*bi*ci
          in3 = n1*ncyz + n2*nclz2 + n3 + 1
          phis=phiv(in3+ifir)
          phisum = phisum + phis
          energy = energy + fi*esvdw*phis
        enddo ! n3
      enddo ! n2
    enddo ! n1
  endif
enddo 
end subroutine

!-----------------------------------------------------------------------
!This subroutine computes only the repulsive potental energy
!for one particle used in subroutine INTERACT in simul.f
subroutine par_vdwgd_spln(parn, energy)
use constamod
use grandmod
use listmod
use gsbpmod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
!local
integer ne,sr
integer i, itype
real  xj, yj, zj
logical*1 ok
!local 
integer ncyz,ncel3,ix,iy,iz,ifir
integer k,l,m,ipx,ipy,ipz
real  xi,yi,zi,ai,bi,ci,fi,esvdw,phisum,phis
real  xc,yc,zc,m3

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
energy = 0.0
ifir = 0
ne = parl(parn)%ne
sr = parl(parn)%sr

do i=sr+1,sr+ne
  itype=et(i)
  xj=r(i)%x
  yj=r(i)%y
  zj=r(i)%z
  ok=xj.le.xbcen2+tranx2.and.xj.ge.xbcen2-tranx2.and. &
     yj.le.ybcen2+trany2.and.yj.ge.ybcen2-trany2.and. &
     zj.le.vzmax.and.zj.ge.vzmin
  if (ok) then
  !  ion cartesian coordinates in the local grid system
    xi = xj + tranx2-xbcen2
    yi = yj + trany2-ybcen2
    zi = zj + tranz2-zbcen2
    esvdw = svdw
    if (Qnmcden) then
      ifir = (itype-1)*ncel3
    else
      if (Qsvdw) esvdw = svdw * scal(itype)
    endif
    ix=nint(xi*idcel2)
    iy=nint(yi*idcel2)
    iz=nint(zi*idcel2)
    if(ix.eq.0) ix=1
    if(iy.eq.0) iy=1
    if(iz.eq.0) iz=1
    if(ix.eq.nclx2-1)ix=nclx2-2
    if(iy.eq.ncly2-1)iy=ncly2-2
    if(iz.eq.nclz2-1)iz=nclz2-2
    phisum=0.0
    do k = ix-1, ix+1
      ipx = k*ncyz
      xc = k*dcel2
      ai = 1.5 - (xi-xc)*idcel2
      ai = m3(ai)
      if (ai.ne.0.0) then 
        do l = iy-1, iy+1
          ipy = l*nclz2 + ipx
          yc = l*dcel2
          bi = 1.5 - (yi-yc)*idcel2
          bi = m3(bi)
          if (bi.ne.0.0) then
            do m = iz-1, iz+1
              ipz = m + ipy + 1
              zc = m*dcel2
              ci = 1.5 - (zi-zc)*idcel2
              fi = ai*bi*m3(ci)
              if (fi.ne.0.0) then
                phis=phiv(ipz+ifir)
                phisum = phisum + phis
                energy = energy + fi*esvdw*phis
              endif
            enddo
          endif
        enddo
      endif
    enddo
  endif
enddo
end subroutine

subroutine par_interact(parn, energy)
use listmod
use grandmod
use constamod
use gsbpmod
implicit none
integer parn
real energy
integer parm,itype,jtype
integer nen,srn,nem,srm
integer i,j,k,l,is
real dist, dist2, dist6, idist, idist2
real pener,qiqj,enonbondiloc
logical*1 Qchr
type(pair) la((nele-parl(parn)%ne)*parl(parn)%ne)

!real emembi,erfpari,estaticfi,evdwgdi,enonbondi
energy = 0.0

! Membrane
if (Qmemb) then 
  call par_membrane(parn, emembi) 
  energy=energy+emembi
endif

!reaction field parameter  
if (Qrfpar) then
  if (parl(parn)%ne.eq.1) then
    call rfparionj(parn, erfpari)
  else
    call par_rfparion(parn, erfpari)
  endif
  energy=energy+erfpari
endif ! Qrfpar

!static external field contribution
if (Qphix) then
  call par_staticf(parn, estaticfi)
  energy=energy+estaticfi
endif ! Qphix

!grid-based repulsive potential
if (Qphiv) then
  if (Qtrln) then
    call par_vdwgd_trln(parn, evdwgdi)
  else 
    call par_vdwgd_spln(parn, evdwgdi)
  endif
  energy=energy+evdwgdi
endif ! Qphiv

! Ignore all internal energies of every particle
! Compute energy between all particles and the test particle (parn)
if (Qnonbond) then
  ! Create List to divide jobs in threads
  l=0
  srn = parl(parn)%sr
  nen = parl(parn)%ne
  do parm=1,npar
    if (parm.ne.parn) then
      nem = parl(parm)%ne
      srm = parl(parm)%sr
      do i = srm+1,srm+nem
        do j=srn+1,srn+nen
          l=l+1
          la(l)%a=i
          la(l)%b=j
        enddo
      enddo
    endif
  enddo
  enonbondi=0.0
  !$omp parallel private (k,i,j,itype,jtype,is,dist2,pener,enonbondiloc,idist2,qiqj,Qchr,idist,dist6,dist)
  enonbondiloc=0.0
  !$omp do
  do k=1,l
    i=la(k)%a
    j=la(k)%b
    itype=et(i)
    jtype=et(j)
    is=etpidx(itype,jtype)
    dist2 = dist2car(r(i),r(j))
    if (Qefpot(is)) then
      call gety(is,dist2,pener,dist)
      enonbondiloc=enonbondiloc+pener
    else
      idist2 = 1.0/dist2
      qiqj=q(i)*q(j)
      Qchr=qiqj.ne.0.0
      if (Qchr.or.Qsrpmfi(is)) idist = sqrt(idist2)
      !electrostatic interaction    
      if (Qchr) enonbondiloc=enonbondiloc+cecd*qiqj*idist
      !Lennard-Jones 6-12 potential 
      if (Qlj(is)) then
        dist6 = (sgp2(is)*idist2)**3
        enonbondiloc=enonbondiloc+epp4(is)*dist6*(dist6-1.0)
      endif
      !water-mediated short-range ion-ion interaction  
      if (Qsrpmfi(is)) then
        if (dist2.le.rth) then
          dist=1.0/idist
          pener = c0(is)*exp((c1(is)-dist)*c2(is))*cos(c3(is)*pi*(c1(is)-dist))+c4(is)*(c1(is)*idist)**6
          if (dist.ge.srpx) pener=pener*exp(-srpk*(dist-srpx))-srpy  ! smoothly fix discontinuity 
          ! Eq. 9 W. Im,and B. Roux J. Mol. Biol. 322:851-869 (2002)
          enonbondiloc=enonbondiloc+pener
        endif
      endif
    endif
  enddo
  !$omp end do
  !$omp critical
  enonbondi=enonbondi+enonbondiloc
  !$omp end critical
  !$omp end parallel 
  energy=energy+enonbondi
endif !Qnonbond
end subroutine
