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
energy = 0.0
ne = parl(parn)%ne
sr = parl(parn)%sr
do iat = 1+sr, ne+sr
  itype=et(iat)
  charge=etypl(itype)%chg
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
integer,intent(in) :: parn
real,intent(out) :: energy
integer ne,sr
real srfe(nelenuc+1:nele)
real reff(nelenuc+1:nele)
integer ncyz,ncel3,i,j,ix,iy,iz,n1,n2,n3,in3,ifir,itype,jtype
real gaux1,gaux2
real aux,tau,de,dist2,rfdn,rfcf,aux0,aux1,aux2,srfeij,reffij
real prefa1,prefa2
real xi,yi,zi,ai,bi,ci,fi
real aisign,bisign,cisign
ncyz=ncly3*nclz3
ncel3=nclx3*ncyz
energy=0.0
srfe=0.0
reff=0.0
ne = parl(parn)%ne
sr = parl(parn)%sr
!     Main loop by atoms
do i=nelenuc+1,nele
  itype = et(i)
  if (etypl(itype)%chg.ne.0.0) then
    if (r(i)%x.le.xbcen3+tranx3.and.r(i)%x.ge.xbcen3-tranx3.and. &
        r(i)%y.le.ybcen3+trany3.and.r(i)%y.ge.ybcen3-trany3.and. &
        r(i)%z.le.zbcen3+tranz3.and.r(i)%z.ge.zbcen3-tranz3) then
      if (Qrfpsin) then
        ifir=0
      else
        ifir=(itype-netnuc-1)*ncel3
      endif
      aux1=0.0e0
      aux2=0.0e0
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
            if(gaux1.eq.0.0) exit
            prefa1=gaux1*idcel3
            prefa2=gaux2*idcel3

    !     Local reaction field parameters     
            aux1=aux1+fi*gaux1
            aux2=aux2+fi*gaux2
          enddo
        enddo
      enddo
      srfe(i)=aux1
      reff(i)=aux2
      tau = celec*etypl(itype)%chg
    ! self reaction field energy minus Born energy
      aux = tau*etypl(itype)%chg*srfe(i)
      ! reaction field energy 
      do j=1+nelenuc,i-1
        if ((i.gt.sr.and.i.le.ne+sr).or.(j.gt.sr.and.j.le.ne+sr)) then
          jtype = et(j)
          if (etypl(jtype)%chg.ne.0.0.and.srfe(j).ne.0.0) then
            dist2 = (r(j)%x-r(i)%x)**2+(r(j)%y-r(i)%y)**2+(r(j)%z-r(i)%z)**2
            srfeij = srfe(j)*srfe(i)
            reffij = reff(j)*reff(i)
            rfdn = 1.0/sqrt(reffij*reffij+dist2)
            rfcf = reffij*rfdn
            aux0 = tau*etypl(jtype)%chg
            aux1 = aux0*rfcf
            ! reaction field energy 
            energy = energy + aux1*srfeij
          endif
        endif
      enddo
    endif
  endif
enddo
end subroutine

subroutine par_interact(parn, energy)
use listmod
use grandmod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
logical*1,intent(in),optional :: alert
logical*1 doForce, Qalert
integer ne,sr
real pener

energy = 0.0

if (parn.gt.npar) return

ne = parl(parn)%ne
sr = parl(parn)%sr

! Membrane
if (Qmemb) then 
  call par_membrane(parn, pener, doForce, Qalert) 
  energy=energy+pener
endif

!reaction field parameter  
if (Qrfpar) then
  call par_rfparion(parn, pener)
  energy=energy+pener
endif ! Qrfpar

!static external field contribution
if (Qphix) then
  call staticf0(xj,yj,zj,jtype)
  dener = dener + egsbpa
endif ! Qphix

!grid-based repulsive potential
if (Qphiv) then
  if (Qtrln) then
     call vdwgd0trln(xj,yj,zj,j,jtype,jtype2,Qalert)
  else
     call vdwgd0spln(xj,yj,zj,j,jtype,jtype2,Qalert)
  endif
  dener = dener + evdwgd
endif ! Qphiv

end subroutine
