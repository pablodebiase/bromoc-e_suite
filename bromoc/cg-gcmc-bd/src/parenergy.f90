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

subroutine par_membrane(parn, energy, force, alert)
use listmod
use grandmod
use constamod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
logical*1,intent(in),optional :: force, alert
logical*1 doForce,Qalert
integer ne,sr
real r1,r2, charge
integer iat, itype
real  ir

doForce = .false.
if (present(force)) doForce = force
Qalert = doForce
if (present(force)) Qalert = alert

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
      if (doForce) f(iat)%z = f(iat)%z- charge*afact*ikappa*exp(ikappa*(r(iat)%z-zmemb1))
    elseif ((r(iat)%z.ge.zmemb1) .and. (r(iat)%z.le.zmemb2)) then ! REGION 2
      energy = energy + charge*afact*(ceps*ikappa*(r(iat)%z-zmemb1) + 1.0) ! Eq. (31) paper
      if (doForce) f(iat)%z = f(iat)%z - charge*afact*ceps*ikappa
    elseif (r(iat)%z.gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
      energy = energy + charge*(voltage-afact*exp(-ikappa*(r(iat)%z-zmemb2))) ! Eq. (31) paper
      if (doForce) f(iat)%z = f(iat)%z - charge*afact*ikappa*exp(-ikappa*(r(iat)%z-zmemb2))
    endif
  endif ! voltage
  ! membrane boundary
  if (ampl1(itype).gt.0.0) then
    call switch1(sw1,dsw1,r(iat)%z,p1(1,itype),p1(2,itype),zmemb1,zmemb2)
    if(sw1.ne.0.0) then ! if particle is not in the bulk
      if (Qpore) then ! cylindrical channel
        r2 = r(iat)%x**2+r(iat)%y**2
        call switch2(sw2,dsw2,r2,rcylinder(itype),p2(1,itype),p2(2,itype),r1)
        if (sw2.ne.0.0) then ! if particle is not inside the pore
          if (sw1.eq.1.0) then ! if particle is inside membrane or in pore wall
            if (sw2.eq.1.0) then ! particle in membrane
              if (Qalert) then 
                warn(itype)=warn(itype)+1
                if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//etypl(itype)%nam,r(iat)%x,r(iat)%y,r(iat)%z
              endif
              energy=energy+ampl2(itype)
            else        ! particle in pore wall
              if (doForce) then
                ir=1.0/r1
                f(iat)%x = f(iat)%x + ampl2(itype)*dsw2*r(iat)%x*ir
                f(iat)%y = f(iat)%y + ampl2(itype)*dsw2*r(iat)%y*ir
              endif
              energy=energy+ampl2(itype)*sw2
            endif
          else ! if particle is in membrane wall or in membrane+pore wall
            if (sw2.eq.1.0) then ! particle in membrane wall
              if (doForce) f(iat)%z = f(iat)%z + ampl1(itype)*dsw1
              energy=energy+ampl1(itype)*sw1
            else ! particle in membrane and pore wall
              if (doForce) then
                ir=1.0/r1
                f(iat)%x = f(iat)%x + 0.5*ampl2(itype)*dsw2*r(iat)%x*ir
                f(iat)%y = f(iat)%y + 0.5*ampl2(itype)*dsw2*r(iat)%y*ir
                f(iat)%z = f(iat)%z + 0.5*ampl1(itype)*dsw1
              endif
              energy = energy + 0.5*(ampl2(itype)*sw2+ampl1(itype)*sw1)
            endif
          endif
        endif
      else ! no cylindrical channel
        if (sw1.eq.1.0) then
          if (Qalert) then 
            warn(itype)=warn(itype)+1
            if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//etypl(itype)%nam,r(iat)%x,r(iat)%y,r(iat)%z
          endif
        endif
        energy = energy + ampl1(itype)*sw1
      endif
    endif
  endif
enddo

return
end subroutine

!     calculate the reaction field energy and forces on each ions
subroutine par_rect_rf(parn, energy, force)
use listmod
use grandmod
use constamod
!use ioxmod
!use gsbpmod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
logical*1,intent(in),optional :: force
logical*1 doForce
integer ne,sr
!local variables
real  rxnbfx,rxnbfy,rxnbfz
real  xg,yg,zg,dx,dy,dz,norm
real  ccc
real  xlpol,ylpol,zlpol,xxs,yys,zs
real  lpolx(nele,xnpol),lpoly(nele,ynpol),lpolz(nele,znpol)
real  dlpolx(xnpol),dlpoly(ynpol),dlpolz(znpol)
real  mq(ntpol**2), charge
integer i,ii,jj,ij,mm,n
integer xpol,ypol,zpol,itype

doForce = .false.
if (present(force)) doForce = force

ne = parl(parn)%ne
sr = parl(parn)%sr

!calculate q_{lm} coefficients
do ii = 1, ntpol
   coef(ii) = 0.0
enddo

do i = 1+sr, ne+sr
  itype = et(i)
  charge = etypl(itype)%chg
  xg = r(i)%x   
  yg = r(i)%y   
  zg = r(i)%z   
  if (xg.ge.xmin.and.xg.le.xmax .and.yg.ge.ymin.and.yg.le.ymax .and.zg.ge.zmin.and.zg.le.zmax) then
    xxs = xscal*xg
    yys = yscal*yg
    zs = zscal*zg
    lpolx(i,1) = 1.0
    lpoly(i,1) = 1.0
    lpolz(i,1) = 1.0
    lpolx(i,2) = xxs
    lpoly(i,2) = yys
    lpolz(i,2) = zs
    lpolx(i,3) = 0.5*(3.0*xxs*xxs-1.0)
    lpoly(i,3) = 0.5*(3.0*yys*yys-1.0)
    lpolz(i,3) = 0.5*(3.0*zs*zs-1.0)
    do n = 3, xnpol-1
      lpolx(i,n+1)=((2.0*n-1)*xxs*lpolx(i,n)-(n-1.0)*lpolx(i,n-1))/n
    enddo
    do n = 3, ynpol-1
      lpoly(i,n+1)=((2.0*n-1)*yys*lpoly(i,n)-(n-1.0)*lpoly(i,n-1))/n
    enddo
    do n = 3, znpol-1
      lpolz(i,n+1)=((2.0*n-1)*zs*lpolz(i,n)-(n-1.0)*lpolz(i,n-1))/n
    enddo
    do ii = 1, ntpol
      xpol = lstpx(ii)
      ypol = lstpy(ii)
      zpol = lstpz(ii)
      norm = bnorm(ii)
      coef(ii) = coef(ii) + charge*norm*lpolx(i,xpol+1)*lpoly(i,ypol+1)*lpolz(i,zpol+1)
    enddo
  endif 
enddo

!construct MQ array to speed up the calculations
do ii = 1, ntpol
   mq(ii) = 0.0
   do jj = 1, ntpol
      ij = (ii-1)*ntpol+jj
      mq(ii) = mq(ii) + mmij(ij)*coef(jj)*rfscal
   enddo
enddo

!reaction field energy calculation    
energy=0.0
do ii = 1, ntpol
   energy = energy + 0.5*coef(ii)*mq(ii)
enddo
energy = energy*celec

!reaction field force calculations     
if (Qforces) then 
  do mm = nelenuc+1, nele
    itype = et(mm)
    charge = etypl(itype)%chg
    xg = r(mm)%x  
    yg = r(mm)%y  
    zg = r(mm)%z  
    if (xg.ge.xmin.and.xg.le.xmax .and.yg.ge.ymin.and.yg.le.ymax .and.zg.ge.zmin.and.zg.le.zmax) then          
      rxnbfx = 0.0
      rxnbfy = 0.0
      rxnbfz = 0.0
      ccc = charge*celec
      xxs = xscal*xg
      yys = yscal*yg
      zs = zscal*zg
      dlpolx(1) = 0.0
      dlpoly(1) = 0.0
      dlpolz(1) = 0.0
      dlpolx(2) = 1.0
      dlpoly(2) = 1.0
      dlpolz(2) = 1.0
      dlpolx(3) = 3.0*xxs
      dlpoly(3) = 3.0*yys
      dlpolz(3) = 3.0*zs
      do n = 3, xnpol-1
        dlpolx(n+1) = xxs*dlpolx(n) + n*lpolx(mm,n)
      enddo
      do n = 3, ynpol-1
        dlpoly(n+1) = yys*dlpoly(n) + n*lpoly(mm,n)
      enddo
      do n = 3, znpol-1
        dlpolz(n+1) = zs*dlpolz(n) + n*lpolz(mm,n)
      enddo
      do n = 1, xnpol
        dlpolx(n) = dlpolx(n)*xscal
      enddo
      do n = 1, ynpol
        dlpoly(n) = dlpoly(n)*yscal
      enddo
      do n = 1, znpol
        dlpolz(n) = dlpolz(n)*zscal
      enddo
      do ii = 1, ntpol
        xpol = lstpx(ii)
        ypol = lstpy(ii)
        zpol = lstpz(ii)
        norm = bnorm(ii)
        xlpol = lpolx(mm,xpol+1)
        ylpol = lpoly(mm,ypol+1)
        zlpol = lpolz(mm,zpol+1)
        dx = norm*dlpolx(xpol+1)*ylpol*zlpol
        dy = norm*xlpol*dlpoly(ypol+1)*zlpol
        dz = norm*xlpol*ylpol*dlpolz(zpol+1)
        rxnbfx = rxnbfx - dx*mq(ii)
        rxnbfy = rxnbfy - dy*mq(ii)
        rxnbfz = rxnbfz - dz*mq(ii)
      enddo
      f(mm)%x = f(mm)%x + rxnbfx*ccc
      f(mm)%y = f(mm)%y + rxnbfy*ccc
      f(mm)%z = f(mm)%z + rxnbfz*ccc
    endif 
  enddo
endif
return
end subroutine

subroutine par_interact(parn, energy, alert)
use listmod
use grandmod
implicit none
integer,intent(in) :: parn
real,intent(out) :: energy
logical*1,intent(in),optional :: alert
logical*1 doForce, Qalert
integer ne,sr
real pener

doForce = .false.
if (present(force)) doForce = force
Qalert = .false.
if (present(force)) Qalert = alert
energy = 0.0

if (parn.gt.npar) return

ne = parl(parn)%ne
sr = parl(parn)%sr

if (Qmemb) then 
  call par_membrane(parn, pener, doForce, Qalert) 
  energy=energy+pener
endif

