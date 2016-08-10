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

subroutine membrane
!Paper => B. Roux, Biophys. J., 73:2980-2989 (1997)
use grandmod
use constamod
use stdiomod
use listmod
use errormod
!LOCAL VARIABLES
implicit none
real sw1, dsw1, sw2, dsw2
real r1,r2, charge
integer iat, itype
real  ir
logical*1 lgiat

ememb = 0.0 ! initilization
lgiat = Qforces

do iat = 1, nele
  itype=et(iat)
  charge=etypl(itype)%chg
  ! Nernst transmembrane potential
  if (voltage.ne.0.0.and.charge.ne.0.0) then
    if (r(iat)%z.lt.zmemb1) then ! REGION 1: z=z(iat)-zmemb1 < 0, lim{z->-inf} pot(1) = 0
      ememb = ememb + charge*afact*exp(ikappa*(r(iat)%z-zmemb1)) ! Eq. (31) paper
      if (lgiat) f(iat)%z = f(iat)%z- charge*afact*ikappa*exp(ikappa*(r(iat)%z-zmemb1))
    elseif ((r(iat)%z.ge.zmemb1) .and. (r(iat)%z.le.zmemb2)) then ! REGION 2
      ememb = ememb + charge*afact*(ceps*ikappa*(r(iat)%z-zmemb1) + 1.0) ! Eq. (31) paper
      if (lgiat) f(iat)%z = f(iat)%z - charge*afact*ceps*ikappa
    elseif (r(iat)%z.gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
      ememb = ememb + charge*(voltage-afact*exp(-ikappa*(r(iat)%z-zmemb2))) ! Eq. (31) paper
      if (lgiat) f(iat)%z = f(iat)%z - charge*afact*ikappa*exp(-ikappa*(r(iat)%z-zmemb2))
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
              warn(itype)=warn(itype)+1
              if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//etypl(itype)%nam,r(iat)%x,r(iat)%y,r(iat)%z
              ememb=ememb+ampl2(itype)
            else        ! particle in pore wall
              if (lgiat) then
                ir=1.0/r1
                f(iat)%x = f(iat)%x + ampl2(itype)*dsw2*r(iat)%x*ir
                f(iat)%y = f(iat)%y + ampl2(itype)*dsw2*r(iat)%y*ir
              endif
              ememb=ememb+ampl2(itype)*sw2
            endif
          else ! if particle is in membrane wall or in membrane+pore wall
            if (sw2.eq.1.0) then ! particle in membrane wall
              if (lgiat) f(iat)%z = f(iat)%z + ampl1(itype)*dsw1
              ememb=ememb+ampl1(itype)*sw1
            else ! particle in membrane and pore wall
              if (lgiat) then
                ir=1.0/r1
                f(iat)%x = f(iat)%x + 0.5*ampl2(itype)*dsw2*r(iat)%x*ir
                f(iat)%y = f(iat)%y + 0.5*ampl2(itype)*dsw2*r(iat)%y*ir
                f(iat)%z = f(iat)%z + 0.5*ampl1(itype)*dsw1
              endif
              ememb = ememb + 0.5*(ampl2(itype)*sw2+ampl1(itype)*sw1)
            endif
          endif
        endif
      else ! no cylindrical channel
        if (sw1.eq.1.0) then
          warn(itype)=warn(itype)+1
          if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//etypl(itype)%nam,r(iat)%x,r(iat)%y,r(iat)%z
        endif
        ememb = ememb + ampl1(itype)*sw1
      endif
    endif
  endif
enddo 

return
end subroutine
