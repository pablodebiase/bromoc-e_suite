!    BROMOCng  -  CG-GCMC-BD
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

subroutine energy()
use grandmod
use gsbpmod
use constamod
use extramod
use efpmod
use listmod
implicit none
integer i, j, itype, jtype, is, q, p
integer neq, srq, nep, srp
real dist, dist2, dist6, idist, idist2
real de, dc
real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
real  eefp,fdf,fdv
real  pener
real  dist12

! Initializations
ener     = 0.0
ememb    = 0.0
estaticf = 0.0
evdwgd   = 0.0
erfpar   = 0.0
eintern  = 0.0
eefpot   = 0.0
eelec    = 0.0
evdw     = 0.0
esrpmf   = 0.0

if (Qenergy) then
  ! Initializations
  if (Qforces) then
     do i = 1, nele
        f(i)%x = 0.0
        f(i)%y = 0.0
        f(i)%z = 0.0
     enddo
  endif

  ! membrane contribution
  if (Qmemb) call membrane
  ! ememb

  ! static external field contribution
  if (Qphix) then
     call staticf1
     ! estaticf
  endif

  ! grid-based repulsive potential
  if (Qphiv) then
     if (Qtrln) then
        call vdwgd1trln
     else
        call vdwgd1spln
     endif
     ! evdwgd
  endif
  
  if (Qrfpar) then
    call rfparion
    ! erfpar
  endif

  ! bonded energy
  if (Qbond) then
    if (nparnuc .gt. 0) call nucenergy(eintern)
    pener=0.0
    do i=1+nparnuc,npar
      if (parl(i)%ne.eq.1) cycle
      ! call charmmenergy(pener)
      eintern=eintern+pener
    enddo
  endif

  ! nonbonded interaction between elements
  if (Qnonbond) then
    if (Qproxdiff) dids(1:5,nelenuc+1:nele)=0.0
    do q = nparnuc+1,npar
      neq=parl(q)%ne
      srq=parl(q)%sr
      do p = 1,q-1
        nep=parl(p)%ne
        srp=parl(p)%sr
        do i = srq+1, srq+neq
          itype  = et(i)
          do j = srp+1, srp+nep
            jtype  = et(j)
            is=etpidx(itype,jtype)
            dist2 = dist2car(r(i),r(j))
            if (Qefpot(is)) then
              if (Qforces) then
                call getyd(is,dist2,eefp,de,dist)
              else
                call gety(is,dist2,eefp,dist)
              endif
              eefpot = eefpot + eefp
              if (Qproxdiff) call proxdiff(i,j,is,dist)
            else 
              idist2 = 1.0/dist2
              if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
              if (Qchr(is)) then
                cofo=fct(is)*idist
                eelec = eelec + cofo
              endif
              if (Qlj(is)) then
                dist6 =(sgp2(is)*idist2)**3
                dist12=dist6**2
                evdw = evdw + epp4(is)*(dist12-dist6) ! van der waals potential
              endif
              if (Qsrpmfi(is)) then 
                if (dist2.le.rth) then
                  dist=1.0/idist
                  esrpmf1 = exp((c1(is)-dist)*c2(is))
                  esrpmf2 = cos(c3(is)*pi*(c1(is)-dist))
                  esrpmf3 = (c1(is)*idist)**6
                  esrpmf0 = c0(is)*esrpmf1*esrpmf2+c4(is)*esrpmf3 
                  if (dist.ge.srpx) then ! smoothly fix discontinuity 
                    fdf=exp(-srpk*(dist-srpx))-srpy
                    fdv=esrpmf0 
                    esrpmf0=fdv*fdf
                  endif
                  esrpmf = esrpmf + esrpmf0
                endif
              endif
            endif  
            if (Qforces) then
              if (.not.Qefpot(is)) then
                de=0.0
                if (Qchr(is)) then
                  de=cofo*idist2
                endif
                if (Qlj(is)) then
                  de=de+epp4(is)*(2.0*dist12-dist6)*6.0*idist2 ! van der waals forces
                endif
                if (Qsrpmfi(is)) then 
                  if (dist2.le.rth) then 
                    dc=(-c2(is)*esrpmf2+c3(is)*pi*sin(c3(is)*pi*(c1(is)-dist)))*c0(is)*esrpmf1-6.0*c4(is)*esrpmf3*idist! forces
                    if (dist.ge.srpx) dc=dc*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
                    de=de-dc*idist
                  endif
                endif
              endif
              if (de.ne.0.0) then 
                f(j)%x = f(j)%x + de*(r(j)%x-r(i)%x)
                f(j)%y = f(j)%y + de*(r(j)%y-r(i)%y)
                f(j)%z = f(j)%z + de*(r(j)%z-r(i)%z)
                f(i)%x = f(i)%x - de*(r(j)%x-r(i)%x)
                f(i)%y = f(i)%y - de*(r(j)%y-r(i)%y)
                f(i)%z = f(i)%z - de*(r(j)%z-r(i)%z)
              endif
            endif
          enddo 
        enddo
      enddo
    enddo
    enonbond = eefpot + eelec + evdw + esrpmf 
  !  write(*,*) eefpot, eelec, evdw, esrpmf
  endif !Qnonbond

  ener = ememb + estaticf + evdwgd + erfpar + eintern + enonbond
!  write(*,*) ener, ememb, estaticf, evdwgd, erfpar, eintern, enonbond
endif !Qenergy

end subroutine

