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

subroutine energy
use grandmod
use gsbpmod
use constamod
use extramod
use efpmod
use listmod
implicit none
integer i, j, itype, jtype, is
integer ne, sr
real dist, dist2, dist6, idist, idist2, idistkd
real de, dc
real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
integer isite1, isite2, isite3, isite4
real  dd,vard,eefp,fdf,fdv,vard2
real  ang, ang0, varang, cte1, pener
real  dihe, dihe0, vardihe,didstmp
real  epsln,sgex2,f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),m1,m2,m3,modval
real  dist12,esolve,n1(3),n2(3),n1v,n2v,v22,v2v,av,bv,iv22,v12,v23
logical*1 ok,ok2

! Initializations
ener     = 0.0
eelec    = 0.0
eefpot   = 0.0
eefpotmx = 0.0
evdw     = 0.0
esrpmf   = 0.0
esrpmfmx = 0.0
ememb    = 0.0
esolve   = 0.0
rgsbpa   = 0.0
egsbpb   = 0.0
evdwgd   = 0.0
ebond    = 0.0
eang     = 0.0
edihe    = 0.0
estack   = 0.0 
ebp      = 0.0
eex      = 0.0
eqq      = 0.0
esolv    = 0.0
eqqmx    = 0.0
evdwmx   = 0.0
econ     = 0.0
eintern  = 0.0


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
  ener = ememb

  ! static external field contribution
  if (Qphix) then
     call staticf1
     ener = ener + egsbpa
  endif

  ! grid-based repulsive potential
  if (Qphiv) then
     if (Qtrln) then
        call vdwgd1trln
     else
        call vdwgd1spln
     endif
     ener = ener + evdwgd
  endif
  
  if (Qrfpar) then
    call rfparion
    ener = ener + erfpar
    egsbpb=erfpar
  endif

  ! bonded energy
  if (Qnobond) then
    do i=1,npar
      if (parl(i)%ne.eq.1) cycle
      if (parl(i)%kind.eq.1) then
        call nucenergy(i,pener)
      else
        !charmmenergy
      endif
      eintern=eintern+pener
    enddo
    ener = ener + eintern
  endif
  ! nonbonded interaction between ions
  if (Qnonbond) then
    do i = 1,npar
      do j=i+1,npar
        itype=et(i)
        jtype=et(j)
        
      enddo
    enddo
    do j = nelenuc+1, nele
      jtype  = et(j)
      do i = nelenuc+1, j-1
        itype  = et(i)
        is=etpidx(itype,jtype)
        dist2 = ((r(i)%x-r(j)%x)**2+(r(i)%y-r(j)%y)**2+(r(i)%z-r(j)%z)**2)
        if (Qefpot(is)) then
          if (Qforces) then
            call getyd(is,dist2,eefp,de,dist)
          else
            call gety(is,dist2,eefp,dist)
          endif
          eefpot = eefpot + eefp
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
            if (i.gt.nelenuc) then
              f(i)%x = f(i)%x + de*(r(i)%x-r(j)%x)
              f(i)%y = f(i)%y + de*(r(i)%y-r(j)%y)
              f(i)%z = f(i)%z + de*(r(i)%z-r(j)%z)
            endif
            f(j)%x = f(j)%x - de*(r(i)%x-r(j)%x)
            f(j)%y = f(j)%y - de*(r(i)%y-r(j)%y)
            f(j)%z = f(j)%z - de*(r(i)%z-r(j)%z)
          endif
        endif
      enddo ! i = nelenuc+1,...,(j-1)
    enddo ! j = nelenuc+1,...,nele 
  endif !Qnonbond
  !ener = ener + eelec + evdw + esrpmf + esrpmfmx + ebond + eang + edihe + estack + ebp + eex + eqq + esolv + eqqmx + evdwmx + eefpot + eefpotmx + econ
  ener = eelec + evdw + esrpmf + esrpmfmx + ebond + eang + edihe + estack + ebp + eex + eqq + esolv + eqqmx + evdwmx + eefpot + eefpotmx + econ
endif !Qenergy

end subroutine

