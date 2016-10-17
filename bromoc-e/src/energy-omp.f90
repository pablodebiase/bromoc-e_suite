!    BROMOC-E
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA with External Force Field for Explicit Particles Program.
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
use efpmod
use listmod
implicit none
integer i, j, k, l, itype, jtype, is 
real dist, dist2, dist6, idist, idist2
real de, dc
real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
real  eefp,fdf,fdv
real  pener, ehcons
real  dist12
real qiqj
real einternloc,eefpotloc,eelecloc,esrpmfloc,evdwloc
logical*1 Qchr
type(car) floc(nele)
type(pair) la(nele*(nele-1)/2)
integer multiat(npar),man

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
ehcons   = 0.0

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
    ! Make list to divide processes
    man=0
    do i=1+nparnuc,npar
      if (parl(i)%ne.gt.1) then 
         man=man+1
         multiat(man)=i
      endif
    enddo
    !$omp parallel private (pener,i,einternloc)
    pener=0.0
    einternloc=0.0    
    !$omp do
    do j=1,man
      i=multiat(j)
      if (ptypl(parl(i)%ptyp)%Qpsf) call charmmenergy(i,pener)
      einternloc=einternloc+pener
    enddo
    !$omp end do
    !$omp critical
    eintern=eintern+einternloc
    !$omp end critical
    !$omp end parallel 
  endif

  ! nonbonded interaction between elements
  if (Qnonbond) then
    if (Qproxdiff) dids(1:5,nelenuc+1:nele)=0.0
    ! Create List to divide jobs in threads
    l=0
    do k=1,nele*(nele-1)/2
      if (pe(lu(k)%a).eq.pe(lu(k)%b)) cycle ! if elements belongs to the same particle skip
      if (pe(lu(k)%a).le.nparnuc.and.pe(lu(k)%b).le.nparnuc) cycle ! if elements belongs to DNA skip, double stranded DNA is handled by nucenergy
      l=l+1
      la(l)=lu(k)
    enddo
    !$omp parallel private (k,i,j,itype,jtype,is,dist2,eefpotloc,eefp,de,idist2,qiqj,Qchr,idist,cofo,eelecloc,dist6,dist12,evdwloc,dist,esrpmf1,esrpmf2,esrpmf3,esrpmf0,fdf,fdv,esrpmfloc,dc,floc)
    eefpotloc=0.0
    eelecloc=0.0
    evdwloc=0.0
    esrpmfloc=0.0
    if (Qforces) then
      do i=1,nele
        call setcarzero(floc(i))
      enddo
    endif
    !$omp do
    do k=1,l
      i=la(k)%a
      j=la(k)%b
      itype=et(i)
      jtype=et(j)
      is=etpidx(itype,jtype)
      dist2 = dist2car(r(i),r(j))
      if (Qefpot(is)) then
        if (Qforces) then
          call getyd(is,dist2,eefp,de,dist)
        else
          call gety(is,dist2,eefp,dist)
        endif
        eefpotloc = eefpotloc + eefp
        if (Qproxdiff) call proxdiff(i,j,is,dist)
      else 
        idist2 = 1.0/dist2
        qiqj=q(i)*q(j)
        Qchr=qiqj.ne.0.0
        if (Qchr.or.Qsrpmfi(is)) idist = sqrt(idist2)
        if (Qchr) then 
          cofo=cecd*qiqj*idist
          eelecloc = eelecloc + cofo
        endif
        if (Qlj(is)) then
          dist6 =(sgp2(is)*idist2)**3
          dist12=dist6**2
          evdwloc = evdwloc + epp4(is)*(dist12-dist6) ! van der waals potential
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
            esrpmfloc = esrpmfloc + esrpmf0
          endif
        endif
      endif  
      if (Qforces) then
        if (.not.Qefpot(is)) then
          de=0.0
          if (Qchr) de=cofo*idist2
          if (Qlj(is)) de=de+epp4(is)*(2.0*dist12-dist6)*6.0*idist2 ! van der waals forces
          if (Qsrpmfi(is)) then 
            if (dist2.le.rth) then 
              dc=(-c2(is)*esrpmf2+c3(is)*pi*sin(c3(is)*pi*(c1(is)-dist)))*c0(is)*esrpmf1-6.0*c4(is)*esrpmf3*idist! forces
              if (dist.ge.srpx) dc=dc*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
              de=de-dc*idist
            endif
          endif
        endif
        if (de.ne.0.0) then 
          floc(j)%x = floc(j)%x + de*(r(j)%x-r(i)%x)
          floc(j)%y = floc(j)%y + de*(r(j)%y-r(i)%y)
          floc(j)%z = floc(j)%z + de*(r(j)%z-r(i)%z)
          floc(i)%x = floc(i)%x - de*(r(j)%x-r(i)%x)
          floc(i)%y = floc(i)%y - de*(r(j)%y-r(i)%y)
          floc(i)%z = floc(i)%z - de*(r(j)%z-r(i)%z)
        endif
      endif
    enddo
    !$omp end do
    !$omp critical
    eefpot = eefpot + eefpotloc
    eelec = eelec + eelecloc
    evdw = evdw + evdwloc
    esrpmf = esrpmf + esrpmfloc
    if (Qforces) then 
      do i=1,nele
        call addcar(f(i),floc(i))
      enddo
    endif
    !$omp end critical
    !$omp end parallel 
    enonbond = eefpot + eelec + evdw + esrpmf 
  endif !Qnonbond

  ! Apply Harmonic Constrains
  do i=1,nefix
    j=efix(i)%fen
    if (efix(i)%fc%x.ne.0.0) then
      if (Qforces) f(j)%x=f(j)%x-efix(i)%fc%x*(r(j)%x-efix(i)%rfx%x)
      ehcons = ehcons + 0.5*efix(i)%fc%x*(r(j)%x-efix(i)%rfx%x)**2
    endif
    if (efix(i)%fc%y.ne.0.0) then
      if (Qforces) f(j)%y=f(j)%y-efix(i)%fc%y*(r(j)%y-efix(i)%rfx%y)
      ehcons = ehcons + 0.5*efix(i)%fc%y*(r(j)%y-efix(i)%rfx%y)**2
    endif
    if (efix(i)%fc%z.ne.0.0) then
      if (Qforces) f(j)%z=f(j)%z-efix(i)%fc%z*(r(j)%z-efix(i)%rfx%z)
      ehcons = ehcons + 0.5*efix(i)%fc%z*(r(j)%z-efix(i)%rfx%z)**2
    endif
  enddo
 
  ! Add all energies
  ener = ememb + estaticf + evdwgd + erfpar + eintern + enonbond + ehcons
endif !Qenergy

end subroutine

