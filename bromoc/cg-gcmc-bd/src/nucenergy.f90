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
subroutine nucenergy(parn, energy)
use apfmod
use grandmod
use constamod
use nucleotmod
use extramod
use listmod
use efpmod
implicit none
integer i, j, itype, jtype, is, parn
real dist, dist2, dist6, idist, idist2, idistkd
real de, dc
real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
real energy
integer isite1, isite2, isite3, isite4
real  dd,vard,eefp,fdf,fdv,vard2
real  ang, ang0, varang, cte1
real  dihe, dihe0, vardihe,didstmp
real  epsln,sgex2,f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),m1,m2,m3,modval
real  dist12,n1(3),n2(3),n1v,n2v,v22,v2v,av,bv,iv22,v12,v23
logical*1 ok,ok2
! Initializations
energy     = 0.0
eefpotmx = 0.0
esrpmfmx = 0.0
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

! interaction between interaction sites (nucleotides)         
!  Bonded interactions  
!  -------------------
!  The strech energy is calculated by means of a 
!  Taylor expansion around natural bond lenght.
!  This contribution contains harmonic and anharmonic 
!  interactions.
do i = 1, nbond
  isite1 = sitebond(i,1)
  isite2 = sitebond(i,2)
  dd = distbond(i) ! natural bond length
  v1(1)=r(isite2)%x-r(isite1)%x
  v1(2)=r(isite2)%y-r(isite1)%y
  v1(3)=r(isite2)%z-r(isite1)%z
  dist = sqrt(dot_product(v1,v1))
  vard = dist - dd
  vard2 = vard**2 
  ebond = ebond + epsnuc*vard2*(1.0+100.0*vard2)
  ok = Qforces .and. stfree(isite1).or.stfree(isite2)
  if (ok) then
    de = epsnuc*vard*(2.0+400.0*vard2)/dist
    f1=de*v1
    if (stfree(isite1)) then   
      f(isite1)%x = f(isite1)%x + f1(1)
      f(isite1)%y = f(isite1)%y + f1(2)
      f(isite1)%z = f(isite1)%z + f1(3)
    endif
    if (stfree(isite2)) then
      f(isite2)%x = f(isite2)%x - f1(1)
      f(isite2)%y = f(isite2)%y - f1(2)
      f(isite2)%z = f(isite2)%z - f1(3)
    endif 
  endif
enddo       
!  The bending energy is calculated 
!  by means of a Taylor expansion around natural bond
!  angle. This contribution contains harmonic interactions. 
do i = 1, nangle
  isite1 = siteangle(i,1)
  isite2 = siteangle(i,2) ! central site
  isite3 = siteangle(i,3) 
  ang0 = valangle(i) ! natural bond angle
  v1(1) = r(isite1)%x - r(isite2)%x
  v2(1) = r(isite3)%x - r(isite2)%x
  v1(2) = r(isite1)%y - r(isite2)%y
  v2(2) = r(isite3)%y - r(isite2)%y
  v1(3) = r(isite1)%z - r(isite2)%z
  v2(3) = r(isite3)%z - r(isite2)%z
  m1=dot_product(v1,v1) 
  m2=dot_product(v2,v2)
  m3=dot_product(v1,v2)
  modval=1.0/sqrt(m1*m2)
  ang = acos(m3*modval)
  varang = ang - ang0
  eang = eang + 700.0*epsnuc*varang**2
  ok = Qforces .and. stfree(isite1).or.stfree(isite2).or.stfree(isite3)
  if (ok) then
    de = 1400.0*epsnuc*varang*modval/sin(ang)
    f1=de*v1
    f2=de*v2
    if (stfree(isite1)) then
      f(isite1)%x = f(isite1)%x + f2(1)
      f(isite1)%y = f(isite1)%y + f2(2)
      f(isite1)%z = f(isite1)%z + f2(3)
    endif
    if (stfree(isite2)) then
      f(isite2)%x = f(isite2)%x - (f1(1)+f2(1))
      f(isite2)%y = f(isite2)%y - (f1(2)+f2(2))
      f(isite2)%z = f(isite2)%z - (f1(3)+f2(3))
    endif
    if (stfree(isite3)) then
      f(isite3)%x = f(isite3)%x + f1(1)
      f(isite3)%y = f(isite3)%y + f1(2)
      f(isite3)%z = f(isite3)%z + f1(3)
    endif
  endif
enddo  
 ! The torsional energy is calculated 
 ! using Fourier series for natural torsional angle. 
do i = 1, ndihe
  isite1 = sitedihe(i,1)
  isite2 = sitedihe(i,2)
  isite3 = sitedihe(i,3)
  isite4 = sitedihe(i,4)
  dihe0 = valdihe(i) ! natural torsion angle
  v1(1) = r(isite1)%x - r(isite2)%x
  v2(1) = r(isite3)%x - r(isite2)%x
  v3(1) = r(isite3)%x - r(isite4)%x
  v1(2) = r(isite1)%y - r(isite2)%y
  v2(2) = r(isite3)%y - r(isite2)%y
  v3(2) = r(isite3)%y - r(isite4)%y
  v1(3) = r(isite1)%z - r(isite2)%z
  v2(3) = r(isite3)%z - r(isite2)%z
  v3(3) = r(isite3)%z - r(isite4)%z
  call cross_product(v1,v2,n1)
  call cross_product(v2,v3,n2)
  v22=dot_product(v2,v2)
  v2v=sqrt(v22)
  av=dot_product(v1,n2)
  bv=dot_product(n1,n2)
  dihe = atan2(v2v*av,bv)
  vardihe = dihe - dihe0
  edihe = edihe + 28.0*epsnuc*(1.0-cos(vardihe))
  ok = Qforces .and. stfree(isite1).or.stfree(isite2).or.stfree(isite3).or.stfree(isite4)
  if (ok) then
    n1v=dot_product(n1,n1)
    n2v=dot_product(n2,n2)
    v12=dot_product(v1,v2)
    v23=dot_product(v2,v3)
    de=28.0*epsnuc*sin(vardihe)*v2v
    f1=-de/n1v*n1          
    f4=de/n2v*n2          
    iv22=1.0/v22
    v12=v12*iv22
    v23=v23*iv22
    f2=-f1+v12*f1-v23*f4
    f3=-f4-v12*f1+v23*f4
    if (stfree(isite1)) then
      f(isite1)%x=f(isite1)%x + f1(1)
      f(isite1)%y=f(isite1)%y + f1(2)
      f(isite1)%z=f(isite1)%z + f1(3)
    endif
    if (stfree(isite2)) then
      f(isite2)%x = f(isite2)%x + f2(1)
      f(isite2)%y = f(isite2)%y + f2(2)
      f(isite2)%z = f(isite2)%z + f2(3)
    endif
    if (stfree(isite3)) then
      f(isite3)%x = f(isite3)%x + f3(1)
      f(isite3)%y = f(isite3)%y + f3(2)
      f(isite3)%z = f(isite3)%z + f3(3)
    endif
    if (stfree(isite4)) then
      f(isite4)%x = f(isite4)%x + f4(1) 
      f(isite4)%y = f(isite4)%y + f4(2)
      f(isite4)%z = f(isite4)%z + f4(3)
    endif
  endif
enddo

!Non-bonded interactions
!-----------------------
!Native contacts
do i = 1, nstack
  isite1 = sitestack(i,1)
  isite2 = sitestack(i,2)
  v1(1)=r(isite2)%x-r(isite1)%x
  v1(2)=r(isite2)%y-r(isite1)%y
  v1(3)=r(isite2)%z-r(isite1)%z
  dist2 = dot_product(v1,v1)
  idist2=1.0/dist2
  cte1 = (sgstack(i)**2*idist2)**3   
  estack = estack + 4.0*epsnuc*cte1*(cte1-1.0)
  ok = Qforces .and. stfree(isite1).or.stfree(isite2)
  if (ok) then
    de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
    f1=de*v1 
    if (stfree(isite1)) then
      f(isite1)%x = f(isite1)%x - f1(1)  
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
    endif                            
    if (stfree(isite2)) then         
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3)
    endif
  endif  
enddo  
! Hydrogen bonding
do i = 1, nbp
  isite1 = sitebp(i,1)
  isite2 = sitebp(i,2)
  if (namsite(isite1).eq.'Gb' .or. namsite(isite1).eq.'Cb') then
    epsln = 2.532*epsnuc*scalepairing
  else  
    epsln = 2.0*epsnuc*scalepairing
  endif
  v1(1)=r(isite2)%x-r(isite1)%x
  v1(2)=r(isite2)%y-r(isite1)%y
  v1(3)=r(isite2)%z-r(isite1)%z
  dist2 = dot_product(v1,v1)
  idist2=1.0/dist2
  cte1 = sgbp(i)*sgbp(i)*idist2
  ebp = ebp + epsln*cte1**5*(20.0*cte1-24.0) 
  ok = Qforces .and. stfree(isite1).or.stfree(isite2)
  if (ok) then
    de = 240.0*epsln*cte1**5*(cte1-1.0)*idist2
    f1=de*v1 
    if (stfree(isite1)) then
      f(isite1)%x = f(isite1)%x - f1(1)
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
    endif                            
    if (stfree(isite2)) then         
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3)   
    endif
  endif                
enddo
! Excluded volume
do i = 1, nex
  isite1 = siteex(i,1)
  isite2 = siteex(i,2)
  v1(1)=r(isite2)%x-r(isite1)%x
  v1(2)=r(isite2)%y-r(isite1)%y
  v1(3)=r(isite2)%z-r(isite1)%z
  dist2 = dot_product(v1,v1)
  idist2=1.0/dist2
  sgex2=sgex(i)**2
  if (dist2.lt.sgex2) then
    cte1 = (sgex2*idist2)**3      
    eex = eex + 4.0*epsnuc*cte1*(cte1-1.0) + epsnuc
    ok = Qforces .and. stfree(isite1).or.stfree(isite2)
    if (ok) then           
      de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
      f1=de*v1 
      if (stfree(isite1)) then
        f(isite1)%x = f(isite1)%x - f1(1)
        f(isite1)%y = f(isite1)%y - f1(2)
        f(isite1)%z = f(isite1)%z - f1(3)
      endif                            
      if (stfree(isite2)) then         
        f(isite2)%x = f(isite2)%x + f1(1)
        f(isite2)%y = f(isite2)%y + f1(2)
        f(isite2)%z = f(isite2)%z + f1(3)
      endif
    endif
  endif
enddo 
! Coulomb interaction
do i = 1, nqq
  isite1 = siteqq(i,1)
  isite2 = siteqq(i,2)
  v1(1)=r(isite2)%x-r(isite1)%x
  v1(2)=r(isite2)%y-r(isite1)%y
  v1(3)=r(isite2)%z-r(isite1)%z
  dist=sqrt(dot_product(v1,v1)) 
  idist=1.0/dist
  if (Qdebyhyb) then
    if(outbox(r(isite1)%x,r(isite1)%y,r(isite1)%z).and.outbox(r(isite2)%x,r(isite2)%y,r(isite2)%z)) Qdeby=.true.
  endif
  if (Qdeby) then
    idistkd=exp(-dist*ikappa)
    eqq = eqq + fctn*idist*idistkd
    ok = Qforces .and. stfree(isite1).or.stfree(isite2)
    if (ok) then
      de = fctn*(dist+kappa)*idist**3*ikappa*idistkd
      f1=de*v1
      if (stfree(isite1)) then
        f(isite1)%x = f(isite1)%x - f1(1)
        f(isite1)%y = f(isite1)%y - f1(2)
        f(isite1)%z = f(isite1)%z - f1(3)
      endif                            
      if (stfree(isite2)) then         
        f(isite2)%x = f(isite2)%x + f1(1)
        f(isite2)%y = f(isite2)%y + f1(2)
        f(isite2)%z = f(isite2)%z + f1(3) 
      endif
    endif      
  else
    eqq = eqq + fctn*idist
    ok = Qforces .and. stfree(isite1).or.stfree(isite2)
    if (ok) then
      de = fctn*idist**3
      f1=de*v1
      if (stfree(isite1)) then
        f(isite1)%x = f(isite1)%x - f1(1)
        f(isite1)%y = f(isite1)%y - f1(2)
        f(isite1)%z = f(isite1)%z - f1(3)
      endif
      if (stfree(isite2)) then
        f(isite2)%x = f(isite2)%x + f1(1)
        f(isite2)%y = f(isite2)%y + f1(2)
        f(isite2)%z = f(isite2)%z + f1(3)
      endif
    endif
  endif
  if (Qdebyhyb) Qdeby=.false.
enddo 
!     Solvent-induced contribution
if (Qsolv) then
  do i = 1, nsolv
    isite1 = siteslv(i,1)
    isite2 = siteslv(i,2)
    v1(1)=r(isite2)%x-r(isite1)%x
    v1(2)=r(isite2)%y-r(isite1)%y
    v1(3)=r(isite2)%z-r(isite1)%z
    dist=sqrt(dot_product(v1,v1)) 
    cte1 = exp((13.38-dist)*0.1875)
    esolv = esolv + epsolv*(1.0-cte1)**2 - epsolv
    ok = Qforces .and. stfree(isite1).or.stfree(isite2)
    if (ok) then
      de = 0.375*epsolv*cte1*(cte1-1.0)/dist
      f1=de*v1
      if (stfree(isite1)) then
        f(isite1)%x = f(isite1)%x - f1(1)
        f(isite1)%y = f(isite1)%y - f1(2)
        f(isite1)%z = f(isite1)%z - f1(3)
      endif                            
      if (stfree(isite2)) then         
        f(isite2)%x = f(isite2)%x + f1(1)
        f(isite2)%y = f(isite2)%y + f1(2)
        f(isite2)%z = f(isite2)%z + f1(3)
      endif
    endif
  enddo
endif 
! APFOR    
if (Qapfor) then
  do i=1,afn
    j=sn(i)
    if (Qforces.and.stfree(j)) then
      f(j)%x=f(j)%x+af(1,i)
      f(j)%y=f(j)%y+af(2,i)
      f(j)%z=f(j)%z+af(3,i)
    endif
  enddo
endif
if (Qcontrans) then
  do i=1,ctn
    j=csn(i)
    if (kx(i).ne.0.0) then
      if (j.eq.0) then
        if (Qunsplit) then
          xcon=sum(r(1:nelenuc1st)%x)*inelenuc*2.0-contrx(i)
          f(1:nelenuc1st)%x=f(1:nelenuc1st)%x-kx(i)*xcon*inelenuc*2.0
          econ = econ + 0.5*kx(i)*xcon**2
          xcon=sum(r(1+nelenuc1st:nelenuc)%x)*inelenuc*2.0-contrx(ctn+1)
          f(1+nelenuc1st:nelenuc)%x=f(1+nelenuc1st:nelenuc)%x-kx(i)*xcon*inelenuc*2.0
          econ = econ + 0.5*kx(i)*xcon**2
        else
          xcon=sum(r(1:nelenuc)%x)*inelenuc-contrx(i)
          f(1:nelenuc)%x=f(1:nelenuc)%x-kx(i)*xcon*inelenuc
          econ = econ + 0.5*kx(i)*xcon**2
        endif
      else
        f(j)%x=f(j)%x-kx(i)*(r(j)%x-contrx(i))
        econ = econ + 0.5*kx(i)*(r(j)%x-contrx(i))**2
      endif 
    endif
    if (ky(i).ne.0.0) then
      if (j.eq.0) then
        if (Qunsplit) then
          ycon=sum(r(1:nelenuc1st)%y)*inelenuc*2.0-contry(i)
          f(1:nelenuc1st)%y=f(1:nelenuc1st)%y-ky(i)*ycon*inelenuc*2.0
          econ = econ + 0.5*ky(i)*ycon**2
          ycon=sum(r(1+nelenuc1st:nelenuc)%y)*inelenuc*2.0-contry(ctn+1)
          f(1+nelenuc1st:nelenuc)%y=f(1+nelenuc1st:nelenuc)%y-ky(i)*ycon*inelenuc*2.0
          econ = econ + 0.5*ky(i)*ycon**2
        else
          ycon=sum(r(1:nelenuc)%y)*inelenuc-contry(i)
          f(1:nelenuc)%y=f(1:nelenuc)%y-ky(i)*ycon*inelenuc
          econ = econ + 0.5*ky(i)*ycon**2
        endif
      else
        f(j)%y=f(j)%y-ky(i)*(r(j)%y-contry(i))
        econ = econ + 0.5*ky(i)*(r(j)%y-contry(i))**2
      endif
    endif
    if (kz(i).ne.0.0) then
      if (j.eq.0) then
        if (Qunsplit) then
          zcon=sum(r(1:nelenuc1st)%z)*inelenuc*2.0-contrz(i)
          f(1:nelenuc1st)%z=f(1:nelenuc1st)%z-kz(i)*zcon*inelenuc*2.0
          econ = econ + 0.5*kz(i)*zcon**2
          zcon=sum(r(1+nelenuc1st:nelenuc)%z)*inelenuc*2.0-contrz(ctn+1)
          f(1+nelenuc1st:nelenuc)%z=f(1+nelenuc1st:nelenuc)%z-kz(i)*zcon*inelenuc*2.0
          econ = econ + 0.5*kz(i)*zcon**2
        else
          zcon=sum(r(1:nelenuc)%z)*inelenuc-contrz(i)
          f(1:nelenuc)%z=f(1:nelenuc)%z-kz(i)*zcon*inelenuc
          econ = econ + 0.5*kz(i)*zcon**2
        endif
      else
        f(j)%z=f(j)%z-kz(i)*(r(j)%z-contrz(i))
        econ = econ + 0.5*kz(i)*(r(j)%z-contrz(i))**2
      endif
    endif
  enddo
endif
dids(1:5,nelenuc+1:nele)=0.0
!  nonbonded interactions between interaction sites and ions
if (Qpar .and. Qnucl .and. Qnonbond) then
  do j = 1, nelenuc
    jtype = et(j)
    do i = nelenuc+1, nele
      itype = et(i)
      is=etpidx(itype,jtype)
! Compute distances dna fragment-ion
      dist2 = (r(i)%x-r(j)%x)**2 + (r(i)%y-r(j)%y)**2 + (r(i)%z-r(j)%z)**2
      ok=.false.
      ok2=Qforces .and.(i.gt.nelenuc.or.stfree(j))
      if (Qefpot(is)) then
        if (ok2) then
          call getyd(is,dist2,eefp,de,dist)
        else
          call gety(is,dist2,eefp,dist)
        endif
        eefpotmx=eefpotmx+eefp
        if (Qproxdiff) then
          if (dist.gt.0.0) then
            if (j.eq.1) then 
              dids(1,i)=dist-efp(is)%xl
              dids(2,i)=dist 
              dids(3,i)=r(i)%x-r(j)%x
              dids(4,i)=r(i)%y-r(j)%y
              dids(5,i)=r(i)%z-r(j)%z
            else
              didstmp=dist-efp(is)%xl
              if (didstmp.lt.dids(1,i)) then
                dids(1,i)=didstmp
                dids(2,i)=dist
                dids(3,i)=r(i)%x-r(j)%x
                dids(4,i)=r(i)%y-r(j)%y
                dids(5,i)=r(i)%z-r(j)%z
              endif
            endif
          endif
        endif
      else
        idist2 = 1.0/dist2
        if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
! Compute Coulomb contribution
        if (Qchr(is)) then
          cofo=fct(is)*idist
          eqqmx = eqqmx + cofo
        endif
! Compute Lennard Jones contribution 
        if (Qlj(is)) then 
          dist6 = (sgp2(is)*idist2)**3
          dist12 = dist6**2
          evdwmx = evdwmx + epp4(is)*(dist12-dist6)
        endif
! Compute Short-Range correction contribution
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
            esrpmfmx = esrpmfmx + esrpmf0
          endif
        endif
      endif
! Compute forces
      if (ok2) then    
        if (.not.Qefpot(is)) then
          de=0.0
          if (Qchr(is)) then
            de = cofo*idist2  ! Coulomb force
          endif
          if (Qlj(is)) then
            de = de+epp4(is)*(2.0*dist12-dist6)*6.0*idist2 ! vdw force
          endif
          if (Qsrpmfi(is)) then
            if (dist2.le.rth) then
              dc=(-c2(is)*esrpmf2 + c3(is)*pi*sin(c3(is)*pi*(c1(is)-dist)))*c0(is)*esrpmf1-6.0*c4(is)*esrpmf3*idist! forces
              if (dist.ge.srpx) dc=dc*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
              de = de - dc*idist
            endif
          endif
        endif
        if (de.ne.0.0) then 
          if (i.gt.nelenuc) then
            f(i)%x = f(i)%x + de*(r(i)%x-r(j)%x)
            f(i)%y = f(i)%y + de*(r(i)%y-r(j)%y)
            f(i)%z = f(i)%z + de*(r(i)%z-r(j)%z)
          endif
          if (stfree(j)) then
            f(j)%x = f(j)%x - de*(r(i)%x-r(j)%x)
            f(j)%y = f(j)%y - de*(r(i)%y-r(j)%y)
            f(j)%z = f(j)%z - de*(r(i)%z-r(j)%z)
          endif
        endif
      endif
    enddo ! i=nelenuc+1,...,nele 
  enddo   ! j=1,...,nelenuc
endif
energy = energy + esrpmfmx + ebond + eang + edihe + estack + ebp + eex + eqq + esolv + eqqmx + evdwmx + eefpotmx + econ
end subroutine
