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
subroutine nucenergy(energy)
use apfmod
use grandmod
use constamod
use nucleotmod
use listmod
use efpmod
implicit none
integer i, j 
real dist, dist2, idist, idist2, idistkd
real de
real energy
integer isite1, isite2, isite3, isite4
real  dd,vard,vard2
real  ang, ang0, varang, cte1
real  dihe, dihe0, vardihe
real  epsln,sgex2,f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),m1,m2,m3,modval
real  n1(3),n2(3),n1v,n2v,v22,v2v,av,bv,iv22,v12,v23
! Initializations
energy   = 0.0
ebond    = 0.0
eang     = 0.0
edihe    = 0.0
estack   = 0.0
ebp      = 0.0
eex      = 0.0
eqq      = 0.0
esolv    = 0.0
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
  if (Qforces) then
    de = epsnuc*vard*(2.0+400.0*vard2)/dist
    f1=de*v1
    f(isite1)%x = f(isite1)%x + f1(1)
    f(isite1)%y = f(isite1)%y + f1(2)
    f(isite1)%z = f(isite1)%z + f1(3)
    f(isite2)%x = f(isite2)%x - f1(1)
    f(isite2)%y = f(isite2)%y - f1(2)
    f(isite2)%z = f(isite2)%z - f1(3)
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
  if (Qforces) then
    de = 1400.0*epsnuc*varang*modval/sin(ang)
    f1=de*v1
    f2=de*v2
    f(isite1)%x = f(isite1)%x + f2(1)
    f(isite1)%y = f(isite1)%y + f2(2)
    f(isite1)%z = f(isite1)%z + f2(3)
    f(isite2)%x = f(isite2)%x - (f1(1)+f2(1))
    f(isite2)%y = f(isite2)%y - (f1(2)+f2(2))
    f(isite2)%z = f(isite2)%z - (f1(3)+f2(3))
    f(isite3)%x = f(isite3)%x + f1(1)
    f(isite3)%y = f(isite3)%y + f1(2)
    f(isite3)%z = f(isite3)%z + f1(3)
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
  if (Qforces) then
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
    f(isite1)%x=f(isite1)%x + f1(1)
    f(isite1)%y=f(isite1)%y + f1(2)
    f(isite1)%z=f(isite1)%z + f1(3)
    f(isite2)%x = f(isite2)%x + f2(1)
    f(isite2)%y = f(isite2)%y + f2(2)
    f(isite2)%z = f(isite2)%z + f2(3)
    f(isite3)%x = f(isite3)%x + f3(1)
    f(isite3)%y = f(isite3)%y + f3(2)
    f(isite3)%z = f(isite3)%z + f3(3)
    f(isite4)%x = f(isite4)%x + f4(1) 
    f(isite4)%y = f(isite4)%y + f4(2)
    f(isite4)%z = f(isite4)%z + f4(3)
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
  if (Qforces) then
    de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
    f1=de*v1 
    f(isite1)%x = f(isite1)%x - f1(1)  
    f(isite1)%y = f(isite1)%y - f1(2)
    f(isite1)%z = f(isite1)%z - f1(3)
    f(isite2)%x = f(isite2)%x + f1(1)
    f(isite2)%y = f(isite2)%y + f1(2)
    f(isite2)%z = f(isite2)%z + f1(3)
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
  if (Qforces) then
    de = 240.0*epsln*cte1**5*(cte1-1.0)*idist2
    f1=de*v1 
    f(isite1)%x = f(isite1)%x - f1(1)
    f(isite1)%y = f(isite1)%y - f1(2)
    f(isite1)%z = f(isite1)%z - f1(3)
    f(isite2)%x = f(isite2)%x + f1(1)
    f(isite2)%y = f(isite2)%y + f1(2)
    f(isite2)%z = f(isite2)%z + f1(3)   
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
    if (Qforces) then           
      de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
      f1=de*v1 
      f(isite1)%x = f(isite1)%x - f1(1)
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3)
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
    if (Qforces) then
      de = fctn*(dist+kappa)*idist**3*ikappa*idistkd
      f1=de*v1
      f(isite1)%x = f(isite1)%x - f1(1)
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3) 
    endif      
  else
    eqq = eqq + fctn*idist
    if (Qforces) then
      de = fctn*idist**3
      f1=de*v1
      f(isite1)%x = f(isite1)%x - f1(1)
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3)
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
    if (Qforces) then
      de = 0.375*epsolv*cte1*(cte1-1.0)/dist
      f1=de*v1
      f(isite1)%x = f(isite1)%x - f1(1)
      f(isite1)%y = f(isite1)%y - f1(2)
      f(isite1)%z = f(isite1)%z - f1(3)
      f(isite2)%x = f(isite2)%x + f1(1)
      f(isite2)%y = f(isite2)%y + f1(2)
      f(isite2)%z = f(isite2)%z + f1(3)
    endif
  enddo
endif 
! APFOR    
if (Qapfor) then
  do i=1,afn
    j=sn(i)
    if (Qforces) then
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
energy = ebond + eang + edihe + estack + ebp + eex + eqq + esolv + econ
end subroutine
