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
subroutine charmmenergy(energy)
use explatmod
implicit none
real energy
energy=0.0
! bonded interaction among explicit atoms
if (Qlbond) call en2cen ! Bond terms
if (Qlang) call en3cen ! Angle terms
if (Qlubs) call ubr ! UB terms
if (Qldih) call en4cen ! Dihedral angle terms
if (Qldef) call improper ! Improper angle terms
if (Qlcmap) call cmapr ! CMAP terms
end subroutine

subroutine en2cen
! Bond terms for explicit atoms
use explatmod
use nucleotmod
use grandmod

implicit none
! local variables
integer kk,ibond,iat,jat,itype
real bondk,rij0,rbij2,rbij,enbond,f1,dij(3),fb(3)
logical*1 ok

!     bond .....
!
!     iat --- jat
!         dij

do ibond = 1,nbonds
  ! force centers
  iat = nsites + bonds(1,ibond)
  jat = nsites + bonds(2,ibond)
  ! bond type
  itype = bonds(3,ibond)
  ! energy constant [Kcal/mole/Angs.**2]
  bondk = stretch(1,itype)
  ! natural bond [Angs.]
  rij0 = stretch(2,itype)
  ! distance
  dij(1) = x(jat) - x(iat)
  dij(2) = y(jat) - y(iat)
  dij(3) = z(jat) - z(iat)
  ! *** Energy contribution
  rbij2 = dot_product(dij,dij)
  rbij = sqrt(rbij2)
  enbond = bondk*(rbij-rij0)*(rbij-rij0)
  ebond = ebond + enbond
  ! *** Forces calculation
  ok = Qforces .and. (.not.fixat(iat).or..not.fixat(jat))
  if (ok) then
    f1 = 2.0*bondk*(rbij-rij0)/rbij
    do kk = 1,3
      fb(kk) = f1*dij(kk)
    end do
    if (.not.fixat(iat)) then
      fx(iat) = fx(iat) + fb(1)
      fy(iat) = fy(iat) + fb(2)
      fz(iat) = fz(iat) + fb(3)
    end if
    if (.not.fixat(jat)) then
      fx(jat) = fx(jat) - fb(1)
      fy(jat) = fy(jat) - fb(2)
      fz(jat) = fz(jat) - fb(3)
    end if
  end if
end do ! next ibond

return
end subroutine

subroutine en3cen
! Angle terms for explicit atoms
use explatmod
use nucleotmod
use grandmod
use constamod

implicit none
! local variables
integer kk,ibend,iat,jat,kat,itype
real bendk,aijk0,r12,r22,r1r2,modval,cst
real bondangle,snt,enbend,force,r12inv,r22inv,fiat(3),fjat(3),fkat(3)
real rji(3),rjk(3)
real pos1(3),pos2(3),pos3(3)
real, parameter :: zero=1.0e-15
logical*1 ok

!     bend angle ......
!
!     iat --- jat --- kat
!     pos1    pos2    pos3

do ibend = 1,nbends
  ! force centers
  iat = nsites + bends(1,ibend)
  jat = nsites + bends(2,ibend) ! central atom
  kat = nsites + bends(3,ibend)
  ! bend angle type
  itype  = bends(4,ibend)
  ! energy constant [Kcal/mole/radians**2]
  bendk  = bend(1,itype)
  ! natural bending angle [radians]
  aijk0  = bend(2,itype)*radians
  ! positions
  pos1(1) = x(iat)
  pos1(2) = y(iat)
  pos1(3) = z(iat)
  pos2(1) = x(jat)
  pos2(2) = y(jat)
  pos2(3) = z(jat)
  pos3(1) = x(kat)
  pos3(2) = y(kat)
  pos3(3) = z(kat)
  ! *** Bond angle calculation
  ! interparticle vectors
  do kk  =  1,3
    rji(kk) =  pos1(kk) - pos2(kk)
    rjk(kk) =  pos3(kk) - pos2(kk)
  end do
  ! bond angle
  r12 = dot_product(rji,rji)
  r22 = dot_product(rjk,rjk)
  r1r2 = 1.0/(r12*r22)
  modval = sqrt(r1r2)
  cst = dot_product(rji,rjk)*modval
  if (cst.lt.-1.0) cst = -1.0
  if (cst.gt.1.0) cst = 1.0
  bondangle = acos(cst)
  ! **** Energy contribution
  enbend = bendk*(bondangle-aijk0)*(bondangle-aijk0)
  eang = eang + enbend
  ok = Qforces .and. (.not.fixat(iat).or..not.fixat(jat).or..not.fixat(kat))
  ! **** Forces calculation
  if (ok) then
    snt = sin(bondangle)
    if (abs(snt).lt.zero) snt = zero
    r12inv = 1.0/r12
    r22inv = 1.0/r22
    force = 2.0*bendk*(bondangle-aijk0)/snt
    do kk = 1,3
      fiat(kk) = force*(modval*rjk(kk)-rji(kk)*cst*r12inv)
      fkat(kk) = force*(modval*rji(kk)-rjk(kk)*cst*r22inv)
      fjat(kk) = - fiat(kk) - fkat(kk)
    end do
    if (.not.fixat(iat)) then
      fx(iat) = fx(iat) + fiat(1)
      fy(iat) = fy(iat) + fiat(2)
      fz(iat) = fz(iat) + fiat(3)
    end if
    if (.not.fixat(jat)) then
      fx(jat) = fx(jat) + fjat(1)
      fy(jat) = fy(jat) + fjat(2)
      fz(jat) = fz(jat) + fjat(3)
    end if
    if (.not.fixat(kat)) then
      fx(kat) = fx(kat) + fkat(1)
      fy(kat) = fy(kat) + fkat(2)
      fz(kat) = fz(kat) + fkat(3)
    end if
  end if
end do ! next ibend

return
end subroutine

subroutine ubr
! UB terms for explicit atoms
use explatmod
use nucleotmod
use grandmod
use explatmod

implicit none
! local variables
integer kk,iub,iat,jat,itype
real ubk,rij0,rbij2,rbij,enub,f1,dij(3),fb(3)
logical*1 ok

!     Urey-Bradley term (1,3 distance) .....
!
!     iat --- jat
!         dij

do iub = 1,nubs
  ! force centers
  iat = nsites + ubs(1,iub)
  jat = nsites + ubs(2,iub)
  ! UB type
  itype = ubs(3,iub)
  ! energy constant [Kcal/mole/Angs.**2]
  ubk = ubt(1,itype)
  ! natural distance [Angs.]
  rij0 = ubt(2,itype)
  ! distance
  dij(1) = x(jat) - x(iat)
  dij(2) = y(jat) - y(iat)
  dij(3) = z(jat) - z(iat)
  ! *** Energy contribution
  rbij2 = dot_product(dij,dij)
  rbij = sqrt(rbij2)
  enub = ubk*(rbij-rij0)*(rbij-rij0)
  eub = eub + enub
  ! *** Forces calculation
  ok = Qforces .and. (.not.fixat(iat).or..not.fixat(jat))
  if (ok) then
    f1 = 2.0*ubk*(rbij-rij0)/rbij
    do kk = 1,3
      fb(kk) = f1*dij(kk)
    end do
    if (.not.fixat(iat)) then
      fx(iat) = fx(iat) + fb(1)
      fy(iat) = fy(iat) + fb(2)
      fz(iat) = fz(iat) + fb(3)
    end if
    if (.not.fixat(jat)) then
      fx(jat) = fx(jat) - fb(1)
      fy(jat) = fy(jat) - fb(2)
      fz(jat) = fz(jat) - fb(3)
    end if
  end if
end do ! next iub

return
end subroutine

subroutine en4cen
! Dihedral angle terms for explicit atoms
use explatmod
use nucleotmod
use grandmod
use constamod

implicit none
! local variables
integer itort,i,j,iat,jat,kat,lat,itype
integer nfolds,kk
integer tcmap,pcmap
real Kdih,delta 
real pos1(3),pos2(3),pos3(3),pos4(3)
real rji(3),rjk(3),rlk(3)
real m(3),n(3)
real m2,n2,im2n2,dotmn,acs,dotjin,phi
real f1,entort,nablai(3),nablal(3),nablaval(3)
real rjk2,rjk1,irjk2,im2,in2,djijk,dlkjk
logical*1 ok1,ok2

!     dihedral angle ......
!
!     iat --- jat --- kat --- lat
!     pos1    pos2    pos3    pos4

do itort = 1,ntorts
  ! force centers
  iat = nsites + torts(1,itort) ! terminal atom
  jat = nsites + torts(2,itort)
  kat = nsites + torts(3,itort)
  lat = nsites + torts(4,itort) ! terminal atom
  ! dihedral angle type
  itype = torts(5,itort)
  ! positions
  pos1(1) = x(iat)
  pos1(2) = y(iat)
  pos1(3) = z(iat)
  pos2(1) = x(jat)
  pos2(2) = y(jat)
  pos2(3) = z(jat)
  pos3(1) = x(kat)
  pos3(2) = y(kat)
  pos3(3) = z(kat)
  pos4(1) = x(lat)
  pos4(2) = y(lat)
  pos4(3) = z(lat)
  ! **** Dihedral angle calculation
  ! interparticle vectors
  do kk  =  1,3
    rji(kk) =  pos1(kk) - pos2(kk)
    rjk(kk) =  pos3(kk) - pos2(kk)
    rlk(kk) =  pos3(kk) - pos4(kk) 
  end do
  ! normal vectors
  call cross_product(rji,rjk,m)
  call cross_product(rjk,rlk,n)
  ! dihedral angle
  m2 = dot_product(m,m)
  n2 = dot_product(n,n)
  im2n2 = 1.0/(m2*n2)
  dotmn = dot_product(m,n)*sqrt(im2n2)
  if (dotmn.lt.-1.0) dotmn = -1.0
  if (dotmn.gt.1.0) dotmn = 1.0
  acs = acos(dotmn)
  dotjin = dot_product(rji,n)
  phi = sign(acs,dotjin) ! IUPAC convention
  ok1 = .false.
  ok2 = .false.
  if (Qlcmap) then
    tcmap = lthetacmap(itort)
    pcmap = lpsicmap(itort)
    ok1 = tcmap.gt.0
    ok2 = pcmap.gt.0
    if (ok1) thetacmap(tcmap) = phi
    if (ok2) psicmap(pcmap) = phi 
  endif ! Qlcmap
  if (Qforces) then
    rjk2 = dot_product(rjk,rjk)
    rjk1 = sqrt(rjk2)
    irjk2 = 1.0/rjk2
    im2 = 1.0/m2
    in2 = 1.0/n2
    djijk = dot_product(rji,rjk)
    dlkjk = dot_product(rlk,rjk)
    do kk = 1,3
      nablai(kk) = rjk1*m(kk)*im2
      nablal(kk) = - rjk1*n(kk)*in2
      nablaval(kk) = (djijk*nablai(kk)-dlkjk*nablal(kk))*irjk2
      if (ok1) then
        nablatcmp(kk,1,tcmap) = nablai(kk)
        nablatcmp(kk,2,tcmap) = - nablai(kk) + nablaval(kk)
        nablatcmp(kk,3,tcmap) = - nablal(kk) - nablaval(kk)
        nablatcmp(kk,4,tcmap) = nablal(kk)
      endif
      if (ok2) then
        nablapcmp(kk,1,pcmap) = nablai(kk)
        nablapcmp(kk,2,pcmap) = - nablai(kk) + nablaval(kk)
        nablapcmp(kk,3,pcmap) = - nablal(kk) - nablaval(kk)
        nablapcmp(kk,4,pcmap) = nablal(kk)
      endif
    enddo
  end if
  do i = 1, nprms(itype)
    j = (i-1)*2
    ! energy constants
    nfolds = ndih(i,itype)  
    Kdih = dih(j+1,itype) ! Kcal/mole
    delta = dih(j+2,itype)*radians ! radians
    ! **** Energy contribution
    entort = Kdih*(1.0+cos(nfolds*phi-delta))
    edihe = edihe + entort 
    ! **** Forces calculation
    if (Qforces) then
      f1 = -Kdih*nfolds*sin(nfolds*phi-delta)
      if (.not.fixat(iat)) then
        fx(iat) = fx(iat) - f1*nablai(1)
        fy(iat) = fy(iat) - f1*nablai(2)
        fz(iat) = fz(iat) - f1*nablai(3)
      end if
      if (.not.fixat(jat)) then
        fx(jat) = fx(jat) + f1*nablai(1) - f1*nablaval(1)
        fy(jat) = fy(jat) + f1*nablai(2) - f1*nablaval(2)
        fz(jat) = fz(jat) + f1*nablai(3) - f1*nablaval(3)
      end if
      if (.not.fixat(kat)) then
        fx(kat) = fx(kat) + f1*nablal(1) + f1*nablaval(1)
        fy(kat) = fy(kat) + f1*nablal(2) + f1*nablaval(2)
        fz(kat) = fz(kat) + f1*nablal(3) + f1*nablaval(3)
      end if
      if (.not.fixat(lat)) then
        fx(lat) = fx(lat) - f1*nablal(1)
        fy(lat) = fy(lat) - f1*nablal(2)
        fz(lat) = fz(lat) - f1*nablal(3)
      end if 
    end if
  end do ! next i
end do ! next itort

return
end subroutine

subroutine improper
! Improper angle terms for explicit atoms
use explatmod
use nucleotmod
use grandmod
use constamod

implicit none
! local variables
integer ideform,iat,jat,kat,lat,itype,kk
real oopsk,omega 
real pos1(3),pos2(3),pos3(3),pos4(3)
real rji(3),rjk(3),rlk(3)
real m(3),n(3)
real m2,n2,im2n2,dotmn,acs,dotjin,phi
real f1,enopbs,nablai(3),nablal(3),nablaval(3)
real rjk2,rjk1,irjk2,im2,in2,djijk,dlkjk
logical*1 ok

!     improper angle ......
!
!     jat --- iat --- kat 
!     pos2    pos1(4) pos3
!              | ---- lat
!                     pos4(1)  

do ideform = 1,ndeforms
  ! force centers
  iat = nsites + deforms(1,ideform) 
  jat = nsites + deforms(2,ideform)
  kat = nsites + deforms(3,ideform)
  lat = nsites + deforms(4,ideform) 
  ! improper angle type
  itype = deforms(5,ideform)
  ! energy constant [Kcal/mole/radians**2]
  oopsk = deform(1,itype)
  ! natural improper angle [radians]
  omega = deform(2,itype)*radians
  ! positions
  pos1(1) = x(iat)
  pos1(2) = y(iat)
  pos1(3) = z(iat)
  pos2(1) = x(jat)
  pos2(2) = y(jat)
  pos2(3) = z(jat)
  pos3(1) = x(kat)
  pos3(2) = y(kat)
  pos3(3) = z(kat)
  pos4(1) = x(lat)
  pos4(2) = y(lat)
  pos4(3) = z(lat)
  ! **** Improper angle calculation
  ! interparticle vectors
  do kk  =  1,3
    rji(kk) =  pos1(kk) - pos2(kk)
    rjk(kk) =  pos3(kk) - pos2(kk)
    rlk(kk) =  pos3(kk) - pos4(kk) 
  end do
  ! normal vectors
  call cross_product(rji,rjk,m)
  call cross_product(rjk,rlk,n)
  ! dihedral angle
  m2 = dot_product(m,m)
  n2 = dot_product(n,n)
  im2n2 = 1.0/(m2*n2)
  dotmn = dot_product(m,n)*sqrt(im2n2)
  if (dotmn.lt.-1.0) dotmn = -1.0
  if (dotmn.gt.1.0) dotmn = 1.0
  acs = acos(dotmn)
  dotjin = dot_product(rji,n)
  phi = sign(acs,dotjin) ! IUPAC convention
  ! *** Energy contribution
  enopbs = oopsk*(phi-omega)*(phi-omega)
  eopbs = eopbs + enopbs
  ! *** Forces calculation
  ok = Qforces .and. (.not.fixat(iat).or..not.fixat(jat).or..not.fixat(kat).or..not.fixat(lat))
  if (ok) then
    rjk2 = dot_product(rjk,rjk)
    rjk1 = sqrt(rjk2)
    irjk2 = 1.0/rjk2
    im2 = 1.0/m2
    in2 = 1.0/n2
    djijk = dot_product(rji,rjk)
    dlkjk = dot_product(rlk,rjk)
    do kk = 1,3
      nablai(kk) = rjk1*m(kk)*im2
      nablal(kk) = - rjk1*n(kk)*in2
      nablaval(kk) = (djijk*nablai(kk)-dlkjk*nablal(kk))*irjk2
    enddo
    f1 = 2.0*oopsk*(phi-omega)
    if (.not.fixat(iat)) then
      fx(iat) = fx(iat) - f1*nablai(1)
      fy(iat) = fy(iat) - f1*nablai(2)
      fz(iat) = fz(iat) - f1*nablai(3)
    end if
    if (.not.fixat(jat)) then
      fx(jat) = fx(jat) + f1*nablai(1) - f1*nablaval(1)
      fy(jat) = fy(jat) + f1*nablai(2) - f1*nablaval(2)
      fz(jat) = fz(jat) + f1*nablai(3) - f1*nablaval(3)
    end if
    if (.not.fixat(kat)) then
      fx(kat) = fx(kat) + f1*nablal(1) + f1*nablaval(1)
      fy(kat) = fy(kat) + f1*nablal(2) + f1*nablaval(2)
      fz(kat) = fz(kat) + f1*nablal(3) + f1*nablaval(3)
    end if
    if (.not.fixat(lat)) then
      fx(lat) = fx(lat) - f1*nablal(1)
      fy(lat) = fy(lat) - f1*nablal(2)
      fz(lat) = fz(lat) - f1*nablal(3)
    end if 
  end if
end do ! next ideform

return
end subroutine

subroutine cmapr
! CMAP terms for explicit atoms
use explatmod
use errormod
use nucleotmod
use grandmod

implicit none
! local variables
integer icmap,i,j,k,itheta,ipsi,n,ic
integer itype,iat,jat,kat,lat,mat,nnat,oat,ppat
real theta,psi,dang,idang,thetag,psig,c(4,4)
real t,u,ansy,ansy1,ansy2

!     CMAP terms 

do icmap = 1,ncmaps
  itype = cmaps(3,icmap)
  ! dihedral angles
  theta = thetacmap(icmap)
  psi = psicmap(icmap)
  ! grid spacing
  dang = gscmap(itype)
  idang = 1.0/dang
  ! grid points
  n = cmap(itype)
  itheta = int((theta+180.0)*idang) + 1
  ipsi = int((psi+180.0)*idang) + 1
  if (itheta.gt.n) itheta = 1
  if (ipsi.gt.n) ipsi = 1
  ic = (itheta-1)*n + ipsi
  thetag = -180.0 + (itheta-1)*dang
  psig = -180.0 + (ipsi-1)*dang
  ! coefficient for the bicubic interpolation
  k = 0
  do i = 1,4  
    do j = 1,4
      k = k + 1
      c(i,j) = ccoef(k,ic,itype)
    end do
  end do
  ! *** Calculate function, gradient and cross derivative at the four grid points
  !     of a rectangular grid cell (numbered counter clockwise from the lower left)
  t = (theta-thetag)*idang
  u = (psi-psig)*idang
  ansy = 0.0
  ansy1 = 0.0
  ansy2 = 0.0
  do i = 4,1,-1
    ansy = t*ansy + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
    ansy2 = t*ansy2 + (3.0*c(i,4)*u+2.0*c(i,3))*u + c(i,2)
    ansy1 = u*ansy1 + (3.0*c(4,i)*t+2.0*c(3,i))*t + c(2,i) 
  end do
  ansy1 = ansy1*idang
  ansy2 = ansy2*idang
  ! *** CMAP energy
  ecmap = ecmap + ansy
  ! *** CMAP forces
  if (Qforces) then
    iat  = nsites + attcmap(1,icmap)
    jat  = nsites + attcmap(2,icmap)
    kat  = nsites + attcmap(3,icmap)
    lat  = nsites + attcmap(4,icmap)
    mat  = nsites + atpcmap(1,icmap)
    nnat = nsites + atpcmap(2,icmap)
    oat  = nsites + atpcmap(3,icmap)
    ppat = nsites + atpcmap(4,icmap)
    if (.not.fixat(iat)) then
      fx(iat) = fx(iat) - ansy1*nablatcmp(1,1,icmap)
      fy(iat) = fy(iat) - ansy1*nablatcmp(2,1,icmap)
      fz(iat) = fz(iat) - ansy1*nablatcmp(3,1,icmap)
    end if
    if (.not.fixat(jat)) then
      fx(jat) = fx(jat) - ansy1*nablatcmp(1,2,icmap)
      fy(jat) = fy(jat) - ansy1*nablatcmp(2,2,icmap)
      fz(jat) = fz(jat) - ansy1*nablatcmp(3,2,icmap)
    end if
    if (.not.fixat(kat)) then
      fx(kat) = fx(kat) - ansy1*nablatcmp(1,3,icmap)
      fy(kat) = fy(kat) - ansy1*nablatcmp(2,3,icmap)
      fz(kat) = fz(kat) - ansy1*nablatcmp(3,3,icmap) 
    end if
    if (.not.fixat(lat)) then
      fx(lat) = fx(lat) - ansy1*nablatcmp(1,4,icmap)
      fy(lat) = fy(lat) - ansy1*nablatcmp(2,4,icmap)
      fz(lat) = fz(lat) - ansy1*nablatcmp(3,4,icmap)
    end if
    if (.not.fixat(mat)) then
      fx(mat) = fx(mat) - ansy2*nablapcmp(1,1,icmap)
      fy(mat) = fy(mat) - ansy2*nablapcmp(2,1,icmap)
      fz(mat) = fz(mat) - ansy2*nablapcmp(3,1,icmap)
    end if
    if (.not.fixat(nnat)) then
      fx(nnat) = fx(nnat) - ansy2*nablapcmp(1,2,icmap)
      fy(nnat) = fy(nnat) - ansy2*nablapcmp(2,2,icmap)
      fz(nnat) = fz(nnat) - ansy2*nablapcmp(3,2,icmap)
    end if
    if (.not.fixat(oat)) then
      fx(oat) = fx(oat) - ansy2*nablapcmp(1,3,icmap)
      fy(oat) = fy(oat) - ansy2*nablapcmp(2,3,icmap)
      fz(oat) = fz(oat) - ansy2*nablapcmp(3,3,icmap)
    end if
    if (.not.fixat(ppat)) then
      fx(ppat) = fx(ppat) - ansy2*nablapcmp(1,4,icmap)
      fy(ppat) = fy(ppat) - ansy2*nablapcmp(2,4,icmap)
      fz(ppat) = fz(ppat) - ansy2*nablapcmp(3,4,icmap)
    end if
  end if
end do ! next icmap 

return
end subroutine
