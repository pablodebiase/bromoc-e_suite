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
subroutine charmmenergy(pari,energy)
use listmod
use grandmod
use constamod
implicit none
integer pari,ptypi,sr,ncmaps
real energy
real,allocatable,dimension(:,:,:) :: nablapcmp, nablatcmp
real,allocatable,dimension(:)     :: thetacmap(:), psicmap(:)
! Get Particle Types from Particle List
ptypi=parl(pari)%ptyp
sr=parl(pari)%sr
! Allocate
if (ptypl(ptypi)%psf(1)%Qlcmap) then
  ncmaps=ptypl(ptypi)%psf(1)%ncmaps
  allocate(nablatcmp(3,4,ncmaps),nablapcmp(3,4,ncmaps))
  allocate(thetacmap(ncmaps),psicmap(ncmaps))
endif
energy=0.0
! bonded interaction among internal particles
! Bond terms
call en2cen(energy)
! Angle terms
call en3cen(energy)
! UB terms
call ubr(energy)
! Dihedral angle terms
call en4cen(energy)
! Improper angle terms
call improper(energy)
! CMAP terms
call cmapr(energy)
! Internal non-bonded interactions term
call nonbonded(energy)
! Deallocate
if (ptypl(ptypi)%psf(1)%Qlcmap) then
  deallocate(nablatcmp,nablapcmp)
  deallocate(thetacmap,psicmap)
endif
contains 
  ! Bond terms for explicit atoms
  subroutine en2cen(ener)
  implicit none
  integer ibond,iat,jat,itype
  real bondk,rij0,rbij,de,dij(3),fb(3),rdif,ener
  !     bond .....
  !
  !     iat --- jat
  !         dij
  do ibond = 1,ptypl(ptypi)%psf(1)%nbonds
    ! force centers
    iat = sr + ptypl(ptypi)%psf(1)%bonds(1,ibond)
    jat = sr + ptypl(ptypi)%psf(1)%bonds(2,ibond)
    ! bond type
    itype = ptypl(ptypi)%psf(1)%bonds(3,ibond)
    ! energy constant [Kcal/mol/Angs.**2]
    bondk = ptypl(ptypi)%psf(1)%stretch(1,itype)
    ! natural bond [Angs.]
    rij0 = ptypl(ptypi)%psf(1)%stretch(2,itype)
    ! distance
    dij(1) = r(jat)%x - r(iat)%x
    dij(2) = r(jat)%y - r(iat)%y
    dij(3) = r(jat)%z - r(iat)%z
    ! *** Energy contribution
    rbij = sqrt(dot_product(dij,dij))
    rdif = rbij-rij0
    ener = ener + bondk*rdif*rdif
    ! *** Forces calculation
    if (Qforces) then
      de = 2.0*bondk*rdif/rbij
      fb = de*dij
      f(iat)%x = f(iat)%x + fb(1)
      f(iat)%y = f(iat)%y + fb(2)
      f(iat)%z = f(iat)%z + fb(3)
      f(jat)%x = f(jat)%x - fb(1)
      f(jat)%y = f(jat)%y - fb(2)
      f(jat)%z = f(jat)%z - fb(3)
    end if
  end do ! next ibond
  end subroutine
  
  ! Angle terms for explicit atoms
  subroutine en3cen(ener)
  implicit none
  integer ibend,iat,jat,kat,itype
  real bendk,aijk0,r11,r22,r12,modval,cst
  real bondangle,force,fiat(3),fjat(3),fkat(3)
  real rji(3),rjk(3),adif,ener
  real pos1(3),pos2(3),pos3(3)
  !     bend angle ......
  !
  !     iat --- jat --- kat
  !     pos1    pos2    pos3
  do ibend = 1,ptypl(ptypi)%psf(1)%nbends
    ! force centers
    iat = sr + ptypl(ptypi)%psf(1)%bends(1,ibend)
    jat = sr + ptypl(ptypi)%psf(1)%bends(2,ibend) ! central atom
    kat = sr + ptypl(ptypi)%psf(1)%bends(3,ibend)
    ! bend angle type
    itype  = ptypl(ptypi)%psf(1)%bends(4,ibend)
    ! energy constant [Kcal/mol/radians**2]
    bendk  = ptypl(ptypi)%psf(1)%bend(1,itype)
    ! natural bending angle [radians]
    aijk0  = ptypl(ptypi)%psf(1)%bend(2,itype)*radians
    ! positions
    pos1(1) = r(iat)%x
    pos1(2) = r(iat)%y
    pos1(3) = r(iat)%z
    pos2(1) = r(jat)%x
    pos2(2) = r(jat)%y
    pos2(3) = r(jat)%z
    pos3(1) = r(kat)%x
    pos3(2) = r(kat)%y
    pos3(3) = r(kat)%z
    ! *** Bond angle calculation
    ! interparticle vectors
    rji =  pos1 - pos2
    rjk =  pos3 - pos2
    ! bond angle
    r11 = dot_product(rji,rji)
    r22 = dot_product(rjk,rjk)
    r12 = dot_product(rji,rjk)
    modval = 1.0/sqrt(r11*r22)
    cst=r12*modval
    if (cst.lt.-1.0) cst = -1.0
    if (cst.gt.1.0) cst = 1.0
    bondangle = acos(cst)
    adif=bondangle-aijk0
    ! **** Energy contribution
    ener = ener + bendk*adif**2
    ! **** Forces calculation
    if (Qforces) then
      force = 2.0*bendk*adif*modval/nonzero(sin(bondangle))
      fiat=force*(rjk-(r12/r11)*rji)
      fkat=force*(rji-(r12/r22)*rjk)
      fjat=-(fiat+fkat)
      f(iat)%x = f(iat)%x + fiat(1)
      f(iat)%y = f(iat)%y + fiat(2)
      f(iat)%z = f(iat)%z + fiat(3)
      f(jat)%x = f(jat)%x + fjat(1)
      f(jat)%y = f(jat)%y + fjat(2)
      f(jat)%z = f(jat)%z + fjat(3)
      f(kat)%x = f(kat)%x + fkat(1)
      f(kat)%y = f(kat)%y + fkat(2)
      f(kat)%z = f(kat)%z + fkat(3)
    end if
  end do ! next ibend
  end subroutine
  
  ! UB terms for explicit atoms
  subroutine ubr(ener)
  implicit none
  integer iub,iat,jat,itype
  real ubk,rij0,rbij,rdif,enub,f1,dij(3),fb(3),ener
  !     Urey-Bradley term (1,3 distance) .....
  !
  !     iat --- jat
  !         dij
  do iub = 1,ptypl(ptypi)%psf(1)%nubs
    ! force centers
    iat = sr + ptypl(ptypi)%psf(1)%ubs(1,iub)
    jat = sr + ptypl(ptypi)%psf(1)%ubs(2,iub)
    ! UB type
    itype = ptypl(ptypi)%psf(1)%ubs(3,iub)
    ! energy constant [Kcal/mol/Angs.**2]
    ubk = ptypl(ptypi)%psf(1)%ubt(1,itype)
    ! natural distance [Angs.]
    rij0 = ptypl(ptypi)%psf(1)%ubt(2,itype)
    ! distance
    dij(1) = r(jat)%x - r(iat)%x
    dij(2) = r(jat)%y - r(iat)%y
    dij(3) = r(jat)%z - r(iat)%z
    ! *** Energy contribution
    rbij = sqrt(dot_product(dij,dij))
    rdif = rbij-rij0
    enub = ubk*rdif*rdif
    ener = ener + enub
    ! *** Forces calculation
    if (Qforces) then
      f1 = 2.0*ubk*rdif/rbij
      fb = f1*dij
      f(iat)%x = f(iat)%x + fb(1)
      f(iat)%y = f(iat)%y + fb(2)
      f(iat)%z = f(iat)%z + fb(3)
      f(jat)%x = f(jat)%x - fb(1)
      f(jat)%y = f(jat)%y - fb(2)
      f(jat)%z = f(jat)%z - fb(3)
    end if
  end do ! next iub
  end subroutine
  
  ! Dihedral angle terms for explicit atoms
  subroutine en4cen(ener)
  implicit none
  integer itort,i,j,iat,jat,kat,lat,itype
  integer nfolds
  integer tcmap,pcmap
  real Kdih,delta 
  real pos1(3),pos2(3),pos3(3),pos4(3)
  real rji(3),rjk(3),rlk(3)
  real m(3),n(3),ener
  real m2,n2,im2n2,dotmn,acs,dotjin,phi
  real f1,entort,nablai(3),nablal(3),nablaval(3)
  real rjk2,rjk1,irjk2,im2,in2,djijk,dlkjk
  logical*1 ok1,ok2
  !     dihedral angle ......
  !
  !     iat --- jat --- kat --- lat
  !     pos1    pos2    pos3    pos4
  do itort = 1,ptypl(ptypi)%psf(1)%ntorts
    ! force centers
    iat = sr + ptypl(ptypi)%psf(1)%torts(1,itort) ! terminal atom
    jat = sr + ptypl(ptypi)%psf(1)%torts(2,itort)
    kat = sr + ptypl(ptypi)%psf(1)%torts(3,itort)
    lat = sr + ptypl(ptypi)%psf(1)%torts(4,itort) ! terminal atom
    ! dihedral angle type
    itype = ptypl(ptypi)%psf(1)%torts(5,itort)
    ! positions
    pos1(1) = r(iat)%x
    pos1(2) = r(iat)%y
    pos1(3) = r(iat)%z
    pos2(1) = r(jat)%x
    pos2(2) = r(jat)%y
    pos2(3) = r(jat)%z
    pos3(1) = r(kat)%x
    pos3(2) = r(kat)%y
    pos3(3) = r(kat)%z
    pos4(1) = r(lat)%x
    pos4(2) = r(lat)%y
    pos4(3) = r(lat)%z
    ! **** Dihedral angle calculation
    ! interparticle vectors
    rji = pos1 - pos2
    rjk = pos3 - pos2
    rlk = pos3 - pos4 
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
    if (ptypl(ptypi)%psf(1)%Qlcmap) then
      tcmap = ptypl(ptypi)%psf(1)%lthetacmap(itort)
      pcmap = ptypl(ptypi)%psf(1)%lpsicmap(itort)
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
      nablai = rjk1*m*im2
      nablal = - rjk1*n*in2
      nablaval = (djijk*nablai-dlkjk*nablal)*irjk2
      if (ok1) then
        nablatcmp(:,1,tcmap) = nablai
        nablatcmp(:,2,tcmap) = - nablai + nablaval
        nablatcmp(:,3,tcmap) = - nablal - nablaval
        nablatcmp(:,4,tcmap) = nablal
      endif
      if (ok2) then
        nablapcmp(:,1,pcmap) = nablai
        nablapcmp(:,2,pcmap) = -nablai + nablaval
        nablapcmp(:,3,pcmap) = -nablal - nablaval
        nablapcmp(:,4,pcmap) = nablal
      endif
    end if
    do i = 1, ptypl(ptypi)%psf(1)%nprms(itype)
      j = (i-1)*2
      ! energy constants
      nfolds = ptypl(ptypi)%psf(1)%ndih(i,itype)
      Kdih = ptypl(ptypi)%psf(1)%dih(j+1,itype) ! Kcal/mol
      delta = ptypl(ptypi)%psf(1)%dih(j+2,itype)*radians ! radians
      ! **** Energy contribution
      entort = Kdih*(1.0+cos(nfolds*phi-delta))
      ener = ener + entort 
      ! **** Forces calculation
      if (Qforces) then
        f1 = -Kdih*nfolds*sin(nfolds*phi-delta)
        f(iat)%x = f(iat)%x - f1*nablai(1)
        f(iat)%y = f(iat)%y - f1*nablai(2)
        f(iat)%z = f(iat)%z - f1*nablai(3)
        f(jat)%x = f(jat)%x + f1*nablai(1) - f1*nablaval(1)
        f(jat)%y = f(jat)%y + f1*nablai(2) - f1*nablaval(2)
        f(jat)%z = f(jat)%z + f1*nablai(3) - f1*nablaval(3)
        f(kat)%x = f(kat)%x + f1*nablal(1) + f1*nablaval(1)
        f(kat)%y = f(kat)%y + f1*nablal(2) + f1*nablaval(2)
        f(kat)%z = f(kat)%z + f1*nablal(3) + f1*nablaval(3)
        f(lat)%x = f(lat)%x - f1*nablal(1)
        f(lat)%y = f(lat)%y - f1*nablal(2)
        f(lat)%z = f(lat)%z - f1*nablal(3)
      end if
    end do ! next i
  end do ! next itort
  
  return
  end subroutine
  
  ! Improper angle terms for explicit atoms
  subroutine improper(ener)
  implicit none
  integer ideform,iat,jat,kat,lat,itype
  real oopsk,omega 
  real pos1(3),pos2(3),pos3(3),pos4(3)
  real rji(3),rjk(3),rlk(3)
  real m(3),n(3),ener
  real m2,n2,im2n2,dotmn,acs,dotjin,phi
  real f1,enopbs,nablai(3),nablal(3),nablaval(3)
  real rjk2,rjk1,irjk2,im2,in2,djijk,dlkjk
  !     improper angle ......
  !     jat --- iat --- kat 
  !     pos2    pos1(4) pos3
  !              | ---- lat
  !                     pos4(1) 
  do ideform = 1,ptypl(ptypi)%psf(1)%ndeforms
    ! force centers
    iat = sr + ptypl(ptypi)%psf(1)%deforms(1,ideform) 
    jat = sr + ptypl(ptypi)%psf(1)%deforms(2,ideform)
    kat = sr + ptypl(ptypi)%psf(1)%deforms(3,ideform)
    lat = sr + ptypl(ptypi)%psf(1)%deforms(4,ideform) 
    ! improper angle type
    itype = ptypl(ptypi)%psf(1)%deforms(5,ideform)
    ! energy constant [Kcal/mol/radians**2]
    oopsk = ptypl(ptypi)%psf(1)%deform(1,itype)
    ! natural improper angle [radians]
    omega = ptypl(ptypi)%psf(1)%deform(2,itype)*radians
    ! positions
    pos1(1) = r(iat)%x
    pos1(2) = r(iat)%y
    pos1(3) = r(iat)%z
    pos2(1) = r(jat)%x
    pos2(2) = r(jat)%y
    pos2(3) = r(jat)%z
    pos3(1) = r(kat)%x
    pos3(2) = r(kat)%y
    pos3(3) = r(kat)%z
    pos4(1) = r(lat)%x
    pos4(2) = r(lat)%y
    pos4(3) = r(lat)%z
    ! **** Improper angle calculation
    ! interparticle vectors
    rji =  pos1 - pos2
    rjk =  pos3 - pos2
    rlk =  pos3 - pos4 
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
    ener = ener + enopbs
    ! *** Forces calculation
    if (Qforces) then
      rjk2 = dot_product(rjk,rjk)
      rjk1 = sqrt(rjk2)
      irjk2 = 1.0/rjk2
      im2 = 1.0/m2
      in2 = 1.0/n2
      djijk = dot_product(rji,rjk)
      dlkjk = dot_product(rlk,rjk)
      nablai = rjk1*m*im2
      nablal = -rjk1*n*in2
      nablaval = (djijk*nablai-dlkjk*nablal)*irjk2
      f1 = 2.0*oopsk*(phi-omega)
      f(iat)%x = f(iat)%x - f1*nablai(1)
      f(iat)%y = f(iat)%y - f1*nablai(2)
      f(iat)%z = f(iat)%z - f1*nablai(3)
      f(jat)%x = f(jat)%x + f1*nablai(1) - f1*nablaval(1)
      f(jat)%y = f(jat)%y + f1*nablai(2) - f1*nablaval(2)
      f(jat)%z = f(jat)%z + f1*nablai(3) - f1*nablaval(3)
      f(kat)%x = f(kat)%x + f1*nablal(1) + f1*nablaval(1)
      f(kat)%y = f(kat)%y + f1*nablal(2) + f1*nablaval(2)
      f(kat)%z = f(kat)%z + f1*nablal(3) + f1*nablaval(3)
      f(lat)%x = f(lat)%x - f1*nablal(1)
      f(lat)%y = f(lat)%y - f1*nablal(2)
      f(lat)%z = f(lat)%z - f1*nablal(3)
    end if
  end do ! next ideform
  end subroutine
  
  ! CMAP terms for explicit atoms
  subroutine cmapr(ener)
  implicit none
  integer icmap,i,j,k,itheta,ipsi,n,ic
  integer itype,iat,jat,kat,lat,mat,nnat,oat,ppat
  real theta,psi,dang,idang,thetag,psig,c(4,4)
  real t,u,ansy,ansy1,ansy2,ener
  !     CMAP terms 
  do icmap = 1,ptypl(ptypi)%psf(1)%ncmaps
    itype = ptypl(ptypi)%psf(1)%cmaps(3,icmap)
    ! dihedral angles
    theta = thetacmap(icmap)
    psi = psicmap(icmap)
    ! grid spacing
    dang = ptypl(ptypi)%psf(1)%gscmap(itype)
    idang = 1.0/dang
    ! grid points
    n = ptypl(ptypi)%psf(1)%cmap(itype)
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
        c(i,j) = ptypl(ptypi)%psf(1)%ccoef(k,ic,itype)
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
    ener = ener + ansy
    ! *** CMAP forces
    if (Qforces) then
      iat  = sr + ptypl(ptypi)%psf(1)%attcmap(1,icmap)
      jat  = sr + ptypl(ptypi)%psf(1)%attcmap(2,icmap)
      kat  = sr + ptypl(ptypi)%psf(1)%attcmap(3,icmap)
      lat  = sr + ptypl(ptypi)%psf(1)%attcmap(4,icmap)
      mat  = sr + ptypl(ptypi)%psf(1)%atpcmap(1,icmap)
      nnat = sr + ptypl(ptypi)%psf(1)%atpcmap(2,icmap)
      oat  = sr + ptypl(ptypi)%psf(1)%atpcmap(3,icmap)
      ppat = sr + ptypl(ptypi)%psf(1)%atpcmap(4,icmap)
      f(iat)%x  = f(iat)%x  - ansy1*nablatcmp(1,1,icmap)
      f(iat)%y  = f(iat)%y  - ansy1*nablatcmp(2,1,icmap)
      f(iat)%z  = f(iat)%z  - ansy1*nablatcmp(3,1,icmap)
      f(jat)%x  = f(jat)%x  - ansy1*nablatcmp(1,2,icmap)
      f(jat)%y  = f(jat)%y  - ansy1*nablatcmp(2,2,icmap)
      f(jat)%z  = f(jat)%z  - ansy1*nablatcmp(3,2,icmap)
      f(kat)%x  = f(kat)%x  - ansy1*nablatcmp(1,3,icmap)
      f(kat)%y  = f(kat)%y  - ansy1*nablatcmp(2,3,icmap)
      f(kat)%z  = f(kat)%z  - ansy1*nablatcmp(3,3,icmap) 
      f(lat)%x  = f(lat)%x  - ansy1*nablatcmp(1,4,icmap)
      f(lat)%y  = f(lat)%y  - ansy1*nablatcmp(2,4,icmap)
      f(lat)%z  = f(lat)%z  - ansy1*nablatcmp(3,4,icmap)
      f(mat)%x  = f(mat)%x  - ansy2*nablapcmp(1,1,icmap)
      f(mat)%y  = f(mat)%y  - ansy2*nablapcmp(2,1,icmap)
      f(mat)%z  = f(mat)%z  - ansy2*nablapcmp(3,1,icmap)
      f(nnat)%x = f(nnat)%x - ansy2*nablapcmp(1,2,icmap)
      f(nnat)%y = f(nnat)%y - ansy2*nablapcmp(2,2,icmap)
      f(nnat)%z = f(nnat)%z - ansy2*nablapcmp(3,2,icmap)
      f(oat)%x  = f(oat)%x  - ansy2*nablapcmp(1,3,icmap)
      f(oat)%y  = f(oat)%y  - ansy2*nablapcmp(2,3,icmap)
      f(oat)%z  = f(oat)%z  - ansy2*nablapcmp(3,3,icmap)
      f(ppat)%x = f(ppat)%x - ansy2*nablapcmp(1,4,icmap)
      f(ppat)%y = f(ppat)%y - ansy2*nablapcmp(2,4,icmap)
      f(ppat)%z = f(ppat)%z - ansy2*nablapcmp(3,4,icmap)
    end if
  end do ! next icmap 
  end subroutine
  
  ! Internal non-bonded interactions term
  subroutine nonbonded(ener)
  implicit none
  integer i,j,k,n,a,b
  real qa,qb,epp4,sgp2,elecpsf,evdwpsf,idist2,idist,dist2,dist6,dist12,de,ener
  ! Compute Nonbonded for 1-4 Pairs
  n=ptypl(ptypi)%psf(1)%np14
  do k=1,n
    a=ptypl(ptypi)%psf(1)%p14(k)%a
    b=ptypl(ptypi)%psf(1)%p14(k)%b
    qa=ptypl(ptypi)%chg(a)
    qb=ptypl(ptypi)%chg(b)
    epp4=ptypl(ptypi)%psf(1)%lj14(k)%epp4
    sgp2=ptypl(ptypi)%psf(1)%lj14(k)%sgp2
    dist2=dist2car(r(sr+a),r(sr+b))
    idist2=1.0/dist2
    idist=sqrt(idist2)
    elecpsf=celec*qa*qb*idist
    dist6=(sgp2*idist2)**3
    dist12=dist6**2
    evdwpsf=epp4*(dist12-dist6) ! van der waals potential
    ener = ener + elecpsf + evdwpsf
    if (Qforces) then
      i=sr+a
      j=sr+b
      de=elecpsf*idist2+epp4*(2.0*dist12-dist6)*6.0*idist2 ! electrostatic & van der waals forces
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
  ! Compute nonbonded for rest of pairs
  n=ptypl(ptypi)%psf(1)%nnbon
  do k=1,n
    a=ptypl(ptypi)%psf(1)%nbon(k)%a
    b=ptypl(ptypi)%psf(1)%nbon(k)%b
    qa=ptypl(ptypi)%chg(a)
    qb=ptypl(ptypi)%chg(b)
    epp4=ptypl(ptypi)%psf(1)%lj(k)%epp4
    sgp2=ptypl(ptypi)%psf(1)%lj(k)%sgp2
    dist2=dist2car(r(sr+a),r(sr+b))
    idist2=1.0/dist2
    idist=sqrt(idist2)
    elecpsf=celec*qa*qb*idist
    dist6=(sgp2*idist2)**3
    dist12=dist6**2
    evdwpsf=epp4*(dist12-dist6) ! van der waals potential
    ener = ener + elecpsf + evdwpsf
    if (Qforces) then
      i=sr+a
      j=sr+b
      de=elecpsf*idist2+epp4*(2.0*dist12-dist6)*6.0*idist2 ! electrostatic & van der waals forces
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
  end subroutine

  real function nonzero(num)
  real num
  real,parameter :: small=1e-15
  if (abs(num).lt.small) then
    nonzero=small
  else
    nonzero=num
  endif
  end function
end subroutine
