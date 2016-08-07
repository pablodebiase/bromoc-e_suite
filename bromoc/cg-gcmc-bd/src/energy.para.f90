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
use apfmod
use ioxmod
use grandmod
use gsbpmod
use constamod
use nucleotmod
use extramod
use efpmod
use listmod
implicit none

integer i, j, itype, jtype, is
real dist, dist2, dist6, idist, idist2, idistkd
real de, dc

real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
integer isite1, isite2, isite3, isite4
real  dd,vard,eefp,fdf,fdv,vard2
real  ang, ang0, varang, cte1
real  dihe, dihe0, vardihe,didstmp
real  epsln,sgex2,f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),m1,m2,m3,modval
real  dist12,esolve,n1(3),n2(3),n1v,n2v,v22,v2v,av,bv,iv22,v12,v23
logical*1 ok,ok2
real bb
integer ini,fin,aa,nth,tid,omp_get_thread_num,omp_get_num_threads
real eefpotloc,eelecloc,evdwloc,esrpmfloc,fxloc(nele),fyloc(nele),fzloc(nele)
real ebondloc,eangloc,ediheloc,estackloc,ebploc,eexloc,eqqloc,esolvloc
real eefpotmxloc,eqqmxloc,evdwmxloc,esrpmfmxloc

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
egsbpa   = 0.0
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


if (Qenergy) then
  ! Initializations
  if (Qforces) then
    f(1:nele)%x = 0.0
    f(1:nele)%y = 0.0
    f(1:nele)%z = 0.0
  endif

  ! membrane contribution
  if (Qmemb) call membrane
  ener = ememb

  ! reaction field contribution
  if (Qmmij) then
    if(shapes.EQ.'RECTBOX ') then
      call rect_rf1
    else if(shapes.eq.'SPHERE  ') then
      call sphe_rf1
    endif
    ener = ener + egsbpb
  endif

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

  dids(1:5,nelenuc+1:nele)=0.0
  aa=(nele-nelenuc)*(nele-nelenuc-1)/2
  !$omp parallel private(ang,ang0,av,bb,bv,cofo,cte1,dc,dd,de,didstmp,dihe,dihe0,dist,dist12,dist2,dist6,eangloc,ebondloc,ebploc,ediheloc,eefp,eefpotloc,eefpotmxloc,eelecloc,eexloc,epsln,eqqloc,eqqmxloc,esolvloc,esrpmf0,esrpmf1,esrpmf2,esrpmf3,esrpmfloc,esrpmfmxloc,estackloc,evdwloc,evdwmxloc,f1,f2,f3,f4,fdf,fdv,fin,fxloc,fyloc,fzloc,i,idist,idist2,idistkd,ini,is,isite1,isite2,isite3,isite4,itype,iv22,j,jtype,m1,m2,m3,modval,n1,n1v,n2,n2v,nth,ok,ok2,Qdeby,sgex2,tid,v1,v12,v2,v22,v23,v2v,v3,varang,vard,vard2,vardihe)
  ebploc=0.0
  eexloc=0.0
  eqqloc=0.0
  esolvloc=0.0
  estackloc=0.0
  ediheloc=0.0
  eangloc=0.0
  ebondloc=0.0
  eefpotmxloc=0.0
  eqqmxloc=0.0
  evdwmxloc=0.0
  esrpmfmxloc=0.0
  eefpotloc=0.0
  eelecloc=0.0
  evdwloc=0.0
  esrpmfloc=0.0
  fxloc(1:nele)=0.0
  fyloc(1:nele)=0.0
  fzloc(1:nele)=0.0
  tid = omp_get_thread_num()
  nth = omp_get_num_threads()
  bb=float(aa)/float(nth)
  if (tid.eq.0) then
    ini=1+nelenuc
  else
    ini=int(sqrt(0.25+2.0*bb*tid)+0.5+1+nelenuc)
  endif
  if (tid+1.eq.nth) then
    fin=nele
  else
    fin=int(sqrt(0.25+2.0*bb*(tid+1))+0.5+nelenuc)
  endif

  ! nonbonded interaction between ions
  if (Qpar) then                 
    if (Qnonbond) then
      do j=ini,fin
!     do j = nelenuc+1, nele
        jtype  = et(j)
        do i = nelenuc+1, j-1
          itype = et(i)
          is=etpidx(itype,jtype)
          dist2 = ((r(i)%x-r(j)%x)**2+(r(i)%y-r(j)%y)**2+(r(i)%z-r(j)%z)**2)
          if (Qefpot(is)) then
            if (Qforces) then
              call getyd(is,dist2,eefp,de,dist)
            else
              call gety(is,dist2,eefp,dist)
            endif
            eefpotloc = eefpotloc + eefp
!            eefpot = eefpot + eefp
          else 
            idist2 = 1.0/dist2
            if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
            if (Qchr(is)) then
              cofo=fct(is)*idist
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
                fxloc(i) = fxloc(i) + de*(r(i)%x-r(j)%x)
                fyloc(i) = fyloc(i) + de*(r(i)%y-r(j)%y)
                fzloc(i) = fzloc(i) + de*(r(i)%z-r(j)%z)
              endif
              fxloc(j) = fxloc(j) - de*(r(i)%x-r(j)%x)
              fyloc(j) = fyloc(j) - de*(r(i)%y-r(j)%y)
              fzloc(j) = fzloc(j) - de*(r(i)%z-r(j)%z)
            endif
          endif
        enddo ! i = nelenuc+1,...,(j-1)
      enddo ! j = nelenuc+1,...,nele
    endif !Qnonbond
  endif  ! Qpar

  ! interaction between interaction sites (nucleotides)         
  if (Qnucl) then
    if (Qnobond) then      
  !  Bonded interactions  
  !  -------------------
  !  The strech energy is calculated by means of a 
  !  Taylor expansion around natural bond lenght.
  !  This contribution contains harmonic and anharmonic 
  !  interactions.
      !$omp do
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
        ebondloc = ebondloc + epsnuc*vard2*(1.0+100.0*vard2)
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = epsnuc*vard*(2.0+400.0*vard2)/dist
          f1=de*v1
          if (stfree(isite1)) then   
            fxloc(isite1) = fxloc(isite1) + f1(1)
            fyloc(isite1) = fyloc(isite1) + f1(2)
            fzloc(isite1) = fzloc(isite1) + f1(3)
          endif
          if (stfree(isite2)) then
            fxloc(isite2) = fxloc(isite2) - f1(1)
            fyloc(isite2) = fyloc(isite2) - f1(2)
            fzloc(isite2) = fzloc(isite2) - f1(3)
          endif 
        endif
      enddo   
      !$omp end do
  !  The bending energy is calculated 
  !  by means of a Taylor expansion around natural bond
  !  angle. This contribution contains harmonic interactions. 

      !$omp do
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
        eangloc = eangloc + 700.0*epsnuc*varang**2
        ok = Qforces .and. stfree(isite1).or.stfree(isite2).or.stfree(isite3)
        if (ok) then
          de = 1400.0*epsnuc*varang*modval/sin(ang)
          f1=de*v1
          f2=de*v2
          if (stfree(isite1)) then
            fxloc(isite1) = fxloc(isite1) + f2(1)
            fyloc(isite1) = fyloc(isite1) + f2(2)
            fzloc(isite1) = fzloc(isite1) + f2(3)
          endif
          if (stfree(isite2)) then
            fxloc(isite2) = fxloc(isite2) - (f1(1)+f2(1))
            fyloc(isite2) = fyloc(isite2) - (f1(2)+f2(2))
            fzloc(isite2) = fzloc(isite2) - (f1(3)+f2(3))
          endif
          if (stfree(isite3)) then
            fxloc(isite3) = fxloc(isite3) + f1(1)
            fyloc(isite3) = fyloc(isite3) + f1(2)
            fzloc(isite3) = fzloc(isite3) + f1(3)
          endif
        endif
      enddo  
      !$omp end do
   ! The torsional energy is calculated 
   ! using Fourier series for natural torsional angle. 

      !$omp do
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
        ediheloc = ediheloc + 28.0*epsnuc*(1.0-cos(vardihe))
!        edihe = edihe + 28.0*epsnuc*(1.0-cos(vardihe))
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
            fxloc(isite1)=fxloc(isite1) + f1(1)
            fyloc(isite1)=fyloc(isite1) + f1(2)
            fzloc(isite1)=fzloc(isite1) + f1(3)
          endif
          if (stfree(isite2)) then
            fxloc(isite2) = fxloc(isite2) + f2(1)
            fyloc(isite2) = fyloc(isite2) + f2(2)
            fzloc(isite2) = fzloc(isite2) + f2(3)
          endif
          if (stfree(isite3)) then
            fxloc(isite3) = fxloc(isite3) + f3(1)
            fyloc(isite3) = fyloc(isite3) + f3(2)
            fzloc(isite3) = fzloc(isite3) + f3(3)
          endif
          if (stfree(isite4)) then
            fxloc(isite4) = fxloc(isite4) + f4(1)
            fyloc(isite4) = fyloc(isite4) + f4(2)
            fzloc(isite4) = fzloc(isite4) + f4(3)
          endif
        endif
      enddo
      !$omp end do
    endif  ! Qnobond 
    if (Qnonbond) then
    !Non-bonded interactions
    !-----------------------
    !Native contacts
      !$omp do
      do i = 1, nstack
        isite1 = sitestack(i,1)
        isite2 = sitestack(i,2)
        v1(1)=r(isite2)%x-r(isite1)%x
        v1(2)=r(isite2)%y-r(isite1)%y
        v1(3)=r(isite2)%z-r(isite1)%z
        dist2 = dot_product(v1,v1)
        idist2=1.0/dist2
        cte1 = (sgstack(i)**2*idist2)**3   
        estackloc = estackloc + 4.0*epsnuc*cte1*(cte1-1.0)
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
          f1=de*v1 
          if (stfree(isite1)) then
            fxloc(isite1) = fxloc(isite1) - f1(1)
            fyloc(isite1) = fyloc(isite1) - f1(2)
            fzloc(isite1) = fzloc(isite1) - f1(3)
          endif                            
          if (stfree(isite2)) then         
            fxloc(isite2) = fxloc(isite2) + f1(1)
            fyloc(isite2) = fyloc(isite2) + f1(2)
            fzloc(isite2) = fzloc(isite2) + f1(3)
          endif
        endif  
      enddo 
      !$omp end do
 
    ! Hydrogen bonding
      !$omp do
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
        ebploc = ebploc + epsln*cte1**5*(20.0*cte1-24.0) 
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = 240.0*epsln*cte1**5*(cte1-1.0)*idist2
          f1=de*v1 
          if (stfree(isite1)) then
            fxloc(isite1) = fxloc(isite1) - f1(1)
            fyloc(isite1) = fyloc(isite1) - f1(2)
            fzloc(isite1) = fzloc(isite1) - f1(3)
          endif                            
          if (stfree(isite2)) then         
            fxloc(isite2) = fxloc(isite2) + f1(1)
            fyloc(isite2) = fyloc(isite2) + f1(2)
            fzloc(isite2) = fzloc(isite2) + f1(3)
          endif
        endif                
      enddo
      !$omp end do

    ! Excluded volume
      !$omp do
      do i = 1, nex
        isite1 = siteex(i,1)
        isite2 = siteex(i,2)
        v1(1)=r(isite2)%x-r(isite1)%x
        v1(2)=r(isite2)%y-r(isite1)%y
        v1(3)=r(isite2)%z-r(isite1)%z
        dist2 = dot_product(v1,v1)
        idist2=1.0/dist2
        sgex2=sger(i)%x**2
        if (dist2.lt.sgex2) then
          cte1 = (sgex2*idist2)**3      
          eexloc = eexloc + 4.0*epsnuc*cte1*(cte1-1.0) + epsnuc
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then           
            de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
            f1=de*v1 
            if (stfree(isite1)) then
              fxloc(isite1) = fxloc(isite1) - f1(1)
              fyloc(isite1) = fyloc(isite1) - f1(2)
              fzloc(isite1) = fzloc(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fxloc(isite2) = fxloc(isite2) + f1(1)
              fyloc(isite2) = fyloc(isite2) + f1(2)
              fzloc(isite2) = fzloc(isite2) + f1(3)
            endif
          endif
        endif
      enddo
      !$omp end do
 
    ! Coulomb interaction
      !$omp do
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
          eqqloc = eqqloc + fctn*idist*idistkd
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = fctn*(dist+kappa)*idist**3*ikappa*idistkd
            f1=de*v1
            if (stfree(isite1)) then
              fxloc(isite1) = fxloc(isite1) - f1(1)
              fyloc(isite1) = fyloc(isite1) - f1(2)
              fzloc(isite1) = fzloc(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fxloc(isite2) = fxloc(isite2) + f1(1)
              fyloc(isite2) = fyloc(isite2) + f1(2)
              fzloc(isite2) = fzloc(isite2) + f1(3)
            endif
          endif      
        else
          eqqloc = eqqloc + fctn*idist
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = fctn*idist**3
            f1=de*v1
            if (stfree(isite1)) then
              fxloc(isite1) = fxloc(isite1) - f1(1)
              fyloc(isite1) = fyloc(isite1) - f1(2)
              fzloc(isite1) = fzloc(isite1) - f1(3)
            endif
            if (stfree(isite2)) then
              fxloc(isite2) = fxloc(isite2) + f1(1)
              fyloc(isite2) = fyloc(isite2) + f1(2)
              fzloc(isite2) = fzloc(isite2) + f1(3)
            endif
          endif
        endif
        if (Qdebyhyb) Qdeby=.false.
      enddo 
      !$omp end do

      !     Solvent-induced contribution
      if (Qsolv) then
        !$omp do
        do i = 1, nsolv
          isite1 = siteslv(i,1)
          isite2 = siteslv(i,2)
          v1(1)=r(isite2)%x-r(isite1)%x
          v1(2)=r(isite2)%y-r(isite1)%y
          v1(3)=r(isite2)%z-r(isite1)%z
          dist=sqrt(dot_product(v1,v1)) 
          cte1 = exp((13.38-dist)*0.1875)
          esolvloc = esolvloc + epsolv*(1.0-cte1)**2 - epsolv
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = 0.375*epsolv*cte1*(cte1-1.0)/dist
            f1=de*v1
            if (stfree(isite1)) then
              fxloc(isite1) = fxloc(isite1) - f1(1)
              fyloc(isite1) = fyloc(isite1) - f1(2)
              fzloc(isite1) = fzloc(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fxloc(isite2) = fxloc(isite2) + f1(1)
              fyloc(isite2) = fyloc(isite2) + f1(2)
              fzloc(isite2) = fzloc(isite2) + f1(3)
            endif
          endif
        enddo
      endif 
    endif ! Qnonbond
    if (tid.eq.0) then 
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
          if (kr(i)%x.ne.0.0) then
            if (j.eq.0) then
              if (Qunsplit) then
                xcon=sum(x(1:nelenuc1st))*inelenuc*2.0-contrr(i)%x
                f(1:nelenuc1st)%x=f(1:nelenuc1st)%x-kr(i)%x*xcon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%x*xcon**2
                xcon=sum(x(1+nelenuc1st:nelenuc))*inelenuc*2.0-contrx(ctn+1)
                f(1+nelenuc1st:nelenuc)%x=f(1+nelenuc1st:nelenuc)%x-kr(i)%x*xcon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%x*xcon**2
              else
                xcon=sum(x(1:nelenuc))*inelenuc-contrr(i)%x
                f(1:nelenuc)%x=f(1:nelenuc)%x-kr(i)%x*xcon*inelenuc
                econ = econ + 0.5*kr(i)%x*xcon**2
              endif
            else
              f(j)%x=f(j)%x-kr(i)%x*(r(j)%x-contrr(i)%x)
              econ = econ + 0.5*kr(i)%x*(r(j)%x-contrr(i)%x)**2
            endif 
          endif
          if (kr(i)%y.ne.0.0) then
            if (j.eq.0) then
              if (Qunsplit) then
                ycon=sum(y(1:nelenuc1st))*inelenuc*2.0-contrr(i)%y
                f(1:nelenuc1st)%y=f(1:nelenuc1st)%y-kr(i)%y*ycon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%y*ycon**2
                ycon=sum(y(1+nelenuc1st:nelenuc))*inelenuc*2.0-contry(ctn+1)
                f(1+nelenuc1st:nelenuc)%y=f(1+nelenuc1st:nelenuc)%y-kr(i)%y*ycon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%y*ycon**2
              else
                ycon=sum(y(1:nelenuc))*inelenuc-contrr(i)%y
                f(1:nelenuc)%y=f(1:nelenuc)%y-kr(i)%y*ycon*inelenuc
                econ = econ + 0.5*kr(i)%y*ycon**2
              endif
            else
              f(j)%y=f(j)%y-kr(i)%y*(r(j)%y-contrr(i)%y)
              econ = econ + 0.5*kr(i)%y*(r(j)%y-contrr(i)%y)**2
            endif
          endif
          if (kr(i)%z.ne.0.0) then
            if (j.eq.0) then
              if (Qunsplit) then
                zcon=sum(z(1:nelenuc1st))*inelenuc*2.0-contrr(i)%z
                f(1:nelenuc1st)%z=f(1:nelenuc1st)%z-kr(i)%z*zcon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%z*zcon**2
                zcon=sum(z(1+nelenuc1st:nelenuc))*inelenuc*2.0-contrz(ctn+1)
                f(1+nelenuc1st:nelenuc)%z=f(1+nelenuc1st:nelenuc)%z-kr(i)%z*zcon*inelenuc*2.0
                econ = econ + 0.5*kr(i)%z*zcon**2
              else
                zcon=sum(z(1:nelenuc))*inelenuc-contrr(i)%z
                f(1:nelenuc)%z=f(1:nelenuc)%z-kr(i)%z*zcon*inelenuc
                econ = econ + 0.5*kr(i)%z*zcon**2
              endif
            else
              f(j)%z=f(j)%z-kr(i)%z*(r(j)%z-contrr(i)%z)
              econ = econ + 0.5*kr(i)%z*(r(j)%z-contrr(i)%z)**2
            endif
          endif
        enddo
      endif
    endif
  endif ! Qnuc
  !  nonbonded interactions between interaction sites and ions
  if (Qpar .and. Qnucl .and. Qnonbond) then
    !$omp do
    do i = nelenuc+1, nele
      itype = et(i)
      do j = 1, nelenuc
        jtype = et(j)
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
          eefpotmxloc=eefpotmxloc+eefp
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
            eqqmxloc = eqqmxloc + cofo
          endif
  ! Compute Lennard Jones contribution 
          if (Qlj(is)) then 
            dist6 = (sgp2(is)*idist2)**3
            dist12 = dist6**2
            evdwmxloc = evdwmxloc + epp4(is)*(dist12-dist6)
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
              esrpmfmxloc = esrpmfmxloc + esrpmf0
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
              fxloc(j) = fxloc(j) - de*(r(i)%x-r(j)%x)
              fyloc(j) = fyloc(j) - de*(r(i)%y-r(j)%y)
              fzloc(j) = fzloc(j) - de*(r(i)%z-r(j)%z)
            endif
          endif
        endif
      enddo   ! j=1,...,nelenuc
    enddo ! i=nelenuc+1,...,nele 
    !$omp end do
  endif
  !$omp critical
  eefpotmx=eefpotmx+eefpotmxloc
  eqqmx=eqqmx+eqqmxloc
  evdwmx=evdwmx+evdwmxloc
  esrpmfmx=esrpmfmx+esrpmfmxloc
  eefpot=eefpot+eefpotloc
  eelec=eelec+eelecloc
  evdw=evdw+evdwloc
  esrpmf=esrpmf+esrpmfloc
  esolv=esolv+esolvloc
  estack=estack+estackloc
  eqq=eqq+eqqloc
  eex=eex+eexloc
  ebp=ebp+ebploc
  edihe=edihe+ediheloc
  ebond=ebond+ebondloc
  eang=eang+eangloc
  f(1:nele)%x=f(1:nele)%x+fxloc(1:nele)
  f(1:nele)%y=f(1:nele)%y+fyloc(1:nele)
  f(1:nele)%z=f(1:nele)%z+fzloc(1:nele)
  !$omp end critical
  !$omp end parallel
  ener = ener + eelec + evdw + esrpmf + esrpmfmx + ebond + eang + edihe + estack + ebp + eex + eqq + esolv + eqqmx + evdwmx + eefpot + eefpotmx + econ
endif                     !Qenergy
return

end subroutine
