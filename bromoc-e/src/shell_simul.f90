!    BROMOC-E
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program with CHARMM Force Field Reader.
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

subroutine shell_simul
use apfmod
use efpmod
use ioxmod
use stdiomod
use grandmod
use gsbpmod
use constamod
use errormod
use nucleotmod
use charfuncmod   !command parser
use splinemod
use charmmmod
use listmod

!Command parser and file name
implicit none
! allocatables
real,allocatable       :: xy(:,:) 
real,allocatable ::    epsLJ(:),sgLJ(:)
real*4,allocatable :: efield(:,:)
integer, allocatable :: nxf(:)  
logical*1, allocatable :: Qefpotread(:)
character*1,allocatable :: secstr(:),sitenam(:)
character com*2048,word*1024,fnam*80,wrd5*5
integer unvec(maxopen), totnumb, nions, kode, maxpart, ncl3, in1, in2
character*80 title
character*4 wrd4
character*4,allocatable,dimension(:) :: nucnam
!character*6 wrd6
real battery
real r1,r2,r6,z1,v1,v2,y1,y2,x1,x2,x3,xm,ym,zm,z2
integer ix1,iy1,iz1,ix2,iy2,iz2 
real*4 idv
real resol,pkind,diffu,totdiffu,pardiffu,inflim,mass
integer ikind
real vc1(3),vc2(3),vc3(3)
logical*1 endlog, logfinal, Qlsprmf, doions, dodna, Qadj, ok 
logical*1 logmemb, logphix, logphiv, logsrpmf,logbuff,Qefpott,Qepwrt,logrfpar,Qnohead
logical*1 Qexpl2nd,Qinputpar,Qonlychden,Qpdb,Qxyz,Qpdbe,Qcrd,Qcrde,Qatexp,Qhomodiffu
real*8 zero
!for time
integer*8       start,finish,timer
character*8  date
character*10 time
character*5  zone
integer   values(8)
integer   is,cnt,mnp,nnp,nn,maxd
integer   wunit,s1,s2,s3,iseed,wallsi4
integer*1 walls,dnaparams
integer nsites,nelem
!for trajectory, fraction of denatued bases files and translocation time
integer    ntras 
real       vfbs,scald,scaldd
!for security outputfile
integer    nsec
type(car) :: rr
!for TEST order
real       maxl,minl,maxlg,minlg
real       dener, rate, totvol
!for MEMBRANE and ION-ION orders
!for effective dielectric constant for DNA
real       conc
!for solvent-induced contribution
real       anumb1, epsn
!for DNA fixed sites
integer      isite
!for DNA matrix rotation
real       sum1, sum2, sum3
real       xold, yold, zold
!forgotten to declare
integer r1i,z1i,itype,jtype,i,j,k,ib,ii,ij,ik
integer iunit,lfnam,ngcmc,nmcm,nbd,ntypold,ntype,ilast
integer*8 ncycle,nsave,nsfbs
integer,allocatable :: iunitv(:)
real cc0,cc1,cc2,cc3,cc4
real xtras, ytras, ztras, rot(3,3)
real,allocatable :: xnat(:), ynat(:), znat(:), rnat(:), phinat(:)
integer parn
! hcons vars
integer hcpar, hcen
type(car) hcr,hck

!Default parameters and options
!------------------------------
unvec        = -1
Qsystem      = .false.
Qatexp       = .false.
Qpar         = .false.
Qbuf         = .false.
Qdeby        = .false.
Qdebyhyb     = .false.
Qsolv        = .false.
Qnucl        = .false.
Qdie         = .false.
Qtraj        = .false.
Qtrajcont    = .false.
Qepwrt       = .false.
Qapfor       = .false.
Qnotrans     = .false.
Qwarn        = .false.
Qcontrans    = .false.
Qcontprint   = .false.
Qunsplit     = .false.
Qmemb        = .false.
Qljpar       = .false.
Qljsin       = .false.
Qsrpmf       = .false.
Qefpott      = .false.
Qcountion    = .false.
Qchden       = .false.
frmt         = ''
Qenergy      = .true.
Qforces      = .false.
Qforceanapot = .false.
Qbond        = .true.
Qecyl        = .false.
Qnonbond     = .true.
Qgr          = .false.
Qphix        = .false.
Qphiv        = .false.
Qlsprmf      = .false.
Qnmcden      = .false.
Qdiffuse     = .false.
Qprofile     = .false.
Qproxdiff    = .false.
Qrfpar       = .false.
Qninfo       = .false.
Qpres        = .false.
Qsvdw        = .false.
Qrfpsin      = .false.
doions       = .false.
dodna        = .false.
Qresintfor   = .false.
iseed        = 3141593
ntype        = 0
nbuffer      = 0
nsites       = 0
istrs        = 0
inuc         = 0
cgnuc        = -1.0
diffnuc      = 0.1
epsnuc       = 0.184
ionstr       = 0.0
temp         = 300.0
lx           = 0.0
ly           = 0.0
lz           = 0.0
cx           = 0.0
cy           = 0.0
cz           = 0.0
cdie         = 80.0
Rsphe        = 0.0
afact        = 0.0
kappa        = 0.0
kbtdna       = 0.0
ikbtdna      = 0.0
setframes    = 0
dnaparams    = 0
scalepairing = 1.0
xs           = 0.0
ys           = 0.0
b            = 0.0
c            = 0.0
d            = 0.0
riffac       = 0.0

! initialize lists
call inivars()
call deltypall()
call delparall()

! allocate space for lists
call resizeetypl(dtype)  ! Resize Element-Type Vectors
call resizeptypl(dtype) ! Resize Particle-Type Vectors
call resizeparl(datom) ! Resize Particle Lists
call resizecvec(datom) ! Resize Element Lists

call header(outu)
start = timer()
call date_and_time(date, time, zone, values)
write(outu,'(6x,a,7(i0,a))') 'Started at (YYYY-MM-DD HH:mm:ss.ms): ',values(1),'-',values(2),'-',values(3),' ',values(5),':',values(6),':',values(7),'.',values(8)
write(outu,*)
write(outu,*)

logfinal = .false. 
do while (.not. logfinal)
  call getlin(com,inpu,outu) ! new commands line
  write(outu,'(/a)')'**********************************************'
  write(outu,'(a)') trim(adjustl(com))
  write(outu,'(a/)')'**********************************************'
  call getfirst(com,wrd5) ! new key word
  wrd5=lcase(wrd5)
  !.....MISCELLANEOUS COMMANDS    
  if  (wrd5.eq.'title') then
  !     ---------------
     call getwrd(com,'name',title)
     write(outu,'(6x,a)') 'TITLE: '//trim(title)
  ! **********************************************************************
  elseif (wrd5.eq.'open')then
  !        ---------------
     ! unit [integer,default=1]
     call gtipar(com,'unit',iunit,1)
     if (iunit.le.0) call error ('shell_simul', 'unit is a zero or negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
     if (unvec(iunit).ne.-1) call error ('shell_simul', 'unit incorrect in OPEN order',faterr)
     call lualloc(unvec(iunit))
     iunit = unvec(iunit)
     ! namefile [character]
     call getwrd(com,'name',fnam)
     lfnam=len_trim(fnam)
     write(outu,105) fnam(1:lfnam),iunit
  105     format(6x,'open file ',a,' as unit ',i3)
     write(outu,*)
     ! write/ read options
     if (check(com,'write')) then
       if (check(com,'file ')) then ! binary file
         open(unit=iunit,file=fnam,form='unformatted')
       else
         open(unit=iunit,file=fnam,form='formatted')
       endif
     elseif (check(com,'read')) then
       if (check(com,'file ')) then
         open(unit=iunit,file=fnam,status='old',form='unformatted')
         rewind(unit=iunit) ! rewinds a file to the beginning
       else
         open(unit=iunit,file=fnam,status='old',form='formatted')
         rewind(unit=iunit) ! rewinds a file to the beginning
       endif
     endif
  ! **********************************************************************
  elseif (wrd5.eq.'close') then
  !        ---------------
     ! unit [integer,default=-1]
     call gtipar(com,'unit',iunit,1)
     if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
     if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in CLOSE order',faterr)
     ij = unvec(iunit) 
     unvec(iunit) = -1
     iunit = ij
     close(unit=iunit)
     call lunalloc (iunit)
     write(outu,'(6x,a,1x,i3)') 'close unit',iunit
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'nucle') then
  !         ---------------
     if (.not.Qsystem) call error('shell_simul','SYSTEM must be defined before NUCLEOTIDES',faterr)
     if (Qatexp) call error ('shell_simul', 'NUCLEOTIDES order is defined after PTYPE order', faterr)
     if (Qpar) call error ('shell_simul', 'NUCLEOTIDES order is defined after PARTICLE order', faterr)
     ! number of strand [default=0]
     call gtipar (com,'strands',istrs,0)
     if (istrs.gt.2 .or. istrs.lt.0) call error ('shell_simul', 'Incorrect numbers of strands in NUCLEOTIDES order', faterr)
     ! number of nucleotides of each strand [default=0]
     call gtipar (com,'nucleot',inuc,0)
     if (inuc.le.0) call error ('shell_simul', 'Incorrect numbers of nucleotides in NUCLEOTIDES order', faterr)
     if (allocated(nucnam)) deallocate(nucnam)
     allocate (nucnam(inuc*2))
     ! allocate variables
     if (Qnucl) then
       deallocate (strand,typenuc,namnucl,namsite,xnat,ynat,znat,rnat,phinat)
       deallocate (sitebond,siteangle,sitedihe,distbond,valangle,valdihe,bond,angle,typbond)
       deallocate (valdihe,bond,angle,typbond,sitestack,sitebp,siteex,siteqq,siteslv)
       deallocate (sgstack,sgbp,sgex,sitenam)
     endif
     maxsite = istrs*inuc*3
     maxbond = maxsite-1
     maxang  = 2*maxsite-2
     maxdihe = 2*maxsite-3
     maxpar  = maxsite*(maxsite-1)/2
     allocate (strand(maxsite),typenuc(maxsite),namnucl(maxsite),namsite(maxsite))
     allocate (xnat(maxsite),ynat(maxsite),znat(maxsite),rnat(maxsite),phinat(maxsite))
     allocate (sitebond(maxbond,2),siteangle(maxang,3),sitedihe(maxdihe,4),distbond(maxbond),valangle(maxang))
     allocate (valdihe(maxdihe),bond(maxpar),angle(maxpar),typbond(maxbond))
     allocate (sitestack(maxpar,2),sitebp(maxpar,2),siteex(maxpar,2),siteqq(maxpar,2),siteslv(maxpar,2))
     allocate (sgstack(maxpar),sgbp(maxpar),sgex(maxpar),sitenam(inuc))
     ! charge nucleotides [default=-1]  
     call gtdpar(com,'charge',cgnuc,-1.0)
     ! diffusion constant for nucleotides [default=0.1]
     call gtdpar(com,'diffusion',diffnuc,0.01)
     if (.not. Qnucl) then
       ! Add Mono-element Particle Types
       call addetyp('S', dif=diffnuc,eps=0.0647,sig=4.16)  ! 1
       call addetyp('P', dif=diffnuc,eps=0.1796,sig=3.14)  ! 2
       call addetyp('Ab',dif=diffnuc,eps=0.4598,sig=4.10)  ! 3
       call addetyp('Tb',dif=diffnuc,eps=0.4598,sig=4.90)  ! 4
       call addetyp('Cb',dif=diffnuc,eps=0.4598,sig=3.40)  ! 5
       call addetyp('Gb',dif=diffnuc,eps=0.4598,sig=3.80)  ! 6
     endif
     if (diffnuc.lt.0.0) call error ('shell_simul', 'Diffusion coefficient for each nucleotide is negative', faterr)
     ! Parameter to calculate bonded and non-bonded potential
     ! energy terms [default=0.26]        
     call gtdpar(com,'dnatemp',dnatemp,temp)
     if (dnatemp.ne.temp) then
       write(outu,'(6x,a)') 'Different Temperature for DNA activated' 
       write(outu,'(6x,a,f10.5)') '   DNA Temperature: ',dnatemp
       write(outu,'(6x,a,f10.5)') '   Defined DNA Diffusivity: ',diffnuc
       write(outu,'(6x,a)') '   Rescaling DNA Diffusivity proportionally to DNA Temperature' 
       diffnuc=dnatemp/temp*diffnuc
       write(outu,'(6x,a,f10.5)') '   New DNA Diffusivity: ',diffnuc
       kbtdna = kboltz*dnatemp/kcalmol
       ikbtdna = 1.0/kbtdna
     endif
     ! Parameter to calculate bonded and non-bonded potential energy terms [real,default=0.184]
     call gtdpar(com,'eps',epsnuc,epsnuc)
     if (epsnuc.lt.0.0) call error ('shell_simul', 'epsnuc for each nucleotide is negative', faterr)
     ! Scale Base pairing energy/forces [real,default=1.0]
     call gtdpar(com,'scalepairing',scalepairing,scalepairing)
     if (scalepairing.lt.0.0) call error ('shell_simul', 'scalepairing is negative', faterr)
     ! Enables solvent-induced contribution
     Qsolv = check(com,'solv') .and. istrs.eq.2
     if (Qsolv) then
       call gtdpar (com,'conc',conc,0.0)
       if (conc.lt.0.0) call error ('shell_simul', 'Molarity is negative in NUCLEOTIDES order', faterr)
       if (conc.eq.0.0) call error ('shell_simul', 'Conc is zero (needed for solv)', warning)
       epsn = 0.504982*epsnuc*(1.0-1.0/(1.40418-0.268231*inuc))
       anumb1 = 0.474876*(1.0 + 1.0/(0.148378+10.9553*conc))
       epsolv = epsn*anumb1
       write(outu,'(6x,a)') 'Solvent-induced contribution enabled (Recommended to improve double stranded DNA interactions)'
       write(outu,'(6x,a,f10.5)') 'Solvent-induced Parameter= ',epsolv
     endif
     ! Enables read and writing sequences in the opposite order
     Qinvstr=check(com,'3-5')
     ! Allows specification of complementary strand to be written explicitly
     Qexpl2nd=check(com,'explicit2nd')
     if (Qexpl2nd.and.istrs.ne.2) call error('shell_simul','explicit2nd ignored, 2nd strand was not specified',warning)
     ! Includes first 5' P in the structure
     QfirstP=check(com,'keepfirstp')
     ! Enables manual input of structural parameters for DNA
     Qinputpar=check(com,'inputparam')
     ! Generic strucutral parameters for DNA
     if (check(com,'depablo')) then
       dnaparams=0
     elseif (check(com,'charmmvac')) then 
       dnaparams=1
     elseif (check(com,'charmmwat')) then 
       dnaparams=2
     elseif (check(com,'charmmwatbc')) then
       dnaparams=3
     else ! charmmwatbc
       dnaparams=3
     endif
     if (Qexpl2nd.and.istrs.ne.2) call error('shell_simul','explicit2nd ignored, 2nd strand was not specified',warning)
     Qninfo=check(com,'printinfo')
     nsites=0
     if (istrs.gt.0) then
       j=inuc*3
       if (.not.QfirstP) j=j-1
       call addptyp(j,'DNA1')
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of date in NUCLEOTIDE order', faterr)
       do j = 1, inuc
         ntype = ntype + 1
         if (ntype.gt.dtype) call error ('shell_simul', 'ntype is greater than dtype', faterr)
         call getfirst(com,word)
         if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of date in NUCLEOTIDE order', faterr)
         ! nucleotide name [character*4]
         wrd4 = ucase(word)
         nucnam(ntype) = wrd4
         ! Phosphate
         if ((QfirstP.and.j.eq.1).or.j.gt.1) then
           nsites = nsites + 1
           if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
           strand(nsites) = 1
           typenuc(nsites) = ntype
           namnucl(nsites) = wrd4(1:1)
           namsite(nsites) = 'P'
           call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
         endif
         ! Base
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 1                 
         typenuc(nsites) = ntype
         namnucl(nsites) = wrd4(1:1)
         if (namnucl(nsites).eq.'A') then
           namsite(nsites) = 'Ab' 
           call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'T') then
           namsite(nsites) = 'Tb'
           call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'C') then
           namsite(nsites) = 'Cb'
           call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'G') then
           namsite(nsites) = 'Gb'
           call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
         else      
           call error ('shell_simul', 'nucleotide name is not correct', faterr)
         endif  
         ! Sugar
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 1
         typenuc(nsites) = ntype
         namnucl(nsites) = wrd4(1:1) 
         namsite(nsites) = 'S'
         call seteleinptyp(nptyp,nsites,getetyp(namsite(nsites)))
       enddo
       call updateptypchg(nptyp)
     endif 
     ! read explicit second strand
     if (Qexpl2nd) then
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of second strand in NUCLEOTIDE order', faterr)
       allocate (secstr(inuc))
       do j=inuc,1,-1
         call getfirst(com,word)
         if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of second strand in NUCLEOTIDE order', faterr)
         ! nucleotide name [character*4]
         secstr(j) = ucase(word)
       enddo
     endif
     if (Qinputpar) then
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       do i=1,6
         do j=1,3
           call getfirst(com,word)
           if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
           !Bases
           cylall(i,j) = chr2real(word)
         enddo
       enddo
       call getfirst(com,word)
       if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       !Internuclear Distance
       din = chr2real(word)
       call getfirst(com,word)
       if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       !Internuclear Angle
       ain = chr2real(word)
       write(outu,'(6x,a)') 'Using User-defined structural parameters for DNA'
     else
       !Cylindrical polar coordinates
       if (dnaparams.eq.0) then 
         !Bases
         cylall(1,1:3) = (/ 0.773, -41.905, -0.051 /)         !Ab
         cylall(2,1:3) = (/ 2.349, -86.119, -0.191 /)         !Tb
         cylall(3,1:3) = (/ 2.296, -85.027, -0.187 /)         !Cb
         cylall(4,1:3) = (/ 0.828, -40.691, -0.053 /)         !Gb
         cylall(5,1:3) = (/ 6.981, -70.197, -1.280 /)         !Sugar
         cylall(6,1:3) = (/ 8.918, -94.038, -2.186 /)         !Phosphate
         din = 3.38
         ain = 36.0
         write(outu,'(6x,a)') 'Using internal De Pablo et al structural parameters for DNA'
       elseif (dnaparams.eq.1) then
         cylall(1,1:3) = (/ 1.5401186532694178, -19.177377752137673, 4.81713596295395771E-002 /)
         cylall(2,1:3) = (/ 2.5874118642077302, -64.556009869361233, 0.22903343233126944      /)
         cylall(3,1:3) = (/ 2.2651275581525985, -69.277179138939957, 0.29825389290205251      /)
         cylall(4,1:3) = (/ 1.6463744611104871, -24.136045962260361, 6.06394873617855691E-003 /)
         cylall(5,1:3) = (/ 7.2988301635144150, -63.884837335190518, 0.45948475145162820      /)
         cylall(6,1:3) = (/ 9.4038525018460799, -89.846061205360229, -0.48592487490127378     /)
         din = 3.3692353762855434 
         ain = 36.790780189800941
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in vacuo'
       elseif (dnaparams.eq.2) then
         cylall(1,1:3) = (/ 0.83216456560990049,-43.802827117829452, 5.84700711558472225E-002 /)
         cylall(2,1:3) = (/ 2.25434746257693770,-87.586682816786521, 7.58267641147968852E-002 /) 
         cylall(3,1:3) = (/ 2.49227490170237240,-87.699281096559588, 5.24429792951401769E-002 /)
         cylall(4,1:3) = (/ 0.89858765931206686,-21.581894968064670, 0.11155621478344679      /)
         cylall(5,1:3) = (/ 7.01490325321993600,-68.281680932114654, 1.34340627519265107E-002 /)
         cylall(6,1:3) = (/ 9.08267487151476870,-92.072561423506826,-1.0389212165439070       /)
         din = 3.4190096123453726 
         ain = 36.291088852071780
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in implicit water'
       elseif (dnaparams.eq.3) then
         cylall(1,1:3) = (/ 2.6019796769397963,-81.191725312579493, 0.15848098213497963      /)
         cylall(2,1:3) = (/ 3.3706977830909244,-96.328909516456363, 0.24317884044268287      /)
         cylall(3,1:3) = (/ 3.4919037468548892,-90.384267464738386, 0.24690507752582472      /)
         cylall(4,1:3) = (/ 2.2051192772253185,-68.561436289845915, 0.10838747455031919      /)
         cylall(5,1:3) = (/ 7.0149106362072891,-68.281465812483361, 1.34359749915350605E-002 /)
         cylall(6,1:3) = (/ 9.0826733119583505,-92.072394528330037,-1.0389241599469166       /)
         din = 3.4189906453957799
         ain = 36.290944969867944
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in implicit water with bases centered at geometric center of bases heavy atoms'
       endif
     endif
     write(outu,'(/6x,a)')       'Particle   xy-dist    xy-ang         z' 
     write(outu,'(6x,a,3f10.4)') 'Ab      ',(cylall(1,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Tb      ',(cylall(2,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Cb      ',(cylall(3,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Gb      ',(cylall(4,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'S       ',(cylall(5,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'P       ',(cylall(6,i),i=1,3)
     write(outu,'(/6x,2(a,f10.4)/)') 'Internucleotide: Distance= ',din,' Angle= ',ain
     if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in NUCLEOTIDE order', faterr)
     if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in NUCLEOTIDE order', faterr)
     if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in NUCLEOTIDE order', faterr)
     ! complementary strand
     nelenuc1st = nsites
     if (istrs.eq.2) then
       j=inuc*3
       if (.not.QfirstP) j=j-1
       call addptyp(j,'DNA2')
       ntypold = ntype      
       do j = 0, inuc-1
         ntype = ntype + 1
         if (ntype.gt.dtype) call error ('shell_simul', 'ntype is greater than dtype', faterr)
         if (Qexpl2nd) then
           wrd4=secstr(j+1)
         else
           ! nucleotide name [character*4]
           if (nucnam(ntypold-j)(1:1).eq.'A') then
             wrd4 = 'T'
           else if (nucnam(ntypold-j)(1:1).eq.'T') then
             wrd4 = 'A'
           else if (nucnam(ntypold-j)(1:1).eq.'C') then
             wrd4 = 'G'
           else
             wrd4 = 'C'
           endif     
         endif 
         nucnam(ntype) = wrd4
         ! Phosphate
         if ((QfirstP.and.j.eq.0).or.j.gt.0) then
           nsites = nsites + 1
           if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
           strand(nsites) = 2
           typenuc(nsites) = ntype
           namnucl(nsites) = wrd4(1:1)
           namsite(nsites) = 'P'
           call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
         endif
         ! Base
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 2
         typenuc(nsites) = ntype
         namnucl(nsites) = wrd4(1:1)           
         if (namnucl(nsites).eq.'A') then
           namsite(nsites) = 'Ab'
           call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'T') then
           namsite(nsites) = 'Tb'
           call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'C') then
           namsite(nsites) = 'Cb'
           call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
         else if (namnucl(nsites).eq.'G') then
           namsite(nsites) = 'Gb'
           call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
         else
           call error ('shell_simul', 'nucleotide name is not correct', faterr)
         endif
         ! Sugar
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 2
         typenuc(nsites) = ntype
         namnucl(nsites) = wrd4(1:1)
         namsite(nsites) = 'S'            
         call seteleinptyp(nptyp,nsites-nelenuc1st,getetyp(namsite(nsites)))
       enddo
       if (Qexpl2nd) deallocate (secstr)
       call updateptypchg(nptyp)
     endif
     fctn=celec*cgnuc**2/cdie
  
     write(outu,*)
     write(outu,'(6x,a,i3,a)') 'There are ',istrs*inuc,' nucleotides'
     if (Qinvstr) then
       write(outu,'(6x,a)') 'DNA strands are written in 3-5 direction'
     else
       write(outu,'(6x,a)') 'DNA strands are written in 5-3 direction'
     endif
     if (QfirstP) then 
       write(outu,'(6x,a)') "First 5' Phosphate included"
     else
       write(outu,'(6x,a)') "First 5' Phosphate not included"
     endif
     write(outu,'(6x,a)') 'EPSILON  -> energy scale [Kcal/mol]'
     write(outu,'(6x,a)') 'CHARGE   -> phosphate sites charge [e]'
     write(outu,'(6x,a)') 'DIFFUSION-> diffusion constant [Ang.**2/ps]'
     write(outu,'(6x,a)') 'STRAND  NAME  EPSILON  CHARGE  DIFFUSION' 
     write(outu,'(6x,a)') '----------------------------------------'
     do i = 1, ntype
       if (i.le.inuc) then
         write(outu,'(6x,i3,4x,a1,2x,f8.3,2x,f8.3,2x,e9.2)') 1, nucnam(i), epsnuc, cgnuc, diffnuc
       else 
         write(outu,'(6x,i3,4x,a1,2x,f8.3,2x,f8.3,2x,e9.2)') 2, nucnam(i), epsnuc, cgnuc, diffnuc        
       endif  
     enddo
     write(outu,*)
     ! Interaction site coordinates for the native structure
     call native_structure(nsites,xnat,ynat,znat,rnat,phinat)
     ! Define number of particle types belonging to dna
     nptnuc=nptyp
     ! Set Charges for P
     if (istrs.ge.1) call setetypchginptyp(1,getetyp('P'),cgnuc)
     if (istrs.eq.2) call setetypchginptyp(2,getetyp('P'),cgnuc)
     ! Add Particle
     if (istrs.gt.0) then
       call addpar(1,1)
       call centerptyp(1)
     endif
     if (istrs.eq.2) then
       call addpar(2,1)
       call centerptyp(2)
     endif
     ! Define number of particles belonging to dna
     nparnuc=npar
     ! Define number of elements belonging to dna
     nelenuc=nele
     ! Define number of elements types belonging to dna
     netnuc=netyp
     ! bond streching terms
     call bbonds
     ! bond angles terms
     call angles
     ! torsinal angles terms
     call dihedral
     ! nonbonded terms
     call go_qq
     Qnucl = .true.
     call allocateqsome()
     call setuplj()
     write(outu,*) 
     ! assert 
     if (nsites.ne.nele) stop 'nsites is not equal to nele'
     if (allocated(nucnam)) deallocate(nucnam)
     deallocate(xnat,ynat,znat,rnat,phinat)
  ! **********************************************************************
  elseif (wrd5.eq.'ptype') then
  !       ---------------
    if (.not.Qsystem) call error ('shell_simul', 'SYSTEM must be defined before PTYPE order', faterr)
    ! unit number for force field [integer,default=0]
    call gtipar (com,'iunprm',iunprm,0)
    if (iunprm.le.0 .or. iunprm.gt.maxopen) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    if (unvec(iunprm).eq.-1) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    iunprm = unvec(iunprm)
    ! unit number for connectivity [integer,default=0]
    call gtipar (com,'iunpsf',iunpsf,0)
    if (iunpsf.le.0 .or. iunpsf.gt.maxopen) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    if (unvec(iunpsf).eq.-1) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    iunpsf = unvec(iunpsf)
    ! Print additional information in outputfile
    ok = check(com,'print')
    ! water viscosity (Kcal ps mole^(-1) Angs.^(-1)) [real,default=0.123]
    call gtdpar(com,'viscwat',viscwat,0.1225)
    if (viscwat.lt.0.0) call error ('shell_simul', 'Wrong value for water viscosity in PTYPE order', faterr)
    ! scale factor for self-diffusion calculation [real,deafult=1.0]
    call gtdpar(com,'scald',scldiff,1.0)
    if (scldiff.lt.0.0) call error ('shell_simul', 'Wrong value for scale factor for self-diffusion calculation in PTYPE order', faterr)
    ! read PRM file
    write(outu,'(6x,a)') 'PARTICLE TYPE: READ PRM FILE'
    call readcharmm(ok)
    ! read PSF file
    write(outu,'(6x,a)') 'PARTICLE TYPE: READ PSF FILE'
    call readpsf(ok)
    ! deallocate PRM routine variables
    deallocate (charmm_label,charmm_mass)
    if (Qchmmbond) deallocate (charmm_btype,charmm_bond)
    if (Qchmmang) deallocate (charmm_atype,charmm_ang,charmm_lub)
    if (Qchmmub) deallocate (charmm_ubtype,charmm_ub)
    if (Qchmmdih) deallocate (charmm_dtype,charmm_dih,charmm_ndih,charmm_nprms)
    if (Qchmmimp) deallocate (charmm_itype,charmm_imp)
    if (Qchmmcmap) deallocate (charmm_icmap,charmm_icmap2,charmm_ncmap,charmm_cmap,charmm_fcmap)
    deallocate (charmm_typen,charmm_nonbonded)
    ! unit number for coordinate [integer,default=0]
    call gtipar (com,'iunpdb',iunpdb,0)
    call gtipar (com,'iuncrd',iuncrd,0)
    if (iunpdb.le.0.and.iuncrd.le.0) call error ('shell_simul', 'Enter a pdb or crd for coordinates in PTYPE order', faterr)
    if (iunpdb.gt.maxopen) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    if (iuncrd.gt.maxopen) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
    if (iunpdb.gt.0) then 
      if (unvec(iunpdb).eq.-1) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
      iunpdb = unvec(iunpdb)
      call loadcoorfrompdbtoptyp(nptyp,iunpdb)
    endif
    if (iuncrd.gt.0) then 
      if (unvec(iuncrd).eq.-1) call error ('shell_simul', 'Incorrect unit in PTYPE order', faterr)
      iuncrd = unvec(iuncrd)
      call loadcoorfromcrdtoptyp(nptyp,iuncrd)
    endif
    call centerptyp(nptyp)
    Qpar=.true.
    Qatexp = .true.
    call allocateqsome()
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'addpa') then ! Add Particle
  if (.not.Qsystem) call error('shell_simul','SYSTEM must be defined previously',faterr)
  ! Print pdb
  Qpdb=check(com,'printpdb')
  ! Print xyz
  Qxyz=check(com,'printxyz')
  ! Print pdbe
  Qpdbe=check(com,'printpdbe')
  ! Print crd
  Qcrd=check(com,'printcrd')
  ! Print crdext
  Qcrde=check(com,'printcrde')
  endlog = .false.
  ok=.true.
  do while (.not.endlog)
    call getlin(com,inpu,outu) ! new commands line 
    endlog = check(com,'end')
    if (.not.endlog) then
      ! name ion type [character*4]
      call getfirst(com,wrd4)
      ! Read Type
      itype=getptyp(wrd4)
      if (itype.eq.0) then
        write(outu,'(6x,a)') 'Particle Type ',wrd4,' not defined yet'
      else
        if (check(com,'solute')) then
          if (.not.ok) call error('shell_simul','ADDPAR: SOLUTE PTYPE cannot be added after SOLVENT PTYPE',faterr)
          call addpar(itype,2) ! As solute
        else
          ok = .false.
          call addpar(itype,3) ! As solvent
        endif
        call gtdpar(com,'x',rr%x,0.0)
        call gtdpar(com,'y',rr%y,0.0)
        call gtdpar(com,'z',rr%z,0.0)
        if (check(com,'randrot')) then ! Random Rotate
          call movepar(npar,rr,norot=.false.)
        else
          call movepar(npar,rr,norot=.true.)
        endif
      endif
    endif
  enddo

  ! print formats 
  if (Qxyz) call printxyz(outu)
  if (Qpdb) call printpdb(outu)
  if (Qcrd) call printcrd(outu)
  if (Qpdbe) call printpdbe(outu)
  if (Qcrde) call printcrde(outu)
  ! **********************************************************************
  elseif (wrd5.eq.'parti') then
  !       ---------------
    if (.not.Qsystem) call error('shell_simul','SYSTEM must be defined previously',faterr)
    if (Qdeby) call error('shell_simul','IMPLICIT IONS defined in system, PARTICLE section is senseless',faterr)
    endlog = .false.
    do while (.not.endlog)
      call getlin(com,inpu,outu) ! new commands line 
      endlog = check(com,'end')
      if (.not.endlog) then
        ! name ion type [character*4]
        call getfirst(com,wrd4)
        ! Add mono particle with element
        itype=getptyp(wrd4)
        if (itype.eq.0) then
          call addmonoptyp(wrd4)
          itype=nptyp
          ! particle charge [real*8,default=0]
          call gtdpar(com,'charge',ptypl(itype)%chg(1),0.0)
          ! Update Particle Charge
          call updateptypchg(itype)
          ! mass
          call gtdpar(com,'mass',mass,0.0)
          if (mass.gt.0.0) then
            ! element mass
            etypl(netyp)%mas=mass
            ! particle mass
            ptypl(nptyp)%mass=mass
          endif
        endif
        ! diffusion constant [real*8,default=0.1]
        call gtdpar(com,'diffusion',diffu,0.0)
        Qhomodiffu=check(com,'homodiffu') 
        if (diffu.eq.0.0) call error ('shell_simul', 'diffusion coefficient omitted or 0', faterr)
        if (diffu.lt.0.0) call error ('shell_simul', 'diffusion coefficient is negative', faterr)
        if (ptypl(itype)%ne.eq.1) then
          etypl(ptypl(itype)%etyp(1))%dif=diffu
        else
          if (Qhomodiffu) then
            do i=1,ptypl(itype)%ne
              etypl(ptypl(itype)%etyp(i))%dif=diffu*(1.0/ptypl(itype)%ne)
            enddo
          else
            totdiffu=0.0
            do i=1,ptypl(itype)%ne
              pardiffu=etypl(ptypl(itype)%etyp(i))%dif
              totdiffu=totdiffu+pardiffu
            enddo
            if (totdiffu.eq.0.0) call error('shell_simul', 'sum of total diffusion coefficients equals zero',faterr)
            do i=1,ptypl(itype)%ne
              ! Do not repeat the type
              ok=.true.
              do j=1,i-1
                 if (ptypl(itype)%etyp(i).eq.ptypl(itype)%etyp(j)) then
                   ok=.false.
                   exit
                 endif
              enddo
              if (ok) then
                pardiffu=etypl(ptypl(itype)%etyp(i))%dif
                etypl(ptypl(itype)%etyp(i))%dif=diffu*pardiffu/totdiffu
              endif
            enddo
          endif
        endif
        ! Add the particle temporarily to the particle list
        ! call addpar(itype,3)
      endif
    enddo
    Qpar = .true.
    call allocateqsome()
    write(outu,*)
    write(outu,'(6x,a,i3,a)') 'There are ',netyp-netnuc,' atom types'
    write(outu,'(6x,a)') 'CHARGE -> ion charge [e]'
    write(outu,'(6x,a)') 'DIFFUSION -> diffusion constant [Ang.**2/ps]'
    write(outu,'(6x,a)') 'NAME---CHARGE---DIFFUSION'
    do i=1,nptyp
      diffu=0.0
      do j=1,ptypl(i)%ne
        diffu=diffu+etypl(ptypl(i)%etyp(j))%dif
      enddo
      write(outu,'(6x,a,2f8.3)') ptypl(i)%nam,ptypl(i)%tchg,diffu
    enddo
    write(outu,*)
  ! **************************************************************************
  elseif (wrd5.eq.'sdiff') then
  !       ---------------
    endlog = .false.
    do while (.not.endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type atnam(itype)
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        ! diffusion constant [real,default=previous value]
        call gtdpar(com,'diffusion',etypl(itype)%dif,0.0)
        if (etypl(itype)%dif.lt.0.0) call error ('shell_simul', 'Diffusion coefficient is negative in SDIFF order', faterr)
      endif
    enddo
    write(outu,*)
    write(word,*) 'DIFFUSION -> diffusion constant [Angs.^2/ps]'
    write(outu,'(6x,a)') trim(adjustl(word))
    write(word,*) 'NAME---DIFFUSION'
    write(outu,'(6x,a)') trim(adjustl(word))
    do i = 1, netyp
      write(word,*) etypl(i)%nam,etypl(i)%dif
      write(outu,'(6x,a)') trim(adjustl(word))
    enddo
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'apfor') then
  !       ----------------
    if (.not. Qnucl) call error ('shell_simul', 'APFOR order is defined before NUCLEOTIDE order', faterr)
    ! number of DNA fixed sites [integer, default=0]
    call gtipar(com,'afn',afn,0)  
    if (afn.gt.nelenuc .or. afn.lt.0) call error ('shell_simul', 'Incorrect number of DNA sites',faterr)
    if (allocated(af)) deallocate (af)
    if (allocated(sn)) deallocate (sn)
    allocate (af(3,afn),sn(afn))
    if (afn.gt.0) then
      do j = 1, afn
        if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in APFOR. New line expected', faterr)
        if (.not.setint(isite,com)) call error ('shell_simul', 'premature end of data in APFOR. Integer Expected.', faterr)
        if (isite.gt.nelenuc .or. isite.le.0) call error ('shell_simul', 'Incorrect DNA site', faterr)
        call gtdpar(com,'fx',af(1,j),0.0)
        call gtdpar(com,'fy',af(2,j),0.0)
        call gtdpar(com,'fz',af(3,j),0.0)
        sn(j)=isite
      enddo
      Qapfor=.true.
    endif  
    if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in APFOR order', faterr)
    if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in APFOR order', faterr)
    if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in APFOR order', faterr)
    if (Qapfor) then
      write(outu,*)
      write(outu,'(6x,a,i3,a)') 'There are ',afn,' DNA sites to be forced:'
      do i=1,afn
        write(outu,'(6x,2i4,x,a,x,3(3x,f10.5))') i,sn(i),namsite(sn(i)),af(1:3,i)
      enddo
    endif           
  ! **********************************************************************
  elseif (wrd5.eq.'notra') then
  !        ---------------
  ! Do not translate, DNA keep the geometric center in a fixed position if off
    if (.not.Qnucl) call error ('shell_simul', 'NOTRANSL must be defined after NUCLEOTIDE', faterr)
    Qnotrans=.not.check(com,'off')
    Qnotrx=check(com,'x')
    Qnotry=check(com,'y')
    Qnotrz=check(com,'z')
    if (.not.(Qnotrx.or.Qnotry.or.Qnotrz)) Qnotrans=.false.
    inelenuc=1.0/nelenuc
    notrx=sum(r(1:nelenuc)%x)*inelenuc
    notry=sum(r(1:nelenuc)%y)*inelenuc
    notrz=sum(r(1:nelenuc)%z)*inelenuc
    if (Qnotrans) then 
      write(outu,'(6x,a)') 'NOTRANSLATION is on'
      if (Qnotrx) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at x =',notrx
      if (Qnotry) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at y =',notry
      if (Qnotrz) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at z =',notrz
    else
      write(outu,'(6x,a)') 'NOTRANSLATION is off'
    endif
  ! **********************************************************************
  elseif (wrd5.eq.'contr') then
  !        ---------------
  ! Constrain translation with an harmonic potential of DNA by the geometric center
    if (.not.Qnucl) call error ('shell_simul', 'CONTRANSL must be defined after NUCLEOTIDE', faterr)
    Qcontprint=check(com,'print')
    call gtipar(com,'ctn',ctn,0)
    if (allocated(kx)) deallocate (kx)
    if (allocated(ky)) deallocate (ky)
    if (allocated(kz)) deallocate (kz)
    if (allocated(csn)) deallocate (csn)
    if (allocated(contrx)) deallocate (contrx)
    if (allocated(contry)) deallocate (contry)
    if (allocated(contrz)) deallocate (contrz)
    allocate (kx(ctn),ky(ctn),kz(ctn),csn(ctn),contrx(ctn+1),contry(ctn+1),contrz(ctn+1))
    Qcontrans=.not.check(com,'off')
    Qunsplit=check(com,'unsplit').and.istrs.eq.2
    if (.not.Qcontrans) ctn=0
    do i = 1, ctn
      if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in CONTRA. New line expected', faterr)
      if (.not.setint(isite,com)) call error ('shell_simul', 'premature end of data in CONTRA. Integer Expected.', faterr)
      if (isite.gt.nelenuc .or. isite.lt.0) call error ('shell_simul', 'Incorrect DNA site', faterr)
      call gtdpar(com,'kx',kx(i),0.0)
      call gtdpar(com,'ky',ky(i),0.0)
      call gtdpar(com,'kz',kz(i),0.0)
      csn(i)=isite
      if (csn(i).eq.0) then
        inelenuc=1.0/nelenuc
        if (Qunsplit) then
          contrx(i)=2.0*sum(r(1:nelenuc1st)%x)*inelenuc
          contry(i)=2.0*sum(r(1:nelenuc1st)%y)*inelenuc
          contrz(i)=2.0*sum(r(1:nelenuc1st)%z)*inelenuc
          contrx(ctn+1)=2.0*sum(r(nelenuc1st+1:nelenuc)%x)*inelenuc
          contry(ctn+1)=2.0*sum(r(nelenuc1st+1:nelenuc)%y)*inelenuc
          contrz(ctn+1)=2.0*sum(r(nelenuc1st+1:nelenuc)%z)*inelenuc
        else
          contrx(i)=sum(r(1:nelenuc)%x)*inelenuc
          contry(i)=sum(r(1:nelenuc)%y)*inelenuc
          contrz(i)=sum(r(1:nelenuc)%z)*inelenuc
        endif
      else
        contrx(i)=r(csn(i))%x
        contry(i)=r(csn(i))%y
        contrz(i)=r(csn(i))%z
      endif
      call gtdpar(com,'x',contrx(i),contrx(i))
      call gtdpar(com,'y',contry(i),contry(i))
      call gtdpar(com,'z',contrz(i),contrz(i))
      Qcontrans=Qcontrans.or.kx(i).ne.0.0.or.ky(i).ne.0.0.or.kz(i).ne.0.0
    enddo
    if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in CONTRA order', faterr)
    if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in CONTRA order', faterr)
    if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in CONTRA order', faterr)
    if (Qcontrans) then
      write(outu,'(/6x,a)') 'CONTRANSLATION is on'
      write(outu,'(6x,a,i3,a)') 'There are ',ctn,' DNA sites/centroid to be constrained:'
      if (Qunsplit) write(outu,'(6x,a)') 'Double stranded DNA will not be splitted'
      do i=1,ctn
        if (csn(i).eq.0) then
          write(outu,'(6x,2i4,2(x,a,x,3(3x,f10.5)))') i,csn(i),'centroid ',kx(i),ky(i),kz(i),'at',contrx(i),contry(i),contrz(i)
        else
          write(outu,'(6x,2i4,2(x,a,x,3(3x,f10.5)))') i,csn(i),namsite(csn(i))//'       ',kx(i),ky(i),kz(i),'at',contrx(i),contry(i),contrz(i)
        endif
      enddo
    else
      write(outu,'(6x,a)') 'CONTRANSLATION is off'
    endif
  ! **********************************************************************
  elseif (wrd5.eq.'hcons') then ! HARMONIC CONSTRAIN
  !        ---------------
  ! Constrain translation with an harmonic potential of DNA by the geometric center
    if (npar.eq.0) call error ('shell_simul', 'No Particles added to List. Add this section after solute particles are added.', faterr)
    if (allocated(efix)) deallocate (efix)
    nefix = 0
    endlog = .false.
    do while (.not.endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = lcase(com(1:3)).eq.'end'
      if (endlog) exit
      call gtipar(com,'par', hcpar, 0) ! Particle Number in List
      if (hcpar.eq.0) call error ('shell_simul', 'particle number not present', faterr)
      if (hcpar.gt.npar) call error ('shell_simul', 'particle selected not yet inserted', faterr)
      if (parl(hcpar)%kind.gt.2) call error ('shell_simul', 'particle selected is not solute', faterr)
      call gtipar(com,'pen', hcen, 0) ! Particle Element Number
      if (hcen.eq.0) call error ('shell_simul', 'particle element number not present', faterr)
      if (hcen.gt.parl(hcpar)%ne) call error ('shell_simul', 'particle element number is greater than the number of elements in particle', faterr)
      call addefix()
      efix(nefix)%fen = parl(hcpar)%sr+hcen
      call gtdpar(com,'k',hck%x,0.0)
      call gtdpar(com,'kx',hck%x,hck%x)
      call gtdpar(com,'ky',hck%y,hck%x)
      call gtdpar(com,'kz',hck%z,hck%x)
      if (hck%x.eq.0.0.and.hck%y.eq.0.0.and.hck%z.eq.0) call error ('shell_simul', 'force constant defined is null; senseless', faterr)
      efix(nefix)%fc=hck
      call gtdpar(com,'rx',hcr%x,r%x)
      call gtdpar(com,'ry',hcr%y,r%x)
      call gtdpar(com,'rz',hcr%z,r%x)
    enddo
    write(outu,'(6x,a)') 'HARMONIC CONSTRAIN for Solutes defined'
    write(outu,'(6x,a)') 'Constrain Number-----Element Number----------kx----ky----kz----------rx----ry----rz'
    do i=1,nefix
      write(outu,'(6x,i5,i5,2(3x,3f10.5))') i, efix(i)%fen, efix(i)%fc%x,efix(i)%fc%y,efix(i)%fc%z,efix(i)%rfx%x,efix(i)%rfx%y,efix(i)%rfx%z
    enddo
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'syste') then
  !        ---------------
    Qwarn = check(com,'showwarn')
    if (Qwarn) write(outu,'(6x,a)') 'Detailed warnings enabled'
    Qdnafree = check(com,'dnafree')
    if (Qdnafree) write(outu,'(6x,a)') 'DNA not forced to be within system boundaries'
    ! random number seed [integer,default=3141593]
    call gtipar(com,'iseed',iseed,iseed)
    call feedseed(iseed) ! if zero or lower will use cpu clock to generate seed
    write(outu,'(6x,a,i0)') 'Using seed: ',iseed
    ! Temperature
    call gtdpar(com,'temp',temp,temp)
    if (temp.lt.0.0) call error ('shell_simul', 'Temperature must be positive',faterr)
    write(outu,'(6x,a,f8.3)') 'Temperature: ',temp
    ! Effective dielectric constant estimated from temperatura and salt concentration
    Qdie = check(com,'calcdie')
    if (Qdie) then
      call gtdpar(com,'conc',conc,0.0)
      if (conc.lt.0.0) then
        call error ('shell_simul', 'Molarity is negative in SIMUL order', faterr)
      elseif (conc.ge.0.0) then
        cdie = (1.0-0.2551*conc+5.151e-2*conc**2-6.889e-3*conc**3)*(249.4-0.788*temp+7.20e-4*temp**2)
        write(outu,'(6x,a,f8.2,a,f8.5,a,f8.4)')'Dielectric Constant estimated from   Temp=',temp,'   Salt Conc=',conc,'  is ',cdie
      endif
    else 
      ! Dielectric solvent [real*8,default=80 (water dielectric
      ! constant)]
      call gtdpar(com,'cdie',cdie,cdie)
      if (cdie.lt.0.0) then
        call error ('shell_simul', 'Dielectric constant is negative in SRPMF order', faterr)
      endif
      write(outu,'(6x,a,f8.3)') 'Dielectric constant: ',cdie
    endif
    ! Logical variable which indicates if Coulomb interactions
    ! are taken into account using the Debye-Hückel approximation
    Qdeby = check(com,'debye')
    Qdebyhyb = check(com,'debyhyb')
    ! Ionic strength [default=0]
    if (Qdeby.and.Qdebyhyb) call error ('shell_simul','Hybrid and Normal Debye-Huckel cannot be combined.',faterr)
    call gtdpar(com,'ionic',ionstr,ionstr)
    if (ionstr.lt.0.0) then
      call error ('shell_simul', 'Ionic strength is negative', faterr)
    elseif (ionstr.gt.0.0) then
       write(outu,'(6x,a,1x,f10.5)') 'Ionic Strength: ',ionstr
    endif
    if ((Qdeby.or.Qdebyhyb).and.ionstr.eq.0.0) then
      if (Qdeby) Qdeby=.false.
      if (Qdebyhyb) Qdebyhyb=.false.
      call error ('shell_simul','Debye-Huckel approximation disabled. Ionic Strength cannot be zero.',warning)
    endif
    if (Qdeby.and.ionstr.gt.0.0) then
       write(outu,'(6x,a)') 'Debye-Huckel approximation enabled'
    endif
    if (Qdebyhyb.and..not.Qdnafree) then
      Qdebyhyb=.false.
      call error ('shell_simul','Hybrid Debye-Huckel approximation disabled. Use it in combination with DNAFREE.',warning)
    endif
    if (.not.Qdebyhyb.and..not.Qdeby.and.Qdnafree) then
      call error ('shell_simul','DNAFREE is recommended to be used with Normal or Hybrid Debye-Huckel (DEBYHYB) approximation.',warning)
    endif
    if (Qdebyhyb.and.Qdnafree.and.ionstr.gt.0.0) then
       write(outu,'(6x,a)') 'Hybrid Debye-Huckel approximation enabled'
    endif

    ! Transmembrane potential [real*8,default=0]
    call gtdpar(com,'volt',voltage,0.0)
    write(outu,'(6x,a,1x,f10.5)') 'Transmembrane Potential (Volts): ',voltage
     
    voltage = voltage*Coulomb/kcalmol
  
    ! Activates calculation of Virial Pressure
    Qpres = check(com,'pres')
    if (Qpres) write(outu,'(6x,a)') 'Virial Pressure enabled'
  
    ! define constant kappa
    if (ionstr.gt.0.0) then 
      kappa = sqrt(const*cdie*temp/ionstr)
      ikappa = 1.0/kappa ! screening factor
    endif
  
    ! define more constants
    kBT      = kboltz*temp/kcalmol
    ikbt = 1.0/kBT
    kbtdna = kBT
    ikbtdna = ikbt
    Qecyl=check(com,'ecyl') 
    Qsphere = check(com,'sphere') ! logical*1 variable
    if (Qecyl.and.Qsphere) call error ('shell_simul','Elliptical Cylinder and Spherical system cannot be used together',faterr)
    if (Qsphere) then
       ! radius of spherical system [real*8,default=0]     
       call gtdpar(com,'radi',rsphe,rsphe)
       if (Rsphe.lt.0.0)  call error ('shell_simul', 'radi is lower than zero in SYSTEM order', faterr)
       tvol=4.0/3.0*pi*rsphe**3
       rsphe2=rsphe**2
       lx = 2.0*rsphe
       ly = 2.0*rsphe
       lz = 2.0*rsphe
       maxl=2.0*rsphe+1.0
    elseif (Qecyl) then
      ! Elliptical Cylinder
       ! minor or mayor axis of ellipse
       call gtdpar(com,'lx',lx,lx)
       if (lx.lt.0.0) call error ('shell_simul', 'LX is lower than zero in SYSTEM order', faterr)
       ! minor or mayor axis of ellipse
       call gtdpar(com,'ly',ly,ly)
       if (ly.lt.0.0) call error ('shell_simul', 'LY is lower than zero in SYSTEM order', faterr)
       ! length of cylinder
       call gtdpar(com,'lz',lz,lz)
       if (lz.lt.0.0) call error ('shell_simul', 'LZ is lower than zero in SYSTEM order', faterr)
       tvol=0.25*pi*lx*ly*lz
       iecx=2.0/lx
       iecy=2.0/ly
       maxl=sqrt(max(lx,ly)**2+lz**2)+1.0
    else 
       ! Orthorhombic box size along the X-axis [real*8,default=0]
       call gtdpar(com,'lx',lx,lx)
       if (lx.lt.0.0) call error ('shell_simul', 'LX is lower than zero in SYSTEM order', faterr)
       ! Orthorhombic box size along the Y-axis [real*8,default=0]
       call gtdpar(com,'ly',ly,ly)
       if (ly.lt.0.0) call error ('shell_simul', 'LY is lower than zero in SYSTEM order', faterr)
       ! Orthorhombic box size along the Z-axis [real*8,default=0]
       call gtdpar(com,'lz',lz,lz)
       if (lz.lt.0.0) call error ('shell_simul', 'LZ is lower than zero in SYSTEM order', faterr)
       tvol=lx*ly*lz
       maxl=sqrt(lx**2+ly**2+lz**2)+1.0
    endif
    call gtdpar(com,'cx',cx,cx)        ! Center of System along X-axis
    call gtdpar(com,'cy',cy,cy)        ! Center of System along Y-axis
    call gtdpar(com,'cz',cz,cz)        ! Center of System along Z-axis
    lx2p = cx+0.5*lx
    ly2p = cy+0.5*ly
    lz2p = cz+0.5*lz
    lx2m = cx-0.5*lx
    ly2m = cy-0.5*ly
    lz2m = cz-0.5*lz
    maxpart=int(avogadro*1e-27*3.0*tvol)+1  ! no more than 3 Molar of particles in the system volume
    if (maxpart.gt.datom) maxpart=datom 
    Qsystem = .true. 
    if (Qsphere) then
      write(outu,'(6x,a,f12.3)') 'Spherical system, Radius ',Rsphe
    elseif (Qecyl) then
      write(outu,'(6x,a,2f12.3,a,f12.3)') 'Elliptical Cylinder system. Ellipse diameters x and y:',lx,ly,' Cylinder Length: ',lz
    else
      write(outu,'(6x,3(a,f12.3))') 'Box:  LX ',lx,'  LY ',ly,'  LZ ',lz
    endif
    write(outu,'(6x,3(a,f12.3))') 'Center:  X ',cx,'  Y ',cy,'  Z ',cz
    write(outu,'(6x,a,f20.3)') 'System Total Volume (Ang**3): ',tvol
  ! **********************************************************************
  elseif (wrd5.eq.'buffe') then
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'BUFFER order is defined before PARTICLE order', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'BUFFER order is defined before SYSTEM order', faterr)
     Qforceanapot=check(com, 'forceanapot')
     endlog = .false.
     nbuffer=0
     do while (.not.endlog)
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         nbuffer = nbuffer + 1
         if (nbuffer.gt.dbuff) call error ('shell_simul', 'nbuffer greater than dbuff', faterr)
         ! Obtention of ion type 
         call getfirst(com,wrd4)
         itype=getptyp(wrd4)
         ibfftyp(nbuffer)=itype
         ! Intrinsic chemical potential [real*8,default=0] 
         call gtdpar(com,'mu',mu(nbuffer),0.0)
         ! To avoid unexpected large deviation from average number 
         ! of ions in the buffer regions, Kb [default=0] 
         Qbufferbias(nbuffer)= check(com,'bufferbias')
         call gtdpar(com,'kb',kb(nbuffer),0.0)
         call gtdpar(com,'lzmin',LZmin(nbuffer),lz2m) 
         call gtdpar(com,'lzmax',LZmax(nbuffer),lz2p) 
         if (LZmin(nbuffer).gt.0.0) then 
           if (LZmin(nbuffer).ge.LZmax(nbuffer)) call error ('shell_simul', 'LZmin => Lzmax in BUFFER order', faterr) 
         else 
           if (LZmax(nbuffer).lt.0.0.and.abs(LZmin(nbuffer)).le.abs(LZmax(nbuffer))) call error('shell_simul', 'LZmin => Lzmax in BUFFER order',faterr) 
         endif 
         if (Qsphere) then
           ! Minimum position of buffer sphere [real*8,default=0]
           call gtdpar(com,'rmin',Rmin(nbuffer),0.0)
           if (Rmin(nbuffer).lt.0.0) call error ('shell_simul', 'Rmin is lower than zero in BUFFER order', faterr)
           if (Rmin(nbuffer).gt.Rsphe) call error ('shell_simul', 'spherical region is too large in BUFFER order', warning)
           ! Maximum position of buffer sphere [real*8,default=r1]
           call gtdpar(com,'rmax',Rmax(nbuffer),Rsphe)
           if (Rmax(nbuffer).lt.0.0) call error ('shell_simul', 'Rmax is lower than zero in BUFFER order', faterr)
           if (Rmax(nbuffer).gt.Rsphe) call error ('shell_simul', 'spherical region is too large in BUFFER order', faterr)
           if (Rmin(nbuffer).ge.Rmax(nbuffer)) call error ('shell_simul', 'Rmin is greater or equal to Rmax in in BUFFER order', faterr)
           if (LZmin(nbuffer).gt.0.0) then 
             if (LZmin(nbuffer).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmin has an incorrect value', faterr) 
             if (LZmax(nbuffer).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmax has an incorrect value', faterr) 
           else 
             if (abs(LZmin(nbuffer)).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmin has an incorrect value', faterr) 
             if (abs(LZmax(nbuffer)).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmax has an incorrect value', faterr) 
           endif 
           ! Obtention of buffer volume (spherical system)
           r1 = Rmin(nbuffer)
           r2 = Rmax(nbuffer)
           logbuff = LZmin(nbuffer).eq.lz2m.and.LZmax(nbuffer).eq.lz2p
           if (logbuff) then 
             v1 = 2.0*twopi*r1**3/3.0 
             v2 = 2.0*twopi*r2**3/3.0 
           else 
             if (LZmin(nbuffer).gt.0.0) z1 = LZmin(nbuffer) 
             if (LZmax(nbuffer).lt.0.0) z1 = abs(LZmax(nbuffer))
             if (r1.gt.z1) then 
               v1 = twopi*(r1**3/3.0-(z1*r1**2*0.5-z1**3/6.0)) 
               v2 = twopi*(r2**3/3.0-(z1*r2**2*0.5-z1**3/6.0)) 
             else 
               v1 = 0.0
               v2 = twopi*(r2**3/3.0-(z1*r2**2*0.5-z1**3/6.0)) 
             endif 
           endif 
           volume(nbuffer) = v2-v1
         elseif (Qecyl) then
           ! Obtention of buffer volume (elliptical cylinder system)      
           volume(nbuffer) = 0.25*pi*lx*ly*(LZmax(nbuffer)-LZmin(nbuffer))
         else
           ! Obtention of buffer volume (orthorombic system)      
           volume(nbuffer) = (LX*LY*(LZmax(nbuffer)-LZmin(nbuffer)))
         endif
         ! Concentration for a specific ion in a buffer
         ! [real*8,default=0]
         call gtdpar(com,'conc',density(nbuffer),0.0)
         if (density(nbuffer).lt.0.0) call error ('shell_simul', 'conc is lower than zero in BUFFER order', faterr)
         ! Obtention of average number of ions in the buffer regions
         density(nbuffer) = density(nbuffer)*(avogadro/liter)*(angstrom**3)
         ! avogadro=6.022045D23; liter=0.001; angstrom=1.0E-10 
         avnum(nbuffer) = density(nbuffer)*volume(nbuffer)
         if (int(avnum(nbuffer)).le.0) call error ('shell_simul','Buffer volume or ion density is too low',warning)
         ! Average number of ions in the buffer regions [real*8,default=0]
         call gtdpar(com,'aver',avnum(nbuffer),avnum(nbuffer))
         if (avnum(nbuffer).lt.0.0) call error ('shell_simul', 'aver is lower than zero in BUFFER order', faterr)
         ! Obtention of concentration for a specific ion in a buffer
         density(nbuffer) = avnum(nbuffer)/volume(nbuffer)
         ! Transmembrane potential for a buffer
         call gtdpar(com,'volt',battery,0.0)
         ! OBTENTION OF THE ELECTROCHEMICAL POTENTIAL  
         mu(nbuffer) = mu(nbuffer)+ptypl(itype)%tchg*battery*Coulomb/kcalmol
       endif
     enddo
     Qbuf = .true.
  
     write(outu,*)
     write(outu,'(6x,a,i3,a)') 'There are ',nbuffer,' buffers'
     write(outu,'(6x,a)') 'MU    -> chemical potential [Kcal/mol]'
     if (Qsphere) then
       if (.not.logbuff) then
         write(outu,'(6x,a)') 'LZmin -> minimum position along Z-axis [Ang.]'
         write(outu,'(6x,a)') 'LZmax -> maximum position along Z-axis [Ang.]'
       endif
       write(outu,'(6x,a)') 'Rmin -> minimum radius of buffer sphere [Ang.]'
       write(outu,'(6x,a)') 'Rmax -> maximum radius of buffer sphere [Ang.]'
     else
       write(outu,'(6x,a)') 'LZmin -> minimum position along Z-axis [Ang.]'
       write(outu,'(6x,a)') 'LZmax -> maximum position along Z-axis [Ang.]'
     endif  
     write(outu,'(6x,a)') 'AVER  -> average number for an ion'
     write(outu,'(6x,a)') 'DENS  -> density for an ion [Ang.**(-3)]'
     write(outu,'(6x,a)') 'VOL   -> volume for a buffer [Ang.**3]'
     write(outu,'(6x,a)') 'KB    -> to avoid unexpected large deviation from an average number of ions'
     if (.not.Qsphere) then
       write(outu,'(6x,a)')   'NAME------MU----------LZmin-------LZmax-------AVER------DENS--------VOL----------KB----'
     else
       if (logbuff) then
         write(outu,'(6x,a)') 'NAME------MU----------Rmin--------Rmax----AVER----DENS----VOL----KB----'
       else
         write(outu,'(6x,a)') 'NAME------MU----------LZmin-------LZmax---Rmin----Rmax----AVER---DENS----VOL----KB----'
       endif
     endif
     do ib = 1, nbuffer
        ! Initializations
        nremove(ib)= 0
        ninsert(ib)= 0
        if (Qsphere) then
          if (logbuff) then
            write(outu,'(6x,a,4f12.5,1x,2e12.4,f10.5)') ptypl(ibfftyp(ib))%nam,mu(ib),Rmin(ib),Rmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
          else
            write(outu,'(6x,a,6f12.5,1x,2e12.4,f10.5)') ptypl(ibfftyp(ib))%nam,mu(ib),LZmin(ib),LZmax(ib),Rmin(ib),Rmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
          endif
        else
          write(outu,'(6x,a,4f12.5,1x,2e12.4,f10.5)') ptypl(ibfftyp(ib))%nam,mu(ib),LZmin(ib),LZmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
        endif
     enddo
     if (Qforceanapot) write (outu,'(6x,a)') 'Forcing Analytical Potential for GCMC. Check that LJ and charges are defined properly.'
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'mass') then
    if (.not.Qatexp.and..not.Qpar.and..not.Qnucl) call error ('shell_simul', 'MASS order is defined before PTYPE, PARTICLE and/or NUCLEOTIDE order', faterr)
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        call gtdpar(com,'mass',etypl(itype)%mas,0.0)
      endif
    enddo

    write(outu,*)
    write(outu,'(6x,a)') 'Mass of element types:'
    write(outu,'(6x,a)') '---------------------'
    write(outu,'(6x,a)') 'type---mass(g/mol)---'
    write(outu,'(6x,a)') '---------------------'
    do  i = 1, netyp
      write(outu,'(6x,a,f12.4)') etypl(i)%nam,etypl(i)%mas
    enddo
    
    write(outu,*)
    write(outu,'(6x,a)') 'Mass of particle types:'
    write(outu,'(6x,a)') '-----------------------'
    write(outu,'(6x,a)') 'ptype----mass(g/mol)---'
    write(outu,'(6x,a)') '-----------------------'
    do i = 1, nptyp
      call updateptypmass(i)
      write(outu,'(6x,a,2f12.4)') ptypl(i)%nam,ptypl(i)%mass
    enddo
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'ljsin') then
    if (.not.Qatexp.and..not.Qpar.and..not.Qnucl) call error ('shell_simul', 'LJSIN order is defined before PTYPE, PARTICLE and/or NUCLEOTIDE order', faterr)
    if (Qljpar) call error ('shell_simul', 'LJSIN is defined after LJPAR', faterr)
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        call gtdpar(com,'epsilon',etypl(itype)%eps,0.0)
        call gtdpar(com,'sigma',etypl(itype)%sig,0.0)
      endif
    enddo
    Qljsin = .true.
  
    write(outu,*)
    write(outu,'(6x,a)') 'LJ Single Parameters:'
    write(outu,'(6x,a)') '---------------------'
    write(outu,'(6x,a)') 'type---epsilon(kcal/mol)---sigma(Ang)'
    write(outu,'(6x,a)') '-------------------------------------'
    do  i = 1, netyp
       write(outu,'(6x,a,2f12.4)') etypl(i)%nam,etypl(i)%eps,etypl(i)%sig
    enddo
    call setuplj()
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'ljpar') then
  !        ---------------     
    if (.not.Qatexp .and. .not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'LJPAR order is defined before PTYPE, PARTICLE and/or NUCLEOTIDE order', faterr)
  !         if (Qionsite) then
  !           call error ('shell_simul', 'Combination rules for LJ'
  !     & //' parameters are desactivated', warning)
  !           Qionsite = .false.
  !         endif
    if (allocated(epsLJ)) deallocate (epsLJ)
    if (allocated(sgLJ)) deallocate (sgLJ)
    allocate (epsLJ(netp),sgLJ(netp))
    epsLJ=0.0
    sgLJ=0.0
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        jtype=getetyp(wrd4)
        is=etpidx(itype,jtype)
        call gtdpar(com,'epsilon',epsLJ(is),0.0)
        call gtdpar(com,'sigma',sgLJ(is),0.0)
      endif
    enddo
    Qljpar = .true.
  
    write(outu,*)
    write(outu,'(6x,a)') 'LJ Pair Parameters:'
    write(outu,'(6x,a)') '-------------------'
    write(outu,'(6x,a)') 'type1---type2---epsilon(kcal/mol)---sigma(Ang)'
    write(outu,'(6x,a)') '----------------------------------------------'
    if (.not.allocated(epp4)) then
      allocate (epp4(netp),sgp2(netp))
      epp4=0.0
      sgp2=0.0
    endif
    call updateuetl()
    do i = 1, netyp
      do j = i, netyp
        is=etpidx(i,j)
        if (epsLJ(is).gt.0.0.and.sgLJ(is).gt.0.0) then
          write(outu,'(6x,a,a,2f12.4,$)') etypl(i)%nam,etypl(j)%nam,epsLJ(is),sgLJ(is) 
          if (Qlj(is)) then
            write(outu,'(a)') ' (replaced)'
          else
            write(outu,*)
            Qlj(is)=.true.
          endif
          epp4(is)=4.0*epsLJ(is)
          sgp2(is)=sgLJ(is)**2
        else
          if (Qlj(is)) then
            write(outu,'(6x,a,a,2f12.4,a)') etypl(i)%nam,etypl(j)%nam,epp4(is)*0.25,sqrt(sgp2(is)),' (from single LJ)'
          else
            if (etul(i).and.etul(j)) write(outu,'(6x,a,x,a,x,a)') 'Warning: Missing LJ parameters to compute pairs for:',etypl(i)%nam,etypl(j)%nam
          endif
        endif
      enddo
    enddo
    deallocate (epsLJ,sgLJ)
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'simul') then
  !        ---------------
     if (.not.Qatexp.and..not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'SIMULATION order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
     if (Qpres.and..not.Qpar) then
       write(outu,'(6x,a)') 'Warning: Pressure cannot be calculated if free particles are absent, so deactivated'
       Qpres=.false.
     endif
     if (Qpar.and..not.(Qefpott.or.Qljsin.or.Qljpar)) call error ('shell_simul', 'No effective potential and no Lennard Jones parameters defined', warning)
     ! number of steps for BD or MC simulations
     ! [integer,default=0]
     call gti8par(com,'ncycle',ncycle,0)
     if (ncycle.le.0) call error ('shell_simul', 'ncycle is equal or lower than zero in SIMUL order', faterr)
     ! number of steps for GCMC [integer,default=0]
     call gtipar(com,'ngcmc',ngcmc,0)
     if (ngcmc.lt.0) call error ('shell_simul','ngcmc is lower than zero in SIMUL order', faterr)
     if (ngcmc.gt.0 .and. .not.Qbuf) then
       call error ('shell_simul', 'GCMC is turned off because buffers have not been defined', warning)
       ngcmc = 0
     endif 
     ! number of steps for MC [integer,default=0]
     call gtipar(com,'nmcm',nmcm,0)
     if (nmcm.lt.0) call error ('shell_simul', 'nmcm is lower than zero in SIMUL order', faterr)
     if (nmcm.gt.0 .and. .not.Qpar) then
       call error ('shell_simul', 'Metropolis MonteCarlo is turned off because ions have not been defined', warning)
       nmcm = 0
     endif 
     if (nmcm.gt.0) then
     ! maximum displacement (mcmax) for MC
     ! [real*8,default=1]
       call gtdpar(com,'mcmax',mcmax,0.5)
       if (mcmax.lt.0.0) call error ('shell_simul', 'mcmax is lower than zero in SIMUL order', faterr)
       mcmax=2.0*mcmax
     endif
     ! number of steps for BD [integer,default=0]  
     call gtipar(com,'nbd',nbd,0)
     if (nbd.lt.0) call error ('shell_simul', 'nbd is lower than zero in SIMUL order', faterr)
     ! maximum displacement (bdmax) for BD
     if (nbd.gt.0) then
       call gtdpar(com,'bdmax',bdmax,1.0)
       if (bdmax.lt.0.0) call error ('shell_simul', 'bdmax is lower than zero in SIMUL order', faterr)
     endif
     ! Logical variable which indicates if a trajectory file
     Qchdencnt = check(com,'chden')
     if (.not.Qchden) Qchdencnt=.false.
     ! Logical variable which indicates if a trajectory file
     ! will be written
     Qtraj = check(com,'traject')
     if (Qtraj) then
       ! unit number of a trajectory file [integer,default=1]
       Qtrajcont=check(com,'trajcont')
       call gtipar(com,'setframes',setframes,0)
       call gtipar(com,'iuntrj',iuntrj,0)
       if (iuntrj.le.0) call error ('shell_simul', 'iuntrj is zero or a negative number', faterr)
       if (iuntrj.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntrj).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntrj = unvec(iuntrj)
     endif
     ! trajectory saving frequency [integer,default=ncycle]
     call gti8par(com,'nsave',nsave,ncycle)
     if (nsave.le.0) call error ('shell_simul', 'nsave is lower than or equal to zero in SIMUL order', faterr)
     if (mod(ncycle,nsave).ne.0) call error ('simul1', 'nsave is not correct. ncycle must be divisible by nsave', faterr)
  
     ! time-step for BD [real*8,default=0.02]
     call gtdpar(com,'dt',dt,0.02)
     if (dt.le.0.0) call error ('shell_simul', 'dt is lower or equal than zero', faterr)
    
     ! Set precalculated diffusion coefficient factors
     kBTdt = dt * ikbt
     do i=1,netnuc
       fact1a(i)=etypl(i)%dif*kbtdna*dt
       fact2a(i)=sqrt(2.0*dt*etypl(i)%dif)
     enddo
     do i=1+netnuc,netyp
       fact1a(i)=etypl(i)%dif*kBTdt
       fact2a(i)=sqrt(2.0*dt*etypl(i)%dif)
     enddo
     write(outu,'(6x,a)') 'Diffusion Factors per Element Type'
     do i=1,netyp
       write(outu,'(6x,i4,x,A4,3(x,f16.8))') i, etypl(i)%nam, etypl(i)%dif, fact1a(i), fact2a(i)
     enddo
     if (Qproxdiff) fact2pd=sqrt(2.0*dt*diff0)
 
     ! print frequency [integer,default=0]
     call gtipar(com,'nprint',nprint,0)
     if (nprint.lt.0) call error ('shell_simul', 'nprint is a negative number',faterr)

     ! Z-position where the number of ion crossing the channel 
     Qcountion=check(com,'countions').and.Qpar
     if (Qcountion) then
       call gtipar(com,'svcntfq',svcntfq,nprint)
       call gtipar(com,'iuncnt',iuncnt,0)
       if (iuncnt.le.0) call error ('shell_simul', 'iuncnt is zero or a negative number', faterr)
       if (iuncnt.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuncnt).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuncnt = unvec(iuncnt)
       call gtcrpar(com,'zcont',cntpts,word)
       if (cntpts.le.0) call error ('shell_simul','zcont cannot be empty or ommited if countions present',faterr)
       allocate (zcont(cntpts),nbackward(nptnuc+1:nptyp,cntpts),nforward(nptnuc+1:nptyp,cntpts))
       call gtcdpar(word,zcont)
     endif

     ! Logical variable which indicates if radial distribution 
     ! function will be printed
     Qgr      = check(com,'rdf')
     if (Qgr .and. .not.Qpar) then
       call error ('shell_simul', 'radial distribution function is turned off because ions have not been defined',warning)
       Qgr = .false.
     endif
     ! Logical variable which indicates if average density profile 
     ! along the Z-axis will be printed
     Qrho     = check(com,'rho')
     if (Qrho .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul', 'average density profile along the Z-axis is turned off because ions and/or system have not been defined', warning)
       Qrho = .false.
     endif
     ! Logical variable which indicates if average density profile 
     ! along the Z-axis will be printed for DNA sites
     Qrdna    = check(com,'rhDNA')
     if (Qrdna .and. (.not.Qnucl.or..not.Qsystem)) then
        call error ('shell_simul', 'average density profile along the Z-axis is turned off because sites and/or system have not been defined', warning)
       Qrdna = .false.
     endif
     ! Logical variable which indicates if various average probability 
     ! distributions will be printed
     Qprob    = check(com,'prob')
     if (Qprob .and. .not.Qpar) then
        call error ('shell_simul', 'average probability distributions are turned off because ions have not been defined', warning)
       Qprob = .false.
     endif 
     ! ion pairing analysis frequency [integer,default=1]
     call gtipar(com,'nanal',nanal,1)
     if (nanal.le.0) call error ('shell_simul', 'nanal is lower or equal than zero in SIMUL order', faterr)
     ! Logical variable which indicates if ion pairing analysis 
     ! (S frequency) will be made
     Qionpair = check(com,'ionpair')
     if (Qionpair .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul','ion pairing analysis is turned off because ions and/or system have not been defined', warning)
       Qionpair = .false.
     endif
     ! Logical variable which indicates if energy profile along 
     ! the Z-axis will be made
     Qenerprofile = check(com,'enerprofile')
     if (Qenerprofile .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul', 'energy profile along the Z-axis is turned off because ions and/or system have not been defined', warning)
       Qenerprofile = .false.
     endif
     ! Logical variable which indicates if fraction of denatured 
     ! bases will be calculated after simulation of DNA
     Qfbases = check(com,'denatured') 
     ! Logical variable which indicates if translocation duration
     ! for DNA through pore will be calculated
     Qfmemb = check(com,'translocation') .and. Qnucl 
     ntras = 1
     if (Qfmemb) then
       ! unit number of a translocation duration file [integer,default=1]
       call gtipar(com,'unitn',iuntfm,1)
       if (iuntfm.le.0) call error ('shell_simul', 'iuntfm is zero or a negative number', faterr)
       if (iuntfm.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntfm).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntfm = unvec(iuntfm)
       call gtipar(com,'ntras',ntras,1)
       ! translocation saving frequency [real*8,deafult=1]
       if (ntras.le.0) call error ('shell_simul', 'ntras < = 0 in SIMUL order',faterr)
       if (mod(ncycle,ntras).ne.0) call error ('simul1', 'ntras is not correct. ncycle must be divisible by ntras', faterr)
     endif
     if (Qprob .or. Qfmemb) then
       ! maximum position of a channel along the 
       ! Z-axis [real*8,default=0]
       call gtdpar(com,'czmax',czmax,0.0)
       ! minimum position of a channel along the 
       ! Z-axis [real*8,default=0]     
       call gtdpar(com,'czmin',czmin,0.0) 
       if (czmin.gt.czmax) call error ('shell_simul', 'minimum position is greater than maximum position of the channel', faterr)
     endif
     if (Qgr) then
       call gtipar(com,'ion',igr,nparnuc+1)
       if (igr.le.0 .or. igr.gt.(nelenuc)) call error ('shell_simul', 'fixed reference ion is not adequate in SIMUL order', faterr)
     endif
     nsfbs = ncycle
     vfbs = 0.0
     if (Qfbases) then 
       if (.not.Qnucl) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because DNA has not been defined', faterr) 
       if (istrs.ne.2) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because DNA has only one strand', faterr) 
       if (inuc.eq.0) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because the number of nucleotides is zero',faterr) 
       ! unit number of a fraction of denatured bases [integer,default=1]
       call gtipar(com,'iunfbs',iunfbs,1) 
       if (iunfbs.le.0) call error ('shell_simul', 'iunfbs is zero or negative number', faterr)
       if (iunfbs.gt.maxopen) call error ('shell_simul', 'iunfbs is greater than maxopen', faterr)
       if (unvec(iunfbs).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order', faterr)
       iunfbs = unvec(iunfbs)
       ! fraction of denatured bases saving frequency
       ! [integer, default=ncycle]
       call gti8par(com,'nsfbs',nsfbs,ncycle)
       if (nsfbs.le.0) call error ('shell_simul', 'nsfbs < = 0 in SIMUL order',faterr)
       if (mod(ncycle,nsfbs).ne.0) call error ('simul1', 'nsfbs is not correct. ncycle must be divisible by nsfbs', faterr)
  
       ! initial time for calculating  fraction of denatured
       ! [real*8,deafult=1]
       call gtdpar(com,'vfbs',vfbs,0.0)
       if (vfbs.lt.0.0) call error ('shell_simul', 'vfbs < 0 n SIMUL order',faterr)
     endif 
     ! Security ouputfile contains coordinates and seed numbers 
     Qsec = check(com,'security') 
     nsec = 1
     if (Qsec) then
       ! unit number of a security file [integer,default=1]
       call gtipar(com,'units',iuntsc,1)
       if (iuntsc.le.0) call error ('shell_simul', 'iuntfm is zero or a negative number', faterr)
       if (iuntsc.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntsc).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntsc = unvec(iuntsc)
       call gtipar(com,'nsec',nsec,1)
       ! saving frequency security [real*8,deafult=1]
       if (nsec.le.0) call error ('shell_simul', 'nsec < = 0 in SIMUL order',faterr)
       if (mod(ncycle,nsec).ne.0) call error ('simul1', 'nsec is not correct. ncycle must be divisible by nsec', faterr)
     endif
     ! turn off total energy
     Qenergy  = .not.check(com,'noenergy')
     ! turn off nonbond energy
     Qnonbond = .not.check(com,'nononbond')
     ! turn off bond energy
     Qbond = .not.check(com,'nobond')
     ! Restrain Internal Forces
     call gtdpar(com,'riffac',riffac,-2.0)
     if (riffac.ge.-1.0) then
         ! Enable it
         Qresintfor = .true.
         ! Auto value if factor is between [-1.0, 0.0)
         if (riffac.lt.0.0) then
             if (dt.le.0.01) then 
                 Qresintfor = .false.  ! Auto value disabled
             else
                 riffac = 0.01 / dt    ! Auto value
             endif
         ! Disable bonded interactions since they will be scaled to 0
         else if (riffac.eq.0.0) then 
             Qbond = .false.
         endif
     ! If lower than -1.0 disable it
     else
         Qresintfor = .false.
     endif
     if (Qresintfor) then
       call checkmass()
       write(outu,'(6x,a,f10.6)') 'Restrained Internal Forces enabled. RIF Factor= ',riffac
     endif

     call simul1(ncycle, ngcmc, nmcm, nbd, nsave, nsfbs, vfbs, ntras, nsec, iseed)
     
     if (sum(warn).gt.0) then
       write(outu,*) '      Warnings summary:'
       do i=1,netyp
         write(outu,*) '               ',etypl(i)%nam,warn(i)
       enddo
     endif
  ! **********************************************************************
  elseif (wrd5.eq.'energ') then
  !        ---------------
     if (.not.Qatexp.and..not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'ENERGY order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
    ! default values 
    ! Qbond = Qnonbond =.true.
    ! Qmemb = Qphix = Qphiv = Qsrpmf = .false.
    Qbond  = .not.check(com,'nobond')
    Qnonbond = .not.check(com,'nononbond')
    if (check(com,'membrane')) then
      if (.not.Qmemb) call error ('shell_simul', 'Planar membrane has not be defined before ENERGY order', faterr)
    else
      logmemb = Qmemb
      Qmemb = .false.
    endif
    if (check(com,'phix')) then
      if (.not.Qphix) call error ('shell_simul', 'Static Field has not be defined before ENERGY order', faterr)
    else
     logphix = Qphix
     Qphix = .false.
    endif
    if (check(com,'phiv')) then
      if (.not.Qphiv) call error ('shell_simul', 'Repulsive term has not be defined before ENERGY order', faterr)
    else
      logphiv = Qphiv
      Qphiv = .false.
    endif
    if (check(com,'srpmf')) then
      if (.not.Qsrpmf) call error ('shell_simul', 'Short-range interaction term has not be defined before ENERGY order', faterr)
    else
      logsrpmf = Qsrpmf
      Qsrpmf = .false.
    endif
    if (check(com,'rfpar')) then
      if (.not.Qrfpar) call error ('shell_simul', 'Reaction Field parameter term has not be defined before ENERGY order', faterr)
    else
      logrfpar = Qrfpar
      Qrfpar = .false.
    endif
    if (check(com,'all')) then 
      Qmemb = logmemb
      Qphix = logphix
      Qphiv = logphiv
      Qsrpmf = logsrpmf
      Qrfpar = logrfpar
      Qbond = .true.
      Qnonbond = .true. 
    endif
 
    write(outu,*)
    write(outu,'(6x,a)') 'ENERGY calculation'
    write(outu,'(6x,a)') '------------------'
    if (Qbond)    write(outu,'(6x,a)') 'Bonding Energy Term'
    if (Qnonbond) write(outu,'(6x,a)') 'Nonbonding Energy Term'
    if (Qmemb)    write(outu,'(6x,a)') 'Planar membrane Term'
    if (Qphix)    write(outu,'(6x,a)') 'External Field Energy Term'
    if (Qphiv)    write(outu,'(6x,a)') 'Repulsive Energy Term'
    if (Qsrpmf)   write(outu,'(6x,a)') 'Short-range Interaction Term'
    if (Qrfpar)   write(outu,'(6x,a)') 'Reaction Field Parameter Term'
  
    Qforces=.true.
    call energy()
    write(outu,'(6x,a,f16.8)') 'Total energy (Forces On) ',ener
    write(outu,'(6x,a)') 'ememb, estaticf, evdwgd, erfpar, eintern, enonbond'
    write(outu,'(6x,6(x,f16.8))') ememb,estaticf,evdwgd,erfpar,eintern,enonbond
    if (check(com,'forces')) then
      write(outu,*)
      write(outu,'(6x,a)') '#, EleNum, Force X, Force Y, Force Z'
      do i=1,nele
        write(outu,'(6x,2i6,3(x,f16.8))') i,et(i),f(i)%x,f(i)%y,f(i)%z
      enddo
    endif
    Qforces=.false.
    call energy()
    write(outu,'(6x,a,f16.8)') 'Total energy (Forces Off) ',ener
    write(outu,'(6x,a)') 'ememb, estaticf, evdwgd, erfpar, eintern, enonbond'
    write(outu,'(6x,6(x,f16.8))') ememb,estaticf,evdwgd,erfpar,eintern,enonbond
    Qmemb = logmemb 
    Qphix = logphix 
    Qphiv = logphiv 
    Qsrpmf = logsrpmf 
    Qrfpar = logrfpar
  ! **********************************************************************
  elseif (wrd5.eq.'inter') then 
  !        ---------------
     if (.not.Qatexp.and..not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'INTERACT order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
    ! atom type [integer,default=1]
    call gtipar(com,'part',parn,1) 
    ! default values
    ! Qnobond = Qnonbond =.true.
    ! Qmemb = Qphix = Qphiv = Qsrpmf = .false.         
    Qbond  = .not.check(com,'nobond')
    Qnonbond = .not.check(com,'nononbond')
    if (check(com,'membrane')) then
      if (.not.Qmemb) call error ('shell_simul', 'Planar membrane has not be defined before INTERACT order', faterr)
    else
      logmemb = Qmemb
      Qmemb = .false.
    endif
    if (check(com,'phix')) then
      if (.not.Qphix) call error ('shell_simul', 'Static Field has not be defined before INTERACT order', faterr)
    else
     logphix = Qphix
     Qphix = .false.
    endif
    if (check(com,'phiv')) then
      if (.not.Qphiv) call error ('shell_simul', 'Repulsive term has not be defined before INTERACT order', faterr)
    else
      logphiv = Qphiv
      Qphiv = .false.
    endif
    if (check(com,'srpmf')) then
      if (.not.Qsrpmf) call error ('shell_simul', 'Short-range interaction term has not be defined before INTERACT order', faterr)
    else
      logsrpmf = Qsrpmf
      Qsrpmf = .false.
    endif
    if (check(com,'rfpar')) then
      if (.not.Qrfpar) call error ('shell_simul', 'Reaction Field Parameter term has not be defined before INTERACT order', faterr)
    else
      logrfpar = Qrfpar
      Qrfpar = .false.
    endif
    if (check(com,'allener')) then
      Qmemb = logmemb
      Qphix = logphix
      Qphiv = logphiv
      Qsrpmf = logsrpmf
      Qrfpar = logrfpar
      Qbond = .true.
      Qnonbond = .true.
    endif
 
    write(outu,*)
    write(outu,'(6x,a)') 'INTERACT calculation'
    write(outu,'(6x,a)') '--------------------'         
    if (Qbond)    write(outu,'(6x,a)') 'Bonding Energy Term' 
    if (Qnonbond) write(outu,'(6x,a)') 'Nonbonding Energy Term'
    if (Qmemb)    write(outu,'(6x,a)') 'Planar membrane Term'
    if (Qphix)    write(outu,'(6x,a)') 'External Field Energy Term'
    if (Qphiv)    write(outu,'(6x,a)') 'Repulsive Energy Term'
    if (Qsrpmf)   write(outu,'(6x,a)') 'Short-range Interaction Term'
    if (Qrfpar)   write(outu,'(6x,a)') 'Reaction Field Parameter Term'
  
    ! Calculate the interaction of particle "parn" with the rest of
    ! the system
    write(outu,'(6x,a)') 'Particle Number, Total Energy, emembi, erfpari, estaticfi, evdwgdi, enonbondi'
    if (check(com,'allpart')) then
      do parn=1,npar
        call par_interact(parn,dener)
        write(outu,'(6x,i6,6(x,f18.6))') parn,dener,emembi,erfpari,estaticfi,evdwgdi,enonbondi
      enddo
    else
      call par_interact(parn,dener)
      write(outu,'(6x,i6,6(x,f18.6))') parn,dener,emembi,erfpari,estaticfi,evdwgdi,enonbondi
    endif
    Qmemb = logmemb 
    Qphix = logphix 
    Qphiv = logphiv 
    Qsrpmf = logsrpmf
    Qrfpar = logrfpar
  ! **********************************************************************
  elseif (wrd5.eq.'membr') then
  !        ---------------
     if (.not.Qatexp.and..not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'MEMBRANE order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'MEMBRANE order is defined before SYSTEM order', faterr)
     ! Turn on/off repulsive membrane potential in cylindrical pore
     ! region
     Qpore = check(com,'pore')
     ! Thickness of a membrane [real*8,default=0]
     call gtdpar(com,'thick',thick2,0.0)
     if (thick2.le.0.0) call error ('shell_simul', 'Thickness of a membrane is negative or zero', faterr)
     ! Position of center of a membrane along the Z-axis
     ! [real*8,default=0]
     call gtdpar(com,'zmemb',zmemb,0.0)
     if (zmemb.gt.lz2p.or.zmemb.lt.lz2m) call error ('shell_simul', 'Position of center of a membrane along the Z-axis is not correct', faterr)
     tmemb  =  thick2 ! thickness of the membrane
     thick2 =  thick2*0.5
     zmemb1 = -thick2+zmemb ! lower limit of the membrane
     zmemb2 =  thick2+zmemb ! upper limit of the membrane
     if (zmemb2.gt.lz2p) call error ('shell_simul', 'Upper limit of the membrane along the Z-axis is not correct', faterr)
     if (zmemb1.lt.lz2m) call error ('shell_simul', 'Lower limit of the membrane along the Z-axis is not correct', faterr)
     ! Dielectric constant of the membrane [real*8,default=2]
     call gtdpar(com,'epsm',epsm,2.0)
     if (epsm.le.0.0) call error ('shell_simul', 'Dielectric constant of the membrane has a null or negative value', faterr)
     if (ionstr.le.0.0) call error ('shell_simul', 'Ion Strength cannot be lower or equal zero. Define it in SYSTEM.', faterr)
     ! Setup constant
     afact = epsm*voltage/(2.0*epsm+cdie*ikappa*tmemb) ! Eq. (32) paper
     ceps=cdie/epsm
     ampl1 = 0.0
     p1 = 1.0
     ampl2 = 0.0
     p2 = 1.0
     rcylinder = 1.0
     endlog = .false.
     do while (.not.endlog)
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type     
         call getfirst(com,wrd4)
         itype=getetyp(wrd4)
         ! Repulsive membrane potential [real*8,default=0]
         call gtdpar(com,'amplmemb',ampl1(itype),0.0)
         if (ampl1(itype).lt.0.0) call error ('shell_simul', 'Repulsive membrane potential has a negative value', faterr)
         ! Switching region between bulk and membrane regions
         ! [real*8,default=1]
         call gtdpar(com,'pmemb',p1(1,itype),1.0)
         if (p1(1,itype).le.0.0) call error ('shell_simul', 'Switching region between bulk and membrane regions has a negative or zero value', faterr)
         p1(2,itype)=1.0/p1(1,itype)
         if (Qpore) then 
           ! Repulsive potential for a cylindrical pore
           ! [real*8,default=0]
           call gtdpar(com,'amplpore',ampl2(itype),0.0)
           if (ampl2(itype).le.0.0) call error ('shell_simul', 'Repulsive cylindrical pore potential has a negative or zero value', faterr)
           ! Switching region between pore and membrane regions
           ! [real*8,default=1]
           call gtdpar(com,'ppore',p2(1,itype),1.0)
           if (p2(1,itype).le.0.0) call error ('shell_simul', 'Switching region between pore and membrane regions has a negative or zero value', faterr)
           p2(2,itype)=1.0/p2(1,itype)
           ! Radius of a cylindrical pore
           call gtdpar(com,'radi',rcylinder(itype),1.0)
           if (rcylinder(itype).le.0.0) call error ('shell_simul', 'Radius of a cylindrical pore has a negative or zero value', faterr)
         endif
       endif
     enddo
     Qmemb = .true.
     do i=1,netyp
       if (ampl1(i).lt.0.0) then 
         write(*,'(6x,a)') 'ERROR: Invalid membrane constant for type '//etypl(i)%nam
         Qmemb=.false.
       endif
       if (ampl2(i).lt.0.0.and.Qpore) then
         write(*,'(6x,a)') 'ERROR: Invalid Pore constant for type '//etypl(i)%nam
         Qmemb=.false.
       endif
     enddo
     if (.not.Qmemb) stop
  
     write(outu,*)
     write(outu,'(6x,a)') 'MEMBRANE parameters:'
     write(outu,'(6x,a)') '--------------------'
     write(outu,'(6x,a,f12.3,a)') 'Membrane thickness ',tmemb,' [Angs]'
     write(outu,'(6x,a,f12.3,a)') 'Membrane center    ',zmemb,' [Angs]'
     write(outu,'(6x,a,f12.3,a)') 'Voltage            ',voltage*kcalmol/Coulomb,' [Volts]'
     write(outu,'(6x,a,f12.3)') 'Dielectric constant of the membrane ',epsm
     write(outu,'(6x,a,f12.3,a)') 'Ionic strength     ',ionstr,' [mol/L]' 
     write(outu,'(6x,a)') 'AMPLMEMB -> repulsive membrane potential [Kcal/mol]'
     write(outu,'(6x,a)') 'PMEMB -> switching region between bulk and membrane regions [Ang.]'
     if (Qpore) then  
       write(outu,'(6x,a)') 'PORE Enabled'
       write(outu,'(6x,a)') 'AMPLPORE -> repulsive potential for a cylindrical pore [Kcal/mol]'
       write(outu,'(6x,a)') 'PPORE -> switching region between pore and membrane regions [Ang.]'
       write(outu,'(6x,a)') 'RCYL -> radius of a cylindrical pore [Ang.]'
       write(outu,'(6x,a)') 'NAME---AMPLMEMB---PMEMB---AMPLPORE---PPORE---RCYL'   
     else
       write(outu,'(6x,a)') 'NAME---AMPLMEMB---PMEMB'   
     endif
     do i = 1, netyp
       if (Qpore) then 
         write(outu,'(6x,a,1x,5f12.3)') etypl(i)%nam,ampl1(i),p1(1,i),ampl2(i),p2(1,i),rcylinder(i)
       else
         write(outu,'(6x,a,1x,5f12.3)') etypl(i)%nam,ampl1(i),p1(1,i)
       endif
     enddo
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'srpmf') then ! short-range ion-ion interactions
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'SRPMF order is defined before PARTICLE order', faterr)
     if (Qefpott) call error ('shell_simul', 'SRPMF must be defined before EFPOT', faterr)
     allocate (c0(netp),c1(is),c2(netp),c3(is),c4(netp))
     c0=0.0
     c1=0.0
     c2=0.0
     c3=0.0
     c4=0.0
     call gtdpar(com,'rth',rth,8.0)
     srpx=rth-0.5
     srpk=6.90775527898/0.5
     srpy=0.001
     rth=rth**2
     endlog = .false.
     do while (.not. endlog)
       call getlin(com,inpu,outu) ! new commands line
       itype = 0 ! counter of ion types
       jtype = 0 ! counter of ion types
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type 
         call getfirst(com,wrd4)
         itype=getetyp(wrd4)
         ! Obtention of ion type
         call getfirst(com,wrd4)
         jtype=getetyp(wrd4)
         is=etpidx(itype,jtype)
         ! c_0 parameter [real*8,default=0]
         call gtdpar(com,'c0',c0(is),0.0)
         ! c_1 parameter [real*8,default=0]
         call gtdpar(com,'c1',c1(is),0.0)
         ! c_2 parameter [real*8,default=0]
         call gtdpar(com,'c2',c2(is),0.0)
         ! c_3 parameter [real*8,default=0]
         call gtdpar(com,'c3',c3(is),0.0)
         ! c_4 parameter [real*8,default=0]
         call gtdpar(com,'c4',c4(is),0.0)
         ! Logical variable which indicates if total ion-ion 
         ! interaction energy file will be written  
         Qlsprmf = check(com,'file')
         if (Qlsprmf) then
           ! Unit number to write the total ion-ion interaction energy
           ! [integer,default=1]             
           call gtipar(com,'unit',iunit,1)
           if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
           if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
           if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in SRPMF order', faterr)
           iunit = unvec(iunit)
           is=etpidx(itype,jtype)
           ! Parameters of the Lennard-Jones 6-12 potential 
           cc0 = c0(is)
           cc1 = c1(is)
           cc2 = c2(is)
           cc3 = c3(is)
           cc4 = c4(is)
           ! r1 -> Interaction distance. This functions is applied
           ! only up to r1=8.0 Angstrom (0.05 increment)
           do r1i = 2, 160
             r1=real(r1i)*0.05
             r2 = r1**2
             r6 = (sgp2(is)/r2)**3
             ! OBTENTION OF LENNARD-JONES 6-12 POTENTIAL
             evdw  = epp4(is)*r6*(r6-1.0)
             ! OBTENTION OF ELECTROSTATIC INTERACTION BETWEEN IONS
             ! OBTENTION OF WATER-MEDIATED SHORT-RANGE ION-ION
             ! INTERACTION  
             esrpmf = cc0*exp((cc1-r1)/cc2)*cos(cc3*pi*(cc1-r1))+cc4*(cc1/r1)**6
             if (evdw+esrpmf.le.50.0) write(iunit,'(4f12.5)') r1,evdw+esrpmf,evdw,esrpmf
           enddo
         endif
       endif
     enddo
     Qsrpmf = .true.
  
     write(outu,'(6x,a)')'Short-range ion-ion interactions are turn on' 
     write(outu,'(6x,a,/,6x,a)') 'Coefficients for the short-range ion-ion interactions','NAME----NAME----C0----C1----C2----C3----C4' 
     do i = 1,netyp
       do j=i,netyp
         is=etpidx(i,j)
         if (c2(is).ne.0.0) then
           write(outu,'(6x,a,1x,a,5f8.3)') etypl(i)%nam,etypl(j)%nam,c0(is),c1(is),c2(is),c3(is),c4(is)
           c2(is)=1.0/c2(is)
           Qsrpmfi(is)=.true.
         else
            if (etul(i).and.etul(j)) write(outu,'(6x,a,x,a,x,a)') 'Warning: Missing SRPMF parameters for:',etypl(i)%nam,etypl(j)%nam
         endif
       enddo
     enddo
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'efpot') then ! effective potential
  !       ---------------
     if (.not.Qatexp.and..not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'EFPOT order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
    cnt=0
    Qepwrt=check(com,'write')
    ! Assume 0 charge for all particles
    if (check(com,'ignoreq')) then
      Qcol=.false.
      write(outu,'(6x,a,i0)') 'EFPOT: Coulombic terms will be assumed 0 in head and tails (except for built potentials)'
    endif
    if (Qepwrt) then 
      call gtipar(com,'unit',wunit,1)
      wunit = unvec(wunit)
    endif
    call gtdpar(com,'scal',scald,kBT)
    call gtdpar(com,'res',res,0.1)
    ires=1.0/res
    if (Qsrpmf) maxl=sqrt(rth)
    call gtdpar(com,'maxdist',maxlg,maxl)
    mnp=int(maxlg*ires)
    call gtdpar(com,'mindist',minlg,2.5)
    nnp=int(minlg*ires)
    minlg=float(nnp)*res
    maxlg=float(mnp)*res
    maxd=int(float(netp)*maxlg*ires)
    call gtipar(com,'maxdata',maxd,maxd)     
    write(outu,'(6x,a,i0)') 'Max number of data points to be stored (temporary storage array): ',maxd 
    if (allocated(efp)) deallocate(efp)
    allocate(efp(netp))
    allocate (xy(2,maxd),Qefpotread(netp),nxf(netp))
    call updateuetl()
    Qefpotread=.false.
    write(outu,'(6x,a)')
    endlog = .false.
    call getetchg() ! Build Element Types Charge List
    allocate (fct(netp))
    fct=0.0
    do i=1,netyp
      do j=i,netyp
        is=etpidx(i,j)
        fct(is)=celec*etchg(i)*etchg(j)/cdie
      enddo
    enddo
    write(outu,'(6x,a)')'WARNING: Do not use EFPOT with useq or build if are using CHARMM Force Field/PSF and a same element type has different charge'
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      itype = 0 ! counter of ion types
      jtype = 0 ! counter of ion types
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        jtype=getetyp(wrd4)
        is=etpidx(itype,jtype)
        if (check(com,'read')) then
          if (check(com,'ignoreq')) then
            Qcol(is)=.false.
            write(outu,'(6x,a,i0)') 'EFPOT: Coulombic terms will be assumed 0 in head and tails for this pair'
          endif
          if (check(com,'useq')) then 
            Qcol(is)=.true.
            write(outu,'(6x,a,i0)') 'EFPOT: Coulombic terms will be used in head and tails for this pair'
          endif
          Qefpot(is) = .true.
          Qefpotread(is) = .true.
          ! Unit number to read the pair effective potential
          ! [integer,default=1]             
          call gtipar(com,'unit',iunit,1)
          if (iunit.le.0) call error('shell_simul','unit is zero or a negative number', faterr)
          if (iunit.gt.maxopen) call error('shell_simul','unit is greater than maxopen', faterr)
          if (unvec(iunit).eq.-1) call error('shell_simul','unit incorrect in EFPOT order', faterr)
          iunit = unvec(iunit)
          call gtdpar(com,'scal',scaldd,scald)
          write(outu,*) '     Scaling ',etypl(itype)%nam,'-',etypl(jtype)%nam,' potential by ',scaldd
          cnt=1
          if (cnt.gt.maxd) call error ('shell_simul','Number of data is greater than expected. Increase maxdata.',faterr) ! Check if maxdata is correct
          read(iunit,*,IOSTAT=kode) xy(1,cnt),xy(2,cnt)
          xy(2,cnt)=xy(2,cnt)*scaldd
          do while (kode.eq.0)
            ! Check if resolution is correct {
            if (cnt.ge.2) then 
              if (xy(1,cnt)-xy(1,cnt-1).lt.0.99*res.or.xy(1,cnt)-xy(1,cnt-1).gt.1.01*res) call error ('shell_simul','Unexpected x spacing. Check that all x data spacing is separated by res.',faterr)
            endif
            ! }
            cnt=cnt+1
            if (cnt.gt.maxd) call error ('shell_simul','Number of data is greater than expected. Increase maxdata.',faterr) ! Check if maxdata is correct
            read(iunit,*,IOSTAT=kode) xy(1,cnt),xy(2,cnt)
            xy(2,cnt)=xy(2,cnt)*scaldd
          enddo
          cnt=cnt-1
          nxf(is)=cnt
          if (allocated(efp(is)%ep)) deallocate (efp(is)%ep)
          allocate (efp(is)%ep(cnt))
          call splinepot(is,cnt,xy(1,1:cnt),xy(2,1:cnt),fct(is))
        elseif (check(com,'build')) then
          if (Qlj(is)) then 
            Qefpot(is)=.true.
            call gtdpar(com,'maxdist',maxl,maxlg)
            mnp=int(maxl*ires)
            call gtdpar(com,'mindist',minl,minlg)
            nnp=int(minl*ires)
            cnt=1+mnp-nnp
            nxf(is)=cnt
            efp(is)%xl=nnp*res
            efp(is)%xl2=(nnp*res)**2
            efp(is)%xu2=(mnp*res)**2
            if (allocated(efp(is)%ep)) deallocate (efp(is)%ep)
            allocate (efp(is)%ep(cnt))
          else
            call error ('shell_simul','Cannot build potential, LJ parameters missing.',faterr)
          endif
        endif
      endif
    enddo
    call updateuetl()
    write(outu,'(6x,a)') 'Effective potential was activated for the following pairs:'
    write(*,*) 'Used Element Type-Element Type pairs'
    do i = netnuc+1,netyp
      do j=1,i
        is=etpidx(i,j)
        if (Qefpot(is).and.Qefpotread(is)) then
          write(outu,'(6x,a,1x,2a,i5,2(a,f8.3))')etypl(i)%nam,etypl(j)%nam,'  Number of Points:',nxf(is),'  From-To: ',efp(is)%xl,' - ',sqrt(efp(is)%xu2)
        elseif (Qefpot(is).and..not.Qefpotread(is)) then
          Qcol(is) = .true. ! if discretize is active, ignore this keyword for that pair 
          nnp=int(efp(is)%xl*ires)
          mnp=int(sqrt(efp(is)%xu2)*ires)
          call discretize(is,nnp,mnp,nxf(is),fct(is))
          write(outu,'(6x,a,1x,2a,i5,2(a,f8.3),a)')etypl(i)%nam,etypl(j)%nam,'  Number of Points:',1+mnp-nnp,'  From-To: ',efp(is)%xl,' - ',mnp*res,'  (built from current parameters)'
        elseif (etul(i).and.etul(j)) then
          write(outu,'(6x,a,x,a,x,a,a)')'WARNING: Missing Effective Potential for:',etypl(i)%nam,etypl(j)%nam,' . Using Coulombic and Lennard Jones and/or SRPMF parameters on the fly'
        endif
      enddo
    enddo
    write(outu,*)'     Resolution: ',res
    if (Qepwrt) then
      do i=1,netyp
        do j=i,netyp
          is=etpidx(i,j)
          if (Qefpot(is)) then
            write(wunit,*) etypl(i)%nam,' - ',etypl(j)%nam
            inflim=efp(is)%xl-1.0
            if (inflim < 0.0) inflim=res*0.1
            do k=int(inflim*ires*10.0),int((sqrt(efp(is)%xu2)+10.0)*ires*10.0)
              x1=k*res*0.1
              call getyd(is,x1**2,x2,x3,r1)
              write(wunit,*) x1,x2,x3
            enddo
          endif
        enddo
      enddo
    endif
    deallocate (xy,nxf,Qefpotread) 
    Qefpott=.true.
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'proxd') then ! user defined proximity diffusion
  !        -----------------     
    if (.not.(Qpar.and.Qnucl)) call error ('shell_simul','PROXDIFF need IONS and DNA',faterr)
    if (.not.(Qefpott)) call error ('shell_simul','PROXDIFF must be defined after EFPOT',faterr)
    call gtdpar(com,'beta',beta,2.93)
    call gtdpar(com,'diff0',diff0,0.101)
    call gtdpar(com,'diffcutoff',diffcutoff,maxl)
    ibeta=1.0/beta
    Qproxdiff=.true.
    write(outu,'(6x,a)') 'PROXDIFF ACTIVATED:'
    write(outu,'(6x,a)') '------------------'
    write(outu,'(6x,a,f8.3,a)') 'beta  =      ',beta,' [Angs]'
    write(outu,'(6x,a,f8.3,a)') 'diff0 =      ',diff0,' [Angs^2/ps]'
    write(outu,'(6x,a,f8.3,a)') 'diffcutoff = ',diffcutoff,' [Angs]'
  elseif (wrd5.eq.'chden') then ! user defined proximity diffusion
  !        -----------------    
    if (.not.Qsystem) call error ('shell_simul', 'CHDEN must be defined after SYSTEM', faterr)
    Qchdenorm = check(com,'norm')
    call gtdpar(com,'dcel',dcel4,0.5)
    call gtdpar(com,'xbcen',xbcen4,cx)
    call gtdpar(com,'ybcen',ybcen4,cy)
    call gtdpar(com,'zbcen',zbcen4,cz)
    nclx4=int(lx/dcel4)+1
    ncly4=int(ly/dcel4)+1
    nclz4=int(lz/dcel4)+1
    call gtipar(com,'nclx',nclx4,nclx4)
    call gtipar(com,'ncly',ncly4,ncly4)
    call gtipar(com,'nclz',nclz4,nclz4)
    tranx4 = 0.5*(nclx4-1)*dcel4
    trany4 = 0.5*(ncly4-1)*dcel4
    tranz4 = 0.5*(nclz4-1)*dcel4
    idcel4=1.0/dcel4
    allocate (chden(nclx4*ncly4*nclz4))
    chden=0.0
    Qchden=.true.
    write(outu,'(6x,a)') 'Charge Density (CHDEN) Activated:'
    write(outu,'(6x,a)') '------------------'
    write(outu,'(6x,a,i6)')  'Number of grid point in X   (nclx) = ',nclx4
    write(outu,'(6x,a,i6)')  'Number of grid point in Y   (ncly) = ',ncly4
    write(outu,'(6x,a,i6)')  'Number of grid point in Z   (nclz) = ',nclz4
    write(outu,'(6x,a,f8.3)') 'Grid spacing                (dcel) = ',dcel4
    write(outu,'(6x,a,f8.3)') 'Center of box in X          (xbcen)= ',xbcen4
    write(outu,'(6x,a,f8.3)') 'Center of box in Y          (ybcen)= ',ybcen4
    write(outu,'(6x,a,f8.3)') 'Center of box in Z          (zbcen)= ',zbcen4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in X from ',xbcen4-tranx4,' to ',xbcen4+tranx4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in Y from ',ybcen4-trany4,' to ',ybcen4+trany4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in Z from ',zbcen4-tranz4,' to ',zbcen4+tranz4
  ! **********************************************************************
  elseif (wrd5.eq.'profi') then ! user defined diffusion constant profile
  !        -----------------     
     if (.not.Qpar) call error ('shell_simul', 'PROFILE order is defined before PARTICLE order', faterr)
     ! Unit for file of diffusion constant profile      
     call gtipar(com,'difunit',iunit,0)
     if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
     if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PROFILE order', faterr)
     iunit = unvec(iunit)
     ! Reading file of diffusion constant profile
     read(iunit,*) nspline
     if (nspline.gt.dspline) call error ('shell_simul', 'nspline is greater than dspline', faterr)
     if (nspline.gt.2) then 
       Qprofile = .true.
       do is = 1, nspline
         read(iunit,*) xs(is),ys(is)
       enddo
       call spline (nspline,xs,ys,b,c,d)
       write(outu,*)'Diffusion coefficient profile activated'
       write(outu,*)'  nspline = ',nspline
       write(outu,*)'  z: defined from ',xs(1),' to ',xs(nspline)
     endif
  ! ****************************************************************
  elseif (wrd5.eq.'diffu') then ! diffusion_constant
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'DIFFUSION order is defined before PARTICLE order', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'DIFFUSION order is defined before SYSTEM order', faterr)
     ! pore length [real*8,default=0]
     call gtdpar(com,'porelength',plength2,0.0)
     if (plength2.lt.0.0) call error ('shell_simul', 'Pore lenght is lower than zero in DIFFUSION order', faterr)
     plength2=plength2*0.5  !store half the pore length
     ! membrane center [real*8,default=0]
     call gtdpar(com,'pcenter',pcenter,0.0)
     if (pcenter.gt.lz2p.or. pcenter.lt.lz2m) call error ('shell_simul', 'Position of center of a pore along the Z-axis is not correct', faterr)
  
     endlog = .false.
     do while (.not. endlog) 
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type 
         call getfirst(com,wrd4)
         itype=getetyp(wrd4)
         call gtdpar(com,'lowfraction',ampl3(itype),0.0)
         call gtdpar(com,'switchlength',p3(itype),1.0)
         ! Writting file of diffusion constant 
         if (check(com,'file')) then
           call gtipar(com,'unit',iunit,1)
           if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
           if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
           if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in DIFFUSION order', faterr)
           iunit = unvec(iunit)
           do z1i = -90, 90, 1
             z1 = real(z1i)
             call switch3(r1,r2,z1,plength2,p3(itype),pcenter)
             write(iunit,'(2f13.5)') z1,etypl(itype)%dif*(ampl3(itype)+(1.0-ampl3(itype))*r1)
           enddo
         endif
       endif  
     enddo
     Qdiffuse = .true.
  
     write(outu,*)
     write(outu,'(6x,a)') 'DIFFUSION parameters:'
     write(outu,'(6x,a)') '---------------------'
     write(outu,'(6x,a,f8.3,a)') 'pore length        ',2*plength2,' [Angs]'
     write(outu,'(6x,a,f8.3,a)') 'membrane center    ',   pcenter,' [Angs]'
  
     do i = netnuc+1, netyp
        write(outu,'(6x,i3,3x,a,1x,5f8.3)') i,etypl(i)%nam,ampl3(i),p3(i)
     enddo
     write(outu,*)
  ! *******************************************************************
  elseif (wrd5.eq.'print') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'PRINT order is defined before PARTICLE and/or NUCLEOTIDE orders', warning)
     if (check(com,'system')) then
       if (Qnucl) write (outu,'(6x,a,i4)') 'Total number of sites ',nelenuc
       if (Qpar) then       
         nions = nele - nelenuc
         write(outu,'(6x,a,i4)') 'Total number of ions  ',nions
         do ib = 1, nbuffer
           write(outu,'(6x,a,2i4)') 'buffer---number of ions ',ib,nat(ib)
         enddo
       endif
     elseif (check(com,'coor')) then
       call gtipar(com,'unit',iunit,0)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)   
       write(outu,*)    
       write(outu,'(6x,a)') 'Configuration has been written'
       if (check(com,'xyz')) then
         call printxyz(iunit)
       else if (check(com,'pdbe')) then ! PDB with Extended Digits
         if (Qsphere) then
           write(iunit,'(A6)') 'CRYST1'
         else
           write(iunit,'(A6,3f9.3,3f7.2)') 'CRYST1',LX,LY,LZ,90.0,90.0,90.0
         endif
         call printpdbe(iunit)
       else if (check(com,'crd')) then
         if (Qsphere) then
           write(iunit,'(a,f12.3)') '* System shape: Spherical, Radius: ',Rsphe
         elseif (Qecyl) then
           write(iunit,'(a,2f12.3,a,f12.3)') '* System shape: Elliptical Cylinder, Ellipse diameters x and y:',lx,ly,' Cylinder Length: ',lz
         else
           write(iunit,'(3(a,f12.3))') '* System shape: Box, Dimensions ->  LX: ',lx,'  LY: ',ly,'  LZ: ',lz
         endif
         write(iunit,'(3(a,f12.3))') '* System Center:  X ',cx,'  Y ',cy,'  Z ',cz
         write(iunit,'(a,f20.3)') '* Total System Volume (Ang**3): ',tvol
         call printcrd(iunit)
       else if (check(com,'crde')) then
         if (Qsphere) then
           write(iunit,'(a,f12.3)') '* System shape: Spherical, Radius: ',Rsphe
         elseif (Qecyl) then
           write(iunit,'(a,2f12.3,a,f12.3)') '* System shape: Elliptical Cylinder, Ellipse diameters x and y:',lx,ly,' Cylinder Length: ',lz
         else
           write(iunit,'(3(a,f12.3))') '* System shape: Box, Dimensions ->  LX: ',lx,'  LY: ',ly,'  LZ: ',lz
         endif
         write(iunit,'(3(a,f12.3))') '* System Center:  X ',cx,'  Y ',cy,'  Z ',cz
         write(iunit,'(a,f20.3)') '* Total System Volume (Ang**3): ',tvol
         call printcrde(iunit)
       else !if (check(com,'pdb')) then
         if (Qsphere) then
           write(iunit,'(A6)') 'CRYST1'
         else
           write(iunit,'(A6,3f9.3,3f7.2)') 'CRYST1',LX,LY,LZ,90.0,90.0,90.0
         endif
         call printpdb(iunit)
       endif  
     elseif (check(com,'title')) then
       write(outu,'(6x,a)') 'TITLE: '//trim(title)
     elseif (check(com,'iseed')) then
       write(outu,'(a,i10)') ' iseed: ',iseed
     elseif (check(com,'dnacenter')) then
       if (Qnucl) then
         xm=0.0
         ym=0.0
         zm=0.0
         do i=1,nelenuc
           xm=xm+r(i)%x
           ym=ym+r(i)%y
           zm=zm+r(i)%z
         enddo
         xm=xm/nelenuc
         ym=ym/nelenuc
         zm=zm/nelenuc
         write(outu,'(6x,a,3(f15.8,1x))') 'DNA GEOMETRIC CENTER: ',xm,ym,zm
       else 
         call error ('shell_simul', 'Cannot PRINT dnacenter, DNA is not defined', warning)
       endif
     elseif (check(com,'static')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0) 
       call gtdpar(com,'y1',y1,0.0) 
       call gtdpar(com,'z1',z1,0.0) 
       call gtdpar(com,'x2',x2,0.0) 
       call gtdpar(com,'y2',y2,0.0) 
       call gtdpar(com,'z2',z2,0.0) 
       call staticplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead)
     elseif (check(com,'repul')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0)
       call gtdpar(com,'y1',y1,0.0)
       call gtdpar(com,'z1',z1,0.0)
       call gtdpar(com,'x2',x2,0.0)
       call gtdpar(com,'y2',y2,0.0)
       call gtdpar(com,'z2',z2,0.0)
       call repulplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead)
     elseif (check(com,'rfpar')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0)
       call gtdpar(com,'y1',y1,0.0)
       call gtdpar(com,'z1',z1,0.0)
       call gtdpar(com,'x2',x2,0.0)
       call gtdpar(com,'y2',y2,0.0)
       call gtdpar(com,'z2',z2,0.0)
       if (Qrfpsin) then 
         call rfparplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead,1)
       else
         call rfparplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead,netyp-netnuc)
       endif
     elseif (check(com,'statxd')) then
!  Prints static field values in x,y,z,pot format in an ASCII file. Allows to generate 3D,2D,1D Plots.
!      [PRINT] statxd unit (integer) {noheader} ix1 (integer) ix2 (integer) iy1 (integer) iy2 (integer) iz1 (integer) iz2 (integer)
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtipar(com,'ix1',ix1,0)
       call gtipar(com,'iy1',iy1,0)
       call gtipar(com,'iz1',iz1,0)
       call gtipar(com,'ix2',ix2,nclx1-1)
       call gtipar(com,'iy2',iy2,ncly1-1)
       call gtipar(com,'iz2',iz2,nclz1-1)
       call statxd(ix1,iy1,iz1,ix2,iy2,iz2,iunit,Qnohead)
     elseif (check(com,'repxd')) then
!  Prints repulsion field values in x,y,z,rep format in an ASCII file. Allows to generate 3D,2D,1D Plots.
!      [PRINT] repxd unit (integer) {noheader} ix1 (integer) ix2 (integer) iy1 (integer) iy2 (integer) iz1 (integer) iz2 (integer)
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtipar(com,'ix1',ix1,0)
       call gtipar(com,'iy1',iy1,0)
       call gtipar(com,'iz1',iz1,0)
       call gtipar(com,'ix2',ix2,nclx2-1)
       call gtipar(com,'iy2',iy2,ncly2-1)
       call gtipar(com,'iz2',iz2,nclz2-1)
       call repxd(ix1,iy1,iz1,ix2,iy2,iz2,iunit,Qnohead)
     elseif (check(com,'chden')) then
       if (Qchden) then
         call gtipar(com,'unit',iunit,1)
         if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
         if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
         if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
         iunit = unvec(iunit)
         zero=0d0
         idv=sng(1.0/dcel4**3)
         ncl3=nclx4*ncly4*nclz4
         write(iunit) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunit) (zero,i=1,6)
         write(iunit) ((idv*chden(i)),i=1,ncl3)
         write(outu,'(6x,A)') 'Charge density map written'
       else
         call error ('shell_simul', 'nothing to print, CHDEN was not called', warning) 
       endif
     elseif (check(com,'efield')) then
       if (allocated(iunitv)) deallocate (iunitv)
       allocate (iunitv(4))
       Qonlychden=check(com,'onlychden')
       call gtipar(com,'xunit',iunitv(1),0)
       if (iunitv(1).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(1).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(1)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(1) = unvec(iunitv(1))
       call gtipar(com,'yunit',iunitv(2),0)
       if (iunitv(2).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(2).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(2)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(2) = unvec(iunitv(2))
       call gtipar(com,'zunit',iunitv(3),0)
       if (iunitv(3).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(3).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(3)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(3) = unvec(iunitv(3))
       call gtipar(com,'munit',iunitv(4),0)
       if (iunitv(4).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(4).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(4)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(4) = unvec(iunitv(4))
       if (Qchden) then
         ncl3=nclx4*ncly4*nclz4
         allocate(efield(4,ncl3))
         efield=0.0
         zero=0d0
         if (Qchdencnt) then 
           do i=0,nclx4-1
             do j=0,ncly4-1
               do k=0,nclz4-1
                 x1=i*dcel4-trany4+xbcen4 
                 y1=j*dcel4-trany4+ybcen4 
                 z1=k*dcel4-tranz4+zbcen4
                 in1=i*ncly4*nclz4+j*nclz4+k+1
                 do ii=0,nclx4-1
                   do ij=0,ncly4-1
                     do ik=0,nclz4-1
                       x2=ii*dcel4-trany4+xbcen4
                       y2=ij*dcel4-trany4+ybcen4 
                       z2=ik*dcel4-tranz4+zbcen4
                       in2=ii*ncly4*nclz4+ij*nclz4+ik+1
                       if (in1.ne.in2.and.chden(in2).ne.0.0) then
                         r2=1.0/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                         r1=sqrt(r2)
                         efield(1,in1)=efield(1,in1)+sng(celec2*chden(in2)*(x2-x1)*r1*r2)
                         efield(2,in1)=efield(2,in1)+sng(celec2*chden(in2)*(y2-y1)*r1*r2)
                         efield(3,in1)=efield(3,in1)+sng(celec2*chden(in2)*(z2-z1)*r1*r2)
                         efield(4,in1)=efield(4,in1)+sng(celec2*chden(in2)*r2)
                       endif
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
         if (.not.Qonlychden) then
           if (Qphix) then
             do i=0,nclx4-1
               do j=0,ncly4-1
                 do k=0,nclz4-1
                   x1=i*dcel4-trany4+xbcen4
                   y1=j*dcel4-trany4+ybcen4
                   z1=k*dcel4-tranz4+zbcen4
                   in1=i*ncly4*nclz4+j*nclz4+k+1
                   call staefield(x1,y1,z1,x2,y2,z2)
                   efield(1,in1)=efield(1,in1)+sng(x2)
                   efield(2,in1)=efield(2,in1)+sng(y2)
                   efield(3,in1)=efield(3,in1)+sng(z2)
                   efield(4,in1)=efield(4,in1)+sng(sqrt(x2**2+y2**2+z2**2))
                 enddo
               enddo
             enddo
           elseif (voltage.ne.0.0.and.Qmemb) then
             x2=afact/Coulomb*kcalmol
             do k=0,nclz4-1
               z1=k*dcel4-tranz4+zbcen4
               if (z1.lt.zmemb1) then ! REGION 1: z=z(iat)-zmemb1 < 0, lim{z->-inf} pot(1) = 0
                 z2 = -x2*ikappa*exp(ikappa*(z1-zmemb1))
               elseif ((z1.ge.zmemb1) .and. (z1.le.zmemb2)) then ! REGION 2
                 z2 = -x2*ceps*ikappa
               elseif (z1.gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
                 z2 = -x2*ikappa*exp(-ikappa*(z1-zmemb2))
               endif
               do i=0,nclx4-1
                 do j=0,ncly4-1
                   x1=i*dcel4-trany4+xbcen4
                   y1=j*dcel4-trany4+ybcen4
                   in1=i*ncly4*nclz4+j*nclz4+k+1
                   efield(3,in1)=efield(3,in1)+sng(z2)
                   efield(4,in1)=efield(4,in1)+sng(z2)
                 enddo
               enddo
             enddo
           endif
         endif
         write(iunitv(1)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(1)) (zero,i=1,6)
         write(iunitv(1)) ((efield(1,i)),i=1,ncl3)
         write(iunitv(2)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(2)) (zero,i=1,6)
         write(iunitv(2)) ((efield(2,i)),i=1,ncl3)
         write(iunitv(3)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(3)) (zero,i=1,6)
         write(iunitv(3)) ((efield(3,i)),i=1,ncl3)
         write(iunitv(4)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(4)) (zero,i=1,6)
         write(iunitv(4)) ((efield(4,i)),i=1,ncl3)
         write(outu,'(6x,A)') 'Electric Field maps written'
         deallocate(efield)
       else
         call error('shell_simul', 'nothing to print, CHDEN was not called', warning)
       endif
       deallocate (iunitv)
     elseif (check(com,'statsup')) then
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'x',x1,0.0)
       call gtdpar(com,'y',y1,0.0)
       call areapot(x1,y1,iunit)
     endif
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'count') then
  !        ---------------
     if (.not.Qbuf) call error ('shell_simul', 'COUNT order is defined before BUFFER order', faterr)
     call count
     ncnt=0
     do i = nparnuc+1, npar
       itype = parl(i)%ptyp
       ncnt(itype) = ncnt(itype) + 1
     enddo
     write(outu,*)
     write(outu,'(6x,a,i4)') 'Total number of ions ',npar-nparnuc
     do ib = 1, nbuffer
       write(outu,'(6x,a,2i4)') 'buffer--number of ions ',ib,nat(ib)
     enddo
     write(outu,'(6x,a)') 'ION----TYPE--BUFFER'
     do i = nparnuc+1, npar
       if (parl(i)%ibuf.ne.0) write(outu,'(6x,i5,1x,a4,1x,i5)') i,ptypl(parl(i)%ptyp)%nam,parl(i)%ibuf
     enddo
  ! **********************************************************************
  elseif (wrd5.eq.'coor') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'COOR order is defined before PARTICLE and/or NUCLEOTIDE orders', faterr)
     if (check(com,'read')) then
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in COOR order', faterr)
       iunit = unvec(iunit)
       !call delparofkind(3)  ! Delete all particles added when ptype were created
       write(outu,'(6x,a,i3)')'Reading coordinates from unit ',iunit
       Qpdb=check(com,'pdb')
       Qpdbe=check(com,'pdbe')
       Qcrd=check(com,'crd')
       Qcrde=check(com,'crde')
       if (Qpdb) then
         write(outu,'(6x,a)')'Format PDB'
       else if (Qpdbe) then
         write(outu,'(6x,a)')'Format PDBE (PDB with extended coordinate digits)'
       else if (Qcrd) then
         write(outu,'(6x,a)')'Format CRD (CHARMM xplor format)'
       else if (Qcrde) then
         write(outu,'(6x,a)')'Format CRDE (CHARM xplor Extended format)'
       else
         call error ('shell_simul', 'File Format not specified in COOR', faterr)
       endif
       do
         read(iunit,'(A)') com
         if (Qpdbe.or.Qpdb) then
           if (com(1:4).eq.'ATOM') exit
         else
           if (com(1:1).ne.'*') then
             read(com, *) nelem
             exit
           endif
         endif
       enddo
       if (Qnucl) then
         ilast=0
         if ((Qcrd.or.Qcrde).and.nelem.lt.nelenuc) call error ('shell_simul', 'Less number of elements than expected in COOR', faterr)
         do i = 1, nelenuc
           if (Qcrd.or.Qcrde) read(iunit,'(a)') com
           if (Qpdbe) then 
             read(com,'(6x,5x,x,5x,5x,I4,4x,3F16.8)') itype,rr%x,rr%y,rr%z
           elseif (Qpdb) then
             read(com,'(6x,5x,x,5x,5x,I4,4x,3F8.3)') itype,rr%x,rr%y,rr%z
           elseif (Qcrd) then
             read(com,'(5x,I5,1x,4x,1x,4x,3F10.5)') itype,rr%x,rr%y,rr%z
           elseif (Qcrde) then
             read(com,'(10x,I10,2x,4x,6x,4x,4x,3F20.10)') itype,rr%x,rr%y,rr%z
           endif
           if (itype.ne.ilast) jtype=1
           call putcoorinpar(itype,jtype,rr%x,rr%y,rr%z)
           ilast=itype
           jtype=jtype+1
           if (Qpdb.or.Qpdbe) read(iunit,'(a)') com
         enddo
       endif
       ilast=0
       if (Qpdbe.or.Qpdb) then
         do while (trim(adjustl(com)).ne.'END')
           if (Qpdbe) then
             read (com,'(6x,5x,x,A5,5x,I4,4x,3F16.8,F6.2)') wrd4,itype,rr%x,rr%y,rr%z,pkind
           else
             read (com,'(6x,5x,x,A5,5x,I4,4x,3F8.3,F6.2)') wrd4,itype,rr%x,rr%y,rr%z,pkind
           endif
           if (ilast.ne.itype) then
             jtype=1
             call addpar(getptyp(wrd4),kind=int(pkind))
           endif
           call putcoorinpar(npar,jtype,rr%x,rr%y,rr%z)
           ilast=itype
           jtype=jtype+1
           read(iunit,'(a)') com
         enddo
       else
         do i=nelenuc+1,nelem
           read(iunit,'(a)') com
           if (Qcrde) then
             read (com,'(10x,I10,2x,A4,6x,4x,4x,3F20.10,2X,I4)') itype,wrd4,rr%x,rr%y,rr%z,ikind
           else
             read (com,'(5x,I5,1x,A4,1x,4x,3F10.5,1X,I4,1X,I4)') itype,wrd4,rr%x,rr%y,rr%z,ikind
           endif
           if (ilast.ne.itype) then
             jtype=1
             call addpar(getptyp(wrd4),kind=ikind)
           endif
           call putcoorinpar(npar,jtype,rr%x,rr%y,rr%z)
           ilast=itype
           jtype=jtype+1
         enddo
       endif
       write(outu,'(6x,a)') 'coordinates have been read'
       call count
       ncnt = 0
       do i = nparnuc+1, npar
         itype = parl(i)%ptyp
         ncnt(itype) = ncnt(itype) + 1
       enddo
     elseif (check(com,'gener')) then
       doions=Qbuf
       if (.not.doions) call error ('shell_simul', 'ions missing after COOR gener', faterr)
       if (doions) then
         !call delparofkind(3)
         do ib = 1, nbuffer
           nat(ib) = nint(avnum(ib))
         enddo
         do ib = 1, nbuffer
           itype = ibfftyp(ib)
           do i = 1, nat(ib)
             do
               call insert(ib,rr%x,rr%y,rr%z) ! find new center of particle
               call addpar(itype,kind=3,ibuf=ib) ! add particle to list
               call movepar(npar,rr) ! locate and rotate
               call par_interact(npar, dener)
               rate = (avnum(ib)/float(nat(ib)+1))*exp(-(dener-mu(ib))*ikBT)
               rate = rate/(1.0+rate)
               if (rndm().le.rate) exit
               call delpar(npar,3)
             enddo
           enddo
         enddo
       endif ! Qbuf
       call count
       ncnt = 0
       do i = nparnuc+1, npar
         itype = parl(i)%ptyp
         ncnt(itype) = ncnt(itype) + 1
       enddo
       write(outu,'(6x,a)') 'coordinates have been generated'
     elseif (check(com,'rot').and.Qnucl) then
       ! rotation for DNA sites coordinates
       Qrot = .true.
       if (check(com,'ref')) then
         call gtipar(com,'s1',s1,1)
         call gtipar(com,'s2',s2,nelenuc)
         call gtipar(com,'s3',s3,2)
         write(outu,'(6x,3(a,x,i0,x))') 'Rotating using references: s1=',s1,'s2=',s2,'s3=',s3
         vc3(1)=r(s1)%x-r(s2)%x
         vc3(2)=r(s1)%y-r(s2)%y
         vc3(3)=r(s1)%z-r(s2)%z
         x3=sqrt(dot_product(vc3,vc3))
         vc3=vc3/x3
         vc1(1)=r(s3)%x-r(s2)%x
         vc1(2)=r(s3)%y-r(s2)%y
         vc1(3)=r(s3)%z-r(s2)%z
         x1=sqrt(dot_product(vc1,vc1))
         vc1=vc1/x1
         call cross_product(vc3,vc1,vc2)
         x2=sqrt(dot_product(vc2,vc2))
         vc2=vc2/x2
         call cross_product(vc2,vc3,vc1)
         x1=sqrt(dot_product(vc1,vc1))
         vc1=vc1/x1
         rot(1,:)=vc1 
         rot(2,:)=vc2
         rot(3,:)=vc3 
       else
       ! rotation matrix elements
         call gtdpar(com,'r11',rot(1,1),1.0)
         call gtdpar(com,'r12',rot(1,2),0.0)
         call gtdpar(com,'r13',rot(1,3),0.0)
         call gtdpar(com,'r21',rot(2,1),0.0)
         call gtdpar(com,'r22',rot(2,2),1.0)
         call gtdpar(com,'r23',rot(2,3),0.0)
         call gtdpar(com,'r31',rot(3,1),0.0)
         call gtdpar(com,'r32',rot(3,2),0.0)
         call gtdpar(com,'r33',rot(3,3),1.0)
       endif
       do i = 1, 3
         sum1 = 0.0
         do j = 1, 3
           sum1 = sum1 + rot(i,j)*rot(i,j)
         enddo
         if (abs(1.0-sum1).ge.1.0e-8) call error ('shel_simul', 'rotation matrix for DNA sites is not orthogonal in COOR order', faterr)
       enddo
       sum1 = 0.0e0
       sum2 = 0.0e0
       sum3 = 0.0e0
       do j = 1, 3
         sum1 = sum1 + rot(1,j)*rot(2,j)
         sum2 = sum2 + rot(2,j)*rot(3,j)
         sum3 = sum3 + rot(1,j)*rot(3,j)
       enddo
       if (abs(sum1).ge.1.0e-8 .or. abs(sum2).ge.1.0e-8 .or.abs(sum3).ge.1.0e-8) call error ('shel_simul', 'rotation matrix for DNA sites is not orthogonal in COOR order', faterr)
       do i = 1, nelenuc
         xold = r(i)%x
         yold = r(i)%y
         zold = r(i)%z
         r(i)%x = xold*rot(1,1) + yold*rot(1,2) + zold*rot(1,3)
         r(i)%y = xold*rot(2,1) + yold*rot(2,2) + zold*rot(2,3)
         r(i)%z = xold*rot(3,1) + yold*rot(3,2) + zold*rot(3,3)
       enddo
       write(outu,'(6x,a)') 'Coordinates for DNA sites have been rotated. Rotation matrix:'
       do i = 1, 3
         write (outu,'(2x,3(1x,f8.3))') (rot(i,j),j=1,3)
       enddo
       ! translocation for DNA sites coordinates
     elseif (check(com,'tras').and.Qnucl) then
       Qtras=.true.
       ! traslation in x direction of DNA B isoform coordinates
       call gtdpar(com,'xtras',xtras,0.0)
       ! traslation in y direction of DNA B isoform coordinate
       call gtdpar(com,'ytras',ytras,0.0)
       ! traslation in z direction of DNA B isoform coordinate 
       call gtdpar(com,'ztras',ztras,0.0)
       do i = 1, nelenuc
         r(i)%x = r(i)%x + xtras
         r(i)%y = r(i)%y + ytras
         r(i)%z = r(i)%z + ztras
       enddo
       write(outu,'(6x,a,1x,a,3(f12.6,a))')'Coordinates for DNA sites have been moved by','(',xtras,',',ytras,',',ztras,')'
     elseif (check(com,'setori').and.Qnucl) then
       xm=0.0
       ym=0.0
       zm=0.0
       do i=1,nelenuc
         xm=xm+r(i)%x
         ym=ym+r(i)%y
         zm=zm+r(i)%z
       enddo
       xm=xm/nelenuc
       ym=ym/nelenuc
       zm=zm/nelenuc
       do i = 1, nelenuc
         r(i)%x = r(i)%x - xm
         r(i)%y = r(i)%y - ym
         r(i)%z = r(i)%z - zm
       enddo
       call gtdpar(com,'x',xm,xm)
       call gtdpar(com,'y',ym,ym)
       call gtdpar(com,'z',zm,zm)
       do i = 1, nelenuc
         r(i)%x = r(i)%x + xm
         r(i)%y = r(i)%y + ym
         r(i)%z = r(i)%z + zm
       enddo
       write(outu,'(6x,a,3(f12.6,a))')'DNA geometric center is now at  (',xm,',',ym,',',zm,')'
     endif
     ! Print out data
     call printcrde(outu)
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'gsbp') then !  generalized solvent boundary potential
  !        ---------------
  !REPULSIVE POTENTIAL
    Qadj = check(com,'adjust')
    Qphiv = check(com,'phiv')
    if (Qphiv) then
      ! Threshold for 27 and 8 cells
      call gtdpar(com,'thold27',thold27,27.0)
      call gtdpar(com,'thold8',thold8,8.0)
      ! Magnitude of grid-based repulsive potential 'phiv'
      ! [real*8,default=0] 
      call gtdpar(com,'svdw',svdw,50.0)
      ! Trilinear interpolation is used for 'phiv'
      ! [default=3rd-order B-spline interpolation]
!      Qtrln = check(com,'trilinear')
      Qtrln = .not.check(com,'bspline')
      if (Qtrln) then
         write(outu,'(6x,a)') 'Trilinear function will be used for repulsive potential'
      else
         write(outu,'(6x,a)') 'B-spline function will be used for repulsive potential'
      endif
      Qnmcden = check(com,'nmcden')
      if (Qnmcden) then
        write(outu,*)
        write(outu,'(6x,a)') 'Different ion-accessible space is used for different ions and sites'
        totnumb=netyp
        if (.not.Qatexp.and..not.Qnucl .and. .not. Qpar) then
          call error ('shell_simul', 'GSBP order is defined before PTYPE, PARTICLE and/or NUCLEOTIDES orders', faterr)
        endif
        if (totnumb.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
        call gtcpar(com,'munit',nn,word)
        if (nn.eq.0) call error ('shell_simul','munit cannot be empty or ommited',faterr)
        if (nn.ne.totnumb) call error ('shell_simul','number of units in munit does not match with number of types',faterr)
        call gtcipar(word,vecphiv)
        do i=1,totnumb
          iunit=vecphiv(i)
          if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
          if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
          if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
          vecphiv(i) = unvec(iunit)
        enddo
        iunit=vecphiv(1)
        write(outu,'(6x,a,10i3)') 'Reading grid-based repulsive potential from unit ',(vecphiv(i),i=1,totnumb)
      else
        call gtipar(com,'phivunit',iunit,1)
        if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
        if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
        if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
        iunit = unvec(iunit)
        write(outu,*)
        write(outu,'(6x,a,i3)')'Reading grid-based repulsive potential from unit ',iunit
      endif
      call readphi(iunit,outu,'PHIV',Qadj)
      ! Maximum position of a grid for 'phiv' along the Z-axis
      ! [real*8,default=max map limit]
      call gtdpar(com,'vzmax',vzmax,zbcen2+tranz2)
      if (vzmax.gt.zbcen2+tranz2) call error('shell_simul', 'vzmax cannot be outside boundaries of repulsion map', faterr)
      ! Minimum position of a grid for 'phiv' along the Z-axis
      ! [real*8,default=-min map limit]         
      call gtdpar(com,'vzmin',vzmin,zbcen2-tranz2)
      if (vzmin.lt.zbcen2-tranz2) call error('shell_simul', 'vzmin cannot be outside boundaries of repulsion map', faterr)
      write(outu,'(6x,A,F8.2,A,F8.2,A,F8.2)') 'Uniformed repusive potential will be scaled by ',svdw,' kcal/mol between ',vzmin,' and ',vzmax 
      call gtipar(com,'repwalls',wallsi4,126)
      if (wallsi4.gt.0) then
        if (wallsi4.gt.126) then
          walls=126
        else
          walls=itoi1(wallsi4)
        endif
        call repwalls(walls)
        write(outu,'(6x,a,i0)') 'REPWALLS activated. Number of walls defined: ',walls
      else
        walls=0
        write(outu,'(6x,a,i0)') 'REPWALLS deactivated. Strongly recommended to be used.'
      endif
      call accessiblevolume(totvol)
      write(outu,'(/6x,a,f20.8)') 'Accessible volume for ions: ', totvol
    endif ! Qphiv

    !STATIC EXTERNAL FIELD
    Qphix = check(com,'phix')
    if (Qphix) then
      call gtipar(com,'phixunit',iunit,1)
      if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
      if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
      if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
      iunit = unvec(iunit)         
      write(outu,*) 
      write(outu,'(6x,a,i3)') 'Reading static external field from unit ',iunit
      call READPHI(iunit,outu,'PHIX',Qadj)
    endif ! Qphix
    Qrfpar = check(com,'rfpar')
    if (Qrfpar) then
      Qrfpsin=check(com,'rfpsingle')
      call gtcpar(com,'rfparunit',nn,word)
      if (nn.eq.0) call error ('shell_simul','rfparunit cannot be empty or ommited',faterr)
      allocate (iunitv(nn))
      call gtcipar(word,iunitv)
      do i=1,nn
        iunit=iunitv(i)
        if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
        if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
        if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
        iunitv(i) = unvec(iunit) 
      enddo
      call gtdpar(com,'rfparfac',reffac,1.5)
      sqrfac=sqrt(reffac)
      write(outu,*) 
      write(outu,'(6x,a$)') 'Reading Reaction field parameters from units:'
      write(outu,*) iunitv
      write(outu,'(6x,a$)') 'Using RFPAR factor for effective radius parameters of:'
      write(outu,*) reffac
      call readrfpar(iunitv,nn,outu,Qadj)
      deallocate (iunitv)
    endif
    write(outu,*)
  ! **************************************************************************
  elseif (wrd5.eq.'svdw') then
    if (.not.Qatexp.and..not.Qpar.and..not.Qnucl) call error ('shell_simul', 'SVDW order is defined before PTYPE, PARTICLE and/or NUCLEOTIDE order', faterr)
    if (.not.Qphiv) call error ('shell_simul', 'Repulsion field has not been read', faterr)
    if (Qnmcden) call error ('shell_simul', 'SVDW not compatible with NMCDEN', faterr)
    allocate (scal(netyp))
    scal=1.0
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type 
        call getfirst(com,wrd4)
        itype=getetyp(wrd4)
        call gtdpar(com,'scale',scal(itype),1.0)
      endif
    enddo
    Qsvdw = .true.
    write(outu,'(/6x,a/6x,a)') 'SVDW Type Scaling Factor enabled','Type   Factor'
    do i=1,netyp
      write(outu,'(6x,a,3x,f10.5)') etypl(i)%nam,scal(i)
    enddo
  elseif (wrd5.eq.'exit') then
  !        ---------------
     write(outu,*)
     finish = timer()-start
     write(outu,'(/6x,a)') 'CPU Time:'
     write(outu,'(8x,a)') 'Format 1:'
     write(com,*) finish,' msec / ',finish/1e3,' sec / ',finish/6e4,' m / ',finish/36e5,' h / ',finish/864e5,' d'
     write(outu,'(8x,a)') trim(com)
     write(outu,'(8x,a)') 'Format 2:'
     values(3)=int(finish/864e5)  ! days
     values(1)=24*values(3) ! hours eq
     values(5)=int(finish/36e5)-values(1) ! hours
     values(1)=60*(values(1)+values(5)) ! min eq
     values(6)=int(finish/6e4)-values(1) ! minutes
     values(1)=60*(values(1)+values(6)) ! sec eq
     values(7)=int(finish/1e3)-values(1) ! seconds
     values(1)=1000*(values(1)+values(7)) ! milisec eq
     values(8)=finish-values(1) ! miliseconds
     write(com,*) values(3),' d, ',values(5),' h, ',values(6),' m, ',values(7),' s, ',values(8),' ms'
     write(outu,'(8x,a)') trim(com)
     call date_and_time(date, time, zone, values)
     write(outu,'(/6x,a,7(i0,a))') 'Finished at (YYYY-MM-DD HH:mm:ss.ms): ',values(1),'-',values(2),'-',values(3),' ',values(5),':',values(6),':',values(7),'.',values(8)
     logfinal = .true.
  else
    write(outu,'(6x,a)') '*ERROR*  Unrecognized command:'
  endif
enddo
contains
  subroutine allocateqsome()
  implicit none
  ! Allocate Qefpot, Qcol, Qlj, Qsrpmfi, warn
  if (allocated(Qefpot)) deallocate(Qefpot)
  allocate (Qefpot(netp))
  Qefpot=.false.
  if (allocated(Qcol)) deallocate(Qcol)
  allocate (Qcol(netp))
  Qcol=.true.
  if (allocated(Qlj)) deallocate(Qlj)
  allocate (Qlj(netp))
  Qlj=.false.
  if (allocated(Qsrpmfi)) deallocate(Qsrpmfi)
  allocate (Qsrpmfi(netp))
  Qsrpmfi=.false.
  if (allocated(warn)) deallocate (warn)
  allocate (warn(netyp))
  warn=0
  ! ASSIGN CHARGES USING CDIE
  cecd=celec/cdie
  end subroutine

  subroutine setuplj()
  implicit none
  if (allocated(epp4)) deallocate(epp4)
  if (allocated(sgp2)) deallocate(sgp2)
  allocate (epp4(netp),sgp2(netp))
  epp4=0.0
  sgp2=0.0
  call updateuetl()
  do i=1,netyp
    do j=i,netyp
      if (etypl(i)%eps.gt.0.0.and.etypl(i)%sig.gt.0.0.and.etypl(j)%eps.gt.0.0.and.etypl(j)%sig.gt.0.0) then
        is=etpidx(i,j)
        epp4(is)=4.0*sqrt(etypl(i)%eps*etypl(j)%eps)
        sgp2(is)=(0.5*(etypl(i)%sig+etypl(j)%sig))**2
        Qlj(is)=.true.
      else
        if (etul(i).and.etul(j)) write(outu,'(6x,a,x,a,x,a)')'Warning: Missing single LJ parameters to compute pairs for:',etypl(i)%nam,etypl(j)%nam
      endif
    enddo
  enddo
  end subroutine
end
