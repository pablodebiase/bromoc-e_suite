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

! Fortran common block for Grand Canonical Monte Carlo
! Dimension of arrays
module listmod
implicit none

! TYPES DEFINITIONS
! Cartesian Type
type :: car
    real :: x
    real :: y
    real :: z
end type car

! Particle Properties Type
type :: parpro
    integer :: ptyp                            !! Identifies the Particle Type with an Index number
    integer :: sr                              !! Starting in the position vector
    integer :: ne                              !! Number of element in the particle
    integer :: kind                            !! Kind of Particle (Nucleotide, Particles, Charmm)
    integer :: ibuf                            !! Buffer Id (only for Particles)
end type parpro

! Particle Type
type :: partype
    integer :: ne                              !! Number of Elements in Particle
    character*4                        :: nam    !! Particle Name
    real                               :: chg    !! Particle Total Charge
    integer,allocatable,dimension(:)   :: etyp   !! Particle Element Types Vector :: Size of ne
    type(car),allocatable,dimension(:) :: r      !! Position Vector of Particle Elements :: Size of ne
    real,allocatable,dimension(:)      :: chgi   !! Charge Vector of Particle Elements :: Size of ne
end type partype

! Element Type
type :: eletype
    character*4 :: nam    ! Element Name
    real        :: chg    ! Element Charge
    real        :: dif    ! Element Diffusivity
    real        :: eps    ! Element Epsilon Lennard Jones
    real        :: sig    ! Element Sigma Lennard Jones
end type eletype

! LIST DEFINITIONS
!! Particles List
integer                                 :: npar       !! Number of Particles 
type(parpro), allocatable, dimension(:) :: parl       !! Particles List Description :: Size of Number of Particles (npar)
real, allocatable, dimension(:)         :: parz       !! z coordinate for particles

!! Element-Type List
integer                                  :: netyp     !! Number of Element Types
type(eletype),allocatable,dimension(:)   :: etypl     !! Element Name List
logical*1,allocatable,dimension(:)       :: etul      !! Element Types Used List

!! Number of Element-Type Pairs
integer                                  :: netp

!! Particle-Type List
integer                                  :: nptyp      !! Number of Particle Types
type(partype), allocatable, dimension(:) :: ptypl    !! Type of Particle list :: Size of nptyp

!! Force and position vectors
integer nele                                !! Total Number of elements
type(car), allocatable, dimension(:) :: r   !! Elements Position Vector
type(car), allocatable, dimension(:) :: f   !! Elements Force Vector
integer, allocatable, dimension(:)   :: et  !! Elements Type List of nele size
integer, allocatable, dimension(:)   :: pt  !! Particle Type List nele size
integer, allocatable, dimension(:)   :: pe  !! Particle List of nele size
real, allocatable, dimension(:)      :: chq !! Charge List of nele size

!! Used Element Types List
integer nuet                      !! Total Number of Used Element Types in et
integer, allocatable, dimension(:) :: uetl

!! Nucleotides
integer nptnuc,netnuc,nelenuc,nparnuc

contains

subroutine merge(a,na,b,nb,c,nc)
integer, intent(in) :: na,nb,nc         ! normal usage: na+nb = nc
integer, intent(in out) :: a(na)        ! b overlays c(na+1:nc)
integer, intent(in)     :: b(nb)
integer, intent(in out) :: c(nc)
integer :: i,j,k
i = 1; j = 1; k = 1;
do while(i <= na .and. j <= nb)
   if (a(i) <= b(j)) then
      c(k) = a(i)
      i = i+1
   else
      c(k) = b(j)
      j = j+1
   endif
   k = k + 1
enddo
do while (i <= na)
   c(k) = a(i)
   i = i + 1
   k = k + 1
enddo
return
end subroutine merge
 
recursive subroutine mergesort(a,n,t)
integer, intent(in) :: n
integer, dimension(n), intent(in out) :: a
integer, dimension((n+1)/2), intent (out) :: t
integer :: na,nb,v
if (n < 2) return
if (n == 2) then
   if (a(1) > a(2)) then
      v = a(1)
      a(1) = a(2)
      a(2) = v
   endif
   return
endif      
na=(n+1)/2
nb=n-na

call mergesort(a,na,t)
call mergesort(a(na+1),nb,t)

if (a(na) > a(na+1)) then
   t(1:na)=a(1:na)
   call merge(t,na,a(na+1),nb,a,n)
endif
return
end subroutine mergesort

subroutine msort(a,n)
integer n
integer, dimension(n) :: a       ! variable type to sort
integer, dimension ((n+1)/2) :: t
call MergeSort(a,n,t) ! order from lower to higher
end subroutine

subroutine updateuetl()
implicit none
integer i,j
logical f
! Allocate if needed
if (allocated(etul)) deallocate(etul)
allocate (etul(netyp))
if (allocated(uetl)) then
   if (size(uetl).lt.netyp) then
     deallocate (uetl)
     allocate (uetl(netyp))
   endif
else
   allocate (uetl(netyp))
endif
! Make a logical list if the etyp is being used 
do i = 1,netyp
   j=1
   etul(i)=.false.
   do while (j.le.nele)
      if (et(j).eq.i) then
         etul(i)=.true.
         exit
      endif
      j=j+1
   enddo
enddo
! Get Unused etyp
nuet=1
uetl(nuet)=et(1)
do i = 2,nele
   j=1
   f=.false.
   do while (j.le.nuet)
     if (et(i).eq.uetl(j)) then
       f=.true.
       exit
     endif
     j=j+1
   enddo
   if (.not.f) then
     nuet=nuet+1
     uetl(nuet)=et(i)
   endif
enddo
! Sort uetl
call msort(uetl(1:nuet),nuet)
end subroutine

subroutine updatetypels()
implicit none
integer newsize,i,j,ne,sr,ptype
if (allocated(et)) deallocate (et)
if (allocated(pt)) deallocate (pt)
if (allocated(pe)) deallocate (pe)
newsize=size(r)
allocate (et(newsize),pt(newsize),pe(newsize))
do i=1,npar
  sr=parl(i)%sr
  ne=parl(i)%ne
  ptype=parl(i)%ptyp
  do j=1,ne
    et(sr+j)=ptypl(ptype)%etyp(j)
    pt(sr+j)=ptype
    pe(sr+j)=i
  enddo
enddo
end subroutine

subroutine getparz()
implicit none
integer i, pls
if (.not.allocated(parz)) then
  allocate(parz(size(parl)))
else
  pls=size(parl)
  if (size(parz).lt.pls) then
    deallocate(parz)
    allocate(parz(pls))
  endif
endif
do i=1+nparnuc,npar
  parz(i)=iparz(i)
enddo
end subroutine

function iparz(i)
implicit none
integer i, j, ne, sr
real iparz
ne=parl(i)%ne
sr=parl(i)%sr
if (ne.eq.1) then
  iparz=r(sr+1)%z
else
  iparz=0.0
  do j=1,ne
    iparz=iparz+r(sr+j)%z
  enddo
  iparz=iparz/ne
endif
end function

subroutine resizeetypl(newsize)
implicit none
integer vdim,newsize
type(eletype),allocatable,dimension(:)     :: etypltmp      !! Element Name List 
! if 0
if (newsize .lt. 1) then
    if (allocated(etypl)) deallocate (etypl)
endif
! Measure Vector Size
vdim=size(etypl)
! Allocate double of original size in temp vectors
allocate (etypltmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(etypl) .and. vdim .gt. 0) etypltmp(1:vdim)=etypl(1:vdim)
! Deallocate original Vector
if (allocated(etypl)) deallocate (etypl)
! Move allocations
call move_alloc(etypltmp,etypl)
end subroutine

subroutine resizeptypl(newsize)
implicit none
integer vdim,newsize
type(partype), allocatable, dimension(:) :: ptypltmp    !! Type of Particle list :: Size of nptyp
! if 0
if (newsize .lt. 1) then
    if (allocated(ptypl)) deallocate (ptypl)
endif
! Measure Vector Size
vdim=size(ptypl)
! Allocate double of original size in temp vectors
allocate (ptypltmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(ptypl) .and. vdim .gt. 0) ptypltmp(1:vdim)=ptypl(1:vdim)
! Deallocate original Vector
if (allocated(ptypl)) deallocate (ptypl)
! Move allocations
call move_alloc(ptypltmp,ptypl)
end subroutine

subroutine resizeparl(newsize)
implicit none
integer vdim,newsize
type(parpro), allocatable, dimension(:) :: parltmp        !! Particles List Description :: Size of Number of Particles (npar)
! if 0
if (newsize .lt. 1) then
    if (allocated(parl)) deallocate (parl)
endif
! Measure Vector Size
vdim=size(parl)
! Allocate double of original size in temp vectors
allocate (parltmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(parl) .and. vdim .gt. 0) parltmp(1:vdim)=parl(1:vdim)
! Deallocate original Vector
if (allocated(parl)) deallocate (parl)
! Move allocations
call move_alloc(parltmp,parl)
end subroutine

subroutine resizecvec(newsize)
implicit none
integer vdim,newsize
type(car), allocatable, dimension(:) :: rtmp   !! Elements Position Vector
type(car), allocatable, dimension(:) :: ftmp   !! Elements Force Vector
integer, allocatable, dimension(:)   :: ettmp !! Elements Type List of nele size
integer, allocatable, dimension(:)   :: pttmp !! Elements Type List of nele size
integer, allocatable, dimension(:)   :: petmp !! Elements Type List of nele size

! if 0
if (newsize .lt. 1) then
    if (allocated(r)) deallocate (r)
    if (allocated(f)) deallocate (f)
    if (allocated(et)) deallocate (et)
    if (allocated(pt)) deallocate (pt)
    if (allocated(pe)) deallocate (pe)
endif
! Measure Vector Size
vdim=size(r)
! Allocate double of original size in temp vectors
allocate (rtmp(newsize),ftmp(newsize),ettmp(newsize),pttmp(newsize),petmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(r) .and. vdim .gt. 0) rtmp(1:vdim)=r(1:vdim)
if (allocated(f) .and. vdim .gt. 0) ftmp(1:vdim)=f(1:vdim)
if (allocated(et) .and. vdim .gt. 0) ettmp(1:vdim)=et(1:vdim)
if (allocated(pt) .and. vdim .gt. 0) pttmp(1:vdim)=pt(1:vdim)
if (allocated(pe) .and. vdim .gt. 0) petmp(1:vdim)=pe(1:vdim)
! Deallocate original Vector
if (allocated(r)) deallocate (r)
if (allocated(f)) deallocate (f)
if (allocated(et)) deallocate (et)
if (allocated(pt)) deallocate (pt)
if (allocated(pe)) deallocate (pe)
! Move allocations
call move_alloc(rtmp,r)
call move_alloc(ftmp,f)
call move_alloc(ettmp,et)
call move_alloc(pttmp,pt)
call move_alloc(petmp,pe)
end subroutine

! Add Mono Particle (Particle with Single Element)
subroutine addmonoptyp(ename,pname,chg,dif,eps,sig)
implicit none
character*(*) ename
character*(*),optional,intent(in) :: pname
real,optional,intent(in)  :: chg    ! Element Charge
real,optional,intent(in)  :: dif    ! Element Diffusivity
real,optional,intent(in)  :: eps    ! Element Epsilon Lennard Jones
real,optional,intent(in)  :: sig    ! Element Sigma Lennard Jones
real            :: lchg,ldif,leps,lsig
lchg=0.0
ldif=0.0
leps=0.0
lsig=0.0
if (present(chg)) lchg=chg
if (present(dif)) ldif=dif
if (present(eps)) leps=eps
if (present(sig)) lsig=sig
call addetyp(ename,lchg,ldif,leps,lsig)
if (present(pname)) then
  call addptyp(1,pname)
else
  call addptyp(1,ename)
endif
ptypl(nptyp)%chg=lchg
ptypl(nptyp)%etyp(1)=netyp
call setcarzero(ptypl(nptyp)%r(1))
end subroutine

! get etyp from ename
function getetyp(ename)
implicit none
integer getetyp
character*(*) ename
getetyp=1
do while (getetyp.le.netyp)
   if (etypl(getetyp)%nam.eq.ename) return
   getetyp=getetyp+1
enddo
getetyp=0
write(*,'(a,a,a)') 'Type name ',trim(adjustl(ename)),' not found'
end function

! get ptyp from pname
function getptyp(pname)
implicit none
integer getptyp
character*(*) pname
getptyp=1
do while (getptyp.le.nptyp)
   if (ptypl(getptyp)%nam.eq.pname) return
   getptyp=getptyp+1
enddo
getptyp=0
write(*,'(a,a,a)') 'Type name ',trim(adjustl(pname)),' not found'
end function

! Edit Element Type
subroutine editetyp(etype,chg,dif,eps,sig)
integer,intent(in)        :: etype  ! Element Type Number
real,optional,intent(in)  :: chg    ! Element Charge
real,optional,intent(in)  :: dif    ! Element Diffusivity
real,optional,intent(in)  :: eps    ! Element Epsilon Lennard Jones
real,optional,intent(in)  :: sig    ! Element Sigma Lennard Jones
if (present(chg)) etypl(etype)%chg=chg
if (present(dif)) etypl(etype)%dif=dif
if (present(eps)) etypl(etype)%eps=eps
if (present(sig)) etypl(etype)%sig=sig
end subroutine

! Add Element Name to list
subroutine addetyp(nam,chg,dif,eps,sig)
implicit none
character*(*)             :: nam
real,optional,intent(in)  :: chg    ! Element Charge
real,optional,intent(in)  :: dif    ! Element Diffusivity
real,optional,intent(in)  :: eps    ! Element Epsilon Lennard Jones
real,optional,intent(in)  :: sig    ! Element Sigma Lennard Jones
netyp=netyp+1
netp=netyp*(netyp+1)/2
if (netyp .gt. size(etypl)) call resizeetypl(netyp)
etypl(netyp)%nam=nam
etypl(netyp)%chg=0.0
etypl(netyp)%dif=0.0
etypl(netyp)%eps=0.0
etypl(netyp)%sig=0.0
if (present(chg)) etypl(netyp)%chg=chg
if (present(dif)) etypl(netyp)%dif=dif
if (present(eps)) etypl(netyp)%eps=eps
if (present(sig)) etypl(netyp)%sig=sig
end subroutine

! Set Element Type in Particle Type
subroutine seteleinptyp(ptypn,nen,etyp)
implicit none
integer ptypn, nen, etyp
ptypl(ptypn)%etyp(nen)=etyp
end subroutine

subroutine updateptypchg(ptypn)
implicit none
integer ptypn,i
real chg
if (.not.allocated(ptypl(ptypn)%etyp)) return
chg=0.0
do i=1,ptypl(ptypn)%ne
  chg=chg+etypl(ptypl(ptypn)%etyp(i))%chg
enddo
ptypl(ptypn)%chg=chg
end subroutine

! Add Particle Type
subroutine addptyp(ne,pnam,chg)
implicit none
integer ne
character*(*) pnam
real,optional,intent(in) :: chg 
if (ne .lt. 1) return
nptyp=nptyp+1
if (nptyp .gt. size(ptypl)) call resizeptypl(nptyp)
ptypl(nptyp)%ne = ne
ptypl(nptyp)%nam = pnam
if (present(chg)) then
  ptypl(nptyp)%chg = chg
else
  ptypl(nptyp)%chg = 0.0
endif
allocate (ptypl(nptyp)%etyp(ne),ptypl(nptyp)%r(ne))
end subroutine

! Add Particle Type To Particle List
subroutine addpar(ptype,kind,ibuf)
implicit none
integer ptype,lastnele,i,j
integer,optional,intent(in) :: kind
integer,optional,intent(in) :: ibuf
if (ptype .gt. 0 .and. ptype .le. nptyp) then
  npar=npar+1
  lastnele=nele
  nele=nele+ptypl(ptype)%ne
  if (nele .gt. size(r)) call resizecvec(nele)
  if (npar .gt. size(parl)) call resizeparl(npar)
  parl(npar)%ptyp=ptype     ! Save Particle Type Idx in Particle List
  parl(npar)%sr=lastnele    ! Save Previous Position of last Added Element in Cartesian Vectors to Particle List
  parl(npar)%ne=ptypl(ptype)%ne  ! Save Number of Elements of this particle in the Particle List
  if (present(kind)) then 
    parl(npar)%kind=kind
  else
    parl(npar)%kind=0
  endif
  if (present(ibuf)) then
    parl(npar)%ibuf=ibuf
  else
    parl(npar)%ibuf=0
  endif
  ! Copy Coordinates from Particle Type Template 
  do i=1,parl(npar)%ne
    j=lastnele+i
    r(j)=ptypl(ptype)%r(i)
    et(j)=ptypl(ptype)%etyp(i)
    pt(j)=ptype
    pe(j)=npar
  enddo
endif
end subroutine

subroutine insertpar(parn,rcent,norot)
implicit none
type(car) :: rcent,arcent  ! Particle Centroid
integer parn ! Particle Number 
logical,optional,intent(in) :: norot
logical rotate
rotate = .true.
if (present(norot)) rotate=.not.norot
if (parl(parn)%ne .eq. 1) then
  call setcarmonopar(parn,rcent)
else
  ! Compute Actual Centroid
  call getcentroid(parn,arcent)
  ! Compute Displacement: Substract Actual Centroid to New Centroid
  call subcar2par(parn,arcent)
  ! Randomly rotate particle
  if (rotate) call uranrot(parn,nocenter=.true.)
  ! Add Displacement to all elements position
  call addcar2par(parn,rcent)
endif
end subroutine

! use center true when rotating a particle that is not at the origin
subroutine movepar(parn,rd,center,norot)
implicit none
type(car) :: rd  ! Vector Shift
integer parn ! Particle Number 
logical,optional,intent(in) :: norot,center
logical rotate, nocenter
rotate = .true.
nocenter = .true.
if (present(norot)) rotate=.not.norot
if (present(center)) nocenter=.not.center
if (parl(parn)%ne .eq. 1) then
  if (nocenter) then 
    call setcarmonopar(parn,rd)
  else
    call setcarmonopar(parn,sumcar(r(parl(parn)%sr+1),rd))
  endif
else
  ! Randomly rotate particle
  if (rotate) call uranrot(parn,nocenter=nocenter)
  ! Add Displacement to all elements position
  call addcar2par(parn,rd)
endif
end subroutine

! Randomly rotates a particle
! only use nocenter true when particle is at the origin
subroutine uranrot(parn,nocenter)
implicit none
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510
real, parameter :: twopi=2.0*pi
integer parn
real,external :: rndm
real phi, theta, psi, cosp, sinp, ocosp
real x,y,z,rot(3,3)
type(car) :: cent
logical,optional,intent(in) :: nocenter
logical center
center = .true.
if (present(nocenter)) center = .not.nocenter
if (parl(parn)%ne .lt. 2) return
theta=acos(2.0*rndm()-1.0) ! from 0 to pi
phi=twopi*rndm()           ! from 0 to 2*pi
psi=twopi*rndm()           ! from 0 to 2*pi
! Make Random Vector
sinp=sin(theta)
x=sinp*cos(phi)
y=sinp*sin(phi)
z=cos(theta)
! Build Rotation Matrix
cosp=cos(psi)
ocosp=1.0-cosp
sinp=sin(psi)
! First Row
rot(1,1)=cosp+x*x*ocosp
rot(1,2)=x*y*ocosp-z*sinp
rot(1,3)=x*z*ocosp+y*sinp
! Second Row
rot(2,1)=x*y*ocosp+z*sinp
rot(2,2)=cosp+y*y*ocosp
rot(2,3)=y*z*ocosp-x*sinp
! Third Row
rot(3,1)=x*z*ocosp-y*sinp
rot(3,2)=y*z*ocosp+x*sinp
rot(3,3)=cosp+z*z*ocosp
if (center) then
  ! Get Centroid
  call getcentroid(parn,cent)
  ! Remove Centroid
  call subcar2par(parn,cent)
endif
! Rotate Particle
call rotatepar(parn,rot)
if (center) then 
  ! Restore Centroid
  call addcar2par(parn,cent)
endif
end subroutine

! Given a rotation matrix, rotates particle
subroutine rotatepar(parn,rot)
implicit none
integer i,ne,sr,parn ! Particle Number 
real rot(3,3)
ne=parl(parn)%ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call rotatecar(r(sr+i),rot)
enddo
end subroutine

subroutine rotatecar(rc,mat)
implicit none
type(car) :: rc  ! Particle Centroid
real mat(3,3),x,y,z
x=mat(1,1)*rc%x+mat(1,2)*rc%y+mat(1,3)*rc%z
y=mat(2,1)*rc%x+mat(2,2)*rc%y+mat(2,3)*rc%z
z=mat(3,1)*rc%x+mat(3,2)*rc%y+mat(3,3)*rc%z
rc%x=x
rc%y=y
rc%z=z
end subroutine

subroutine getcentroid(parn,rc)
implicit none
integer i,ne,sr
real ine
integer,intent(in)    :: parn ! Particle Number 
type(car),intent(out) :: rc  ! Particle Centroid
ne=parl(parn)%ne
sr=parl(parn)%sr
if (ne.eq.1) then
  rc=r(sr+1)
else
  ine=1.0/ne
  ! Compute Actual Centroid
  ! Zero vector
  call setcarzero(rc)
  ! Add each position vector to rc
  do i=1,ne
     call addcar(rc,r(sr+i))
  enddo
  ! divide rc by number of elements
  call mulcar(rc,ine)
endif
end subroutine

subroutine putcoorinptyp(ptypn,elen,x,y,z)
implicit none
integer ptypn,elen
real x,y,z
call setcar(ptypl(ptypn)%r(elen),x,y,z)
end subroutine

subroutine getptypcentroid(ptypn,rc)
implicit none
integer,intent(in)    :: ptypn ! Particle Type Number 
type(car),intent(out) :: rc  ! Particle Type Centroid
integer i,ne
real ine
ne=ptypl(ptypn)%ne
if (ne.eq.1) then
  rc=ptypl(ptypn)%r(1)
else
  ine=1.0/ne
  ! Compute Actual Centroid
  ! Zero vector
  call setcarzero(rc)
  ! Add each position vector to rc
  do i=1,ne
     call addcar(rc,ptypl(ptypn)%r(i))
  enddo
  ! divide rc by number of elements
  call mulcar(rc,ine)
endif
end subroutine

subroutine centerptyp(ptypn)
implicit none
integer ptypn,i
type(car) :: rc
call getptypcentroid(ptypn,rc)
do i=1,ptypl(ptypn)%ne
  call subcar(ptypl(ptypn)%r(i),rc)
enddo
end subroutine

subroutine putcoorinpar(parn,elen,x,y,z)
implicit none
integer parn,elen
real x,y,z
call setcar(r(parl(parn)%sr+elen),x,y,z)
end subroutine

subroutine putmonopar(parn,x,y,z)
implicit none
integer parn
real x,y,z
call setcar(r(parl(parn)%sr+1),x,y,z)
end subroutine

subroutine setcarmonopar(parn,rcor)
implicit none
integer parn
type(car) rcor
r(parl(parn)%sr+1)=rcor
end subroutine

subroutine addcar2par(parn,rc)
implicit none
type(car) :: rc  ! Cartesian coordinates
integer i,ne,sr,parn ! Particle Number 
ne=parl(parn)%ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call addcar(r(sr+i),rc)
enddo
end subroutine 

subroutine subcar2par(parn,rc)
implicit none
type(car) :: rc  ! Cartesian coordinates
integer i,ne,sr,parn ! Particle Number 
ne=parl(parn)%ne
sr=parl(parn)%sr
! Add each position vector to rc
do i=1,ne
   call subcar(r(sr+i),rc)
enddo
end subroutine

function getcar(x,y,z)
implicit none
type(car) :: getcar
real :: x,y,z
getcar%x=x
getcar%y=y
getcar%z=z
end function

subroutine setcar(r,x,y,z)
implicit none
type(car) :: r
real :: x,y,z
r%x=x
r%y=y
r%z=z
end subroutine

subroutine setcarzero(r)
implicit none
type(car) :: r
call setcar(r,0.0,0.0,0.0)
end subroutine

!sum car
function sumcar(r1,r2)
implicit none
type(car) :: sumcar,r1,r2
sumcar%x=r1%x+r2%x
sumcar%y=r1%y+r2%y
sumcar%z=r1%z+r2%z
end function

! Add Second Cartesian Type to first
subroutine addcar(r1,r2)
implicit none
type(car) :: r1,r2
r1%x=r1%x+r2%x
r1%y=r1%y+r2%y
r1%z=r1%z+r2%z
end subroutine

! Subtract Cartesian Type
subroutine subcar(r1,r2)
implicit none
type(car) :: r1,r2
r1%x=r1%x-r2%x
r1%y=r1%y-r2%y
r1%z=r1%z-r2%z
end subroutine

! Square distance between two vectors in car type
function dist2car(r1,r2)
implicit none
real dist2car
type(car) :: r1,r2
dist2car=(r1%x-r2%x)**2+(r1%y-r2%y)**2+(r1%z-r2%z)**2
end function

! Dot Product of Cartesian Types
function dotcar(r1,r2)
implicit none
real dotcar
type(car) :: r1,r2
dotcar=r1%x*r2%x+r1%y*r2%y+r1%z*r2%z
end function

function dotcarvec(r,v)
implicit none
real dotcarvec,v(3)
type(car) :: r
dotcarvec=v(1)*r%x+v(2)*r%y+v(3)*r%z
end function

! Scalar Multiplication of Cartesian Types
subroutine mulcar(r,sc)
implicit none
type(car) :: r
real sc
r%x=sc*r%x
r%y=sc*r%y
r%z=sc*r%z
end subroutine

subroutine delparofkind(kind)
implicit none
integer :: i,kind
i=npar
do while (i.ge.1)
   if (parl(i)%kind.eq.kind) call delpar(i)
   i=i-1
enddo
end subroutine

subroutine delpar(parn)
implicit none
integer i,j,ii,jj,n,ne,sr,parn
! If parn is out of range exit
!if (parn.gt.npar.or.parn.lt.1) return
! if particle not of kind 3 return
!if (parl(parn)%kind.ne.3) return
! if last particle or no more particles
if (npar .le. 1) then
    nele=0
    npar=0
    return
endif
! If parn is the last particle just remove
if (parn .eq. npar) then
    nele=nele-parl(parn)%ne
    npar=npar-1
    return
endif
! Look for last particle of same number of elements
ne=parl(parn)%ne
n=npar
do while (parl(n)%ne .ne. ne .and. n.gt.parn)
    n=n-1
enddo
! If there is another of particle with the same ne closer to the end copy that one to parn's place
    ! Copy coordinates
if (n .gt. parn) then
  i  = parl(parn)%sr + 1
  ii = parl(parn)%sr + parl(parn)%ne
  j  = parl(n)%sr    + 1
  jj = parl(n)%sr    + parl(n)%ne
   r(i:ii)= r(j:jj)
  et(i:ii)=et(j:jj)
  pt(i:ii)=pt(j:jj)
  pe(i:ii)=parn
endif
! Move all coordinates backwards and remove that n particle from the coordinate vector
sr=parl(n)%sr
do i=sr+1,nele-ne
    j=i+ne
    r(i)=r(j)
    et(i)=et(j)
    pt(i)=pt(j)
    pe(i)=pe(j)-1
enddo
nele=nele-ne

! If there is another of particle with the same ne closer to the end copy that one to parn's place
! Copy particle 
if (n .gt. parn) parl(parn)%ptyp=parl(n)%ptyp

! Move all particles backwards and remove one particle from particle list
do i=n,npar-1
   parl(i)%ptyp=parl(i+1)%ptyp
   parl(i)%sr=parl(i+1)%sr-ne
   parl(i)%ne=parl(i+1)%ne
   parl(i)%kind=parl(i+1)%kind
   parl(i)%ibuf=parl(i+1)%ibuf
enddo
npar=npar-1
end subroutine

subroutine delparall()
implicit none
npar=0
nele=0
end subroutine

subroutine deltypall()
implicit none
nptyp=0 ! number of particle types
netyp=0 ! number of element types
end subroutine

function etpidx(itype,jtype)
implicit none
integer itype,jtype,etpidx,maxi,mini

maxi=MAX(itype,jtype)
mini=MIN(itype,jtype)
etpidx=maxi*(maxi-1)/2+mini
end function

subroutine printpdb(nunit)
implicit none
integer i
integer :: nunit
real ibuff

!write (nunit,'(A6,I5)') 'REMARK',nele
do i=1,nele
  ibuff=parl(pe(i))%ibuf
  write (nunit,'(A6,I5,x,A5,A5,I4,4x,3F8.3,2F6.2,6x,A4)') 'ATOM  ',i,etypl(et(i))%nam,ptypl(pt(i))%nam,pe(i),r(i)%x,r(i)%y,r(i)%z,ibuff,etypl(et(i))%chg,''
enddo
write (nunit,'(A)') 'END'
end subroutine

! Extended PDB
subroutine printpdbe(nunit)
implicit none
integer i
integer :: nunit
real ibuff
!write (nunit,'(A6,I5)') 'REMARK',nele
do i=1,nele
  ibuff=parl(pe(i))%ibuf
  write (nunit,'(A6,I5,x,A5,A5,I4,4x,3F16.8,2F6.2,6x,A4)') 'ATOM  ',i,etypl(et(i))%nam,ptypl(pt(i))%nam,pe(i),r(i)%x,r(i)%y,r(i)%z,ibuff,etypl(et(i))%chg,''
enddo
write (nunit,'(A)') 'END'
end subroutine

subroutine printxyz(nunit)
implicit none
integer i
integer :: nunit

! Print Particle Data
write(nunit,*) nele
write(nunit,*)
do i=1,nele
    write(nunit,*) etypl(et(i))%nam,r(i)%x,r(i)%y,r(i)%z
enddo
end subroutine

subroutine printlists(nunit)
implicit none
integer i,j,k
integer :: nunit

! Print Particle Data
do i=1,npar
   write(nunit,*) i,parl(i)%ptyp,parl(i)%sr,parl(i)%ne
   do j=1,parl(i)%ne
      k=j+parl(i)%sr
      write(nunit,*) '      ',ptypl(parl(i)%ptyp)%etyp(j),r(k)%x,r(k)%y,r(k)%z
   enddo
enddo
end subroutine

end module
