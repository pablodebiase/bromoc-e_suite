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
end type parpro

! Particle Type
type :: partype
    integer :: ne                              !! Number of Elements in Particle
    character*4                      :: pnam   !! Particle Name
    integer,allocatable,dimension(:) :: etyp   !! Particle Element Types Vector :: Size of ne
    type(car),allocatable,dimension(:) :: r    !! Position Vector of Particle Elements :: Size of ne
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

!! Element-Type List
integer                                  :: netyp     !! Number of Element Types
type(eletype),allocatable,dimension(:)   :: etypl     !! Element Name List 

!! Element-Type Pairs List
integer                                  :: netp

!! Particle-Type List
integer                                  :: nptyp      !! Number of Particle Types
type(partype), allocatable, dimension(:) :: ptypl    !! Type of Particle list :: Size of nptyp

!! Force and position vectors
integer nele                                !! Total Number of elements
type(car), allocatable, dimension(:) :: r  !! Elements Position Vector
type(car), allocatable, dimension(:) :: f  !! Elements Force Vector
integer, allocatable, dimension(:) :: etypls !! Elements Type List of nele size
integer, allocatable, dimension(:) :: ptypls !! Particle Type List nele size
integer, allocatable, dimension(:) :: pels !! Particle List of nele size

!! Used Element Types List
integer nuet                      !! Total Number of Used Element Types in etypls
integer, allocatable, dimension(:) :: uetl

contains

subroutine Merge(A,NA,B,NB,C,NC)
integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
integer, intent(in)     :: B(NB)
integer, intent(in out) :: C(NC)
integer :: I,J,K
I = 1; J = 1; K = 1;
do while(I <= NA .and. J <= NB)
   if (A(I) <= B(J)) then
      C(K) = A(I)
      I = I+1
   else
      C(K) = B(J)
      J = J+1
   endif
   K = K + 1
enddo
do while (I <= NA)
   C(K) = A(I)
   I = I + 1
   K = K + 1
enddo
return
end subroutine merge
 
recursive subroutine MergeSort(A,N,T)
integer, intent(in) :: N
integer, dimension(N), intent(in out) :: A
integer, dimension((N+1)/2), intent (out) :: T
integer :: NA,NB,V
if (N < 2) return
if (N == 2) then
   if (A(1) > A(2)) then
      V = A(1)
      A(1) = A(2)
      A(2) = V
   endif
   return
endif      
NA=(N+1)/2
NB=N-NA

call MergeSort(A,NA,T)
call MergeSort(A(NA+1),NB,T)

if (A(NA) > A(NA+1)) then
   T(1:NA)=A(1:NA)
   call Merge(T,NA,A(NA+1),NB,A,N)
endif
return
end subroutine MergeSort

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
if (allocated(uetl)) then
   if (size(uetl).lt.netyp) then
     deallocate (uetl)
     allocate (uetl(netyp))
   endif
else
   allocate (uetl(netyp))
endif
! Update Database
call updatetypels()
! Get Unused etyp
nuet=1
uetl(nuet)=etypls(1)
do i = 2,nele
   j=1
   f=.false.
   do while (j.le.nuet)
     if (etypls(i).eq.uetl(j)) then
       f=.true.
       exit
     endif
     j=j+1
   enddo
   if (.not.f) then
     nuet=nuet+1
     uetl(nuet)=etypls(i)
   endif
enddo
! Sort uetl
call msort(uetl(1:nuet),nuet)
end subroutine

subroutine updatetypels()
implicit none
integer newsize,i,j,ne,sr,ptype
if (allocated(etypls)) deallocate (etypls)
if (allocated(ptypls)) deallocate (ptypls)
if (allocated(pels)) deallocate (pels)
newsize=size(r)
allocate (etypls(newsize),ptypls(newsize),pels(newsize))
do i=1,npar
  sr=parl(i)%sr
  ne=parl(i)%ne
  ptype=parl(i)%ptyp
  do j=1,ne
    ptypls(sr+j)=ptype
    etypls(sr+j)=ptypl(ptype)%etyp(j)
    pels(sr+j)=i
  enddo
enddo
end subroutine

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
type(car), allocatable, dimension(:) :: rtmp  !! Elements Position Vector
type(car), allocatable, dimension(:) :: ftmp  !! Elements Force Vector
! if 0
if (newsize .lt. 1) then
    if (allocated(r)) deallocate (r)
    if (allocated(f)) deallocate (f)
endif
! Measure Vector Size
vdim=size(r)
! Allocate double of original size in temp vectors
allocate (rtmp(newsize),ftmp(newsize))
! Copy Vectors
if (vdim .gt. newsize) vdim = newsize
if (allocated(r) .and. vdim .gt. 0) rtmp(1:vdim)=r(1:vdim)
if (allocated(f) .and. vdim .gt. 0) ftmp(1:vdim)=f(1:vdim)
! Deallocate original Vector
if (allocated(r)) deallocate (r)
if (allocated(f)) deallocate (f)
! Move allocations
call move_alloc(rtmp,r)
call move_alloc(ftmp,f)
end subroutine

! Add Mono Particle (Particle with Single Element)
subroutine addmonopar(pename,pname)
implicit none
character*(*) pename
character*(*),optional,intent(in) :: pname
call addetyp(pename)
if (present(pname)) then
  call addptyp(1,pname)
else
  call addptyp(1,pename)
endif
ptypl(nptyp)%etyp(1)=netyp
call setcarzero(ptypl(nptyp)%r(1))
end subroutine

! get etyp from ename
function getetyp(ename)
implicit none
integer getetyp,i
character*(*) ename
logical doloop
doloop=.true.
i=1
do while (i.le.netyp)
   if (etypl(i)%nam.eq.ename) exit
   i=i+1
enddo
getetyp=i
end function

! Add Element Name to list
subroutine addetyp(nam,chg,dif,eps,sig)
implicit none
character*(*)             :: nam
real,optional,intent(in)  :: chg    ! Element Charge
real,optional,intent(in)  :: dif    ! Element Diffusivity
real,optional,intent(in)  :: eps    ! Element Epsilon Lennard Jones
real,optional,intent(in)  :: sig    ! Element Sigma Lennard Jones
netyp=netyp+1
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

! Add Particle Type
subroutine addptyp(ne,pnam)
implicit none
integer ne
character*(*) pnam
if (ne .lt. 1) return
nptyp=nptyp+1
if (nptyp .gt. size(ptypl)) call resizeptypl(nptyp)
ptypl(nptyp)%ne = ne
ptypl(nptyp)%pnam = pnam
allocate (ptypl(nptyp)%etyp(ne),ptypl(nptyp)%r(ne))
end subroutine

! Add Particle Type To Particle List
subroutine addpar(ptype)
implicit none
integer ptype,lastnele,i
if (ptype .gt. 0 .and. ptype .le. nptyp) then
  npar=npar+1
  lastnele=nele
  nele=nele+ptypl(ptype)%ne
  if (nele .gt. size(r)) call resizecvec(nele)
  if (npar .gt. size(parl)) call resizeparl(npar)
  parl(npar)%ptyp=ptype     ! Save Particle Type Idx in Particle List
  parl(npar)%sr=lastnele    ! Save Previous Position of last Added Element in Cartesian Vectors to Particle List
  parl(npar)%ne=ptypl(ptype)%ne  ! Save Number of Elements of this particle in the Particle List
  ! Copy Coordinates from Particle Type Template 
  do i=1,parl(npar)%ne
    r(lastnele+i)=ptypl(ptype)%r(i)
  enddo
endif
end subroutine

subroutine movepar(parn,rcent)
implicit none
type(car) :: rcent,arcent  ! Particle Centroid
integer parn ! Particle Number 
if (parl(parn)%ne .eq. 1) then
  call setcarmonopar(parn,rcent)
else
  ! Compute Actual Centroid
  call getcentroid(arcent,parn)
  ! Compute Displacement: Substract Actual Centroid to New Centroid
  call subcar2par(parn,arcent)
  ! Add Displacement to all elements position
  call addcar2par(parn,rcent)
endif
end subroutine

! Randomly rotates a particle
subroutine uranrot(parn)
implicit none
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510
real, parameter :: twopi=2.0*pi
integer parn
real,external :: rndm
real phi, theta, psi, cosp, sinp, ocosp
real x,y,z,rot(3,3)
type(car) :: cent
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
! Get Centroid
call getcentroid(cent,parn)
! Remove Centroid
call subcar2par(parn,cent)
! Rotate Particle
call rotatepar(parn,rot)
! Restore Centroid
call addcar2par(parn,cent)
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

subroutine getcentroid(rc,parn)
implicit none
type(car) :: rc  ! Particle Centroid
integer i,ne,sr,parn ! Particle Number 
real ine
ne=parl(parn)%ne
ine=1.0/ne
sr=parl(parn)%sr
! Compute Actual Centroid
! Zero vector
call setcarzero(rc)
! Add each position vector to rc
do i=1,ne
   call addcar(rc,r(sr+i))
enddo
! divide rc by number of elements
call mulcar(rc,ine)
end subroutine

subroutine putcoorinptyp(ptypn,elen,x,y,z)
implicit none
integer ptypn,elen
real x,y,z
call setcar(ptypl(ptypn)%r(elen),x,y,z)
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

subroutine delpar(parn)
implicit none
integer i,n,ne,sr,parn
if (npar .le. 1) then
    nele=0
    npar=0
    return
endif
! If no more particles or parn is out of range exit
if (parn.gt.npar.or.parn.lt.1) return
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
if (n .gt. parn) r(parl(parn)%sr+1:parl(parn)%sr+parl(parn)%ne)=r(parl(n)%sr+1:parl(n)%sr+parl(n)%ne)
! Move all coordinates backwards and remove that n particle from the coordinate vector
sr=parl(n)%sr
do i=sr+1,nele-ne
    r(i)=r(i+ne)
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

call updatetypels()
do i=1,nele
  write (nunit,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,etypl(etypls(i))%nam,ptypl(ptypls(i))%pnam,pels(i),r(i)%x,r(i)%y,r(i)%z
enddo
write (nunit,'(A)') 'END'
end subroutine

subroutine printxyz(nunit)
implicit none
integer i
integer :: nunit

! Print Particle Data
call updatetypels()
write(nunit,*) nele
write(nunit,*)
do i=1,nele
    write(nunit,*) etypl(etypls(i))%nam,r(i)%x,r(i)%y,r(i)%z
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
