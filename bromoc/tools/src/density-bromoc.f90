!    DENSITY - Computes the ions density along z-axis using BROMOC trajectory
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

module comun
implicit none
integer*4,allocatable :: gfi(:),gff(:),fg(:),fp(:),ifp(:),typ(:),gft(:,:)
integer*4 gfn,fn,nsc,ion,ntop,itype,nframe,nn,maxntop,snd,fag,dna,pn,ni,bins
real*8 runtime
real*8 dt,uli,lli,cr,ori(3),res
logical inpopen
character*4, allocatable :: atnam2(:)
character bs*8
integer*4,allocatable :: fti(:),ftf(:),ftl(:),pi(:),pni(:),pnf(:)
character*4,allocatable :: fl(:),gfl(:)
real*8,allocatable :: rt(:,:),q(:),pz(:)
integer*8,allocatable :: densi(:,:)
logical termon
end module

program density_bromoc
use comun
implicit none
real*8 aa,bb
integer*4 h,i,j,k,l,m,a,b,c,d,nii,narg,arg
character inpfile*256,line*256
integer*4 ul(256),ll(256),num,fac

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

a=index(inpfile,'.',back=.true.)

! Read monk input file
call readmonkhead(inpfile,.true.)

! Build residue type list
call readarg('Input Sequence of Species Numbers: ',narg,arg,line)
call findparm(line,256,num,ll,ul)
fn=num
gfn=1
allocate (gfi(gfn),gff(gfn),gfl(gfn),fl(fn),fg(fn),fp(fn),ifp(itype))
do i=1,fn
  read(line(ll(i):ul(i)),*) fp(i)
  ifp(fp(i))=i  
  fl(i)=atnam2(fp(i))
enddo
gfl(1)='IONS'
gfi(1)=1
gff(1)=fn
fg(1:fn)=1

ion=1    ! the group fragment for ions
ni=gff(ion)-gfi(ion)+1

allocate (fti(ni),ftf(ni),q(ni),pni(ni),pnf(ni))

call readarg('Lower Limit, Upper Limit, Cylinder Radius, Origin (x,y), Resolution: ',narg,arg,line)
read(line,*) lli,uli,cr,ori(1),ori(2),res
cr=cr**2

bins=int((uli-lli)/res)+1
allocate (densi(bins,ni+1))
densi = 0

!Read first frame
call readmonkbody()

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

do while (inpopen)                        ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc              ! print frame number
  call fragilist()                        ! build fragment elements list for ions
  call calcdens()                         ! compute density
  call readmonkbody()                     ! read next frame
enddo
inpfile=inpfile(1:a)
inpfile(a+1:)='dens'
! Print out result
call printout(inpfile)

write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine calcdens()
use comun
implicit none
integer i,j,k,l,m,c,d,a,b
real*8 v,dv
logical go
do j=gfi(ion),gff(ion) ! for each ion type
  m=j-gfi(ion)+1
  c=fti(m) ! first atom of the selected ion type
  d=ftf(m) ! last atom of the selected ion type
  do i=c,d
    a=ftl(i)
    dv=(rt(1,a)-ori(1))**2+(rt(2,a)-ori(2))**2
    v=rt(3,a)
    go=(v.le.uli.and.v.ge.lli.and.dv.le.cr)
    if (go) then
      b=int((v-lli)/res)+1
      densi(b,m)=densi(b,m)+1
      densi(b,ni+1)=densi(b,ni+1)+1
    endif
  enddo
enddo
end subroutine

subroutine fragilist()
use comun
implicit none
integer i,j,k,m

if (allocated(ftl)) deallocate (ftl)
allocate (ftl(ntop))

k=0
m=0
do i=gfi(ion),gff(ion)
  m=m+1
  fti(m)=k+1       ! initial number in ftl list for fragment type i
  do j=1,ntop
    if (fp(i).eq.typ(j)) then
      k=k+1
      ftl(k)=j     ! list that contain the position for each fragment type sorted by fragment type
    endif
  enddo
  ftf(m)=k         ! final number in ftl list for fragment type i
enddo
end subroutine

subroutine readmonkhead(inpfile,yes)
use comun
implicit none
integer i, j, kode
character*256 inpfile
character*2048 line
logical yes

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind')
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) itype                  ! number of ions and nucleotides
if (.not.allocated(atnam2)) allocate (atnam2(itype))
read(14,iostat=kode) (atnam2(j),j=1,itype)  ! ion and nucleotides types in char
if (yes) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(line,*)'Fragment types: ',(j,' '//atnam2(j),j=1,itype)
  write(*,'(A)') trim(line)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(14)
endif
end subroutine

subroutine readmonkbody()
use comun
implicit none
integer kode,i,k
real*4, allocatable :: rtt(:,:)

if (allocated(typ)) deallocate (typ)
if (allocated(rt)) deallocate (rt)
if (allocated(rtt)) deallocate (rtt)
read(14,iostat=kode) runtime !simulation time for each step (ns)
read(14,iostat=kode) ntop ! number of particles
k=kode
if (k.eq.0) then 
  allocate (rtt(3,ntop),rt(3,ntop),typ(ntop))
  read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(14,iostat=kode) (rtt(1,i),i=1,ntop)
  read(14,iostat=kode) (rtt(2,i),i=1,ntop)
  read(14,iostat=kode) (rtt(3,i),i=1,ntop)
endif
if (kode.eq.0) then
  nsc=nsc+1
!  write(*,*) runtime,ntop,(typ(i),i=1,ntop)
  rt(1:3,1:ntop)=dble(rtt(1:3,1:ntop))
else
  close(14)
  inpopen=.false.
endif
if (k.eq.0) deallocate (rtt)
end subroutine

subroutine checkmaxntop()
use comun
implicit none
integer kode,i
real*4, allocatable :: rtt(:,:)

maxntop=0
if (inpopen) kode=0
do while (kode.eq.0)
  if (allocated(typ)) deallocate (typ)
  if (allocated(rtt)) deallocate (rtt)
  read(14,iostat=kode) runtime !simulation time for each step (ns)
  read(14,iostat=kode) ntop ! number of particles
  if (kode.eq.0) then
    if (ntop.gt.maxntop) maxntop=ntop
    allocate (rtt(3,ntop),typ(ntop))
    read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
    read(14,iostat=kode) (rtt(1,i),i=1,ntop)
    read(14,iostat=kode) (rtt(2,i),i=1,ntop)
    read(14,iostat=kode) (rtt(3,i),i=1,ntop)
  endif
  if (kode.ne.0) close(14)
enddo
end subroutine

subroutine printout(filename)
use comun
implicit none
integer i,j,k,h,m
real*8 x,itot
character filename*256,line*2048
itot=1.0/sum(densi(1:bins,ni+1))

! Print out result
open(unit=1,file=trim(filename))

write(line,*) '# z   ',(fl(i),i=gfi(ion),gff(ion)),'  Total'
  write(1,'(A)') trim(line)
do i=1,bins
  x=(i-1)*res+lli
  write(line,*) x,(densi(i,j)*itot,j=1,ni+1)
  write(1,'(A)') trim(line)
enddo
close(1)
end subroutine

subroutine findparm(str,dimn,num,llim,ulim)
implicit none
integer dimn,i,length,num
integer ulim(dimn),llim(dimn)
character ( len = dimn ) str
logical chng
length=len_trim(str)
chng=.false.
num=0
ulim(1:dimn)=0
llim(1:dimn)=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
end subroutine

! read arguments
subroutine readarg(ques,narg,num,text)
implicit none
integer*4 narg,num
character text*(*),ques*(*)

num=num+1
write(*,'(/A$)') ques
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
  write(*,'(A)') trim(text)
else
  text=''
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='DENSITY'
prver='version 1.0'
prdesc='Computes the ions density along z-axis using BROMOC trajectory (.btr)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='23 Apr 2016'
lastdate='23 Apr 2016'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

