!    CRD2PDB - Converts CRD (CHARMM Coordinates) to PDB
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

program crd2pdb
implicit none
! pdb and crd
real*8 :: w,e
integer*4 ::  nop,ires
character*4 :: typ,res,segid,resid
real*8 :: rt(3)

integer*4 i,j,k,narg,arg,kode,ll(256),ul(256),num
character inpfile*256,outfile*256,ln*256
logical once,readexton

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('Input CHARMM CRD filename: ',narg,arg,inpfile)
call readarg('Output PDB filename: ',narg,arg,outfile)

open(unit=1,file=inpfile,IOSTAT=kode)
open(unit=2,file=outfile)
once=.true.
do while (once)
  read(1,'(A)') ln
  ln=adjustl(ln)
  if (ln(1:1).ne.'*') then
    call findparm(ln,num,ll,ul)
    if (num.eq.0) stop 'Error reading CRD'
    read(ln(ll(1):ul(1)),*) nop
    readexton=.false.
    if (num.eq.2) then
      if (ln(ll(2):ul(2)).eq.'EXT') readexton=.true.
    endif
    once=.false.
  endif
enddo
if (readexton) then
  write(*,'(/A)') 'Input CRD: EXTended format detected.'
else
  write(*,'(/A)') 'Input CRD: Regular format detected.'
endif

write(*,'(A,I0)') 'Number of atoms: ',nop
e=1.0
do j=1,nop
  if (readexton) then
    read(1,'(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)') k,ires,res,typ,(rt(i),i=1,3),segid,resid,w
  else
    read(1,'(I5,I5,1x,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)') k,ires,res,typ,(rt(i),i=1,3),segid,resid,w
  endif
  write (2,'(A6,I5,3(x,A4),4x,3F8.3,2F6.2,6x,A4)') 'ATOM  ',k,typ,res,resid,rt(1),rt(2),rt(3),e,w,segid
enddo
write (2,'(A)') 'END'
close(1)
close(2)
write(*,'(/A)') 'Normal termination of DNACDF'
end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*32

prname='CRD2PDB'
prver='version 1.0'
prdesc='Converts CHARMM CRD to PDB'
author='Pablo M. De Biase'
startdate='21 Aug 2012'
lastdate='21 Aug 2012'

write(*,'(/A)') trim(prname)//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
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
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! identify words from line
subroutine findparm(str,num,llim,ulim)
implicit none
integer i,length,num
integer ulim(*),llim(*)
character str*(*)
logical chng
length=len_trim(str)
chng=.false.
num=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
      ulim(num)=0
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
end subroutine

