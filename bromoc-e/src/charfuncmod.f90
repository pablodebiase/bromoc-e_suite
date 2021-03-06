!    BROMOC-E
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

module charfuncmod
implicit none
contains
  ! check if a single char is a number 
  function isnum(chr)
  implicit none
  integer i
  character*1 chr
  character*10,parameter :: nums='0123456789'
  logical*1 isnum
  
  isnum=.false.
  
  do i=1,10
    if (chr.eq.nums(i:i)) then
      isnum=.true.
      return
    endif
  enddo
  end function

  function isrnum(chr)
  implicit none
  integer i
  character*1 chr
  character*17,parameter :: nums='0123456789.-+deED'
  logical*1 isrnum

  isrnum=.false.

  do i=1,17
    if (chr.eq.nums(i:i)) then
      isrnum=.true.
      return
    endif
  enddo
  end function
  
  ! convert integer to character
  function num2str(n,mx)
  implicit none
  integer n,num,i,a,mx
  character num2str*(mx),numero*10
  numero='0123456789'
  num=n
  do i=mx,1,-1
    a=int(num/(10**(i-1)))
    num2str(mx-i+1:mx-i+1)=numero(a+1:a+1)
    num=num-a*10**(i-1)
  enddo
  end function

  ! convert to lowercase
  function lcase(inchar)
  implicit none
  integer i
  integer*1 s
  character ( len = * ) inchar
  character lcase*(len_trim(inchar))
  lcase=trim(inchar)
  do i=1,len_trim(inchar)
    s=itoi1(iachar(lcase(i:i)))
    if (s.ge.65.and.s.le.90) lcase(i:i)=achar(s+32)
  enddo
  end function
  
  ! convert default real into real*4 
  function sng(realn)
  implicit none
  real realn
  real*4 sng
  sng=realn
  end function
 
  ! convert to uppercase
  function ucase(inchar)
  implicit none
  integer i
  integer*1 s
  character ( len = * ) inchar
  character ucase*(len_trim(inchar))
  ucase=trim(inchar)
  do i=1,len_trim(inchar)
    s=itoi1(iachar(ucase(i:i)))
    if (s.ge.97.and.s.le.122) ucase(i:i)=achar(s-32)
  enddo
  end function

  ! convert default integer to integer1  
  function itoi1(i)
  implicit none
  integer i
  integer*1 itoi1
  itoi1=i
  end function

  ! convert integer8 to default integer
  function i8toi(i8)
  implicit none
  integer*8 i8
  integer i8toi
  i8toi=i8
  end function

  ! convert character to integer
  function chr2int(str)
  implicit none
  integer chr2int,kode
  character str*(*)
  read(str,*,iostat=kode) chr2int
  if (kode.ne.0) stop 'Not an integer'
  end function

  ! convert character to integer
  function chr2int8(str)
  implicit none
  integer*8 chr2int8
  integer kode
  character str*(*)
  read(str,*,iostat=kode) chr2int8
  if (kode.ne.0) stop 'Not an integer'
  end function

  ! convert character to real8
  function chr2real(str)
  implicit none
  integer kode
  real chr2real
  character str*(*)
  read(str,*,iostat=kode) chr2real
  if (kode.ne.0) stop 'Not a real'
  end function
  
  ! get parameter pn from line
  function getparm(str,num,llim,ulim,pn)
  implicit none
  integer num,pn
  integer llim(*),ulim(*)
  character getparm*(ulim(pn)-llim(pn)+1),str*(*)
  if (pn.gt.num.or.pn.lt.1.or.num.lt.1) then
    getparm=''
  else
    getparm=str(llim(pn):ulim(pn))
  endif
  end function

  function isint(str)
  implicit none
  character str*(*)
  integer test,kode
  logical*1 isint
  read(str,*,iostat=kode) test
  if (kode.eq.0) then
    isint=.true.
  else
    isint=.false.
  endif
  end function
  
  function isreal(str)
  implicit none
  character str*(*)
  integer kode
  real test
  logical*1 isreal
  read(str,*,iostat=kode) test
  if (kode.eq.0) then
    isreal=.true.
  else
    isreal=.false.
  endif
  end function

  function isword(str)
  implicit none
  integer i,s
  character str*(*)
  logical*1 isword
  isword=.true.
  if (len_trim(str).eq.0) isword=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (.not.((s.ge.65.and.s.le.90).or.(s.ge.97.and.s.le.122))) then
      isword=.false.
      return
    endif
  enddo
  end function

  function iswordnum(str)
  implicit none
  integer i,s
  character str*(*)
  logical*1 iswordnum
  iswordnum=.true.
  if (len_trim(str).eq.0) iswordnum=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (.not.((s.ge.65.and.s.le.90).or.(s.ge.97.and.s.le.122).or.(s.ge.48.and.s.le.57))) then
      iswordnum=.false.
      return
    endif
  enddo
  end function
  
  function isalfa(str)
  implicit none
  integer i,s
  character str*(*)
  logical*1 isalfa
  isalfa=.true.
  if (len_trim(str).eq.0) isalfa=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (s.le.32.or.s.ge.127) then
      isalfa=.false.
      return
    endif
  enddo
  end function

  function check(stfin,stin)
  implicit none
  character*(*) stfin,stin
  character stfout*(len(stfin))
  integer, dimension(len_trim(stfin)) :: ll,ul
  integer n,i,j
  logical*1 check
  check=.false.
  stfout=''
  call findparm(stfin,n,ll,ul)
  do i=1,n
    if (lcase(getparm(stfin,n,ll,ul,i)).eq.lcase(adjustl(stin))) then
      check=.true.
      do j=1,n
        if (j.ne.i) stfout=trim(stfout)//' '//getparm(stfin,n,ll,ul,j)
      enddo
      stfin=adjustl(stfout)
      return 
    endif
  enddo
  end function

  ! reads next line skipping blank lines. True if end is not reached
  function setline(u,line)
  implicit none
  integer u,kode
  character*(*) line
  logical*1 setline
  read(u,'(a)',iostat=kode) line
  do while (len_trim(line).eq.0.and.kode.eq.0)
    read(u,'(a)',iostat=kode) line
  enddo
  setline=kode.eq.0
  end function

  ! reads the first word. If no word returns setword false
  function setword(word, line)
  implicit none
  character*(*) word,line
  logical*1 setword
  call getfirst(line,word)
  setword=len_trim(word).ne.0
  end function
  
  function setint(intvar,line)
  implicit none
  character*(*) line
  character*(len_trim(line)) word
  integer intvar
  logical*1 setint
  call getfirst(line,word)
  setint=isint(word)
  if (setint) intvar=chr2int(word)
  end function

  ! count the number of parameters in line
  integer function countparm(str)
  implicit none
  integer i,length,num
  character str*(*)
  logical chng
  length=len_trim(str)
  chng=.false.
  num=0
  do i=1,length
    if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
      chng=.false.
    else
      if (.not.chng) num=num+1
      chng=.true.
    endif
  enddo
  countparm=num
  return
  end function

  ! get parameter pn from line
  integer function getiprm(str,pn)
  implicit none
  integer num,pn,kode
  character str*(*)
  integer llim(len_trim(str)),ulim(len_trim(str))
  call findparm(str,num,llim,ulim)
  if (pn.gt.num.or.pn.lt.1.or.num.lt.1) then
    getiprm=0
  else
    read(str(llim(pn):ulim(pn)),*,iostat=kode) getiprm
    if (kode.ne.0) getiprm=0
  endif
  return
  end function

end module
