program readprmtest
use charfuncmod
implicit none
character com*2048,word*4
logical ok
open(unit=1,file='../test/par_all36_pnclm-newions-namd-newMg-newLi.prm')

call getlin(com,1,6)
do while (lcase(com(1:3)).ne.'end')
  if (checkiflabel(com(1:4))) write(*,'(A)') '*******************'
  write(*,'(A)') trim(com)
  call getlin(com,1,6)
end do
write(*,'(A)') trim(com)
write(*,'(A,I0,A)') '|',getiprm('hola como 8',3),'|'
close(1)
contains
  logical*1 function checkiflabel(word)
  character word*(*)
  character wrd*(len_trim(word))
  checkiflabel = .false.
  if (len_trim(lcase(word)).lt.3) return
  wrd=lcase(word)
  if (wrd.eq.'atom') then 
    checkiflabel = .true.
  elseif (wrd.eq.'bond') then
    checkiflabel = .true.
  elseif (wrd.eq.'angl') then
    checkiflabel = .true.
  elseif (wrd.eq.'thet') then
    checkiflabel = .true.
  elseif (wrd.eq.'dihe') then
    checkiflabel = .true.
  elseif (wrd.eq.'impr') then
    checkiflabel = .true.
  elseif (wrd.eq.'imph') then
    checkiflabel = .true.
  elseif (wrd.eq.'nonb') then
    checkiflabel = .true.
  elseif (wrd.eq.'nbon') then
    checkiflabel = .true.
  elseif (wrd.eq.'cmap') then
    checkiflabel = .true.
  elseif (wrd.eq.'nbfi') then
    checkiflabel = .true.
  elseif (wrd.eq.'hbon') then
    checkiflabel = .true.
  elseif (wrd(1:3).eq.'phi') then
    checkiflabel = .true.
  endif
  return
  end function

  ! get parameter pn from line
  integer function getiprm(str,pn)
  implicit none
  integer num,pn,kode
  character str*(*)
  character getprm*(len(trim(adjustl(str))))
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


end program
