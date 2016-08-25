subroutine psf_imp(ideform,iat1,iat2,iat3,iat4,val,psf_btype)

use errormod
use explatmod
use charmmmod

implicit none
integer ideform,iat1,iat2,iat3,iat4
integer val(natt),psf_btype(4,ndeforms)
! local variables
integer ib1,ib2,ibonds,itype,jtype,ktype,ltype,j
logical*1 ok,ok1,ok2

if (iat1.le.0 .or. iat2.le.0 .or. iat3.le.0 .or. iat4.le.0 .or. iat1.gt.natt .or. iat2.gt.natt .or. iat3.gt.natt .or.iat4.gt.natt & 
    .or. iat1.eq.iat2 .or. iat1.eq.iat3 .or. iat1.eq.iat4 .or. iat2.eq.iat3 .or. iat2.eq.iat4 .or. iat3.eq.iat4) call error ('readpsf', 'Wrong improper angle', faterr)
! check there is not duplicate improper angles
ok = .false.
j = 0
do while (j.lt.ideform .and. .not.ok)
  j = j + 1
  ok = (iat1.eq.deforms(1,j).and.iat2.eq.deforms(2,j).and.iat3.eq.deforms(3,j).and.iat4.eq.deforms(4,j)) .or. &
       (iat1.eq.deforms(4,j).and.iat2.eq.deforms(3,j).and.iat3.eq.deforms(2,j).and.iat4.eq.deforms(1,j))
enddo
if (ok) call error ('psf_imp', 'Improper angles are duplicated', faterr)
ideform = ideform + 1
deforms(1,ideform) = iat1   
deforms(2,ideform) = iat2 
deforms(3,ideform) = iat3
deforms(4,ideform) = iat4 
! determine which bonds make up each deformation  
! a) iat1 is the central atom
ib1 = deforms(1,ideform)
ib2 = deforms(2,ideform)
if (ib2.lt.ib1) then
  j = ib1
  ib1 = ib2
  ib2 = j
endif
ok1 = .false.
ibonds = 0
do while (ibonds.lt.nbonds .and. .not.ok1)
  ibonds = ibonds + 1
  ok1 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
enddo
if (ok1) then
  ib1 = deforms(1,ideform)
  ib2 = deforms(3,ideform)
  if (ib2.lt.ib1) then
    j = ib1
    ib1 = ib2
    ib2 = j
  endif
  ok1 = .false.
  ibonds = 0
  do while (ibonds.lt.nbonds .and. .not.ok1)
    ibonds = ibonds + 1
    ok1 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
  enddo
endif
if (ok1) then
  ib1 = deforms(1,ideform)
  ib2 = deforms(4,ideform)
  if (ib2.lt.ib1) then
    j = ib1
    ib1 = ib2
    ib2 = j
  endif
  ok1 = .false.
  ibonds = 0
  do while (ibonds.lt.nbonds .and. .not.ok1)
    ibonds = ibonds + 1
    ok1 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
  enddo
endif
! b) iat4 is the central atom
ib1 = deforms(4,ideform)
ib2 = deforms(2,ideform)
if (ib2.lt.ib1) then
  j = ib1
  ib1 = ib2
  ib2 = j
endif
ok2 = .false.
ibonds = 0
do while (ibonds.lt.nbonds .and. .not.ok2)
  ibonds = ibonds + 1
  ok2 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
enddo
if (ok2) then
  ib1 = deforms(4,ideform)
  ib2 = deforms(3,ideform)
  if (ib2.lt.ib1) then
    j = ib1
    ib1 = ib2
    ib2 = j
  endif
  ok2 = .false.
  ibonds = 0
  do while (ibonds.lt.nbonds .and. .not.ok2)
    ibonds = ibonds + 1
    ok2 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
  enddo
endif
if (ok2) then
  ib1 = deforms(4,ideform)
  ib2 = deforms(1,ideform)
  if (ib2.lt.ib1) then
    j = ib1
    ib1 = ib2
    ib2 = j
  endif
  ok2 = .false.
  ibonds = 0
  do while (ibonds.lt.nbonds .and. .not.ok2)
    ibonds = ibonds + 1
    ok2 = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
  enddo
endif
ok = ok1 .or. ok2
if (.not.ok) call error ('psf_imp', 'Central atom is not found', faterr)
iat1 = psf_atomtype(iat1)
iat2 = psf_atomtype(iat2)
iat3 = psf_atomtype(iat3)
iat4 = psf_atomtype(iat4)
itype = val(iat1)
jtype = val(iat2)
ktype = val(iat3)
ltype = val(iat4)
j = 0
ok = .false.
do while (j.lt.ndeformt .and. .not.ok)
  j = j + 1
  ok = (itype.eq.psf_btype(1,j).and.jtype.eq.psf_btype(2,j).and.ktype.eq.psf_btype(3,j).and.ltype.eq.psf_btype(4,j)) .or. &
       (itype.eq.psf_btype(4,j).and.jtype.eq.psf_btype(3,j).and.ktype.eq.psf_btype(2,j).and.ltype.eq.psf_btype(1,j))
enddo
if (ok) then
  deforms(5,ideform) = j
else
  ndeformt = ndeformt + 1
  psf_btype(1,ndeformt) = itype 
  psf_btype(2,ndeformt) = jtype 
  psf_btype(3,ndeformt) = ktype
  psf_btype(4,ndeformt) = ltype
  deforms(5,ideform) = ndeformt
endif

return
end subroutine
