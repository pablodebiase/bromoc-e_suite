subroutine psf_bond(ibonds,iat1,iat2,val,psf_btype)

use errormod
use explatmod
use charmmmod

implicit none
integer ibonds,iat1,iat2
integer val(natt),psf_btype(2,nbonds)
! local variables
integer itype,jtype,j
logical*1 ok


if (iat1.le.0 .or. iat2.le.0 .or. iat1.gt.natt .or. iat2.gt.natt .or. iat1.eq.iat2) call error ('psf_bond', 'Wrong bond atoms', faterr)
ibonds = ibonds + 1
if (iat1.lt.iat2) then
  bonds(1,ibonds) = iat1
  bonds(2,ibonds) = iat2
else
  bonds(1,ibonds) = iat2
  bonds(2,ibonds) = iat1
endif
iat1 = psf_atomtype(iat1)
iat2 = psf_atomtype(iat2)
itype = val(iat1)
jtype = val(iat2)
if (itype.gt.jtype) then
  j = jtype
  jtype = itype
  itype = j
endif
j = 0
ok = .false.
do while (j.lt.nbondt .and. .not.ok)
  j = j + 1
  ok = itype.eq.psf_btype(1,j) .and. jtype.eq.psf_btype(2,j)
enddo
if (ok) then
  bonds(3,ibonds) = j
else
  nbondt = nbondt + 1
  psf_btype(1,nbondt) = itype
  psf_btype(2,nbondt) = jtype
  bonds(3,ibonds) = nbondt
endif

return
end subroutine
