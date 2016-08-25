subroutine psf_bend(ibends,iat1,iat2,iat3,val,psf_btype)

use errormod
use explatmod
use charmmmod

implicit none
integer ibends,iat1,iat2,iat3
integer val(natt),psf_btype(3,nbends)
! local variables
integer ib1,ib2,ibonds,itype,jtype,ktype,j
logical*1 ok

if (iat1.le.0 .or. iat2.le.0 .or. iat3.le.0 .or. iat1.gt.natt .or. iat2.gt.natt .or. iat3.gt.natt .or. & 
    iat1.eq.iat2 .or. iat1.eq.iat3 .or. iat2.eq.iat3) call error ('psf_bend', 'Wrong bond angle', faterr)
ibends = ibends + 1
bends(1,ibends) = iat1
bends(2,ibends) = iat2 ! central atom
bends(3,ibends) = iat3
! determine which bonds make up each bend 
ib1 = iat1
ib2 = iat2
if (ib2.lt.ib1) then
  j = ib1
  ib1 = ib2
  ib2 = j
endif
ok = .false.
ibonds = 0
do while (ibonds.lt.nbonds .and. .not.ok) 
  ibonds = ibonds + 1
  ok = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
enddo
if (.not.ok) call error ('psf_bend', 'bend: first bond not found', faterr)
ib1 = iat2
ib2 = iat3
if (ib2.lt.ib1) then
  j = ib1
  ib1 = ib2
  ib2 = j
endif
ok = .false.
ibonds = 0
do while (ibonds.lt.nbonds .and. .not.ok) 
  ibonds = ibonds + 1
  ok = ib1.eq.bonds(1,ibonds) .and. ib2.eq.bonds(2,ibonds)
enddo
if (.not.ok) call error ('psf_bend', 'bend: second bond not found', faterr)
iat1 = psf_atomtype(iat1)
iat2 = psf_atomtype(iat2)
iat3 = psf_atomtype(iat3)
itype = val(iat1)
jtype = val(iat2)
ktype = val(iat3)
j = 0
ok = .false.
do while (j.lt.nbendt .and. .not.ok)
  j = j + 1
  ok = (itype.eq.psf_btype(1,j).and.jtype.eq.psf_btype(2,j).and.ktype.eq.psf_btype(3,j)) .or. &
       (itype.eq.psf_btype(3,j).and.jtype.eq.psf_btype(2,j).and.ktype.eq.psf_btype(1,j))
enddo
if (ok) then
  bends(4,ibends) = j
else
  nbendt = nbendt + 1
  psf_btype(1,nbendt) = itype
  psf_btype(2,nbendt) = jtype
  psf_btype(3,nbendt) = ktype
  bends(4,ibends) = nbendt
endif

return
end subroutine
