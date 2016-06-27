subroutine psf_dih(itort,iat1,iat2,iat3,iat4,val,psf_btype)

use errormod
use explatmod
use charmmmod

implicit none
integer itort,iat1,iat2,iat3,iat4
integer val(natt),psf_btype(4,ntorts)
! local variables
integer ib1,ib2,ibonds,itype,jtype,ktype,ltype,j
logical*1 ok

if (iat1.le.0 .or. iat2.le.0 .or. iat3.le.0 .or. iat4.le.0 .or. iat1.gt.natt .or. iat2.gt.natt .or. iat3.gt.natt .or. iat4.gt.natt &
    .or. iat1.eq.iat2 .or. iat1.eq.iat3 .or. iat1.eq.iat4 .or. iat2.eq.iat3 .or. iat2.eq.iat4 .or. iat3.eq.iat4) call error ('psf_dih', 'Wrong dihedral angle', faterr)
! check there is not duplicate dihedral angles
ok = .false.
j = 0
do while (j.lt.itort .and. .not.ok)
  j = j + 1
  ok = (iat1.eq.torts(1,j).and.iat2.eq.torts(2,j).and.iat3.eq.torts(3,j).and.iat4.eq.torts(4,j)) .or. &
       (iat1.eq.torts(4,j).and.iat2.eq.torts(3,j).and.iat3.eq.torts(2,j).and.iat4.eq.torts(1,j))
enddo
if (ok) call error ('psf_dih', 'Dihedral angles are duplicated', faterr)
itort = itort + 1
torts(1,itort) = iat1 ! terminal atom
torts(2,itort) = iat2 
torts(3,itort) = iat3
torts(4,itort) = iat4 ! terminal atom
! determine which bonds make up each torsion
ib1 = torts(1,itort)
ib2 = torts(2,itort)
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
if (.not.ok) call error ('psf_dih', 'torsion: first bond not found', faterr)
ib1 = torts(2,itort)
ib2 = torts(3,itort)
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
if (.not.ok) call error ('psf_dih', 'torsion: second bond not found', faterr)
ib1 = torts(3,itort)
ib2 = torts(4,itort)
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
if (.not.ok) call error ('psf_dih', 'torsion: third bond not found', faterr)
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
do while (j.lt.ntortt .and. .not.ok)
  j = j + 1
  ok = (itype.eq.psf_btype(1,j).and.jtype.eq.psf_btype(2,j).and.ktype.eq.psf_btype(3,j).and.ltype.eq.psf_btype(4,j)) .or. &
       (itype.eq.psf_btype(4,j).and.jtype.eq.psf_btype(3,j).and.ktype.eq.psf_btype(2,j).and.ltype.eq.psf_btype(1,j))
enddo
if (ok) then
  torts(5,itort) = j
else
  ntortt = ntortt + 1
  psf_btype(1,ntortt) = itype
  psf_btype(2,ntortt) = jtype
  psf_btype(3,ntortt) = ktype
  psf_btype(4,ntortt) = ltype
  torts(5,itort) = ntortt
endif

return
end subroutine
