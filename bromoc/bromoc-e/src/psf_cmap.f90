subroutine psf_cmap(i,j,iat1,iat2,iat3,iat4,psf_mass)

use errormod
use explatmod
use charmmmod

implicit none
integer i,j,iat1,iat2,iat3,iat4,itype,jtype,ktype,ltype
real psf_mass(natt)
! local variables
integer itort
real, parameter :: tol=1.0e-2
logical*1 ok

itype = psf_atomtype(iat1)
jtype = psf_atomtype(iat2)
ktype = psf_atomtype(iat3)
ltype = psf_atomtype(iat4)
if (abs(psf_mass(itype)-12.011).lt.tol .and. abs(psf_mass(ltype)-12.011).lt.tol) then ! theta angle (C-N-C_alpha-C)
  if ((abs(psf_mass(jtype)-14.007).lt.tol.and.abs(psf_mass(ktype)-12.011).lt.tol) .or. &
      (abs(psf_mass(ktype)-12.011).lt.tol.and.abs(psf_mass(jtype)-14.007).lt.tol)) then
    j = 1
  else
    call error ('psf_cmap', 'Wrong atoms in cross-term section', faterr)
  endif
else if (abs(psf_mass(itype)-14.007).lt.tol .and. abs(psf_mass(ltype)-14.007).lt.tol) then ! psi angle (N-C_alpha-C-N)
  if (abs(psf_mass(jtype)-12.011).lt.tol.and.abs(psf_mass(ktype)-12.011).lt.tol) then
    j = 2
  else
    call error ('psf_cmap', 'Wrong atoms in cross-term section', faterr)
  endif
else
  call error ('psf_cmap', 'Wrong atoms in cross-term section', faterr)
endif
! check dihedral angle
itort = 0
ok = .false.
do while (itort.lt.ntorts .and. .not.ok)
  itort = itort + 1
  ok = (iat1.eq.torts(1,itort).and.iat2.eq.torts(2,itort).and.iat3.eq.torts(3,itort).and.iat4.eq.torts(4,itort)) .or. &
       (iat4.eq.torts(1,itort).and.iat3.eq.torts(2,itort).and.iat2.eq.torts(3,itort).and.iat1.eq.torts(4,itort))
enddo
if (.not.ok) call error ('psf_cmap', 'CMAP: dihedral angle not found', faterr) 
! assignament
cmaps(j,i) = itort
if (j.eq.1) then
  lthetacmap(itort) = i
else
  lpsicmap(itort) = i
endif

return
end subroutine
