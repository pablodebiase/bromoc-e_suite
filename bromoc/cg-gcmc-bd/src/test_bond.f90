logical*1 function test_bond(itype,jtype)

use charmmmod
implicit none
integer itype, jtype
! logical variables
integer i1, i2, k

if (itype.gt.jtype) then
  i1 = jtype
  i2 = itype 
else
  i1 = itype
  i2 = jtype
endif
k = 0
test_bond = .false.
do while (k.lt.chmmbond .and. .not.test_bond)
  k = k + 1
  test_bond = i1.eq.charmm_btype(1,k) .and. i2.eq.charmm_btype(2,k)
enddo

end function
