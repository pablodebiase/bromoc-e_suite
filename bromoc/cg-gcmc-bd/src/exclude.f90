subroutine exclude

use explatmod
implicit none
! local variables
integer iat,nex,ibond,ix,iex
integer jat,kat,ibend,itmp,i14a,i14b,ibond1,ibond2
logical*1 ok,HasBeenUsed,Found14

allocate (listmex(natt),listex(natt,natt),listm14(natt),list14(natt,natt))
listmex = 0
listm14 = 0

! MAKE A LIST FOR 1-2 AND 1-3 TERMS
do iat = 1, natt - 1  ! search tables for each atom
  nex = 0  ! initialized number excluded to 0
  ! determine all bond interactions with iat
  ! NOTE: bonds(2,) > bonds(1,)
  do ibond = 1, nbonds
    if (bonds(1,ibond) .eq. iat ) then ! exclude
      ix = bonds(2,ibond)
      ok = .false.
      iex = 0
      do while (iex.lt.nex .and. .not.ok)
        iex = iex + 1
        ok = listex(iex,iat).eq.ix
      end do
      if (.not.ok) then
        nex = nex + 1
        listex(nex,iat) = ix
      end if
    end if
  end do ! next ibond 
  ! determine all bend interactions with iat as first atom or last atom
  do ibend = 1, nbends
    if (bends(1,ibend) .eq. iat ) then ! exclude
      ix = bends(3,ibend)
      if (ix .gt. iat) then
        ok = .false.
        iex = 0
        do while (iex.lt.nex .and. .not.ok)
          iex = iex + 1
          ok = listex(iex,iat).eq.ix
        end do
        if (.not.ok) then
          nex = nex + 1
          listex(nex,iat) = ix
        end if
      end if
    else if (bends(3,ibend) .eq. iat) then
      ix = bends(1,ibend)
      if (ix .gt. iat) then
        ok = .false.
        iex = 0
        do while (iex.lt.nex .and. .not.ok)
          iex = iex + 1
          ok = listex(iex,iat).eq.ix
        end do
        if (.not.ok) then
          nex = nex + 1
          listex(nex,iat) = ix
        end if
      end if
    end if
  end do ! next ibend
  listmex(iat) = nex
end do ! next iat

! MAKE A LIST FOR 1-4 TERMS
do ibend = 1, nbends
  iat = bends(1,ibend)
  jat = bends(2,ibend)
  kat = bends(3,ibend)
  do ibond = 1, nbonds
    Found14 = .false.
    ibond1 = bonds(1,ibond)
    ibond2 = bonds(2,ibond)
    if (ibond1.eq.kat .and. ibond2.ne.jat) then
      i14a = iat
      i14b = ibond2
      Found14 = .true.
    end if
    if (ibond2.eq.kat .and. ibond1.ne.jat) then
      i14a = iat
      i14b = ibond1
      Found14 = .true.
    end if
    if (ibond1.eq.iat .and. ibond2.ne.jat) then
      i14a = kat
      i14b = ibond2
      Found14 = .true.
    end if
    if (ibond2.eq.iat .and. ibond1.ne.jat) then
      i14a = kat
      i14b = ibond1
      Found14 = .true.
    end if
    if (Found14 .and. i14a.eq.i14b) Found14 = .false.
    if (Found14) then
      if (i14a.gt.i14b) then
        itmp = i14a
        i14a = i14b
        i14b = itmp
      end if
      HasBeenUsed = .false. 
      iex = 0
      do while (iex.lt.listm14(i14a) .and. .not.HasBeenUsed)
        iex = iex + 1
        HasBeenUsed = list14(iex,i14a).eq.i14b
      end do
      ! matters for cyclic molecules 
      iex = 0
      do while (iex.lt.listmex(i14a) .and. .not.HasBeenUsed)
        iex = iex + 1
        HasBeenUsed = listex(iex,i14a).eq.i14b
      end do

      if (.not.HasBeenUsed) then
        listm14(i14a) = listm14(i14a) + 1
        list14(listm14(i14a),i14a) = i14b
      end if
    end if 
  end do ! next ibond
end do ! next ibend

return
end subroutine
