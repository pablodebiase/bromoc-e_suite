subroutine readpsf(Qprint)

use errormod
use stdiomod
use charmmmod
use grandmod
use charfuncmod
use explatmod
use constamod

implicit none
logical*1 Qprint
! local variables
integer :: dummyi, i, j, k, l, m, IDat, intert, intert2, ncharge, ilines, irest, ibonds, ibends, iubs, itort, itort1, itort2, ideform
integer :: iat1, iat2, iat3, iat4, iat5, iat6, iat7, iat8, iat9
integer :: itype, jtype, ktype, ltype, mtype, nntype, otype, ptype, icmap
integer :: n, nmax
integer, allocatable :: psf_nq(:), val(:), psf_btype(:,:), psf_bendtyp(:), psf_ubtp(:,:), cmaptype(:)
real, parameter :: tol=1.0e-3, cte=2**(-5.0/6.0)
real :: atomchar, atommass, Qnet
real, allocatable :: psf_qat(:,:), psf_mass(:)
character*7, allocatable :: psf_non_labels(:)
character :: dummyc*144, atname*7, intg*2, word*5, wrtline*2048
logical*1 :: dobond, doang, dodih, dodef, docmap
logical*1 :: ok, ok1, ok2

! *** OBTENTION OF DIMENSIONS
nbonds   = 0
nbends   = 0
ntorts   = 0
ndeforms = 0
ncmaps   = 0
Qlbond = .false.
Qlang  = .false.
Qldih  = .false.
Qldef  = .false.
Qlcmap = .false.
dobond = .false.
doang = .false.
dodih = .false.
dodef = .false.
docmap = .false.
! *** read header
read(iunpsf,*)

! *** read title section
read(iunpsf,*)
read(iunpsf,*) dummyi,dummyc
do i = 1,dummyi
  read(iunpsf,*)
enddo

! *** read atoms section
read(iunpsf,*)
read(iunpsf,*) natt,dummyc
if (natt.le.0) call error ('readpsg', 'Wrong psf format', faterr)
do i = 1,natt
  read(iunpsf,*)
enddo

! *** rest of sections
ok = .true.
ok1 = .true.
do while (ok.and.ok1) 
  read(iunpsf,*,iostat=i)
  if (i.eq.0) then
    read(iunpsf,*,iostat=i) dummyi,dummyc
    if (dummyi.lt.0) call error ('readpsf', 'Wrong psf format', faterr)
    ok2 = dummyi.gt.0
  else if (i.gt.0) then
    call error ('readpsf', 'Wrong psf format', faterr)
  endif
  if (i.eq.0) then
    word = trim(adjustl(dummyc)) 
    if (word.eq.'!NBON') then
      if (.not.Qlbond) then 
        Qlbond = ok2
      else
        call error ('readpsf', 'Wrong psf format', faterr)
      endif
    else if (word.eq.'!NTHE') then
      if (.not.Qlang) then
        Qlang  = ok2
      else
        call error ('readpsf', 'Wrong psf format', faterr)
      endif
    else if (word.eq.'!NPHI') then
      if (.not.Qldih) then
        Qldih  = ok2
      else
        call error ('readpsf', 'Wrong psf format', faterr)
      endif
    else if (word.eq.'!NIMP') then
      if (.not.Qldef)  then
        Qldef  = ok2
      else
        call error ('readpsf', 'Wrong psf format', faterr)
      endif
    else if (word.eq.'!NCRT') then
      if (.not.Qlcmap) then
        Qlcmap = ok2
      else
        call error ('readpsf', 'Wrong psf format', faterr)
      endif
    else 
      ok1 = .false.
    endif
  else if (i.gt.0) then
   call error ('readpsf', 'Wrong psf format', faterr)
  else
   ok = .false. ! EOF
  endif
  if (Qlbond .and. .not.dobond) then
    nbonds = dummyi
    ilines = nbonds/4 
    if (4*ilines.lt.nbonds) ilines = ilines + 1
    do i = 1,ilines
      read(iunpsf,*)
    enddo
    dobond = .true.
  endif
  if (Qlang .and. .not.doang) then
    nbends = dummyi
    ilines = nbends/3 
    if (3*ilines.lt.nbends) ilines = ilines + 1
    do i = 1,ilines
      read(iunpsf,*)
    enddo
    doang = .true.
  endif
  if (Qldih .and. .not.dodih) then
    ntorts = dummyi
    ilines = ntorts/2 
    if (2*ilines.lt.ntorts) ilines = ilines + 1
    do i = 1,ilines
      read(iunpsf,*)
    enddo
    dodih = .true.
  endif
  if (Qldef .and. .not.dodef) then
    ndeforms = dummyi
    ilines = ndeforms/2 
    if (2*ilines.lt.ndeforms) ilines = ilines + 1
    do i = 1,ilines
      read(iunpsf,*)
    enddo
    dodef = .true.
  endif
  if (Qlcmap .and. .not.docmap) then
    ncmaps = dummyi
    do i = 1,ncmaps
      read(iunpsf,*)
    enddo
    docmap = .true.
  endif
enddo
! *** checking format
if (.not.Qlbond .and. (Qlang.or.Qldih.or.Qldef.or.Qlcmap)) call error ('readpsf', 'Wrong psf format', faterr)
if (.not.Qlang .and. (Qldih.or.Qldef.or.Qlcmap)) call error ('readpsf', 'Wrong psf format', faterr)
if (.not.Qldih .and. Qlcmap) call error ('readpsf', 'Wrong psf format', faterr)
! *** checking compatibilty with CHARMM paprameters file
if (Qlbond .and. .not.Qchmmbond) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr)
if (Qlang .and. .not.Qchmmang) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr)
if (Qldih .and. .not.Qchmmdih) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr)
if (Qldef .and. .not.Qchmmimp) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr)
if (Qlcmap .and. .not.Qchmmcmap) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr) 


! *** OBTENTION OF PARAMETERS
rewind(unit=iunpsf)

! *** read header
read(iunpsf,*)

! *** read title section
read(iunpsf,*)
read(iunpsf,*) dummyi,dummyc
do i = 1,dummyi
  read(iunpsf,*)
enddo

! *** read atoms section
read(iunpsf,*)
read(iunpsf,*) natt,dummyc
allocate(psf_nq(natt),psf_qat(natt,natt),psf_mass(natt),val(natt),psf_atomtype(natt),psf_atomtype2(natt),psf_non_labels(natt))
! obtention of maxtypes and some arrays for more assignaments
Qnet = 0.0
maxtypes = 0
psf_nq = 0
do i = 1,natt
  read(iunpsf,*) IDat,dummyc,dummyi,dummyc,dummyc,atname,atomchar,atommass,dummyi
  if (IDat.le.0 .or. IDat.gt.natt) call error ('readpsf', 'Wrong atom ID', faterr)
  Qnet = Qnet + atomchar
  j = 0
  ok = .false.
  ok2 = .false.
  do while (j.lt.maxtypes .and. .not.ok)
    j = j + 1
    ok = lcase(atname).eq.lcase(psf_non_labels(j)) 
    if (ok) then
      if (abs(atommass-psf_mass(j)).ge.tol) call error ('readpsf', 'Wrong atom mass', faterr)
      k = 0
      do while (k.lt.psf_nq(j) .and..not.ok2)
        k = k + 1
        ok2 = abs(atomchar-psf_qat(j,k)).lt.tol
      enddo
    endif
  enddo
  if (.not.ok) then
    maxtypes = maxtypes + 1
    psf_non_labels(maxtypes) = atname
    psf_mass(maxtypes) = atommass
    call fatnam(charmm_label,1,chmmntype,atname,itype) ! atom type
    val(maxtypes) = itype
    psf_atomtype(i) = maxtypes
    psf_nq(maxtypes) = psf_nq(maxtypes) + 1
    psf_atomtype2(i) = psf_nq(maxtypes)
    psf_qat(maxtypes,psf_nq(maxtypes)) = atomchar
  else if (ok .and. .not.ok2) then
    psf_atomtype(i) = j
    psf_nq(j) = psf_nq(j) + 1
    psf_atomtype2(i) = psf_nq(j)
    psf_qat(j,psf_nq(j)) = atomchar
  else
    psf_atomtype(i) = j
    psf_atomtype2(i) = k
  endif
enddo ! next i
! obtention of maxcharges
maxcharges = 0
do itype = 1, maxtypes
  do i = 1, psf_nq(itype)
    maxcharges = maxcharges + 1
  enddo
enddo
! allocations
allocate(typen(maxtypes,maxtypes),non_labels(maxtypes),sdat(maxtypes),qat(maxcharges), &
         atom_labels(maxcharges),non_of_charge(maxcharges),psf_charge(2,maxcharges))
! assign typen initially
intert = 0 
do itype = 1, maxtypes
  do jtype = itype, maxtypes
    intert = intert + 1
    typen(itype,jtype) = intert
    typen(jtype,itype) = intert
  enddo
enddo
allocate(nonbonded(4,intert))
! non-bonded forces field
ncharge = 0
do itype = 1, maxtypes
  i = val(itype) 
  non_labels(itype) = charmm_label(i)
  intert = typen(itype,itype)
  j = charmm_typen(i,i) ! nonbonded pair
  do k = 1, 4
    nonbonded(k,intert) = charmm_nonbonded(k,j)
  enddo
  sdat(itype) = scldiff*kbt/(6.0*pi*viscwat*cte*nonbonded(2,intert)) 
  do jtype = 1, psf_nq(itype)
    write(intg,'(i2)') jtype
    ncharge = ncharge + 1
    atom_labels(ncharge) = trim(adjustl(non_labels(itype)))//trim(adjustl(intg))
    non_of_charge(ncharge) = itype
    qat(ncharge) = psf_qat(itype,jtype)
    psf_charge(1,ncharge) = itype
    psf_charge(2,ncharge) = jtype 
  enddo
enddo 
do itype = 1, maxtypes-1
  i = val(itype)
  do jtype = itype+1, maxtypes
    j = val(jtype)
    intert = typen(itype,jtype)
    intert2  = charmm_typen(i,j)
    do k = 1, 4
      nonbonded(k,intert) = charmm_nonbonded(k,intert2)
    enddo
  enddo
enddo
deallocate(psf_non_labels,psf_nq,psf_qat)

! *** read bond section
if (Qlbond) then
  read(iunpsf,*)
  read(iunpsf,*) nbonds,dummyc
  ! allocations
  allocate (bonds(3,nbonds),psf_btype(2,nbonds))
  ! obtention of nbondt, bonds and an array for more assignaments
  ilines = nbonds/4 ! lines to read
  nbondt = 0
  ibonds = 0
  do i = 1,ilines
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8 
    ! first bond
    call psf_bond(ibonds,iat1,iat2,val,psf_btype)
    ! second bond
    call psf_bond(ibonds,iat3,iat4,val,psf_btype)
    ! third bond
    call psf_bond(ibonds,iat5,iat6,val,psf_btype)
    ! fourth bond
    call psf_bond(ibonds,iat7,iat8,val,psf_btype)
  enddo ! next i
  irest = nbonds - ibonds ! rest of bonds to read in the last line
  if (irest.eq.1) then
    read(iunpsf,*) iat1,iat2
    call psf_bond(ibonds,iat1,iat2,val,psf_btype)
  else if (irest.eq.2) then
    read(iunpsf,*) iat1,iat2,iat3,iat4
    ! first bond
    call psf_bond(ibonds,iat1,iat2,val,psf_btype) 
    ! second bond
    call psf_bond(ibonds,iat3,iat4,val,psf_btype)
  else if (irest.eq.3) then
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6
    ! first bond
    call psf_bond(ibonds,iat1,iat2,val,psf_btype)
    ! second bond
    call psf_bond(ibonds,iat3,iat4,val,psf_btype)
    ! third bond
    call psf_bond(ibonds,iat5,iat6,val,psf_btype)
  endif
  ! allocations
  if (nbondt.gt.chmmbond) call error ('readpsf', 'Wrong number of bond types', faterr)
  allocate(stretch(3,nbondt))
  ! bond terms
  do i = 1, nbondt
    j = 0
    ok = .false.
    do while (j.lt.chmmbond .and. .not.ok)
      j = j + 1
      ok = psf_btype(1,i).eq.charmm_btype(1,j) .and. psf_btype(2,i).eq.charmm_btype(2,j)
    enddo
    if (ok) then
      stretch(1,i) = charmm_bond(1,j) ! Kb
      stretch(2,i) = charmm_bond(2,j) ! b0
      stretch(3,i) = charmm_bond(2,j) ! b0 (const.)
    else
      call error ('readpsf', 'Wrong bond term', faterr)
    endif
  enddo
  ! deallocations
  deallocate(psf_btype)
endif ! Qlbond

! *** read bond angle section
if (Qlang) then
  read(iunpsf,*)
  read(iunpsf,*) nbends,dummyc
  ! allocations
  allocate (bends(4,nbends),psf_btype(3,nbends))
  ! obtention of nbendt, bends and an array for more assignaments
  ilines = nbends/3 ! lines to read
  nbendt = 0
  ibends = 0
  do i = 1,ilines
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8,iat9 
    ! first bond angle
    call psf_bend(ibends,iat1,iat2,iat3,val,psf_btype) 
    ! second bond angle
    call psf_bend(ibends,iat4,iat5,iat6,val,psf_btype)
    ! third bond angle
    call psf_bend(ibends,iat7,iat8,iat9,val,psf_btype)
  enddo ! next i
  irest = nbends - ibends ! rest of bond angles to read in the last line
  if (irest.eq.1) then
    read(iunpsf,*) iat1,iat2,iat3
    call psf_bend(ibends,iat1,iat2,iat3,val,psf_btype)
  else if (irest.eq.2) then
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6
    ! first bond angle
    call psf_bend(ibends,iat1,iat2,iat3,val,psf_btype)
    ! second bond angle 
    call psf_bend(ibends,iat4,iat5,iat6,val,psf_btype)
  endif
  ! allocations
  if (nbendt.gt.chmmang) call error ('readpsf', 'Wrong number of angle types', faterr)
  allocate(bend(2,nbendt),psf_bendtyp(nbendt))
  ! bond angle terms
  do i = 1, nbendt
    j = 0
    ok = .false.
    do while (j.lt.chmmang .and. .not.ok)
      j = j + 1
      ok = (psf_btype(1,i).eq.charmm_atype(1,j).and.psf_btype(2,i).eq.charmm_atype(2,j).and.psf_btype(3,i).eq.charmm_atype(3,j)) .or. &
           (psf_btype(1,i).eq.charmm_atype(3,j).and.psf_btype(2,i).eq.charmm_atype(2,j).and.psf_btype(3,i).eq.charmm_atype(1,j))
    enddo
    if (ok) then
      bend(1,i) = charmm_ang(1,j) ! Ktheta
      bend(2,i) = charmm_ang(2,j) ! Theta0
      psf_bendtyp(i) = j  
    else
      call error ('readpsf', 'Wrong bond angle term', faterr)
    endif
  enddo
  ! obtention of UB terms
  nubs = 0
  do ibends = 1, nbends
    i = bends(4,ibends)
    j = psf_bendtyp(i)
    if (charmm_lub(j).gt.0) nubs = nubs + 1
  enddo
  Qlubs = nubs.gt.0
  if (Qlubs .and. .not.Qchmmub) call error ('readpsf', 'PSF and CHARMM parameter files are not compatibles', faterr)
  if (Qlubs) then
    ! allocations
    allocate(ubs(3,nubs),psf_ubtp(3,nubs))
    ! obtention of nubt, ubs and an array for more assignaments
    i = 0
    nubt = 0
    do ibends = 1, nbends
      iat1 = bends(1,ibends)
      iat3 = bends(3,ibends)
      itype = bends(4,ibends)
      j = psf_bendtyp(itype)  
      if (charmm_lub(j).gt.0) then
        i = i + 1
        ubs(1,i) = iat1
        ubs(2,i) = iat3
        iat1 = psf_atomtype(iat1)
        iat3 = psf_atomtype(iat3)
        itype = val(iat1)
        jtype = val(iat3)
        k = 0
        ok = .false.
        do while (k.lt.nubt .and. .not.ok)
          k = k + 1
          ok = (itype.eq.psf_ubtp(1,k).and.jtype.eq.psf_ubtp(2,k)) .or. &
               (itype.eq.psf_ubtp(2,k).and.jtype.eq.psf_ubtp(1,k)) 
        enddo
        if (ok) then  
          ubs(3,i) = k
        else
          nubt = nubt + 1
          ubs(3,i) = nubt
          psf_ubtp(1,nubt) = itype
          psf_ubtp(2,nubt) = jtype
          psf_ubtp(3,nubt) = charmm_lub(j)
        endif
      endif
    enddo
    ! allocations
    if (nubt.gt.chmmub) call error ('readpsf', 'Wrong number of UB types', faterr)
    allocate(ubt(2,nubt))
    ! UB terms
    do iubs = 1, nubt
      i = psf_ubtp(3,iubs)
      ubt(1,iubs) = charmm_ub(1,i) ! Kub
      ubt(2,iubs) = charmm_ub(2,i) ! S0
    enddo
    ! deallocations
    deallocate(psf_ubtp)
  endif ! Qlubs
  ! deallocations
  deallocate(psf_btype,psf_bendtyp)
endif ! Qlang

! *** read dihedral angle section
if (Qldih) then
  read(iunpsf,*)
  read(iunpsf,*) ntorts,dummyc
  ! allocations
  allocate (torts(5,ntorts),psf_btype(4,ntorts))
  ! obtention of ntortt, torts and an array for more assignaments
  ilines = ntorts/2 ! lines to read
  ntortt = 0
  itort = 0
  do i = 1,ilines
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8
    ! first dihedral angle
    call psf_dih(itort,iat1,iat2,iat3,iat4,val,psf_btype)
    ! second dihedral angle
    call psf_dih(itort,iat5,iat6,iat7,iat8,val,psf_btype)
  enddo ! next i
  irest = ntorts - itort ! rest of dihedral angles to read in the last line
  if (irest.eq.1) then
    read(iunpsf,*) iat1,iat2,iat3,iat4
    call psf_dih(itort,iat1,iat2,iat3,iat4,val,psf_btype)
  endif
  ! allocations
  allocate(dih(2*charmm_nmax,ntortt),ndih(charmm_nmax,ntortt),nprms(ntortt))
  ! dihedral angle terms
  do i = 1, ntortt
    j = 0
    ok = .false.
    do while (j.lt.chmmdih .and. .not.ok)
      j = j + 1
      if (charmm_dtype(1,j).eq.0 .and. charmm_dtype(4,j).eq.0) then ! X-A-B-X
        ok = (psf_btype(2,i).eq.charmm_dtype(2,j).and.psf_btype(3,i).eq.charmm_dtype(3,j)) .or. & 
             (psf_btype(3,i).eq.charmm_dtype(2,j).and.psf_btype(2,i).eq.charmm_dtype(3,j))
      else if (charmm_dtype(1,j).eq.0) then ! X-A-B-C
        ok = (psf_btype(2,i).eq.charmm_dtype(2,j).and.psf_btype(3,i).eq.charmm_dtype(3,j).and.psf_btype(4,i).eq.charmm_dtype(4,j)) .or. &
             (psf_btype(3,i).eq.charmm_dtype(2,j).and.psf_btype(2,i).eq.charmm_dtype(3,j).and.psf_btype(1,i).eq.charmm_dtype(4,j))
      else if (charmm_dtype(4,j).eq.0) then ! A-B-C-X
        ok = (psf_btype(1,i).eq.charmm_dtype(1,j).and.psf_btype(2,i).eq.charmm_dtype(2,j).and.psf_btype(3,i).eq.charmm_dtype(3,j)) .or. &
             (psf_btype(4,i).eq.charmm_dtype(1,j).and.psf_btype(3,i).eq.charmm_dtype(2,j).and.psf_btype(2,i).eq.charmm_dtype(3,j))
      else ! A-B-C-D
        ok = (psf_btype(1,i).eq.charmm_dtype(1,j).and.psf_btype(2,i).eq.charmm_dtype(2,j).and.psf_btype(3,i).eq.charmm_dtype(3,j).and.psf_btype(4,i).eq.charmm_dtype(4,j)) .or. &
             (psf_btype(4,i).eq.charmm_dtype(1,j).and.psf_btype(3,i).eq.charmm_dtype(2,j).and.psf_btype(2,i).eq.charmm_dtype(3,j).and.psf_btype(1,i).eq.charmm_dtype(4,j))
      endif
    enddo
    if (ok) then
      k = charmm_nprms(j) 
      nprms(i) = k ! number of terms (multiple dihedral angles)
      do l = 1, k
        m = (l-1)*2
        dih(m+1,i) = charmm_dih(m+1,j) ! Kchi
        dih(m+2,i) = charmm_dih(m+2,j) ! delta
        ndih(l,i) = charmm_ndih(l,j) ! n
      enddo  
    else
      call error ('readpsf', 'Wrong dihedral angle term', faterr)
    endif
  enddo
  ! deallocations
  deallocate(psf_btype)
endif ! Qldih

! *** read improper angle section
if (Qldef) then
  read(iunpsf,*)
  read(iunpsf,*) ndeforms,dummyc
  ! allocations
  allocate (deforms(5,ndeforms),psf_btype(4,ndeforms))
  ! obtention of ndeformt, deforms and an array for more assignaments
  ilines = ndeforms/2 ! lines to read
  ndeformt = 0
  ideform = 0
  do i = 1,ilines
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8
    ! first improper angle
    call psf_imp(ideform,iat1,iat2,iat3,iat4,val,psf_btype)
    ! second improper angle
    call psf_imp(ideform,iat5,iat6,iat7,iat8,val,psf_btype)
  enddo ! next i
  irest = ndeforms - ideform ! rest of improper angles to read in the last line
  if (irest.eq.1) then
    read(iunpsf,*) iat1,iat2,iat3,iat4
    call psf_imp(ideform,iat1,iat2,iat3,iat4,val,psf_btype)
  endif
  ! allocations
  allocate(deform(2,ndeformt))
  ! improper angle terms
  do i = 1, ndeformt
    j = 0
    ok = .false.
    do while (j.lt.chmmimp .and. .not.ok)
      j = j + 1
      if (charmm_itype(1,j).ne.0 .and. charmm_itype(4,j).ne.0) then ! A-B-C-D; A-X-X-B
        if (charmm_itype(2,j).eq.0 .and. charmm_itype(3,j).eq.0) then ! A-X-X-B
          ok = (psf_btype(1,i).eq.charmm_itype(1,j).and.psf_btype(4,i).eq.charmm_itype(4,j)) .or. &
               (psf_btype(4,i).eq.charmm_itype(1,j).and.psf_btype(1,i).eq.charmm_itype(4,j))
        else ! A-B-C-D
          ok = (psf_btype(1,i).eq.charmm_itype(1,j).and.psf_btype(2,i).eq.charmm_itype(2,j).and.psf_btype(3,i).eq.charmm_itype(3,j).and.psf_btype(4,i).eq.charmm_itype(4,j)) .or. &
               (psf_btype(4,i).eq.charmm_itype(1,j).and.psf_btype(3,i).eq.charmm_itype(2,j).and.psf_btype(2,i).eq.charmm_itype(3,j).and.psf_btype(1,i).eq.charmm_itype(4,j))
        endif
      else if (charmm_itype(2,j).ne.0 .and. charmm_itype(3,j).ne.0) then ! X-A-B-C; X-A-B-X
        if (charmm_itype(1,j).eq.0 .and. charmm_itype(4,j).ne.0) then ! X-A-B-C
          ok = (psf_btype(2,i).eq.charmm_itype(2,j).and.psf_btype(3,i).eq.charmm_itype(3,j).and.psf_btype(4,i).eq.charmm_itype(4,j)) .or. &
               (psf_btype(3,i).eq.charmm_itype(2,j).and.psf_btype(2,i).eq.charmm_itype(3,j).and.psf_btype(1,i).eq.charmm_itype(4,j))
        else if (charmm_itype(1,j).eq.0 .and. charmm_itype(4,j).eq.0) then ! X-A-B-X
          ok = (psf_btype(2,i).eq.charmm_itype(2,j).and.psf_btype(3,i).eq.charmm_itype(3,j)) .or. &
               (psf_btype(3,i).eq.charmm_itype(2,j).and.psf_btype(2,i).eq.charmm_itype(3,j))
        endif
      else ! X-X-A-B
        ok = (psf_btype(3,i).eq.charmm_itype(3,j).and.psf_btype(4,i).eq.charmm_itype(4,j)) .or. &
             (psf_btype(2,i).eq.charmm_itype(3,j).and.psf_btype(1,i).eq.charmm_itype(4,j))
      endif
    enddo
    if (ok) then
      deform(1,i) = charmm_imp(1,j) ! Kpsi
      deform(2,i) = charmm_imp(2,j) ! psi0
    else
      call error ('readpsf', 'Wrong improper angle term', faterr)
    endif
  enddo
  ! deallocations
  deallocate(psf_btype)
endif ! Qldef

! *** read cross-term section
if (Qlcmap) then
  read(iunpsf,*)
  read(iunpsf,*) ncmaps,dummyc
  ! allocations
  allocate(cmaps(3,ncmaps),lthetacmap(ntorts),lpsicmap(ntorts),psf_btype(8,ncmaps))
  allocate(thetacmap(ncmaps),psicmap(ncmaps),attcmap(4,ncmaps),atpcmap(4,ncmaps),nablatcmp(3,4,ncmaps),nablapcmp(3,4,ncmaps))
  ! obtention of ncmap, cmaps, lthetacmap, lpsicmap and an array for more assignaments
  lthetacmap = 0
  lpsicmap = 0
  ncmap = 0
  do i = 1,ncmaps
    read(iunpsf,*) iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8
    if (iat1.le.0 .or. iat2.le.0 .or. iat3.le.0 .or. iat4.le.0 .or. iat1.gt.natt .or. iat2.gt.natt .or. iat3.gt.natt .or.iat4.gt.natt) & 
      call error ('readpsf', 'Wrong first dihedral angle in cross-term section', faterr)
    if (iat5.le.0 .or. iat6.le.0 .or. iat7.le.0 .or. iat8.le.0 .or. iat5.gt.natt .or. iat6.gt.natt .or. iat7.gt.natt .or.iat8.gt.natt) & 
      call error ('readpsf', 'Wrong second dihedral angle in cross-term section', faterr)
    ok = (iat2.eq.iat5.and.iat3.eq.iat6.and.iat4.eq.iat7) .or. &
         (iat1.eq.iat6.and.iat2.eq.iat7.and.iat3.eq.iat8) .or. &
         (iat2.eq.iat8.and.iat3.eq.iat7.and.iat4.eq.iat6) .or. &
         (iat1.eq.iat7.and.iat2.eq.iat6.and.iat3.eq.iat5)
    if (.not.ok) call error ('readpsf', 'Wrong relation between dihedral angles in cross-term section', faterr)
    ! get angles   
    ! first angle
    call psf_cmap(i,j,iat1,iat2,iat3,iat4,psf_mass)
    iat1 = psf_atomtype(iat1)
    iat2 = psf_atomtype(iat2)
    iat3 = psf_atomtype(iat3)
    iat4 = psf_atomtype(iat4)
    itype = val(iat1)
    jtype = val(iat2)
    ktype = val(iat3)
    ltype = val(iat4)
    ! second angle
    call psf_cmap(i,l,iat5,iat6,iat7,iat8,psf_mass)
    if (j.eq.l) call error ('readpsf', 'Wrong relation between dihedral angles in cross-term section', faterr)
    iat5 = psf_atomtype(iat5)
    iat6 = psf_atomtype(iat6)
    iat7 = psf_atomtype(iat7)
    iat8 = psf_atomtype(iat8)
    mtype = val(iat5)
    nntype = val(iat6)
    otype = val(iat7)
    ptype = val(iat8)
    k = 0
    ok = .false.
    do while (k.lt.ncmap .and. .not.ok)
      k = k + 1
      if (j.eq.1) then
        ! theta angle (C-N-C_alpha-C)
        ok1 = (itype.eq.psf_btype(1,k).and.jtype.eq.psf_btype(2,k).and.ktype.eq.psf_btype(3,k).and.ltype.eq.psf_btype(4,k)) .or. &
              (itype.eq.psf_btype(4,k).and.jtype.eq.psf_btype(3,k).and.ktype.eq.psf_btype(2,k).and.ltype.eq.psf_btype(1,k)) 
        ! psi angle (N-C_alpha-C-N)
        ok2 = (mtype.eq.psf_btype(5,k).and.nntype.eq.psf_btype(6,k).and.otype.eq.psf_btype(7,k).and.ptype.eq.psf_btype(8,k)) .or. &
              (mtype.eq.psf_btype(8,k).and.nntype.eq.psf_btype(7,k).and.otype.eq.psf_btype(6,k).and.ptype.eq.psf_btype(5,k))
      else
        ! theta angle (C-N-C_alpha-C)
        ok1 = (mtype.eq.psf_btype(1,k).and.nntype.eq.psf_btype(2,k).and.otype.eq.psf_btype(3,k).and.ptype.eq.psf_btype(4,k)) .or. &
              (mtype.eq.psf_btype(4,k).and.nntype.eq.psf_btype(3,k).and.otype.eq.psf_btype(2,k).and.ptype.eq.psf_btype(1,k))
        ! psi angle (N-C_alpha-C-N)
        ok2 = (itype.eq.psf_btype(5,k).and.jtype.eq.psf_btype(6,k).and.ktype.eq.psf_btype(7,k).and.ltype.eq.psf_btype(8,k)) .or. &
              (itype.eq.psf_btype(8,k).and.jtype.eq.psf_btype(7,k).and.ktype.eq.psf_btype(6,k).and.ltype.eq.psf_btype(5,k))
      endif 
      ok = ok1.and.ok2
    enddo
    if (ok) then
      cmaps(3,i) = k
    else
      ncmap = ncmap + 1
      if (j.eq.1) then
        ! theta angle (C-N-C_alpha-C)
        psf_btype(1,ncmap) = itype
        psf_btype(2,ncmap) = jtype
        psf_btype(3,ncmap) = ktype
        psf_btype(4,ncmap) = ltype
        ! psi angle (N-C_alpha-C-N)
        psf_btype(5,ncmap) = mtype
        psf_btype(6,ncmap) = nntype
        psf_btype(7,ncmap) = otype
        psf_btype(8,ncmap) = ptype
      else
        ! theta angle (C-N-C_alpha-C)
        psf_btype(1,ncmap) = mtype
        psf_btype(2,ncmap) = nntype
        psf_btype(3,ncmap) = otype
        psf_btype(4,ncmap) = ptype
        ! psi angle (N-C_alpha-C-N)
        psf_btype(5,ncmap) = itype
        psf_btype(6,ncmap) = jtype
        psf_btype(7,ncmap) = ktype
        psf_btype(8,ncmap) = ltype
      endif
      cmaps(3,i) = ncmap
    endif
  enddo ! next i
  ! allocations
  if (ncmap.gt.chmmcmap) call error ('readpsf', 'Wrong number of cross-term types', faterr)
  allocate(cmap(ncmap),gscmap(ncmap),cmaptype(ncmap))
  ! CMAP terms
  nmax = 0
  do i = 1, ncmap
    j = 0
    ok = .false.
    do while (j.lt.chmmcmap .and. .not.ok)
      j = j + 1
      ! theta angle (C-N-C_alpha-C)
      ok1 = (psf_btype(1,i).eq.charmm_icmap(1,j).and.psf_btype(2,i).eq.charmm_icmap(2,j).and.psf_btype(3,i).eq.charmm_icmap(3,j).and.psf_btype(4,i).eq.charmm_icmap(4,j)) .or. &
            (psf_btype(1,i).eq.charmm_icmap(4,j).and.psf_btype(2,i).eq.charmm_icmap(3,j).and.psf_btype(3,i).eq.charmm_icmap(2,j).and.psf_btype(4,i).eq.charmm_icmap(1,j))
      ! psi angle (N-C_alpha-C-N)
      ok2 = (psf_btype(5,i).eq.charmm_icmap(5,j).and.psf_btype(6,i).eq.charmm_icmap(6,j).and.psf_btype(7,i).eq.charmm_icmap(7,j).and.psf_btype(8,i).eq.charmm_icmap(8,j)) .or. &
            (psf_btype(5,i).eq.charmm_icmap(8,j).and.psf_btype(6,i).eq.charmm_icmap(7,j).and.psf_btype(7,i).eq.charmm_icmap(6,j).and.psf_btype(8,i).eq.charmm_icmap(5,j))
      ok = ok1.and.ok2   
    enddo
    if (ok) then
      ! grid points and grid spacing
      cmap(i) = charmm_ncmap(j)
      gscmap(i) = charmm_cmap(j)
      n = charmm_ncmap(j)*charmm_ncmap(j)
      if (n.gt.nmax) nmax = n 
      cmaptype(i) = j
    else
      call error ('readpsf', 'Wrong CMAP term', faterr)
    endif
  enddo
  ! allocations
  allocate(fcmap(nmax,ncmap),ftcmap(nmax,ncmap),fpcmap(nmax,ncmap),ftpcmap(nmax,ncmap),ccoef(16,nmax,ncmap))
  do i = 1, ncmap
    n = cmap(i)
    l = cmaptype(i)
    do j = 1, n
      itype = n*(j-1)
      do k = 1, n
        fcmap(itype+k,i) = charmm_fcmap(itype+k,l) 
      enddo
    enddo
  enddo 
  call cmapg() ! Obtain gradients and cross derivatives at the grid points
  call cmapc() ! Obtain the CMAP coefficient for the cubic interpolation
  ! Obtain attcmap and atpcmap
  do itort = 1, ntorts
    i = lthetacmap(itort)
    j = lpsicmap(itort)
    ok1 = i.gt.0
    ok2 = j.gt.0
    if (ok1) then
      attcmap(1,i) = torts(1,itort)
      attcmap(2,i) = torts(2,itort)
      attcmap(3,i) = torts(3,itort)
      attcmap(4,i) = torts(4,itort)
    endif
    if (ok2) then
      atpcmap(1,j) = torts(1,itort)
      atpcmap(2,j) = torts(2,itort)
      atpcmap(3,j) = torts(3,itort)
      atpcmap(4,j) = torts(4,itort)
    endif
  enddo
  deallocate(psf_btype,cmaptype)
endif ! Qlcmap
deallocate(psf_mass,val)

! *** assign additional variables
nch = 1 
allocate(chain(natt),nbondsch(1),fixed(natt),ghost(natt))
do i = 1, natt
  chain(i) = 1
enddo
nbondsch(1)   = nbonds
fixed         = .false.
ghost         = .false.

if (Qprint) then
  ! Write outputfile
  write(outu,'(a)') "Nonboned self-terms"
  write(outu,'(a)') "Label--epsilon[Kcal/mole]--sigma[Ang.]--epsilon1,4[Kcal/mole]--sigma1,4[Ang.]--D[Ang.**2/ps]--q[e]"
  do i = 1,maxcharges
    itype = non_of_charge(i)
    intert =  typen(itype,itype)
    write(wrtline,*) atom_labels(i),nonbonded(1,intert),nonbonded(2,intert),nonbonded(3,intert),nonbonded(4,intert),sdat(itype),qat(i)
    write(outu,'(a)') trim(adjustl(wrtline))
  enddo
  write(wrtline,*) 'TOTAL CHARGE = ',Qnet,' [e]'
  write(outu,'(a)') trim(adjustl(wrtline))
  if (Qlbond) then
    write(outu,'(a)') "Bonds"
    write(outu,'(a)') "atom1--atom2--Label1--Label2--Kbond[Kcal/mole/Ang.**2]--R0[Ang.]--Rconst[Ang.]"
    do i = 1,nbonds
      iat1 = bonds(1,i) 
      iat2 = bonds(2,i) 
      itype = psf_atomtype(iat1)
      jtype = psf_atomtype(iat2)
      ktype = bonds(3,i)
      write(wrtline,*) iat1,' ',iat2,' ',non_labels(itype),non_labels(jtype),stretch(1,ktype),stretch(2,ktype),stretch(3,ktype)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif
  if (Qlang) then
    write(outu,'(a)') "Bond angles"
    write(outu,'(a)') "atom1--atom2--atom3--Label1--Label2--Label3--Kbend[Kcal/mole/rad**2]--Theta0[degress]"
    do i = 1,nbends
      iat1 = bends(1,i) 
      iat2 = bends(2,i) ! central atom
      iat3 = bends(3,i) 
      itype = psf_atomtype(iat1)
      jtype = psf_atomtype(iat2)
      ktype = psf_atomtype(iat3)
      ltype = bends(4,i)
      write(wrtline,*) iat1,' ',iat2,' ',iat3,' ',non_labels(itype),non_labels(jtype),non_labels(ktype),bend(1,ltype),bend(2,ltype)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
    if (Qlubs) then
      write(outu,'(a)') "Urey-Bradley"
      write(outu,'(a)') "atom1--atom2--Label1--Label2----KUB[Kcal/mole/Ang.**2]--S0[Ang.]"
      do i = 1,nubs
        iat1 = ubs(1,i) 
        iat2 = ubs(2,i)
        itype = psf_atomtype(iat1)
        jtype = psf_atomtype(iat2)
        ktype = ubs(3,i)
        write(wrtline,*) iat1,' ',iat2,' ',non_labels(itype),non_labels(jtype),ubt(1,ktype),ubt(2,ktype)   
        write(outu,'(a)') trim(adjustl(wrtline))
      enddo
    endif
  endif
  if (Qldih) then
    write(outu,'(a)') "Dihedral angles"
    write(outu,'(a)') "atom1--atom2--atom3--atom4--Label1--Label2--Label3--Label4--(n--Kdih[Kcal/mole]--delta[degrees])..."
    do i = 1,ntorts
      iat1 = torts(1,i) ! terminal atom
      iat2 = torts(2,i)  
      iat3 = torts(3,i) 
      iat4 = torts(4,i) ! terminal atom
      itype = psf_atomtype(iat1)
      jtype = psf_atomtype(iat2)
      ktype = psf_atomtype(iat3)
      ltype = psf_atomtype(iat4)
      mtype = torts(5,i)
      write(wrtline,*) iat1,' ',iat2,' ',iat3,' ',iat4,' ',non_labels(itype),non_labels(jtype),non_labels(ktype),non_labels(ltype), & 
                      (ndih(j,mtype),dih((j-1)*2+1,mtype),dih((j-1)*2+2,mtype),j=1,nprms(mtype))
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif
  if (Qldef) then
    write(outu,'(a)') "Improper angles"
    write(outu,'(a)') "atom1--atom2--atom3--atom4--Label1--Label2--Label3--Label4--K[Kcal/mole/rad**2]--Omega0[degrees]"
    do i = 1,ndeforms
      iat1 = deforms(1,i)  
      iat2 = deforms(2,i) 
      iat3 = deforms(3,i) 
      iat4 = deforms(4,i) 
      itype = psf_atomtype(iat1)
      jtype = psf_atomtype(iat2)
      ktype = psf_atomtype(iat3)
      ltype = psf_atomtype(iat4)
      mtype = deforms(5,i)
      write(wrtline,*) iat1,' ',iat2,' ',iat3,' ',iat4,' ',non_labels(itype),non_labels(jtype),non_labels(ktype),non_labels(ltype),deform(1,mtype),deform(2,mtype)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif
  if (Qlcmap) then
    write(outu,'(a)') "CMAP terms"
    write(outu,'(a)') "atom1--atom2--atom3-atom4--atom5--atom6--atom7--atom8--Label1--Label2--Label3--Label4--Label5--Label6--Label7--Label8----grid points--grid spacing[degrees]"
    do i = 1,ncmaps
      itort1 = cmaps(1,i)
      itort2 = cmaps(2,i)
      ! theta angle (C-N-C_alpha-C)
      iat1 = torts(1,itort1) ! terminal atom
      iat2 = torts(2,itort1)  
      iat3 = torts(3,itort1) 
      iat4 = torts(4,itort1) ! terminal atom
      itype = psf_atomtype(iat1)
      jtype = psf_atomtype(iat2)
      ktype = psf_atomtype(iat3)
      ltype = psf_atomtype(iat4)
      ! psi angle (N-C_alpha-C-N)
      iat5 = torts(1,itort2) ! terminal atom
      iat6 = torts(2,itort2)  
      iat7 = torts(3,itort2) 
      iat8 = torts(4,itort2) ! terminal atom
      mtype = psf_atomtype(iat5)
      nntype = psf_atomtype(iat6)
      otype = psf_atomtype(iat7)
      ptype = psf_atomtype(iat8)
      ! CMAP type
      icmap = cmaps(3,i) 
      write(wrtline,*) iat1,' ',iat2,' ',iat3,' ',iat4,' ',iat5,' ',iat6,' ',iat7,' ',iat8,' ',non_labels(itype),non_labels(jtype),non_labels(ktype),non_labels(ltype), &
                       non_labels(mtype),non_labels(nntype),non_labels(otype),non_labels(ptype),cmap(icmap),gscmap(icmap)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif
endif

return
end subroutine
