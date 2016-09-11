subroutine readpsf(Qprint)
use errormod
use stdiomod
use charmmmod
use grandmod
use charfuncmod
use constamod
use listmod

implicit none

! ex explatmod
! New variables related to explicit atoms in GCMC/BD code
! ======================================================
!

! Non-bonded terms
! ----------------
!
! maxtypes              -> Number of nonbonded types 
! maxcharges            -> Number of atom types
! maxnnb                -> Number of atom type pairs
! typen                 -> Pointer for non-bonded pairs
! non_of_charge         -> Pointer for non-bonded type
! nonbonded             -> LJ parameters
! non_labels            -> Name for nonbonded types
! atom_labels           -> Name for atom types
! sdat                  -> Self-diffusion constants for nonbonded types
! qat                   -> Atom type charges

integer     maxtypes, maxcharges!, maxnnb
integer,allocatable :: typen(:,:), non_of_charge(:)
real,allocatable :: nonbonded(:,:), sdat(:), qat(:)
character*7,allocatable :: non_labels(:) 
character*9,allocatable :: atom_labels(:)
! Bonded terms
! ------------
!
! nbondt                -> Number of bond types
! nbendt                -> Number of bond angle types
! nubt                  -> Number of Urey-Bradley term types
! ntortt                -> Number of dihedral angle types
! ndeformt              -> Number of improper angle types
! ncmap                 -> Number of CMAP terms
! stretch               -> Bond parameters
! bend                  -> Bond angle parameters
! ubt                   -> Urey-Bradley parameters
! dih, ndih, nprms      -> Dihedral angle parameters
! deform                -> Improper angle parameters
! cmap                  -> CMAP grid points
! gscmap                -> CMAP grid spacing
! fcmap                 -> CMAP interpolating function
! ftcmap, fpcmap        -> CMAP gradients with respect to theta and psi angles
! ftpcmap               -> CMAP cross derivatives
! ccoef                 -> CMAP coefficient for the bicubic interpolation 
! wt                    -> Matrix necessary for the CMAP coefficient calculations

integer   nbondt, nbendt, nubt, ntortt, ndeformt, ncmap
integer,allocatable :: ndih(:,:), nprms(:), cmap(:)
real,allocatable :: stretch(:,:), bend(:,:), ubt(:,:), dih(:,:), deform(:,:), gscmap(:), fcmap(:,:), ftcmap(:,:), fpcmap(:,:), ftpcmap(:,:), ccoef(:,:,:)
real wt(16,16)
data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4, &
         10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4, &
         4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2, &
         10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2, & 
         0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2, &
         10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2, &
         5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1, &
         10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/         

! Connectivity terms
! ------------------
!
! nch                   -> Number of chains
! natt                  -> Number of explicit atoms
! Qlbond                -> Logical variable which indicates if CHARMM FF includes bonds
! nbondsch              -> Number of bonds for each chain 
! nbonds                -> Number of bonds
! Qlang                 -> Logical variable which indicates if CHARMM FF includes bond angles
! nbends                -> Number of bond angles
! Qlubs                 -> Logical variable which indicates if CHARMM FF includes UB terms
! nubs                  -> Number of Urey-Bradley terms
! Qldih                 -> Logical variable which indicates if CHARMM FF includes dihedral angles
! ntorts                -> Number of dihedral angles
! Qldef                 -> Logical variable which indicates if CHARMM FF includes improper angles
! ndeforms              -> Number of improper angles
! Qlcmap                -> Logical variable which indicates if CHARMM FF includes CMAP terms
! ncmaps                -> Number of CMAP terms
! natfx                 -> Number of fixed atoms
! nghosts               -> Number of ghost atoms
! chain                 -> Pointer which indicates the chain of each atom
! bonds                 -> Connectivity bonds indices
! bends                 -> Connectivity bond angle indices
! ubs                   -> Connectivity Urey-Bradley indices
! torts                 -> Connectivity dihedral angle indices
! deforms               -> Connectivity improper angle indices
! cmaps                 -> Connectivity CMAP indices
! lthetacmap            -> Integer variable which indicates if a dihedral angle is the theta angle for CMAP
! lpsicmap              -> Integer variable which indicates if a dihedral angle is the psi angle for CMAP
! thetacmap             -> Theta angle for CMAP
! psicmap               -> Psi angle for CMAP
! attcmap               -> Atoms which form the theta angle for CMAP
! atpcmap               -> Atoms which form the psi angle for CMAP
! nablatcmp             -> Gradient for atoms which form the theta angle for CMAP
! nablapcmp             -> Gradient for atoms which form the psi angle for CMAP
! fixed                 -> Logical variable which indicates if an atom is a fixed atom
! ghost                 -> Logical variable which indicates if an atom is a ghost atom

integer nch, natt, nbonds, nbends, nubs, ntorts, ndeforms, ncmaps 
!integer natfx, nghosts
integer,allocatable :: chain(:), nbondsch(:), bonds(:,:), bends(:,:)
integer, allocatable :: ubs(:,:), torts(:,:), deforms(:,:), cmaps(:,:)
integer, allocatable :: lthetacmap(:), lpsicmap(:), attcmap(:,:), atpcmap(:,:)
logical,allocatable :: fixed(:), ghost(:)
real, allocatable ::  thetacmap(:), psicmap(:), nablatcmp(:,:,:), nablapcmp(:,:,:)
logical*1 :: Qlbond, Qlang, Qlubs, Qldih, Qldef, Qlcmap

! Energy terms
!------------------
!
! eub                   -> Urey-Bradley energy
! eopbs                 -> improper energy
! ecmap                 -> CMAP energy
!real eub, eopbs, ecmap

! 1-2, 1-3 and 1-4 terms lists
!-----------------------------
!
! listmex               -> Number of atoms for 1-2 and 1-3 terms list
! listex                -> Pointer for atoms belonging to 1-2 and 1-3 terms list
! listm14               -> Number of atoms for 1-4 terms list
! list14                ->  Pointer for atoms belonging to 1-4 terms list

integer, allocatable :: listmex(:), listm14(:), listex(:,:), list14(:,:)

! ex explatmod

logical*1 Qprint
! local variables
integer :: dummyi, i, j, k, l, m, IDat, intert, intert2, ncharge, ilines, irest, ibonds, ibends, iubs, itort, itort1, itort2, ideform
integer :: iat1, iat2, iat3, iat4, iat5, iat6, iat7, iat8, iat9
integer :: itype, jtype, ktype, ltype, mtype, nntype, otype, ptype, icmap
integer :: n, nmax
integer, allocatable :: psf_nq(:), val(:), psf_btype(:,:), psf_bendtyp(:), psf_ubtp(:,:), cmaptype(:)
real, parameter :: tol=1.0e-3
real :: atomchar, atommass, Qnet
real, allocatable :: psf_qat(:,:), psf_mass(:)
character*7, allocatable :: psf_non_labels(:)
character :: dummyc*144, atname*7, intg*2, word*5, wrtline*2048
character*4 :: ptname
logical*1 :: dobond, doang, dodih, dodef, docmap
logical*1 :: ok, ok1, ok2
integer etypn
real, allocatable :: lj14(:,:)

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
  read(iunpsf,*) IDat,ptname,dummyi,dummyc,dummyc,atname,atomchar,atommass,dummyi
  if (i.eq.1) call addptyp(natt, ptname) !ListMod
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
  ! ListMod {
  call addetyp(atname(1:4))
  etypn=getetyp(atname(1:4))
  call seteleinptyp(nptyp,i,etypn)
  ptypl(nptyp)%chg(i)=atomchar
  ! } ListMod 
enddo ! next i
! ListMod {
call updateptypchg(nptyp)
call setpsfinptyp(nptyp)
! } ListMod 
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
allocate (lj14(2,netyp))
ncharge = 0
do itype = 1, maxtypes
  i = val(itype) 
  non_labels(itype) = charmm_label(i)
  intert = typen(itype,itype)
  j = charmm_typen(i,i) ! nonbonded pair
  do k = 1, 4
    nonbonded(k,intert) = charmm_nonbonded(k,j)
  enddo
  sdat(itype) = scldiff*kbt/(6.0*pi*viscwat*nonbonded(2,intert)) ! Stokes-Einstein equation 
  ! ListMod {
  ! Set Diff, Eps, Sigma in ListMod
  m=getetyp(psf_non_labels(itype))
  call editetyp(m,sdat(itype),nonbonded(1,intert),nonbonded(2,intert))
  lj14(1,m)=nonbonded(3,intert)
  lj14(2,m)=nonbonded(4,intert)
  ! } ListMod
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

! *** assign additional variables
nch = 1 
allocate(chain(natt),nbondsch(1),fixed(natt),ghost(natt))
chain(1:natt)=1
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

call exclude()
! 1-2 (bonds) and 1-3 (angles) pairs
!listmex(atomtypenumber)-> number of entries in second dimension in listex
!listex(entry number, atomtypenumber)->atomtypenumber-pair
! 1-4 (dih) pairs
!listm14(atomtypenumber)-> number of entries in second dimension in list14
!list14(entry number, atomtypenumber)->atomtypenumber-pair
! Sort Lists
do i=1,natt
  m=listmex(i)
  if (m.gt.0) call msort(listex(1:m,i),m)
  m=listm14(i)
  if (m.gt.0) call msort(list14(1:m,i),m)
enddo
! New 1-4 Pair List
! ListMod {
m=sum(listm14(1:natt))
ptypl(nptyp)%psf(1)%np14=m
allocate(ptypl(nptyp)%psf(1)%p14(m))
l=0
do i=1,natt
  k=listm14(i)
  do j=1,k
    l=l+1
    ptypl(nptyp)%psf(1)%p14(l)%a=i
    ptypl(nptyp)%psf(1)%p14(l)%b=list14(j,i)
  enddo
enddo
if (l.ne.m) call error ('psf_p14', 'Numbers do not match', faterr)
! Add list14 to listex
do i=1,natt
  m=listm14(i)+listmex(i)
  if (m.gt.natt) call error ('psf_p14', 'ListEx+List14 greater than num of atoms', faterr)
  listex(listmex(i)+1:m,i)=list14(1:listm14(i),i)
  listmex(i)=m
enddo
! Sort new listex
do i=1,natt
  m=listmex(i)
  if (m.gt.0) call msort(listex(1:m,i),m)
enddo
! Make list of complement of listex
ptypl(nptyp)%psf(1)%nnbon=natt*(natt-1)/2-sum(listmex(1:natt))
if (ptypl(nptyp)%psf(1)%nnbon.gt.0) then
  allocate(ptypl(nptyp)%psf(1)%nbon(ptypl(nptyp)%psf(1)%nnbon))
  k=0
  do i=1,natt-1
    do j=i+1,natt
      ok=.true.
      do m=1,listmex(i)
        if (j.eq.listex(m,i)) then
          ok=.false.
          exit
        endif
      enddo
      if (ok) then
        k=k+1
        ptypl(nptyp)%psf(1)%nbon(k)%a=i    
        ptypl(nptyp)%psf(1)%nbon(k)%b=j
      endif
    enddo
  enddo
  if (k.ne.ptypl(nptyp)%psf(1)%nnbon) call error ('psf_nbon', 'k and nnbon do not match', faterr)
endif
! Add the 1,4 eps,sig in a new internal list into psf
m=ptypl(nptyp)%psf(1)%np14
allocate(ptypl(nptyp)%psf(1)%lj(m))
do i=1,m
  j=ptypl(nptyp)%etyp(ptypl(nptyp)%psf(1)%p14(i)%a)
  k=ptypl(nptyp)%etyp(ptypl(nptyp)%psf(1)%p14(i)%b)
  ptypl(nptyp)%psf(1)%lj(i)%epp4=4.0*sqrt(lj14(1,j)*lj14(1,k))
  ptypl(nptyp)%psf(1)%lj(i)%sgp2=(0.5*(lj14(2,j)+lj14(2,k)))**2
enddo
! } ListMod
deallocate(lj14)
deallocate(list14,listm14,listex,listmex)
deallocate(psf_non_labels,psf_nq,psf_qat)
deallocate(psf_mass,val)
! ListMod {
if (Qlbond) then
  ! copy nbonds
  ptypl(nptyp)%psf%nbonds=nbonds
  ! transfer bonds
  if (allocated(ptypl(nptyp)%psf(1)%bonds)) deallocate (ptypl(nptyp)%psf(1)%bonds)
  call move_alloc(bonds,ptypl(nptyp)%psf(1)%bonds)
  ! transfer stretch
  if (allocated(ptypl(nptyp)%psf(1)%stretch)) deallocate (ptypl(nptyp)%psf(1)%stretch)
  call move_alloc(stretch,ptypl(nptyp)%psf(1)%stretch)
endif
if (Qlang) then
  if (Qlubs) then
    ! copy nubs
    ptypl(nptyp)%psf%nubs=nubs
    ! transfer ubs
    if (allocated(ptypl(nptyp)%psf(1)%ubs)) deallocate (ptypl(nptyp)%psf(1)%ubs)
    call move_alloc(ubs,ptypl(nptyp)%psf(1)%ubs)
    ! transfer ubt
    if (allocated(ptypl(nptyp)%psf(1)%ubt)) deallocate (ptypl(nptyp)%psf(1)%ubt)
    call move_alloc(ubt,ptypl(nptyp)%psf(1)%ubt)
  endif ! Qlubs
  ! copy nbends
  ptypl(nptyp)%psf%nbends=nbends
  ! transfer bends
  if (allocated(ptypl(nptyp)%psf(1)%bends)) deallocate (ptypl(nptyp)%psf(1)%bends)
  call move_alloc(bends,ptypl(nptyp)%psf(1)%bends)
  ! transfer bend
  if (allocated(ptypl(nptyp)%psf(1)%bend)) deallocate (ptypl(nptyp)%psf(1)%bend)
  call move_alloc(bend,ptypl(nptyp)%psf(1)%bend)
endif ! Qlang
if (Qldih) then
  ! copy ntorts
  ptypl(nptyp)%psf%ntorts=ntorts
  ! transfer torts
  if (allocated(ptypl(nptyp)%psf(1)%torts)) deallocate (ptypl(nptyp)%psf(1)%torts)
  call move_alloc(torts,ptypl(nptyp)%psf(1)%torts)
  ! transfer ndih
  if (allocated(ptypl(nptyp)%psf(1)%ndih)) deallocate (ptypl(nptyp)%psf(1)%ndih)
  call move_alloc(ndih,ptypl(nptyp)%psf(1)%ndih)
  ! transfer dih
  if (allocated(ptypl(nptyp)%psf(1)%dih)) deallocate (ptypl(nptyp)%psf(1)%dih)
  call move_alloc(dih,ptypl(nptyp)%psf(1)%dih)
endif ! Qldih
if (Qldef) then
  ! copy ndeforms
  ptypl(nptyp)%psf%ndeforms=ndeforms
  ! transfer deforms
  if (allocated(ptypl(nptyp)%psf(1)%deforms)) deallocate (ptypl(nptyp)%psf(1)%deforms)
  call move_alloc(deforms,ptypl(nptyp)%psf(1)%deforms)
  ! transfer deform 
  if (allocated(ptypl(nptyp)%psf(1)%deform)) deallocate (ptypl(nptyp)%psf(1)%deform)
  call move_alloc(deform,ptypl(nptyp)%psf(1)%deform)
endif ! Qldef
if (Qlcmap) then
  ! copy ncmaps
  ptypl(nptyp)%psf%ncmaps=ncmaps
  ! transfer cmaps
  if (allocated(ptypl(nptyp)%psf(1)%cmaps)) deallocate (ptypl(nptyp)%psf(1)%cmaps)
  call move_alloc(cmaps,ptypl(nptyp)%psf(1)%cmaps)
  ! transfer attcmap
  if (allocated(ptypl(nptyp)%psf(1)%attcmap)) deallocate (ptypl(nptyp)%psf(1)%attcmap)
  call move_alloc(attcmap,ptypl(nptyp)%psf(1)%attcmap)
  ! transfer gscmap
  if (allocated(ptypl(nptyp)%psf(1)%gscmap)) deallocate (ptypl(nptyp)%psf(1)%gscmap)
  call move_alloc(gscmap,ptypl(nptyp)%psf(1)%gscmap)
  ! transfer fcmap
  if (allocated(ptypl(nptyp)%psf(1)%fcmap)) deallocate (ptypl(nptyp)%psf(1)%fcmap)
  call move_alloc(fcmap,ptypl(nptyp)%psf(1)%fcmap)
  ! transfer ftcmap
  if (allocated(ptypl(nptyp)%psf(1)%ftcmap)) deallocate (ptypl(nptyp)%psf(1)%ftcmap)
  call move_alloc(ftcmap,ptypl(nptyp)%psf(1)%ftcmap)
  ! transfer fpcmap
  if (allocated(ptypl(nptyp)%psf(1)%fpcmap)) deallocate (ptypl(nptyp)%psf(1)%fpcmap)
  call move_alloc(fpcmap,ptypl(nptyp)%psf(1)%fpcmap)
  ! transfer ftpcmap
  if (allocated(ptypl(nptyp)%psf(1)%ftpcmap)) deallocate (ptypl(nptyp)%psf(1)%ftpcmap)
  call move_alloc(ftpcmap,ptypl(nptyp)%psf(1)%ftpcmap)
  ! transfer ccoef
  if (allocated(ptypl(nptyp)%psf(1)%ccoef)) deallocate (ptypl(nptyp)%psf(1)%ccoef)
  call move_alloc(ccoef,ptypl(nptyp)%psf(1)%ccoef)
endif ! Qlcmap
! } ListMod 

deallocate (psf_atomtype,psf_atomtype2)
deallocate (typen,non_of_charge,nonbonded,sdat,qat,non_labels,atom_labels,psf_charge)
if (allocated(bonds)) deallocate (bonds)
if (allocated(stretch)) deallocate (stretch)
if (allocated(bends)) deallocate (bends)
if (allocated(bend)) deallocate (bend)
if (allocated(ubs)) deallocate (ubs)
if (allocated(ubt)) deallocate (ubt)
if (allocated(torts)) deallocate (torts)
if (allocated(dih)) deallocate (dih)
if (allocated(ndih)) deallocate (ndih)
if (allocated(deforms)) deallocate (deforms)
if (allocated(deform)) deallocate (deform)
if (allocated(cmaps)) deallocate (cmaps,lthetacmap,lpsicmap,thetacmap,psicmap,attcmap,atpcmap,nablatcmp,nablapcmp,cmap,gscmap,fcmap,ftcmap,fpcmap,ftpcmap,ccoef)
deallocate (nprms)
deallocate (chain,nbondsch,fixed,ghost) 
return
contains
 
  subroutine psf_bond(ibonds,iat1,iat2,val,psf_btype)
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
  
  subroutine psf_bend(ibends,iat1,iat2,iat3,val,psf_btype)
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
  
  subroutine psf_dih(itort,iat1,iat2,iat3,iat4,val,psf_btype)
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
  
  subroutine psf_imp(ideform,iat1,iat2,iat3,iat4,val,psf_btype)
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
  
  subroutine psf_cmap(i,j,iat1,iat2,iat3,iat4,psf_mass)
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

  subroutine cmapc()
  implicit none
  ! local variables
  integer icmap,n,n2,i,j,j1,j2,j3,j4,k,l,m
  real y(4),y1(4),y2(4),y12(4),d,c(16)
  
  do icmap = 1,ncmap
    n = cmap(icmap)
    n2 = n*n
    j1 = 0
    do i = 1,n
      do j = 1,n
        ! Grid points for a rectangular grid cell (numbered counter clockwise from the lower left)
        j1 = j1 + 1
        j2 = j1 + n
        j3 = j2 + 1
        j4 = j1 + 1
        ! Symmetry considerations for grid points
        if (j.eq.n) then
          j3 = j3 - n
          j4 = j4 - n
        end if
        if (i.eq.n) then
          j2 = j2 - n2
          j3 = j3 - n2
        end if
        ! Function at the four grid points
        y(1) = fcmap(j1,icmap)
        y(2) = fcmap(j2,icmap)   
        y(3) = fcmap(j3,icmap)
        y(4) = fcmap(j4,icmap) 
        ! Gradients at the four grid points 
        y1(1) = ftcmap(j1,icmap)
        y1(2) = ftcmap(j2,icmap)
        y1(3) = ftcmap(j3,icmap)
        y1(4) = ftcmap(j4,icmap)
        y2(1) = fpcmap(j1,icmap)
        y2(2) = fpcmap(j2,icmap)
        y2(3) = fpcmap(j3,icmap)
        y2(4) = fpcmap(j4,icmap)
        ! Cross derivative at the four grid points
        y12(1) = ftpcmap(j1,icmap)
        y12(2) = ftpcmap(j2,icmap)
        y12(3) = ftpcmap(j3,icmap)
        y12(4) = ftpcmap(j4,icmap)
        ! Obtain the CMAP coefficient for the bicubic interpolation
        d = gscmap(icmap)
        call bcucof(y,y1,y2,y12,d,c)
        m = 0
        do k = 1,4
          do l = 1,4
            m = m + 1
            ccoef(m,j1,icmap) = c(m)
          end do
        end do
      end do ! next j
    end do ! next i
  end do ! next icmap
  
  return
  end subroutine

  subroutine cmapg()
  implicit none
  ! local variables
  integer icmap,n,n2,i,i1,j,j1
  real dang,idang
  real, parameter :: cte1 = 1.0/3.0, cte2 = 1.0/6.0
  real, allocatable :: xp(:),yp(:),y2p(:)
  real, allocatable :: xt(:),yt(:),y2t(:)
  real, allocatable :: csstmp1(:),csstmp2(:)
  
  do icmap = 1,ncmap
    ! total number of grid points in each direction
    n = cmap(icmap)
    n2 = n*n
    allocate (xp(n+1),yp(n+1),y2p(n+1))
    allocate (xt(n+1),yt(n+1),y2t(n+1))
    allocate (csstmp1(n2),csstmp2(n2))
    ! grid spacing
    dang = gscmap(icmap)
    idang = 1.0/dang
    ! *** Gradient calculation
    ! Psi angle
    do i = 1,n
      i1 = n*(i-1)
      do j = 1,n
        j1 = i1 + j
        ! Grid points
        xp(j) = -180.0 + (j-1)*dang
        ! Interpolating function at the grid points 
        yp(j) = fcmap(j1,icmap)
      end do
      xp(n+1) = 180.0
      yp(n+1) = yp(1)
      ! Second derivative of the interpolating function at the grid points
      call cmapspline(xp,yp,n+1,y2p)
      ! Psi gradient
      do j = 1,n
        j1 = i1 + j
        fpcmap(j1,icmap) = (yp(j+1)-yp(j))*idang - y2p(j)*dang*cte1 - y2p(j+1)*dang*cte2
      end do
    end do ! next i
    ! Theta angle
    do j = 1,n
      do i = 1,n
        i1 = (i-1)*n + j
        ! Grid points
        xt(i) = -180.0 + (i-1)*dang
        ! Interpolating function at the grid points
        yt(i) = fcmap(i1,icmap)
      end do
      xt(n+1) = 180.0
      yt(n+1) = yt(1)
      ! Second derivative of the interpolating function at the grid points
      call cmapspline(xt,yt,n+1,y2t)
      ! Theta gradient
      do i = 1,n 
        i1 = (i-1)*n + j
        ftcmap(i1,icmap) = (yt(i+1)-yt(i))*idang - y2t(i)*dang*cte1 - y2t(i+1)*dang*cte2
      end do
    end do ! next j
    ! *** Cross derivative calculation
    ! Cross derivative (theta-psi)
    do i = 1,n
      i1 = n*(i-1)
      do j = 1,n
        j1 = i1 + j
        ! Interpolating function at the grid points 
        yp(j) = ftcmap(j1,icmap) 
      end do
      yp(n+1) = yp(1)
      ! Second derivative of the interpolating function at the grid points
      call cmapspline(xp,yp,n+1,y2p)
      ! Cross derivative 
      do j = 1,n
        j1 = i1 + j
        csstmp1(j1) = (yp(j+1)-yp(j))*idang - y2p(j)*dang*cte1 - y2p(j+1)*dang*cte2
      end do
    end do ! next i
    ! Cross derivative (psi-theta)
    do j = 1,n
      do i = 1,n
        i1 = (i-1)*n + j
        ! Interpolating function at the grid points
        yt(i) = fpcmap(i1,icmap)
      end do
      yt(n+1) = yt(1)
      ! Second derivative of the interpolating function at the grid points
      call cmapspline(xt,yt,n+1,y2t)
      ! Cross derivative
      do i = 1,n 
        i1 = (i-1)*n + j
        csstmp1(i1) = (yt(i+1)-yt(i))*idang - y2t(i)*dang*cte1 - y2t(i+1)*dang*cte2
      end do
    end do ! next j
    ! Cross derivative (average)
    j1 = 0
    do i = 1,n
      do j = 1,n
        j1 = j1 + 1
        ftpcmap(j1,icmap) = (csstmp1(j1)+csstmp2(j1))*0.5
      end do
    end do   
    deallocate (xp,yp,y2p)
    deallocate (xt,yt,y2t)
    deallocate (csstmp1,csstmp2)
  end do ! next icmap
  
  return
  end subroutine

  subroutine cmapspline(x,y,n,y2)
  implicit none
  ! INPUT:
  ! x(1)<x(2)<...<x(n) -> grid points
  ! y(1:n)             -> tabulated function
  ! OUPUT:
  ! y2(1:n)            -> second derivative of the interpolating function at the grid points
  integer n
  real x(n),y(n),y2(n)
  ! local variables
  integer i
  real p,sig,u(n-1)
  
  ! The lower boundary condition is set to be "natural"
  y2(1) = 0.0
  u(1) = 0.0
  ! Descomposition loop of the tridiagonal algortihm. y2 and u are used 
  ! for temporary storage of the descomposed factors
  do i = 2,n-1
    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p = sig*y2(i-1) + 2.0
    y2(i) = (sig-1.0)/p
    u(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) &
         - sig*u(i-1))/p 
  end do
  ! The upper boundary condition is set to be "natural
  y2(n) = 0.0
  ! Backsubstitution loop of the tridiagonal algorithm
  do i = n-1,1,-1
    y2(i) = y2(i)*y2(i+1) + u(i)
  end do
  return
  end subroutine

  subroutine bcucof(y,y1,y2,y12,d,c)
  implicit none
  ! INPUT:
  ! d           -> length of the grid cell 
  ! y,y1,y2,y12 -> function, gradients, and cross derivatives at the four grid points
  !                of a rectangular grid cell (numbered counterclockwise from the lower left)
  ! OUPUT:
  ! c           -> coefficients for the bicubic interpolation
  real d,c(16),y(4),y1(4),y2(4),y12(4)
  ! local variables
  integer i,j
  real d2,xx,x(16)
  
  d2 = d*d
  ! Pack a temporary vector x
  do i = 1,4
    x(i) = y(i)
    x(i+4) = y1(i)*d
    x(i+8) = y2(i)*d
    x(i+12) = y12(i)*d2
  end do
  ! Matrix multiply by the stored table
  do i = 1,16
    xx = 0.0
    do j = 1,16
      xx = xx + wt(i,j)*x(j)
    end do
    c(i) = xx
  end do
  return
  end subroutine

  subroutine exclude()
  implicit none
  ! local variables
  integer iat,nex,ibond,ix,iex
  integer jat,kat,ibend,itmp,i14a,i14b,ibond1,ibond2
  logical*1 ok,HasBeenUsed,Found14
  allocate (listmex(natt),listex(natt,natt),listm14(natt),list14(natt,natt))
  listmex = 0
  listm14 = 0
  list14 = 0
  listex = 0
  ! MAKE A LIST FOR 1-2 AND 1-3 TERMS
  do iat = 1, natt - 1  ! search tables for each atom
    nex = 0  ! initialized number excluded to 0
    ! determine all bond interactions with iat
    ! NOTE: bonds(2,) > bonds(1,)
    do ibond = 1, nbonds
      if (bonds(1,ibond) .eq. iat) then ! exclude
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
end subroutine
