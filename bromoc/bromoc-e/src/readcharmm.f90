subroutine readcharmm(Qprint)

use errormod
use stdiomod
use charmmmod
use charfuncmod

implicit none
logical*1 Qprint
! local variables
integer i,j,k,l,m,n,o,p,imax,itype,jtype,ktype,ltype,mtype,ntype,otype,ptype,dih1,dih2,diht
integer, allocatable :: tmp_dtype(:,:),tmp_ndih(:),tmp_nterms(:)
real, parameter :: tol=1.0e-2,cte=2**(-1.0/6.0)
real, allocatable :: tmp_dih(:,:)
character com*2048,word*1024,wrtline*2048 
character*7 atnam1,atnam2,atnam3,atnam4
logical*1 ok,test_bond
logical*1 endlog,bondlog,anglog,dihlog,implog,cmaplog,nonblog,nbflog,hbonlog

! *** OBTENTION OF DIMENSIONS
chmmntype = 0 ! Number of CHARMM atom types
chmmbond  = 0 ! Number of CHARMM bond types
chmmang   = 0 ! Number of CHARMM bond angle types
chmmub    = 0 ! Number of CHARMM Urey-Bradley types
chmmdih   = 0 ! Number of CHARMM dihedral angle types
chmmimp   = 0 ! Number of CHARMM improper angle types
chmmcmap  = 0 ! Number of Cross-term energy correction map types
chmmnonb  = 0 ! Number of nonbonded pair types
chmmnbfix = 0 ! Number of VDW interactions between specific atom pair types to be modified
chmmhbond = 0 ! Number of hydrogen bond types
! SECTION B: ATOMS
call getlin(com,iunff,outu)
call getfirst(com,word)
word = lcase(word)
if (word(1:4).ne.'atom') call error ('readcharmm', 'ATOMS section is not found', faterr)
ok = .true.
do while (ok)
  call getlin(com,iunff,outu)  
  call getfirst(com,word)
  word = lcase(word)
  endlog  = word(1:3).eq.'end' 
  bondlog = word(1:4).eq.'bond'
  anglog  = word(1:4).eq.'angl'.or.word(1:4).eq.'thet'
  dihlog  = word(1:4).eq.'dihe'.or.word(1:3).eq.'phi'
  implog  = word(1:4).eq.'impr'.or.word(1:4).eq.'imph'
  cmaplog = word(1:4).eq.'cmap'
  nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
  nbflog  = word(1:4).eq.'nbfi'
  hbonlog = word(1:4).eq.'hbon'
  ok = .not.endlog .and. .not.bondlog .and. .not.anglog .and. .not.dihlog .and. .not.implog .and. .not.cmaplog .and. .not.nonblog .and. .not.nbflog .and. &
       .not.hbonlog
  if (ok) then
    if (word(1:4).eq.'mass') then
      chmmntype = chmmntype + 1
    else
      call error ('readcharmm', 'Unknown command', faterr)
    endif
  endif
enddo
if (chmmntype.eq.0) call error ('readcharmm', 'Wrong number of atom types', faterr)
! SECTION C: BONDS
if (bondlog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)  
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    anglog  = word(1:4).eq.'angl'.or.word(1:4).eq.'thet'
    dihlog  = word(1:4).eq.'dihe'.or.word(1:3).eq.'phi'
    implog  = word(1:4).eq.'impr'.or.word(1:4).eq.'imph'
    cmaplog = word(1:4).eq.'cmap'
    nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.anglog .and. .not.dihlog .and. .not.implog .and. .not.cmaplog .and. .not.nonblog .and. .not.nbflog .and. &
         .not.hbonlog
    if (ok) chmmbond = chmmbond + 1
  enddo
endif
Qchmmbond = chmmbond.ne.0
! SECTION D: ANGLES AND UREY-BRADLEY
if (anglog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)  
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    dihlog  = word(1:4).eq.'dihe'.or.word(1:3).eq.'phi'
    implog  = word(1:4).eq.'impr'.or.word(1:4).eq.'imph'
    cmaplog = word(1:4).eq.'cmap'
    nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.dihlog .and. .not.implog .and. .not.cmaplog .and. .not.nonblog .and. .not.nbflog .and. .not.hbonlog
    if (ok) then 
      chmmang = chmmang + 1
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word) 
      call getfirst(com,word)
      if (len_trim(word).ne.0) chmmub = chmmub + 1
    endif
  enddo
endif
Qchmmang = chmmang.ne.0
Qchmmub = chmmub.ne.0
! SECTION E: DIHEDRALS
if (dihlog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    implog  = word(1:4).eq.'impr'.or.word(1:4).eq.'imph'
    cmaplog = word(1:4).eq.'cmap'
    nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.implog .and. .not.cmaplog .and. .not.nonblog .and. .not.nbflog .and. .not.hbonlog
    if (ok) chmmdih = chmmdih +1
  enddo
endif
Qchmmdih = chmmdih.ne.0
! SECTION F: IMPROPER
if (implog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    cmaplog = word(1:4).eq.'cmap'
    nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.cmaplog .and. .not.nonblog .and. .not.nbflog .and. .not.hbonlog
    if (ok) chmmimp = chmmimp +1
  enddo
endif
Qchmmimp = chmmimp.ne.0
! SECTION G: CMAP
if (cmaplog) then
  imax = 0
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    nonblog = word(1:4).eq.'nonb'.or.word(1:4).eq.'nbon'
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.nonblog .and. .not.nbflog .and. .not.hbonlog    
    if (ok) then
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      call getfirst(com,word)
      i = chr2int(word)
      if (mod(360,i).ne.0) call error ('readcharmm', 'Wrong number of CMAP grid points', faterr)
      if (i.gt.imax) imax = i
      k = i/5
      if (k*5.lt.i) k = k + 1
      chmmcmap = chmmcmap +1
      do j = 1, i
        do l = 1, k
          call getlin(com,iunff,outu)
        enddo
      enddo  
    endif
  enddo
endif 
Qchmmcmap = chmmcmap.ne.0   
! SECTION H: NONBONDED
if (nonblog) then
  i = 0
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    nbflog  = word(1:4).eq.'nbfi'
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.nbflog .and. .not.hbonlog   
    if (ok) i = i + 1 
  enddo
  if (i.ne.chmmntype) call error ('readcharmm', 'Wrong number of atom types in NONBONDED section', faterr)
  do itype = 1, chmmntype
    do jtype = itype, chmmntype
      chmmnonb = chmmnonb + 1
    enddo
  enddo
else
  call error ('readcharmm', 'NONBONDED section is not found', faterr)
endif
! SECTION I: NBFIX
if (nbflog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    hbonlog = word(1:4).eq.'hbon'
    ok = .not.endlog .and. .not.hbonlog    
    if (ok) chmmnbfix = chmmnbfix + 1
  enddo 
  if (chmmnbfix.gt.chmmnonb) call error ('readcharmm', 'Wrong number of VDW interactions between specific atom pair types to be modified', faterr)
endif
Qchmmnbfix = chmmnbfix.ne.0 
! SECTION J: HBONDS
if (hbonlog) then
  ok = .true.
  do while (ok)
    call getlin(com,iunff,outu) 
    call getfirst(com,word)
    word = lcase(word)
    endlog  = word(1:3).eq.'end' 
    ok = .not.endlog     
    if (ok) chmmhbond = chmmhbond + 1
  enddo
endif
Qchmmhbond = chmmhbond.ne.0
! SECTION K: END
if (.not.endlog) call error ('readcharmm', 'END command not found', faterr)
! CHECKING THE FILE FORMAT
if (.not.Qchmmbond .and. bondlog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmang .and. anglog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmdih .and. dihlog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmimp .and. implog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmcmap .and. cmaplog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmnbfix .and. nbflog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
! if (.not.Qchmmhbond .and. hbonlog) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmbond .and. (Qchmmang.or.Qchmmub.or.Qchmmdih.or.Qchmmimp.or.Qchmmcmap)) &
  call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmang .and. (Qchmmdih.or.Qchmmimp.or.Qchmmcmap)) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (.not.Qchmmdih .and. Qchmmcmap) call error ('readcharmm', 'Wrong format for CHARMM forces field parameter file', faterr)
if (Qchmmhbond) call error ('readcharmm', 'Hydrogen bond types are not supported in this implementation', warning)
! ALLOCATIONS
allocate(charmm_label(chmmntype),charmm_mass(chmmntype))
if (Qchmmbond) allocate(charmm_btype(2,chmmbond),charmm_bond(2,chmmbond))
if (Qchmmang) then
  allocate(charmm_atype(3,chmmang),charmm_ang(2,chmmang),charmm_lub(chmmang))
  if (Qchmmub) allocate(charmm_ubtype(2,chmmub),charmm_ub(2,chmmub))
  charmm_lub = 0 
endif
if (Qchmmimp) allocate(charmm_itype(4,chmmimp),charmm_imp(2,chmmimp))
if (Qchmmcmap) then
  imax=imax*imax
  allocate(charmm_icmap(8,chmmcmap),charmm_icmap2(2,chmmcmap),charmm_ncmap(chmmcmap),charmm_cmap(chmmcmap),charmm_fcmap(imax,chmmcmap))
endif
allocate(charmm_typen(chmmntype,chmmntype),charmm_nonbonded(4,chmmnonb))

! *** OBTENTION OF PARAMETERS
rewind(unit=iunff)
! SECTION B: ATOMS
call getlin(com,iunff,outu)
do i = 1, chmmntype
  call getlin(com,iunff,outu)
  call getfirst(com,word) 
  call getfirst(com,word)
  call getfirst(com,charmm_label(i)) ! atom type name
  call getfirst(com,word)
  charmm_mass(i) = chr2real(word) ! atom type mass
enddo ! next i
! SECTION C: BONDS
! V(bond) = Kb(b-b0)**2
! Kb: Kcal/mole/A**2
! b0 : A
if (Qchmmbond) then
  call getlin(com,iunff,outu)
  do i = 1, chmmbond    
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,itype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,jtype)
    if (itype.gt.jtype) then
      charmm_btype(1,i) = jtype
      charmm_btype(2,i) = itype
    else
      charmm_btype(1,i) = itype
      charmm_btype(2,i) = jtype
    endif
    call getfirst(com,word)
    charmm_bond(1,i) = chr2real(word) ! Kb
    if (charmm_bond(1,i).lt.0.0) call error ('readcharmm', 'Wrong Kb value in BONDS section', faterr)
    call getfirst(com,word)
    charmm_bond(2,i) = chr2real(word) ! b0
    if (charmm_bond(2,i).le.0.0) call error ('readcharmm', 'Wrong b0 value in BONDS section', faterr)
  enddo ! next i
endif ! Qchmmbond
! SECTION D: ANGLES AND UREY-BRADLEY
! V(angle) = Ktheta(Theta-Theta0)**2
! Ktheta: Kcal/mole/rad**2
! Theta0: degrees
! V(Urey-Bradley) = Kub(S-S0)**2
! Kub: Kcal/mole/A**2
! S0: A
if (Qchmmang) then
  call getlin(com,iunff,outu)
  j = 0
  do i = 1, chmmang
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,itype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,jtype) 
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,ktype)
    ! check bond 1
    if (.not.test_bond(itype,jtype)) then 
      call error ('readcharmm', 'bond is not found in ANGLES section', warning)      
      write(wrtline,*) 'ANGLE: ',charmm_label(itype),charmm_label(jtype),charmm_label(ktype)
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'BOND MISSING: ',charmm_label(itype),charmm_label(jtype)
      write(outu,'(a)') trim(adjustl(wrtline))
    endif
    ! check bond2
    if (.not.test_bond(jtype,ktype)) then 
      call error ('readcharmm', 'bond is not found in ANGLES section', warning)
      write(wrtline,*) 'ANGLE: ',charmm_label(itype),charmm_label(jtype),charmm_label(ktype)
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'BOND MISSING: ',charmm_label(jtype),charmm_label(ktype)
      write(outu,'(a)') trim(adjustl(wrtline))
    endif
    charmm_atype(1,i) = itype
    charmm_atype(2,i) = jtype ! central atom
    charmm_atype(3,i) = ktype
    call getfirst(com,word)
    charmm_ang(1,i) = chr2real(word) ! Ktheta
    if (charmm_ang(1,i).lt.0.0) call error ('readcharmm', 'Wrong value for Ktheta in ANGLES section', faterr)
    call getfirst(com,word)
    charmm_ang(2,i) = chr2real(word) ! Theta0
    if (charmm_ang(2,i).lt.0.0 .or. charmm_ang(2,i).gt.180.0) call error ('readcharmm', 'Wrong value for Theta0 in ANGLES section', faterr)
    call getfirst(com,word)
    if (len_trim(word).ne.0) then
      j = j + 1
      charmm_lub(i) = j
      if (itype.gt.ktype) then
        charmm_ubtype(1,j) = ktype
        charmm_ubtype(2,j) = itype
      else
        charmm_ubtype(1,j) = itype
        charmm_ubtype(2,j) = ktype
      endif
      charmm_ub(1,j) = chr2real(word) ! Kub
      if (charmm_ub(1,j).lt.0.0) call error ('readcharmm', 'Wrong value for Kub in ANGLES section', faterr)
      call getfirst(com,word)
      charmm_ub(2,j) = chr2real(word) ! S0
      if (charmm_ub(2,j).le.0.0) call error ('readcharmm', 'Wrong value for S0 in ANGLES section', faterr)
    endif
  enddo ! next i
endif ! Qchmmang
! SECTION E: DIHEDRALS
! V(dihedral) = Kchi(1+cos(n(chi)-delta))
! Kchi: Kcal/mole
! n: multiplicity
! delta : degrees
if (Qchmmdih) then
  call getlin(com,iunff,outu)
  allocate(tmp_dtype(4,chmmdih),tmp_dih(2,chmmdih),tmp_ndih(chmmdih),tmp_nterms(chmmdih))
  tmp_nterms = 1
  charmm_nmax = 1
  diht = 0
  do i = 1, chmmdih
    call getlin(com,iunff,outu)
    call getfirst(com,word)
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,itype)
    else
      itype = 0
    endif
    call getfirst(com,word)
    if (lcase(word).eq.'x') call error ('readcharmm', 'Wildcars are only valid for terminal atoms in DIHEDRALS section', faterr)
    call fatnam(charmm_label,1,chmmntype,word,jtype)
    call getfirst(com,word)
    if (lcase(word).eq.'x') call error ('readcharmm', 'Wildcars are only valid for terminal atoms in DIHEDRALS section', faterr)
    call fatnam(charmm_label,1,chmmntype,word,ktype)
    call getfirst(com,word)
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,ltype)
    else
      ltype = 0
    endif
    ! check bond 1
    if (itype.gt.0) then
      if (.not.test_bond(itype,jtype)) then 
        call error ('readcharmm', 'bond is not found in DIHEDRALS section', warning)
        if (ltype.eq.0) then
          atnam4 = 'X'
        else
          atnam4 = charmm_label(ltype)
        endif
        write(wrtline,*) 'DIHEDRAL: ',charmm_label(itype),charmm_label(jtype),charmm_label(ktype),atnam4
        write(outu,'(a)') trim(adjustl(wrtline))
        write(wrtline,*) 'BOND MISSING: ',charmm_label(itype),charmm_label(jtype)
        write(outu,'(a)') trim(adjustl(wrtline))
      endif
    endif
    ! check bond 2
    if (.not.test_bond(jtype,ktype)) then 
      call error ('readcharmm', 'bond is not found in DIHEDRALS section', warning)
     if (itype.eq.0) then
       atnam1 = 'X'
      else
       atnam1 = charmm_label(itype)
      endif
     if (ltype.eq.0) then
       atnam4 = 'X'
      else
       atnam4 = charmm_label(ltype)
      endif
      write(wrtline,*) 'DIHEDRAL: ',atnam1,charmm_label(jtype),charmm_label(ktype),atnam4
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'BOND MISSING: ',charmm_label(jtype),charmm_label(ktype)
      write(outu,'(a)') trim(adjustl(wrtline))
    endif
    ! check bond 3
    if (ltype.gt.0) then
      if (.not.test_bond(ktype,ltype)) then 
        call error ('readcharmm', 'bond is not found in DIHEDRALS section', warning)
        if (itype.eq.0) then
          atnam1 = 'X'
        else
          atnam1 = charmm_label(itype)
        endif
        write(wrtline,*) 'DIHEDRAL: ',atnam1,charmm_label(jtype),charmm_label(ktype),charmm_label(ltype)
        write(outu,'(a)') trim(adjustl(wrtline))
        write(wrtline,*) 'BOND MISSING: ',charmm_label(ktype),charmm_label(ltype)
        write(outu,'(a)') trim(adjustl(wrtline))
      endif
    endif
    tmp_dtype(1,i) = itype ! terminal atom
    tmp_dtype(2,i) = jtype 
    tmp_dtype(3,i) = ktype   
    tmp_dtype(4,i) = ltype ! terminal atom
    call getfirst(com,word)
    tmp_dih(1,i) = chr2real(word) ! Kchi
    if (tmp_dih(1,i).lt.0.0) then 
      call error ('readcharmm', 'Unusual value for Kchi in DIHEDRALS section', warning)
     if (itype.eq.0) then
       atnam1 = 'X'
      else
       atnam1 = charmm_label(itype)
      endif
     if (ltype.eq.0) then
       atnam4 = 'X'
      else
       atnam4 = charmm_label(ltype)
      endif
      write(wrtline,*) 'DIHEDRAL: ',atnam1,charmm_label(jtype),charmm_label(ktype),atnam4
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'Kchi = ',tmp_dih(1,i)
      write(outu,'(a)') trim(adjustl(wrtline))      
    endif
    call getfirst(com,word)
    tmp_ndih(i) = chr2int(word) ! n
    if (tmp_ndih(i).lt.1 .or. tmp_ndih(i).gt.6) call error ('readcharmm', 'Wrong value for periodicity in DIHEDRALS section', faterr)
    call getfirst(com,word)
    tmp_dih(2,i) = chr2real(word) ! delta
    if (tmp_dih(2,i).ne.0.0 .and. tmp_dih(2,i).ne.180.0) then 
      call error ('readcharmm', 'Unusual value for phase in DIHEDRALS section', warning)
     if (itype.eq.0) then
       atnam1 = 'X'
      else
       atnam1 = charmm_label(itype)
      endif
     if (ltype.eq.0) then
       atnam4 = 'X'
      else
       atnam4 = charmm_label(ltype)
      endif
      write(wrtline,*) 'DIHEDRAL: ',atnam1,charmm_label(jtype),charmm_label(ktype),atnam4
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'PHASE VALUE = ',tmp_dih(2,i)
      write(outu,'(a)') trim(adjustl(wrtline))
    endif
    ! multiple dihedral angles
    j = 0
    ok = .false.
    do while (j.lt.(i-1) .and. .not.ok)
      j = j + 1 
      ok = .true.
      do k = 1, 4
        ok = ok .and. tmp_dtype(k,i).eq.tmp_dtype(k,j)
      enddo
      if (.not.ok) then
        ok = .true.
        do k = 1, 4 
          ok = ok .and. tmp_dtype(k,i).eq.tmp_dtype(5-k,j)
        enddo
      endif
    enddo
    if (ok) then
      tmp_nterms(j) = tmp_nterms(j) + 1
      if (tmp_nterms(j).gt.charmm_nmax) charmm_nmax = tmp_nterms(j)
    else
      diht = diht + 1
    endif
  enddo ! next i
  ! REORGANIZATION
  allocate (charmm_dtype(4,diht),charmm_dih(2*charmm_nmax,diht),charmm_ndih(charmm_nmax,diht),charmm_nprms(diht))
  charmm_nprms = 1
  j = 0
  ! 1) Dihedral types without wildcards (A-B-C-D)
  do i = 1, chmmdih
    if (tmp_dtype(1,i).gt.0.and.tmp_dtype(4,i).gt.0) then
      ! multiple dihedrals 
      k = 0
      ok = .false.
      do while (k.lt.j .and. .not.ok) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(l,k)
        end do
        if (.not.ok) then
          ok = .true.
          do l = 1, 4
            ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(5-l,k)
          enddo
        endif
      enddo
      if (.not.ok) then
        j = j + 1
        do l = 1, 4
          charmm_dtype(l,j) = tmp_dtype(l,i)
        enddo
        do l = 1, 2
          charmm_dih(l,j) = tmp_dih(l,i)
        enddo
        charmm_ndih(1,j) = tmp_ndih(i)
      else ! multiple
        charmm_nprms(k) = charmm_nprms(k) + 1
        m = charmm_nprms(k)
        n = (m-1)*2
        do l = 1, 2
          charmm_dih(n+l,k) = tmp_dih(l,i)
        enddo
        charmm_ndih(m,k) = tmp_ndih(i)
      endif
    endif 
  enddo ! next i
  ! 2) Dihedral types with 1 wildcard (X-A-B-C;A-B-C-X)
  do i = 1, chmmdih
    if ((tmp_dtype(1,i).eq.0.and.tmp_dtype(4,i).ne.0) .or. (tmp_dtype(1,i).ne.0.and.tmp_dtype(4,i).eq.0)) then
      ! multiple dihedrals 
      k = 0
      ok = .false.
      do while (k.lt.j .and. .not.ok) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(l,k)
        end do
        if (.not.ok) then
          ok = .true.
          do l = 1, 4
            ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(5-l,k)
          enddo
        endif
      enddo
      if (.not.ok) then
        j = j + 1
        do l = 1, 4
          charmm_dtype(l,j) = tmp_dtype(l,i)
        enddo
        do l = 1, 2
          charmm_dih(l,j) = tmp_dih(l,i)
        enddo
        charmm_ndih(1,j) = tmp_ndih(i)
      else ! multiple
        charmm_nprms(k) = charmm_nprms(k) + 1
        m = charmm_nprms(k)
        n = (m-1)*2
        do l = 1, 2
          charmm_dih(n+l,k) = tmp_dih(l,i)
        enddo
        charmm_ndih(m,k) = tmp_ndih(i)
      endif
    endif
  enddo ! next i 
  ! 3) Dihedral types with 2 wildcards (X-A-B-X)
  do i = 1, chmmdih
    if (tmp_dtype(1,i).eq.0.and.tmp_dtype(4,i).eq.0) then
      ! multiple dihedrals 
      k = 0
      ok = .false.
      do while (k.lt.j .and. .not.ok) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(l,k)
        end do
        if (.not.ok) then
          ok = .true.
          do l = 1, 4
            ok = ok .and. tmp_dtype(l,i).eq.charmm_dtype(5-l,k)
          enddo
        endif
      enddo
      if (.not.ok) then
        j = j + 1
        do l = 1, 4
          charmm_dtype(l,j) = tmp_dtype(l,i)
        enddo
        do l = 1, 2
          charmm_dih(l,j) = tmp_dih(l,i)
        enddo
        charmm_ndih(1,j) = tmp_ndih(i)
      else ! multiple
        charmm_nprms(k) = charmm_nprms(k) + 1
        m = charmm_nprms(k)
        n = (m-1)*2
        do l = 1, 2
          charmm_dih(n+l,k) = tmp_dih(l,i)
        enddo
        charmm_ndih(m,k) = tmp_ndih(i)
      endif
    endif       
  enddo ! next i 
  deallocate(tmp_dtype,tmp_dih,tmp_ndih,tmp_nterms)
  chmmdih = diht
endif ! Qchmmdih
! SECTION F: IMPROPER
! V(improper) = Kpsi(psi-psi0)**2
! Kpsi: Kcal/mole/rad**2
! psi0: degrees
! Ordinarily, improper dihedrals are given a multiplicity of 0, which imposes a harmonic restoring potential
! instead of a cosine function.  In this case, the central atom must be either the first or the last atom 
if (Qchmmimp) then
  call getlin(com,iunff,outu)
  allocate(tmp_dtype(4,chmmimp),tmp_dih(2,chmmimp))
  do i = 1, chmmimp
    j = 0
    call getlin(com,iunff,outu)
    call getfirst(com,word)     
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,itype)
    else
      itype = 0
      j = j + 1
    endif
    call getfirst(com,word)     
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,jtype)
    else
      jtype = 0
      j = j + 1
    endif
    call getfirst(com,word)     
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,ktype)
    else
      ktype = 0
      j = j + 1
    endif
    call getfirst(com,word)     
    if (lcase(word).ne.'x') then 
      call fatnam(charmm_label,1,chmmntype,word,ltype)
    else
      ltype = 0
      j = j + 1
    endif
    if (j.gt.2) call error ('readcharmm', 'Wrong wildcars in IMPROPER section', faterr)
    tmp_dtype(1,i) = itype 
    tmp_dtype(2,i) = jtype 
    tmp_dtype(3,i) = ktype   
    tmp_dtype(4,i) = ltype 
    call getfirst(com,word)
    tmp_dih(1,i) = chr2real(word) ! Kpsi
    if (tmp_dih(1,i).lt.0.0) then 
      call error ('readcharmm', 'Unusual Kpsi in IMPROPER section', warning)
      if (itype.gt.0) then
        atnam1 = charmm_label(itype)
      else
        atnam1 = 'X'
      endif
      if (jtype.gt.0) then
        atnam2 = charmm_label(jtype)
      else
        atnam2 = 'X'
      endif
      if (ktype.gt.0) then
        atnam3 = charmm_label(ktype)
      else
        atnam3 = 'X'
      endif
      if (ltype.gt.0) then
        atnam4 = charmm_label(ltype)
      else
        atnam4 = 'X'
      endif
      write(wrtline,*) 'IMPROPER: ',atnam1,atnam2,atnam3,atnam4
      write(outu,'(a)') trim(adjustl(wrtline))
      write(wrtline,*) 'Kpsi = ',tmp_dih(1,i)
      write(outu,'(a)') trim(adjustl(wrtline))
    endif
    call getfirst(com,word)
    j = chr2int(word)
    if (j.ne.0) call error ('readcharmm', 'Wrong multiplicity in IMPROPER section', faterr)
    call getfirst(com,word)
    tmp_dih(2,i) = chr2real(word) ! psi0
    if (tmp_dih(2,i).lt.-180.0 .or. tmp_dih(2,i).gt.180.0) call error ('readcharmm', 'Wrong psi0 in IMPROPER section', faterr) 
  enddo ! next i
  ! REORGANIZATION
  j = 0
  ! 1) Improper types without wildcards (A-B-C-D)
  do i = 1, chmmimp
    if (tmp_dtype(1,i).gt.0.and.tmp_dtype(2,i).gt.0..and.tmp_dtype(3,i).gt.0..and.tmp_dtype(4,i).gt.0) then
      ! multiple impropers 
      k = 0
      do while (k.lt.j) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_itype(l,k)
        end do
        if (ok) call error ('readcharmm', 'Multiple impropers is not implemented', faterr)
      enddo
      j = j + 1
      do k = 1, 4
        charmm_itype(k,j) = tmp_dtype(k,i)
      enddo
      do k = 1, 2
        charmm_imp(k,j) = tmp_dih(k,i)
      enddo
    endif 
  enddo ! next i
  ! 2) Improper types with 2 wildcards (A-X-X-B)
  do i = 1, chmmimp
    if (tmp_dtype(2,i).eq.0.and.tmp_dtype(3,i).eq.0) then
      ! multiple impropers 
      k = 0
      do while (k.lt.j) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_itype(l,k)
        end do
        if (ok) call error ('readcharmm', 'Multiple impropers is not implemented', faterr)
      enddo
      j = j + 1
      do k = 1, 4
        charmm_itype(k,j) = tmp_dtype(k,i)
      enddo
      do k = 1, 2
        charmm_imp(k,j) = tmp_dih(k,i)
      enddo
    endif 
  enddo ! next i
  ! 3) Improper types with 1 wildcard (X-A-B-C)
  do i = 1,  chmmimp
    if (tmp_dtype(2,i).eq.0.and.tmp_dtype(2,i).gt.0..and.tmp_dtype(3,i).gt.0..and.tmp_dtype(4,i).gt.0) then
      ! multiple impropers 
      k = 0
      do while (k.lt.j) 
        k = k + 1
        ok = .true.
        do l = 1, k
          ok = ok .and. tmp_dtype(l,i).eq.charmm_itype(l,k)
        end do
        if (ok) call error ('readcharmm', 'Multiple impropers is not implemented', faterr)
      enddo
      j = j + 1
      do k = 1, 4
        charmm_itype(k,j) = tmp_dtype(k,i)
      enddo
      do k = 1, 2
        charmm_imp(k,j) = tmp_dih(k,i)
      enddo
    endif 
  enddo ! next i
  ! 4) Improper types with 2 wildcards (X-A-B-X)
  do i = 1, chmmimp
    if (tmp_dtype(1,i).eq.0.and.tmp_dtype(4,i).eq.0) then
      ! multiple impropers 
      k = 0
      do while (k.lt.j) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_itype(l,k)
        end do
        if (ok) call error ('readcharmm', 'Multiple impropers is not implemented', faterr)
      enddo
      j = j + 1
      do k = 1, 4
        charmm_itype(k,j) = tmp_dtype(k,i)
      enddo
      do k = 1, 2
        charmm_imp(k,j) = tmp_dih(k,i)
      enddo
    endif 
  enddo ! next i
  ! 5) Improper types with 2 wildcards (X-X-A-B)
  do i = 1, chmmimp
    if (tmp_dtype(1,i).eq.0.and.tmp_dtype(2,i).eq.0) then
      ! multiple impropers
      k = 0
      do while (k.lt.j) 
        k = k + 1
        ok = .true.
        do l = 1, 4
          ok = ok .and. tmp_dtype(l,i).eq.charmm_itype(l,k)
        end do
        if (ok) call error ('readcharmm', 'Multiple impropers is not implemented', faterr)
      enddo
      j = j + 1
      do k = 1, 4
        charmm_itype(k,j) = tmp_dtype(k,i)
      enddo
      do k = 1, 2
        charmm_imp(k,j) = tmp_dih(k,i)
      enddo
    endif 
  enddo ! next i
  deallocate(tmp_dtype,tmp_dih)
endif ! Qchmmimp
! SECTION G: CMAP
if (Qchmmcmap) then
  call getlin(com,iunff,outu)
  do i = 1, chmmcmap
    call getlin(com,iunff,outu)
    ! get angles   
    ! first angle: theta angle (C-N-C_alpha-C)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,itype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,jtype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,ktype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,ltype)
    if (abs(charmm_mass(itype)-12.011).lt.tol .and. abs(charmm_mass(ltype)-12.011).lt.tol) then 
      ok = (abs(charmm_mass(jtype)-14.007).lt.tol.and.abs(charmm_mass(ktype)-12.011).lt.tol) .or. &
           (abs(charmm_mass(jtype)-12.011).lt.tol.and.abs(charmm_mass(ktype)-14.007).lt.tol)
      if (.not.ok) call error ('readcharmm', 'Wrong atom types in CMAP section', faterr)
    else
      call error ('readcharmm', 'Wrong atom types in CMAP section', faterr)
    endif   
    charmm_icmap(1,i) = itype 
    charmm_icmap(2,i) = jtype 
    charmm_icmap(3,i) = ktype   
    charmm_icmap(4,i) = ltype
    ! check dihedral angle 
    ok = .false. 
    k = 0
    do while (k.lt.chmmdih .and. .not.ok)
      k = k + 1
      ok = (charmm_icmap(1,i).eq.charmm_dtype(1,k).and.charmm_icmap(2,i).eq.charmm_dtype(2,k).and.charmm_icmap(3,i).eq.charmm_dtype(3,k).and. &
            charmm_icmap(4,i).eq.charmm_dtype(4,k)) .or. (charmm_icmap(4,i).eq.charmm_dtype(1,k).and.charmm_icmap(3,i).eq.charmm_dtype(2,k).and. &
            charmm_icmap(2,i).eq.charmm_dtype(3,k).and.charmm_icmap(1,i).eq.charmm_dtype(4,k))
    enddo
    if (.not.ok) call error ('readcharmm', 'dihedral angle is not found in CMAP section', faterr)
    if (ok) charmm_icmap2(1,i)  = k
    ! second angle: psi angle (N-C_alpha-C-N)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,itype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,jtype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,ktype)
    call getfirst(com,word)
    call fatnam(charmm_label,1,chmmntype,word,ltype)
    if (abs(charmm_mass(itype)-14.007).lt.tol .and. abs(charmm_mass(ltype)-14.007).lt.tol) then 
      ok = abs(charmm_mass(jtype)-12.011).lt.tol.and.abs(charmm_mass(ktype)-12.011).lt.tol
      if (.not.ok) call error ('readcharmm', 'Wrong atom types in CMAP section', faterr)
    else
      call error ('readcharmm', 'Wrong atom types in CMAP section', faterr)
    endif
    charmm_icmap(5,i) = itype
    charmm_icmap(6,i) = jtype
    charmm_icmap(7,i) = ktype
    charmm_icmap(8,i) = ltype
    ! check dihedral angle 
    ok = .false.
    k = 0
    do while (k.lt.chmmdih .and. .not.ok)
      k = k + 1
      ok = (charmm_icmap(5,i).eq.charmm_dtype(1,k).and.charmm_icmap(6,i).eq.charmm_dtype(2,k).and.charmm_icmap(7,i).eq.charmm_dtype(3,k).and. &
            charmm_icmap(8,i).eq.charmm_dtype(4,k)) .or. (charmm_icmap(8,i).eq.charmm_dtype(1,k).and.charmm_icmap(7,i).eq.charmm_dtype(2,k).and. &
            charmm_icmap(6,i).eq.charmm_dtype(3,k).and.charmm_icmap(5,i).eq.charmm_dtype(4,k))
    enddo
    if (.not.ok) call error ('readcharmm', 'dihedral angle is not found in CMAP section', faterr)
    if (ok) charmm_icmap2(2,i)  = k
    ! grid points and grid spacing
    call getfirst(com,word)
    charmm_ncmap(i) = chr2int(word) ! grid points
    charmm_cmap(i) = 360.0/charmm_ncmap(i) ! grid spacing
    k = charmm_ncmap(i)/5
    if (k*5.lt.charmm_ncmap(i)) k = k + 1
    do j = 1, charmm_ncmap(i)
      m = charmm_ncmap(i)*(j-1)
      n = 0
      do l = 1, k
         call getlin(com,iunff,outu)
         if (l.lt.k) then
           do p = 1, 5
             n = n + 1
             call getfirst(com,word)
             charmm_fcmap(m+n,i) = chr2real(word)
           enddo
         else
           o = charmm_ncmap(i) - n
           do p = 1, o
             n = n + 1
             call getfirst(com,word)
             charmm_fcmap(m+n,i) = chr2real(word)
           enddo
         endif
      enddo
    enddo
  enddo  ! next i
endif ! Qchmmcmap
! SECTION H: NONBONDED
! The first number is ignored, the second term is the well-depth (epsilon) and the third
! term is the (minimum radius)/2. A second set of 3 numbers may be specified to
! indicate the VDW parameters to be used for the calculation of 1-4 nonbonded interactions
i = 0
do itype = 1, chmmntype
  do jtype = itype, chmmntype
    i = i + 1
    charmm_typen(itype,jtype) = i
    charmm_typen(jtype,itype) = i
  enddo
enddo
call getlin(com,iunff,outu)
do i = 1, chmmntype
  call getlin(com,iunff,outu)
  call getfirst(com,word) 
  call fatnam(charmm_label,1,chmmntype,word,itype) ! atom type
  j = charmm_typen(itype,itype) ! nonbonded pair
  call getfirst(com,word) ! ignored
  call getfirst(com,word)
  charmm_nonbonded(1,j) = chr2real(word) ! well-depth (epsilon)
  if (charmm_nonbonded(1,j).gt.0.0) call error ('readcharmm', 'Wrong value for well-depth (epsilon) in NONBONDED section', faterr)
  charmm_nonbonded(1,j) = abs(charmm_nonbonded(1,j))
  call getfirst(com,word)
  charmm_nonbonded(2,j) = 2.0*cte*chr2real(word) ! minimum radius
  if (charmm_nonbonded(2,j).lt.0.0) call error ('readcharmm', 'Wrong value for minimum radius in NONBONDED section', faterr)
  call getfirst(com,word) ! ignored (if 1-4 nonbonded interactions)
  if (len_trim(word).ne.0) then ! 1-4 nonbonded interactions
    call getfirst(com,word)
    charmm_nonbonded(3,j) = chr2real(word) ! well-depth (epsilon)
    if (charmm_nonbonded(3,j).gt.0.0) call error ('readcharmm', 'Wrong value for well-depth (epsilon) in NONBONDED section', faterr)
    charmm_nonbonded(3,j) = abs(charmm_nonbonded(3,j))
    call getfirst(com,word)
    charmm_nonbonded(4,j) = 2.0*cte*chr2real(word) ! minimum radius
    if (charmm_nonbonded(4,j).lt.0.0) call error ('readcharmm', 'Wrong value for minimum radius in NONBONDED section', faterr)
  else
    charmm_nonbonded(3,j) = charmm_nonbonded(1,j) 
    charmm_nonbonded(4,j) = charmm_nonbonded(2,j)
  endif
enddo ! next i
! assign default for off diagonal
do itype = 1, chmmntype-1
  i = charmm_typen(itype,itype)
  do jtype = itype+1, chmmntype
    j = charmm_typen(jtype,jtype)
    k = charmm_typen(itype,jtype)
    charmm_nonbonded(1,k) = sqrt(charmm_nonbonded(1,i)*charmm_nonbonded(1,j))
    charmm_nonbonded(2,k) = 0.5*(charmm_nonbonded(2,i)+charmm_nonbonded(2,j))
    charmm_nonbonded(3,k) = sqrt(charmm_nonbonded(3,i)*charmm_nonbonded(3,j))
    charmm_nonbonded(4,k) = 0.5*(charmm_nonbonded(4,i)+charmm_nonbonded(4,j))
  enddo
enddo
! SECTION I: NBFIX
if (Qchmmnbfix) then
  call getlin(com,iunff,outu)
  do i = 1, chmmnbfix
    call getlin(com,iunff,outu)
    call getfirst(com,word) 
    call fatnam(charmm_label,1,chmmntype,word,itype) ! first atom type
    call getfirst(com,word) 
    call fatnam(charmm_label,1,chmmntype,word,jtype) ! second atom type
    if (itype.eq.jtype) call error ('readcharmm', 'Same types in NBFIX section', warning)
    j = charmm_typen(itype,jtype)
    call getfirst(com,word)
    charmm_nonbonded(1,j) = chr2real(word) ! well-depth (epsilon)
    if (charmm_nonbonded(1,j).gt.0.0) call error ('readcharmm', 'Wrong value for well-depth (epsilon) in NBFIX section', faterr)
    charmm_nonbonded(1,j) = abs(charmm_nonbonded(1,j))
    call getfirst(com,word)
    charmm_nonbonded(2,j) = cte*chr2real(word) ! minimum radius
    if (charmm_nonbonded(2,j).lt.0.0) call error ('readcharmm', 'Wrong value for minimum radius in NBFIX section', faterr)
    call getfirst(com,word)
    if (len_trim(word).ne.0) then ! 1-4 nonbonded interactions
      charmm_nonbonded(3,j) = chr2real(word) ! well-depth (epsilon)
      if (charmm_nonbonded(3,j).gt.0.0) call error ('readcharmm', 'Wrong value for well-depth (epsilon) in NBFIX section', faterr)
      charmm_nonbonded(3,j) = abs(charmm_nonbonded(3,j)) 
      call getfirst(com,word)
      charmm_nonbonded(4,j) = cte*chr2real(word) ! minimum radius
      if (charmm_nonbonded(4,j).lt.0.0) call error ('readcharmm', 'Wrong value for minimum radius in NBFIX section', faterr)    
    else
      charmm_nonbonded(3,j) = charmm_nonbonded(1,j) 
      charmm_nonbonded(4,j) = charmm_nonbonded(2,j)
    endif 
  enddo ! next i
endif ! Qchmmnbfix

if (Qprint) then
  ! Write outputfile
  write(outu,'(a)') "Atom types"
  write(outu,'(a)') "Label---mass"
  do i = 1, chmmntype
    write(wrtline,*) charmm_label(i),charmm_mass(i)
    write(outu,'(a)') trim(adjustl(wrtline))
  enddo
  if (Qchmmbond) then
    write(outu,'(a)') "Bond types"
    write(outu,'(a)') "Label1---Label2---Kb[Kcal/mole/A**2]---b0[A]"
    do i = 1, chmmbond
      itype = charmm_btype(1,i)
      jtype = charmm_btype(2,i)
      write(wrtline,*) charmm_label(itype),charmm_label(jtype),charmm_bond(1,i),charmm_bond(2,i)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif ! Qchmmbond
  if (Qchmmang) then
    write(outu,'(a)') "Bond angle and UB types"
    write(outu,'(a)') "Label1---Label2---Label3---Ktheta[Kcal/mole/rad**2]---Theta0[degrees]---Kub[Kcal/mole/A**2]---S0[A]"
    do i = 1, chmmang
      itype = charmm_atype(1,i) 
      jtype = charmm_atype(2,i) ! central atom
      ktype = charmm_atype(3,i) 
      j = charmm_lub(i)
      if (j.gt.0) then 
        write(wrtline,*) charmm_label(itype),charmm_label(jtype),charmm_label(ktype),charmm_ang(1,i),charmm_ang(2,i),charmm_ub(1,j),charmm_ub(2,j)
      else
        write(wrtline,*) charmm_label(itype),charmm_label(jtype),charmm_label(ktype),charmm_ang(1,i),charmm_ang(2,i)  
      endif
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif ! Qchmmang
  if (Qchmmdih) then
    write(outu,'(a)') "Dihedral angle types"
    write(outu,'(a)') "Label1---Label2---Label3---Label4---(Kchi[Kcal/mole]---n(multiplicity)---delta[degrees])..."
    do i = 1, chmmdih
      itype = charmm_dtype(1,i) ! terminal atom
      jtype = charmm_dtype(2,i)
      ktype = charmm_dtype(3,i)
      ltype = charmm_dtype(4,i) ! terminal atom
      if (itype.gt.0) then
        atnam1 = charmm_label(itype)
      else
        atnam1 = 'X'
      endif
      if (ltype.gt.0) then
        atnam2 = charmm_label(ltype)
      else
        atnam2 = 'X'
      endif
      write(wrtline,*) atnam1,charmm_label(jtype),charmm_label(ktype),atnam2, & 
                       (charmm_dih((j-1)*2+1,i),charmm_ndih(j,i),charmm_dih((j-1)*2+2,i),j=1,charmm_nprms(i))
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif ! Qchmmdih
  if (Qchmmimp) then
    write(outu,'(a)') "Improper angle types"
    write(outu,'(a)') "Label1---Label2---Label3---Label4---Kpsi[Kcal/mole/rad**2]---psi0[degrees]"
    do i = 1, chmmimp
      itype = charmm_itype(1,i)
      jtype = charmm_itype(2,i) 
      ktype = charmm_itype(3,i)    
      ltype = charmm_itype(4,i)  
      if (itype.gt.0) then
        atnam1 = charmm_label(itype)
      else
        atnam1 = 'X'
      endif
      if (jtype.gt.0) then
        atnam2 = charmm_label(jtype)
      else
        atnam2 = 'X'
      endif
      if (ktype.gt.0) then
        atnam3 = charmm_label(ktype)
      else
        atnam3 = 'X'
      endif
      if (ltype.gt.0) then
        atnam4 = charmm_label(ltype)
      else
        atnam4 = 'X'
      endif
      write(wrtline,*) atnam1,atnam2,atnam3,atnam4,charmm_imp(1,i),charmm_imp(2,i)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif ! Qchmmimp
  if (Qchmmcmap) then
    write(outu,'(a)') "CMAP types"
    write(outu,'(a)') "Label1---Label2---Label3---Label4---dihedral type---Label5---Label6---Label7---Label8---dihedral type---Number of grid points---Grid spacing"
    do i = 1, chmmcmap
      ! theta angle (C-N-C_alpha-C)
      itype = charmm_icmap(1,i)
      jtype = charmm_icmap(2,i)
      ktype = charmm_icmap(3,i)
      ltype = charmm_icmap(4,i)
      dih1  = charmm_icmap2(1,i) ! type dihedral angle
      ! psi angle (N-C_alpha-C-N)
      mtype = charmm_icmap(5,i)
      ntype = charmm_icmap(6,i)
      otype = charmm_icmap(7,i)
      ptype = charmm_icmap(8,i)
      dih2  = charmm_icmap2(2,i) ! type dihedral angle
      write(wrtline,*) charmm_label(itype),charmm_label(jtype),charmm_label(ktype),charmm_label(ltype),dih1,' ', &
                       charmm_label(mtype),charmm_label(ntype),charmm_label(otype),charmm_label(ptype),dih2,' ', &
                       charmm_ncmap(i),charmm_cmap(i)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  endif ! Qchmmcmap
  write(outu,'(a)') "NONBONDED terms"
  write(outu,'(a)') "Label1---Label2---epsilon[Kcal/mole]---sigma[A]---epsilon1,4[Kcal/mole]---sigma1,4[A]"
  do itype = 1, chmmntype
    do jtype = itype, chmmntype
      i = charmm_typen(itype,jtype)
      write(wrtline,*) charmm_label(itype),charmm_label(jtype),charmm_nonbonded(1,i),charmm_nonbonded(2,i),charmm_nonbonded(3,i),charmm_nonbonded(4,i)
      write(outu,'(a)') trim(adjustl(wrtline))
    enddo
  enddo
endif ! Qprint
 
return
end subroutine
