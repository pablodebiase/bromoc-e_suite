! New variables related to PRM, PSF and DCD reading routines
! ==========================================================
!
module charmmmod
implicit none

! PRM routine
! -----------
!
! SECTION B: ATOMS
! chmmntype             -> Number of CHARMM atom types
! charmm_label          -> CHARMM atom type names
! charmm_mass           -> CHARMM atom type masses
!
! SECTION C: BONDS
! V(bond) = Kb(b-b0)**2
! Kb: Kcal/mole/A**2
! b0 : A
! Qchmmbond             -> Logical variable which indicates if there are CHARMM bond types
! chmmbond              -> Number of CHARMM bond types
! charmm_btype          -> CHARMM atom types related to each CHARMM bond type
! charmm_bond           -> Kb and b0 parameters
!
! SECTION D: ANGLES AND UREY-BRADLEY
! V(angle) = Ktheta(Theta-Theta0)**2
! Ktheta: Kcal/mole/rad**2
! Theta0: degrees
! V(Urey-Bradley) = Kub(S-S0)**2
! Kub: Kcal/mole/A**2
! S0: A
! Qchmmang              -> Logical variable which indicates if there are CHARMM bond angle types
! Qchmmub               -> Logical variable which indicates if there are CHARMM Urey-Bradley types
! chmmang               -> Number of CHARMM bond angle types
! chmmub                -> Number of CHARMM Urey-Bradley types
! charmm_atype          -> CHARMM atom types related to each CHARMM bond angle type
! charmm_ubtype         -> CHARMM atom types related to each CHARMM Urey-Bradley type
! charmm_ang            -> Ktheta and Theta0 parameters
! charmm_ub             -> Kub and S0 parameters
! charmm_lub            -> Integer variable which indicates the CHARMM Urey-Bradley type
!                          associated with a given CHARMM bond angle type 
!                          (0 value means that there is not a CHARMM Urey-Bradley type)
!
! SECTION E: DIHEDRALS
! V(dihedral) = Kchi(1+cos(n(chi)-delta))
! Kchi: Kcal/mole
! n: multiplicity
! delta : degrees
! Qchmmdih              -> Logical variable which indicates if there are CHARMM dihedral angle types
! chmmdih               -> Number of CHARMM dihedral angle types
! charmm_dtype          -> CHARMM atom types related to each CHARMM dihedral angle type
! charmm_dih            -> Kchi and delta parameters
! charmm_ndih           -> n parameters
! charmm_nprms          -> Number of terms (multiple dihedral angles)a
! charmm_nmax           -> Maximum number of terms
!
! SECTION F: IMPROPER
! V(improper) = Kpsi(psi-psi0)**2
! Kpsi: Kcal/mole/rad**2
! psi0: degrees
! Ordinarily, improper dihedrals are given a multiplicity of 0, which imposes a harmonic restoring potential
! instead of a cosine function.  In this case, the central atom must be either the first or the last atom
! Qchmmimp              -> Logical variable which indicates if there are CHARMM improper angle types
! chmmimp               -> Number of CHARMM improper angle types
! charmm_itype          -> CHARMM atom types related to each CHARMM improper angle type
! charmm_imp            -> Kpsi and psi0 parameters
!
! SECTION G: CMAP
! Qchmmcmap             -> Logical variable which indicates if there are Cross-term energy correction map types
! chmmcmap              -> Number of Cross-term energy correction map types
! charmm_icmap          -> CHARMM atom types related to each CMAP type
! charmm_icmap2         -> CHARMM dihedral angle types related to each CMAP type
! charmm_ncmap          -> Number of grid points for each CMAP type
! charmm_cmap           -> Grid spacing for each CMAP type 
! charmm_fcmap          -> Energy correction map for each CMAP type
!
! SECTION H: NONBONDED
! chmmnonb              -> Number of nonbonded pair types
! charmm_typen          -> Pointer for nonbonded pairs
! charmm_nonbonded      -> Nonbonded parameters
!
! SECTION I: NBFIX
! Qchmmnbfix            -> Logical variable which indicates if there are VDW interactions between specific atom pair types to be modified
! chmmnbfix             -> Number of VDW interactions between specific atom pair types to be modified
!
! SECTION G: HBOND
! Qchmmhbond            -> Logical variable which indicates if there are hydrogen bond types
! chmmhbond             -> Number of hydrogen bond types

integer   chmmntype,chmmbond,chmmang,chmmub,chmmdih,chmmimp,chmmcmap,chmmnonb,chmmnbfix,chmmhbond,charmm_nmax  
integer,allocatable :: charmm_btype(:,:),charmm_atype(:,:),charmm_ubtype(:,:),charmm_dtype(:,:),charmm_itype(:,:),charmm_icmap(:,:),charmm_icmap2(:,:)
integer, allocatable :: charmm_lub(:),charmm_ndih(:,:),charmm_nprms(:),charmm_ncmap(:),charmm_typen(:,:)
real,allocatable :: charmm_mass(:),charmm_bond(:,:),charmm_ang(:,:),charmm_ub(:,:),charmm_dih(:,:),charmm_imp(:,:),charmm_cmap(:),charmm_fcmap(:,:)
real, allocatable :: charmm_nonbonded(:,:)
character*7,allocatable :: charmm_label(:)
logical*1 Qchmmbond,Qchmmang,Qchmmub,Qchmmdih,Qchmmimp,Qchmmcmap,Qchmmnbfix,Qchmmhbond

! PSF routine
! -----------
!
! psf_atomtype, psf_atomtype2 -> Integer pointers for atomtypes
! psf_charge                  -> Integer pointer for atomtypes
! viscwat                     -> Water viscosity [Pa s]
! scldiff                     -> Scale factor for self-diffusion calculation

integer, allocatable :: psf_atomtype(:), psf_atomtype2(:), psf_charge(:,:)
real viscwat, scldiff

! DCD routine
! -----------
!
! DCDUnitCell           -> Logical variable which indicates if unit cell should be read in DCD file

logical*1 DCDUnitCell

end module 
