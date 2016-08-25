! New variables related to explicit atoms in GCMC/BD code
! ======================================================
!
module explatmod
implicit none

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

integer     maxtypes, maxcharges, maxnnb
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
integer natfx, nghosts
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

real eub, eopbs, ecmap

! 1-2, 1-3 and 1-4 terms lists
!-----------------------------
!
! listmex               -> Number of atoms for 1-2 and 1-3 terms list
! listex                -> Pointer for atoms belonging to 1-2 and 1-3 terms list
! listm14               -> Number of atoms for 1-4 terms list
! list14                ->  Pointer for atoms belonging to 1-4 terms list

integer, allocatable :: listmex(:), listm14(:), listex(:,:), list14(:,:)

end module 
