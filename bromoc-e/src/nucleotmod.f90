!    BROMOC-E
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! New variables related to DNA insertion in GCMC/BD code
! ======================================================
!
! Strands, nucleotides and interaction sites
! ------------------------------------------
!
! istrs                 -> Number of strands
! inuc                  -> Number of nucleotides in each strand
! cgnuc                 -> Nucleotide charge
! diffnuc               -> Nucleotide diffusion constant
! epsnuc                -> Nucleotide Lennard-Jones potential parameter 
! Qnucl                 -> Logical variable which indicates if nucleotide order
!                          have been called
! Qdie                  -> Logical variable which indicates if effective dielectric 
!                          constant for DNA have been called
! Qsolv                 -> Logical variable which indicates if solvent-induced 
!                          contributions have been called
! Qatexp                -> Logical variable which indicates if atoms order
!                          have been called
! Qassign               -> Logical variable which indicates if assign order
!                          have been called
! Qpar                  -> Logical variable which indicates if particles order
!                          have been called
! Qsystem               -> Logical variable which indicates if system order
!                          have been called
! Qbuf                  -> Logical variable which indicates if buffers order
!                          have been called
! Qtraj                 -> Logical variable which indicates if DNA and/or ions 
!                          trajectories have to be written in a outputfile
! Qfmemb                -> Logical variable which indicates if traslocation duration 
!                          for cylindrical pores has to be calculated
! eps(P,S,Ab,Tb,Cb,Gb)  -> Lennard-Jones parameter for interaction between ions 
!                          and sites (combination rules)
! sg(P,S,Ab,Tb,Cb,Gb)   -> Lennard-Jones parameter for interaction between ions
!                          and sites (combination rules)
! epsolv                -> Energy scale for solvnet-induced contribution
! maxsite               -> Maximum number of interaction sites
! maxpar                -> Maximum number of interaction pairs
! strand(maxsite)       -> This integer vector indicates the strand for the 
!                          interaction site
! typenuc(maxsite)       -> This integer vector indicates the nucleotide type 
!                          for the interaction site
! namnucl(maxsite)      -> This character vector store the name nucleotides
!                          (A,T,C,G) 
! namsite(maxsite)      -> This character vector stores the name interaction 
!                          sites P=Phosphate, S=Sugar, (Ab,Tb,Cb,Gb)=Bases
! (x,y,z,r,phi)nat(maxsite) -> Interaction site coordinates for the native 
!                              structure
! Qtras                 -> Logical variable which indicates if there is a traslation of DNA 
!                          sites positions
! Qrot                  -> Logical variable which indicates if there is a rotation of DNA
!                          sites positions
!                        -> These integer vectors are used in wrtrraj routine
! nparnuc               -> nucleotides in npar 
! nelenuc               -> nucleotides in nele
! netnuc                -> nucleotides in netyp
! nptnuc                -> nucleotides in nptyp

module nucleotmod
implicit none
integer     istrs, inuc, maxsite, extraP,nelenuc1st
integer,allocatable     ::  strand(:), typenuc(:)
character*1,allocatable :: namnucl(:)
character*2,allocatable :: namsite(:)
real        cgnuc, diffnuc, epsnuc, epsolv, fctn, scalepairing
real        notrx,notry,notrz,inelenuc  ! NOTRAN
integer     setframes
!CONTRA
real        xcon,ycon,zcon
real,allocatable :: kx(:),ky(:),kz(:),contrx(:),contry(:),contrz(:) 
integer     ctn
integer,allocatable ::  csn(:)
logical*1   Qcontrans,Qcontprint,Qunsplit
logical*1   Qnucl, Qassign, Qpar, Qsystem, Qbuf, Qtraj, Qtrajcont, Qdie, Qsolv, Qfmemb
logical*1   Qtras, Qrot, Qnotrans, Qnotrx, Qnotry, Qnotrz, Qdnafree, Qinvstr, QfirstP
real cylall(6,3),din,ain

! Bonded terms
! ------------
!
! maxbond               -> Maximum number of bonds
! nbond                 -> Number of bonds
! sitebond(maxbond,2)   -> This integer matrix indicates the bonded 
!                          interaction sites
! distbond(maxbond)     -> This real vector stores the equilibrium 
!                          bond lenghts
! bond(maxpar)          -> Logical supervector which indicates if two sites 
!                          are forming a bond
! typbond(maxbond)      -> This integer vector indicates if the bond 
!                          is either intranucleotide (0) or 
!                          internucleotide (1)
!
! maxang                -> Maximum number of bond angles
! nangle                -> Number of bond angles
! siteangle(maxang,3)   -> This integer matrix indicates the interaction 
!                          sites which are forming a bond angle
! valangle(maxang)      -> This real vector stores the equilibrium bond
!                          angles
! angle(maxpar)         -> Logical supervector which indicates if two sites 
!                          are forming abond angle
!
! maxdihe               -> Maximum number of dihedral angles
! ndihe                 -> Number of dihedral angles
! sitedihe(maxdihe,4)   -> This integer matrix indicates the interaction
!                          sites which are forming a dihedral angle
! valdihe(maxdihe)      -> This real vector stores the equilibrium dihedral
!                          angles
! (dSAb,dSTb,dSCb,dSGb,dPS5,dPS3) -> Natural bond lenghts
! (phPSAb,phPSTb,phPSCb,phPSGb,phPSAb2,phPSTb2,phPSCb2,phPSGb2,phSPS,phPSP) ->
!                                    Natural bond angles
! (dhAbSPS,dhTbSPS,dhCbSPS,dhGbSPS,dhSPSAb,dhSPSTb,dhSPSCb,dhSPSGb,dhSPSP,dhPSPS) -> 
!                                    Natural dihedral angles

integer   maxbond, maxang, maxdihe  
integer   nbond, nangle, ndihe
integer,allocatable ::   sitebond(:,:),siteangle(:,:),sitedihe(:,:)
real,allocatable ::      distbond(:),valangle(:),valdihe(:)
logical*1,allocatable :: bond(:),angle(:)
integer,allocatable ::   typbond(:)
real      dSAb, dSTb, dSCb, dSGb, dPS5, dPS3
real      phPSAb, phPSTb, phPSCb, phPSGb, phPSAb2, phPSTb2, phPSCb2, phPSGb2, phSPS, phPSP
real      dhAbSPS, dhTbSPS, dhCbSPS, dhGbSPS, dhSPSAb, dhSPSTb, dhSPSCb, dhSPSGb, dhSPSP, dhPSPS 

! Non-bonded terms
! ----------------
!
! nstack                -> Number of native contacts
! nbp                   -> Number of hydrogen bondings
! nex                   -> Number of excluded volume interactions
! nqq                   -> Number of Coulomb interactions
! nsolv                 -> Number of solvent-induced contributions
! sitestack(maxpar,2)   -> This integer matrix indicates the interaction
!                          sites which are forming a native contact
! sitebp(maxpar,2)      -> This integer matrix indicates the interaction
!                          sites which are forming a hydrogen bonding
! siteex(maxpar,2)      -> This integer matrix indicates the interaction
!                          sites which are forming a excluded volume interaction
! siteqq(maxpar,2)      -> This integer matrix indicates the interaction
!                          sites which are forming a Coulomb interaction
! siteslv(maxpar,2)     -> This integer matrix indicates the interaction
!                          sites which are forming a solvent-induced contributions
! sgstack(maxpar)       -> sigma parameter in intra-strand native contacts term
! sgbp(maxpar)          -> sigma parameter in hydrogen bonding terms
! sgex(maxpar)          -> sigma parameter in excluded volume

integer           maxpar
integer           nstack, nbp, nex, nqq, nsolv
integer,allocatable ::           sitestack(:,:),sitebp(:,:),siteex(:,:),siteqq(:,:),siteslv(:,:)
real,allocatable ::              sgstack(:),sgbp(:),sgex(:)

! Potential energy contributions
! ------------------------------
! ebond    -> strech energy
! eang     -> bending energy
! edihe    -> torsional energy
! estack   -> native contacts
! ebp      -> hydrogen bonding
! eex      -> excluded volume
! eqq      -> coulomb interaction
! esolv    -> solvent-induced contribution
! eqqmx    -> coulomb interaction between ions and sites
! evdwmx   -> LJ potential between ions and sites
! Qdeby    -> Logical variable which indicates if Coulomb
!             interactions are taken into account using the
!             Debye-Hückel approximation 
! Qionsite -> Logical variable which indicates if there are 
!             ions-sites interactions using combination rules 
!             for obtaining LJ parameters
! Qljpar   -> Logical variable which indicates if there are
!             ions-sites interactions using LJ parameters
!             directly (i.e, don't using combination rules)
! ionstr   -> Ionic strength [Mol/L]

real ebond, eang, edihe
real estack, ebp, eex, eqq, esolv, econ
real eqqmx, evdwmx
logical*1 Qdeby, Qljsin, Qljpar, Qninfo, Qdebyhyb
real  ionstr

! Fraction of denatured bases
! ---------------------------

! Qfbases -> Logical variable which indicates if fraction of denatured bases will be
!            calculated after simulation of DNA
! fbases  -> Fraction of denatured bases

real fbases
logical*1 Qfbases

end module 
