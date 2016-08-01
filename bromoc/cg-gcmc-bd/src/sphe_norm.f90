!    BROMOC  -  CG-GCMC-BD
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

subroutine sphe_norm(nmpol,bnorm,srdist)
!-----------------------------------------------------------------------
!calculate the normalization constants of the spherical harmonics basis functions
!input: SRDIST -> radius sphere
!       NMPOL  -> number of multipoles
!output: BNORM -> Normalization constant

use constamod
implicit none
integer nmpol
real  bnorm(*),srdist
!local
integer l,m,norder
real  sr2,lpart,rpart,upfacto,dnfacto,factori,ipi,itwopi

!always the same order in spherical harmonics
norder=1
sr2=srdist*srdist
rpart=sr2*srdist

itwopi=1.0/twopi
lpart=1.5*itwopi
ipi=0.5*itwopi
bnorm(norder)=sqrt(lpart/rpart)
do l=1,nmpol-1
   lpart=(2.0*l+1.0)*(2.0*l+3.0)*ipi
   rpart=rpart*sr2
   norder=norder+1
   bnorm(norder)=sqrt(lpart/rpart)
   lpart=(2.0*l+1.0)*(2.0*l+3.0)*itwopi       ! change for m > 0
   do m=1,l
      norder=norder+1 
      upfacto=factori(l-m)
      dnfacto=factori(l+m)
      bnorm(norder)=sqrt(lpart*upfacto/(rpart*dnfacto))
      norder=norder+1 
      bnorm(norder)=bnorm(norder-1)
   enddo
enddo

return
end subroutine
