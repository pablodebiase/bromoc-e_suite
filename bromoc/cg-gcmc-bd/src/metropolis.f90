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

subroutine metropolis(nmcm)
! Perform nmcm of metropolis monte carlo with constant number of particle
use grandmod
use constamod
use listmod
implicit none
integer nmcm
integer parn, imove, ne, sr, pne
!real*16 bltz,eold, enew
real bltz,eold, enew
type(car) :: rr
type(car),allocatable,dimension(:) :: rrori

pne = 0
do imove = 1, nmcm
   !pick one atom and new shift to position
   call move(parn,rr)
   ! backup previous position
   ne = parl(parn)%ne
   sr = parl(parn)%sr
   if (allocated(rrori).and.ne.ne.pne) deallocate (rrori)
   if (.not.allocated(rrori)) allocate (rrori(ne))
   rrori = r(sr+1:sr+ne)
   !calculate old energy
   call par_interact(parn, eold)
   !attempt move and rotate
   call movepar(parn,rr,center=.true.)
   !calculate new energy
   call par_interact(parn, enew)
   if (enew.le.eold) then
      !accept the move
      !ener = ener + (enew-eold)
   else                            
      bltz = exp(-(enew-eold)*ikbt)
      if (rndm().lt.bltz) then
         !accept the move
         !ener = ener + (enew-eold)
      else
         ! restore original position
         r(sr+1:sr+ne) = rrori
      endif
   endif 
   pne = ne
enddo
end subroutine
