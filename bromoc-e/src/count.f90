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

subroutine count()
!     Count all the ions and assign them to their appropriate buffer
!     nbuffer: the number of buffers in the system
!     nele: number of elements in the system
!     buffer number ib concerns ions of type ibfftyp(ib)
!     particle i belongs to the buffer number ibuf 
!     (note that if ibuf is zero, then particle does not belong to any buffer
!      or is a fix ion)
use grandmod
use constamod
use stdiomod
use errormod
use nucleotmod
use listmod
implicit none

integer ib,itype,i
real  radius2
logical*1 endlog
type(car) :: rr

!     Initializations
do ib = 1, nbuffer
   nat(ib) = 0 
enddo
do i = nparnuc+1, npar
  if (parl(i)%kind.ne.3) cycle
! To which buffer does a particle of type (i) at this location belongs?
  parl(i)%ibuf = 0
  itype = parl(i)%ptyp
  endlog = .false.
  ib = 1
  do while (ib.le.nbuffer .and. .not.endlog)
    if (itype.eq.ibfftyp(ib)) then
      call getcentroid(i,rr)
      if (rr%z.ge.LZmin(ib) .and. rr%z.le.LZmax(ib)) then
        if (Qsphere) then
          radius2 = (rr%x-cx)**2+(rr%y-cy)**2+(rr%z-cz)**2
          if (radius2.gt.Rsphe2) call error ('count', 'particles outside main system', faterr)      
          if ((radius2.ge.Rmin(ib)**2).and.(radius2.le.Rmax(ib)**2)) then
            parl(i)%ibuf = ib
            nat(ib) = nat(ib) + 1
            endlog = .true.
          endif
        else
          if (Qecyl) then
            radius2=((rr%x-cx)*iecx)**2+((rr%y-cy)*iecy)**2
            if (radius2.gt.1.0) call error ('count', 'particle/s outside main system', faterr)
          else
            if (rr%x.lt.lx2m.or.rr%x.gt.lx2p) call error ('count', 'particles outside main system', faterr)
            if (rr%y.lt.ly2m.or.rr%y.gt.ly2p) call error ('count', 'particles outside main system', faterr)
          endif
          if (rr%z.lt.lz2m.or.rr%z.gt.lz2p) call error ('count', 'particles outside main system', faterr)
          parl(i)%ibuf = ib
          nat(ib) = nat(ib) + 1
          endlog = .true.
        endif
      endif
    endif
    ib = ib + 1
  enddo
enddo

do ib = 1, nbuffer
  ntotat(ib) = ntotat(ib) + nat(ib)
enddo

return
end
