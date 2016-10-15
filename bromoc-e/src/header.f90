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

!-----------------------------------------------------------------------
subroutine header (lus)
implicit none
integer lus
!.....header - prints out a presentation of the program in the lus
!logical unit.
write (lus,500)
!.....Formats:
500  format ( &
 1x,'*************************************************************************************'/ &
 1x,'*******                              BROMOC-E                                 *******'/ &
 1x,'*******        Simulation Program for Electrodiffusion, Gran Canonical        *******'/ &
 1x,'*******        Monte Carlo, Brownian Dynamics, Coarse Grain Model DNA         *******'/ &
 1x,'*******        and CHARMM force field reader for Solute and Solvent           *******'/ &
 1x,'*******        (c) 2016 Pablo M. De Biase,   Superluminal Technologies Ltd.   *******'/ &
 1x,'*******                       Version 5.00 (2016/10/11)                       *******'/ &
 1x,'*************************************************************************************'/ &
 1x,'*******                                                                       *******'/ & 
 1x,'******* Copyright (C) 2016 Pablo M. De Biase (pablodebiase@gmail.com)         *******'/ &
 1x,'*******                                                                       *******'/ &
 1x,'******* This program is free software: you can redistribute it and/or modify  *******'/ &
 1x,'******* it under the terms of the GNU General Public License as published by  *******'/ &
 1x,'******* the Free Software Foundation, either version 3 of the License, or     *******'/ &
 1x,'******* (at your option) any later version.                                   *******'/ &
 1x,'*******                                                                       *******'/ &
 1x,'******* This program is distributed in the hope that it will be useful,       *******'/ &
 1x,'******* but WITHOUT ANY WARRANTY; without even the implied warranty of        *******'/ &
 1x,'******* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *******'/ &
 1x,'******* GNU General Public License for more details.                          *******'/ &
 1x,'*******                                                                       *******'/ &
 1x,'******* You should have received a copy of the GNU General Public License     *******'/ &
 1x,'******* along with this program.  If not, see <http://www.gnu.org/licenses/>. *******'/ &
 1x,'*******                                                                       *******'/ &
 1x,'************************************************************************************'//)
end

