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

subroutine help (lus)
!.....help - prints out a quick reference guide in the lus logical unit.
implicit none
integer lus
write (lus,500)
write (lus,510)
!.....Formats:
500 format (/1x,'To run:    bromoc [input [output]]'/1x,'bromoc -h gives help'//) 
510 format (1x, 'A manual is provided with source code.'/1x,'For more information:'/ &
1x,'cjfqct@hotmail.com'/1x,'snoskov@ucalgary.ca'/ 1x,'pablodebiase@gmail.com')
end subroutine
