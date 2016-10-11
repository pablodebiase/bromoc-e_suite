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

subroutine splinepot(is,nn,xx,yy,kf)
use efpmod
use grandmod
implicit none
integer nn,is,i
real xx(nn),yy(nn),d1,d2,dd1,dd2,kf
real p(3,nn)

efp(is)%xl=xx(1)
efp(is)%xl2=xx(1)**2
efp(is)%xu2=xx(nn)**2
if (.not.Qcol(is)) kf=0.0
! Shifting potential
if(xx(nn).eq.0.0.or.xx(1).eq.0.0) stop 'Division by zero'
efp(is)%sc%c=(xx(nn)*yy(nn)-kf)*xx(nn)**5
d2=-kf/xx(nn)**2-6.0*efp(is)%sc%c/xx(nn)**7
dd2=2.0*kf/xx(nn)**3+42.0*efp(is)%sc%c/xx(nn)**8
dd1=1.0/(xx(1)**12-xx(2)**12)
efp(is)%sc%a=(yy(2)-yy(1))*xx(2)**12*xx(1)**12*dd1
efp(is)%sc%b=(yy(1)*xx(1)**12-yy(2)*xx(2)**12)*dd1
d1=-12.0*efp(is)%sc%a/xx(1)**13
dd1=156.0*efp(is)%sc%a/xx(1)**14
call squarespline(nn,xx,yy,p,d1,d2,dd2)
do i=1,nn
   efp(is)%ep(i)%a = p(1,i) 
   efp(is)%ep(i)%b = p(2,i) 
   efp(is)%ep(i)%c = p(3,i) 
enddo
end subroutine
