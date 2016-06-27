!    PBEQ2VTK - Converts CHARMM PBEQ Binary Potential Map file format 
!               into VTK (Visualization Toolkit) format 
!    Copyright (C) 2015 Pablo M. De Biase (pablodebiase@gmail.com)
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

module gsbp
implicit none
integer*4 nclx, ncly, nclz, ncyz, mini
real*8  dcel, xbcen, ybcen, zbcen
real*8  epsw, epspp, conc, tmemb, zmemb1, epsm
real*8 tranx, trany, tranz, idcel
integer*4 ncel3
real*4,allocatable ::  phi(:)
end module

program pbeq2vtk
implicit none
integer*4 nsc
integer narg,arg
character inpfile*256,outfile*256,line*256,bs*8
real*8 low,up
logical*1 inpopen

bs=repeat(achar(8),len(bs))

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('CHARMM PBEQ Binary Potential Map (.pbeq) filename: ',narg,arg,inpfile)
open(1,file=inpfile,form='unformatted')
call readphi(1,6)
close(1)

! 5 stdin, 6 stdout, 0 stderr
call readarg('Virtualization Toolkit format (Binary VTK) Output filename: ',narg,arg,outfile)

open (unit=2,file=outfile,form='unformatted',access='stream',action='write',position='rewind')
call writevtk(2,6)
close(2)

write(*,'(/A)') 'Normal termination'
end program

! read arguments
subroutine readarg(ques,narg,num,text)
implicit none
integer*4 narg,num
character text*(*),ques*(*)

num=num+1
write(*,'(/A$)') ques
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
  write(*,'(A)') trim(text)
else
  text=''
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='PBEQ2VTK'
prver='version 1.0'
prdesc='Converts CHARMM PBEQ Binary Potential Map file format into VTK (Visualization Toolkit) format'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='25 Feb 2015'
lastdate='25 Feb 2015'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

subroutine readphi(unit,outu)
!-----------------------------------------------------------------------
! read INPUT static external field PHIX or grid-based repulsion potential PHIV
!            and miscelaneous parameters
use gsbp
!Input
implicit none
integer*4 unit, outu
integer*4 i
integer*4 ifir, ilas

read(unit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
read(unit) epsw,epspp,conc,tmemb,zmemb1,epsm
tranx = 0.5*(nclx-1)*dcel
trany = 0.5*(ncly-1)*dcel
tranz = 0.5*(nclz-1)*dcel
ncel3  = nclx*ncly*nclz
idcel = 1.0/dcel
ncyz=ncly*nclz

!Writting in output file            
write(outu,*) 'Reading PBEQ ...'
write(outu,*) 'Number of grid point in X   (nclx) = ',nclx 
write(outu,*) 'Number of grid point in Y   (ncly) = ',ncly 
write(outu,*) 'Number of grid point in Z   (nclz) = ',nclz 
write(outu,*) 'Grid spacing                (dcel) = ',dcel
write(outu,*) 'Center of box in X          (xbcen)= ',xbcen
write(outu,*) 'Center of box in Y          (ybcen)= ',ybcen
write(outu,*) 'Center of box in Z          (zbcen)= ',zbcen
write(outu,*) 
write(outu,*) 'Solvent dielectric constant (epsw) = ',epsw
write(outu,*) 'Protein dielectric constant (epsp) = ',epspp
write(outu,*) 'Salt concentration          (conc) = ',conc
if (tmemb.gt.0.0) then
  write(outu,*)
  write(outu,*) 'Membrane thickness along Z  (tmemb)= ',tmemb
  write(outu,*) 'Membrane position along Z   (zmemb)= ',zmemb1
  write(outu,*) 'Membrane dielectric constant(epsm) = ',epsm
endif
write(outu,*)
write(outu,*) 'Box in X from ',xbcen-tranx,' to ',xbcen+tranx
write(outu,*) 'Box in Y from ',ybcen-trany,' to ',ybcen+trany
write(outu,*) 'Box in Z from ',zbcen-tranz,' to ',zbcen+tranz

!Grid-based repulsion potential        
ifir = 1
ilas = ncel3
if (allocated(phi)) deallocate (phi)
allocate(phi(ifir:ilas))
read(unit) (phi(i),i=ifir,ilas)
return
end subroutine

subroutine writevtk(unit,outu)
use gsbp
implicit none
integer*4 unit,outu,i
integer*1 newline
character line*256

newline=10
write(outu,*) 'Writing Binary VTK ...'
write(unit) '# vtk DataFile Version 3.0'
write(unit) newline
write(unit) 'File created using PBEQ2VTK from BROMOC Tools'
write(unit) newline
!write(unit) 'ASCII' 
write(unit) 'BINARY'
write(unit) newline
write(unit) 'DATASET STRUCTURED_POINTS'
write(unit) newline
write(line,'(A,I0,x,I0,x,I0)') 'DIMENSIONS ',nclx,ncly,nclz
write(unit) trim(line)
write(unit) newline
write(line,*) 'SPACING ',dcel,dcel,dcel
write(unit) adjustl(trim(line))
write(unit) newline
write(line,*) 'ORIGIN ',xbcen,ybcen,zbcen
write(unit) adjustl(trim(line))
write(unit) newline
write(line,'(A,I0)') 'POINT_DATA ',ncel3
write(unit) trim(line)
write(unit) newline
write(unit) 'SCALARS density float 1'  ! check this
write(unit) newline
write(unit) 'LOOKUP_TABLE default'     ! check this
write(unit) newline
! for ASCII
!do i=1,ncel3 
!  write(line,*) phi(i)
!  write(unit) adjustl(trim(line))
!  write(unit) newline
!enddo
write(unit) phi(1:ncel3)         ! check this

end subroutine
