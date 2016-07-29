program listmodtester
use listmod
implicit none
integer i,k,npnp
real,external :: rndm

! initialize lists
call deltypall()
call delparall()

! empty lists
call resizeetypl(0)
call resizeptypl(0)
call resizeparl(0)
call resizecvec(0)

! Add element and particle types
! Define Particle 1 (POT) and Element 1 (K)
call addmonopar('K','POT')  ! 1
call addmonopar('Cl','CLA') ! 2
call addmonopar('Na','SOD') ! 3
! Add element types
call addetyp('H')   ! 4
call addetyp('C')   ! 5
call addetyp('O')   ! 6
call addetyp('N')   ! 7
call addetyp('S')   ! 8 
call addetyp('P')   ! 9
! Add particle types
! Define Particle 2: GLU
call addptyp(4,'GLU')
call seteleinptyp(nptyp,1,4)
call seteleinptyp(nptyp,2,5)
call seteleinptyp(nptyp,3,6)
call seteleinptyp(nptyp,4,7)
call putcoorinptyp(nptyp,1,0.0,0.0,0.0)
call putcoorinptyp(nptyp,2,1.0,0.0,0.0)
call putcoorinptyp(nptyp,3,0.0,1.0,0.0)
call putcoorinptyp(nptyp,4,0.0,0.0,1.0)
! Define Particle 4: CCC
call addptyp(4,'CSCP')
call seteleinptyp(nptyp,1,5)
call seteleinptyp(nptyp,2,8)
call seteleinptyp(nptyp,3,5)
call seteleinptyp(nptyp,4,9)
call putcoorinptyp(nptyp,1,0.0,0.0,0.0)
call putcoorinptyp(nptyp,2,1.0,0.0,0.0)
call putcoorinptyp(nptyp,3,0.0,1.0,0.0)
call putcoorinptyp(nptyp,4,0.0,0.0,1.0)

! Add the following particle types to list
do i=1,20
    call addpar(getrand(nptyp))
enddo

! Randomly move particles and rotate
do i=1,npar
   call movepar(i,getcar((rand()-0.5)*10.0,(rand()-0.5)*10.0,(rand()-0.5)*10.0))
   call uranrot(i)
enddo

write(*,*) '>>>>>>>>>> LISTS'
call printlists(6)

write(*,*) '>>>>>>>>>> XYZ'
call printxyz(6)

write(*,*) '>>>>>>>>>> PDB'
call printpdb(6)

write(*,*) '>>>>>>>>>>>>',0,3,npar
call delpar(3)
call printpdb(6)

npnp=npar+1
do i=1,npnp
  k=getrand(npar)
  write(*,*) '>>>>>>>>>>>>',i,k,npar
  call delpar(k)
  call printpdb(6)
  !call printlists(6)
enddo

contains

function getrand(imax)
implicit none
integer getrand,imax
getrand=int(imax*rndm()+1)
end function

end program
