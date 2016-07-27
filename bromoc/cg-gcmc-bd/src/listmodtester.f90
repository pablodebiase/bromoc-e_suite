program listmodtester
use listmod
implicit none
integer i,k,npnp
real,external :: rndm

! initialize lists
call deltypall()
call delparall()

! empty lists
call resizeenam(0)
call resizeptypl(0)
call resizeparl(0)
call resizecvec(0)

! Add element types
call addenam('K')   ! 1
call addenam('H')   ! 2
call addenam('C')   ! 3
call addenam('O')   ! 4
call addenam('N')   ! 5
call addenam('Cl')  ! 6
call addenam('S')   ! 7 
call addenam('P')   ! 8
call addenam('Na')  ! 9
! Add particle ypes
! Define Particle 1: POT
call addptyp(1,'POT')
ptypl(nptyp)%etyp(1) = 1
call setcarzero(ptypl(nptyp)%r(1))
! Define Particle 2: GLU
call addptyp(4,'GLU')
ptypl(nptyp)%etyp(1) = 2
ptypl(nptyp)%etyp(2) = 3
ptypl(nptyp)%etyp(3) = 4
ptypl(nptyp)%etyp(4) = 5
call setcarzero(ptypl(nptyp)%r(1))
call setcar(ptypl(nptyp)%r(2),1.0,0.0,0.0)
call setcar(ptypl(nptyp)%r(3),0.0,1.0,0.0)
call setcar(ptypl(nptyp)%r(4),0.0,0.0,1.0)
! Define Particle 3: CLA
call addptyp(1,'CLA')
ptypl(nptyp)%etyp(1) = 6
call setcarzero(ptypl(nptyp)%r(1))
! Define Particle 4: CCC
call addptyp(4,'CSCP')
ptypl(nptyp)%etyp(1) = 3
ptypl(nptyp)%etyp(2) = 7
ptypl(nptyp)%etyp(3) = 3
ptypl(nptyp)%etyp(4) = 8
call setcarzero(ptypl(nptyp)%r(1))
call setcar(ptypl(nptyp)%r(2),1.0,0.0,0.0)
call setcar(ptypl(nptyp)%r(3),0.0,1.0,0.0)
call setcar(ptypl(nptyp)%r(4),0.0,0.0,1.0)
! Define Particle 5: SOD
call addptyp(1,'SOD')
ptypl(nptyp)%etyp(1) = 9
call setcarzero(ptypl(nptyp)%r(1))

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
