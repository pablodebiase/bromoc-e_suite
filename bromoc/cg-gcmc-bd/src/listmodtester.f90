program tt2
use listmod
implicit none
integer i,j,k,npnp
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
call printlists()

write(*,*) '>>>>>>>>>> XYZ'
call printxyz()

write(*,*) '>>>>>>>>>> PDB'
call printpdb()

write(*,*) '>>>>>>>>>>>>',0,3,npar
call delpar(3)
call printpdb()

npnp=npar+1
do i=1,npnp
  k=getrand(npar)
  write(*,*) '>>>>>>>>>>>>',i,k,npar
  call delpar(k)
  call printpdb()
  !call printlists()
enddo

contains

function getrand(imax)
implicit none
integer getrand,imax
getrand=int(imax*rndm()+1)
end function

subroutine printpdb()
implicit none
integer i
call updatetypel()
do i=1,nele
  write (*,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,enam(etypl(i)),ptypl(petypl(i))%pnam,pel(i),r(i)%x,r(i)%y,r(i)%z
enddo
write (*,'(A)') 'END'
end subroutine

subroutine printxyz()
implicit none
integer i
! Print Particle Data
call updatetypel()
write(*,*) nele
write(*,*)
do i=1,nele
    write(*,*) enam(etypl(i)),r(i)%x,r(i)%y,r(i)%z
enddo
end subroutine

subroutine printlists()
implicit none
integer i
! Print Particle Data
do i=1,npar
   write(*,*) i,parl(i)%ptyp,parl(i)%sr,parl(i)%ne
   do j=1,parl(i)%ne
      k=j+parl(i)%sr
      write(*,*) '      ',ptypl(parl(i)%ptyp)%etyp(j),r(k)%x,r(k)%y,r(k)%z
   enddo
enddo
end subroutine

end program
