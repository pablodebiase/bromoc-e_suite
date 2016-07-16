program tt2
use listmod
implicit none
integer n,i,j,k,ne
character*4,allocatable,dimension(:) :: enam,pnam

! Allocate Particle Type List for 3 types
ntyp=3
call resizeptypl(ntyp)
allocate (enam(9),pnam(3))
! etyp
enam(1)='K'
enam(2)='H'
enam(3)='C'
enam(4)='O'
enam(5)='N'
enam(6)='Cl'
enam(7)='S'
enam(8)='P'
enam(9)='Na'
pnam(1)='POT'
pnam(2)='GLU'
pnam(3)='CLA'
! Define Particle 1: POT
n=1
ne=1 ! number of elements per particle
ptypl(n)%ne = ne 
allocate (ptypl(n)%etyp(ne),ptypl(n)%r(ne))
ptypl(n)%etyp(1) = 1
call setcar(ptypl(n)%r(1),0.0,0.0,0.0)
! Define Particle 2: GLU
n=2
ne=4
ptypl(n)%ne = ne 
allocate (ptypl(n)%etyp(ne),ptypl(n)%r(ne))
ptypl(n)%etyp(1) = 2
ptypl(n)%etyp(2) = 3
ptypl(n)%etyp(3) = 4
ptypl(n)%etyp(4) = 5
call setcar(ptypl(n)%r(1),0.0,0.0,0.0)
call setcar(ptypl(n)%r(2),1.0,0.0,0.0)
call setcar(ptypl(n)%r(3),0.0,1.0,0.0)
call setcar(ptypl(n)%r(4),0.0,0.0,1.0)
! Define Particle 3: CLA
n=3
ne=1 ! number of elements per particle
ptypl(n)%ne = ne
allocate (ptypl(n)%etyp(ne),ptypl(n)%r(ne))
ptypl(n)%etyp(1) = 6
call setcar(ptypl(n)%r(1),0.0,0.0,0.0)

call resizeparl(100)
call resizecvec(1000)

! Add the following particle types to list
call addpar(1) 
call addpar(2) 
call addpar(1) 
call addpar(3) 
call addpar(2) 
call addpar(1) 
call addpar(3) 
call addpar(3) 
call addpar(1) 
call addpar(1) 
call addpar(2) 

! Randomly move particles and rotate
do i=1,npar
   call movepar(i,getcar((rand()-0.5)*10.0,(rand()-0.5)*10.0,(rand()-0.5)*10.0))
   call uranrot(i)
enddo

call printlists()

contains

subroutine printxyz()
! Print Particle Data
call updatetypel()
write(*,*) nele
write(*,*)
do i=1,nele
    write(*,*) enam(etypl(i)),r(i)%x,r(i)%y,r(i)%z
enddo
end subroutine

subroutine printlists()
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
