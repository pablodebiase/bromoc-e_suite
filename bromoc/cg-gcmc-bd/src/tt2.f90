program tt2
use listmod
implicit none
integer n,i,j,k,ne
! Allocate Particle Type List for 3 types
ntyp=3
call resizeptypl(ntyp)
! etyp
! 1: K
! 2: H
! 3: C
! 4: O
! 5: N
! 6: Cl
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

! Print Particle Data
do i=1,npar
   write(*,*) i,parl(i)%ptyp,parl(i)%sr,parl(i)%ne
   do j=1,parl(i)%ne
      k=j+parl(i)%sr
      write(*,*) '      ',ptypl(parl(i)%ptyp)%etyp(j),r(k)%x,r(k)%y,r(k)%z
   enddo
enddo
end program
