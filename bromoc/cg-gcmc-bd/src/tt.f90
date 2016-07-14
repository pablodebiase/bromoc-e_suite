program tt
type :: car
    real :: x
    real :: y
    real :: z
end type car

! Particle Type
type :: partype
    integer :: ne                              !! Number of Elements in Particle
    integer,allocatable,dimension(:) :: etyp   !! Particle Element Types Vector :: Size of ne
    type(car),allocatable,dimension(:) :: r    !! Position Vector of Particle Elements :: Size of ne
end type partype

type(car),allocatable,dimension(:) :: r
type(partype),allocatable,dimension(:) :: ptypl !! part type list
integer i,j,k,np
allocate (r(100))
np=3
allocate (ptypl(np))
do i=1,np
   ptypl(i)%ne=i*10
  allocate (ptypl(i)%r(i*10),ptypl(i)%etyp(i*10))
enddo

k=0
do i=1,np
  do j=1,ptypl(i)%ne
    ptypl(i)%etyp(j)=k
    ptypl(i)%r(j)%x=k*0.1
    ptypl(i)%r(j)%y=k*0.2
    ptypl(i)%r(j)%z=k*0.3
    k=k+1
  enddo
enddo
do i=1,100
   r(i)%x = i*0.1
   r(i)%y = i*0.2
   r(i)%z = i*0.3
enddo

do i=1,100
   write(*,*) r(i)%x,r(i)%y,r(i)%z
enddo

do i=1,np
  do j=1,ptypl(i)%ne
    write(*,*) i,np,j,ptypl(i)%ne,ptypl(i)%etyp(j),ptypl(i)%r(j)%x,ptypl(i)%r(j)%y,ptypl(i)%r(j)%z
  enddo
enddo

np=np+1
call increaseptypl(np)
i=35
ptypl(np)%ne=i
allocate (ptypl(np)%r(i),ptypl(np)%etyp(i))

do j=1,ptypl(np)%ne
  ptypl(np)%etyp(j)=k
  ptypl(np)%r(j)%x=0.1
  ptypl(np)%r(j)%y=0.2
  ptypl(np)%r(j)%z=0.3
enddo

do i=1,np
  do j=1,ptypl(i)%ne
    write(*,*) i,np,j,ptypl(i)%ne,ptypl(i)%etyp(j),ptypl(i)%r(j)%x,ptypl(i)%r(j)%y,ptypl(i)%r(j)%z
  enddo
enddo

contains 

subroutine increaseptypl(newsize)
implicit none
integer vdim,newsize
type(partype), allocatable, dimension(:) :: ptypltmp    !! Type of Particle list :: Size of ntyp

! Measure Vector Size
vdim=size(ptypl)

! Allocate double of original size in temp vectors
allocate (ptypltmp(newsize))
! Copy Vectors
ptypltmp(1:vdim)=ptypl
! Deallocate original Vector
deallocate (ptypl)
! Move allocations
call move_alloc(ptypltmp,ptypl)
end subroutine

end

