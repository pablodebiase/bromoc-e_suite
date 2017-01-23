subroutine resintfor()
use grandmod
use listmod
implicit none
integer i, j, ne, sr, pne
real ine
type(car) :: cm, netfor, avfor, nettrq, avtrq, trq, rt
type(car),allocatable :: intfor(:),trqfor(:)

pne=0
do j=1,npar
  if (parl(j)%ne.eq.1) cycle
  ne=parl(j)%ne
  ine=1.0/ne
  sr=parl(j)%sr
  ! Compute centroid of particle
  call getcentmass(j,cm)
  ! Compute Net Force
  call setcarzero(netfor)
  do i=1+sr,ne+sr
      call addcar(netfor,f(i))
  enddo
  avfor=netfor
  call mulcar(avfor,ine)
  ! Copy forces
  if (allocated(intfor)) then
      if (pne.ne.ne) then
        deallocate(intfor,trqfor)
        allocate(intfor(ne),trqfor(ne))
      endif
  else
    allocate(intfor(ne),trqfor(ne))
  endif
  do i=1,ne
    intfor(i)=f(i+sr)
  enddo
  ! Substract Average Net Force to each element
  do i=1,ne
    call subcar(intfor(i),avfor)
  enddo
  ! Compute Net Torque
  call setcarzero(nettrq)
  do i=1,ne
    rt=r(i+sr)
    call subcar(rt,cm)
    call carcrossproduct(rt,intfor(i),trq)
    call addcar(nettrq,trq)
  enddo
  ! Compute Average Torque
  avtrq=nettrq
  call mulcar(avtrq,ine)
  ! Get Torque forces
  do i=1,ne
    rt=r(i+sr)
    call subcar(rt,cm)
    call unicar(rt) ! Expensive
    call carcrossproduct(avtrq,rt,trqfor(i))
  enddo
  ! Substract Torque forces
  do i=1,ne
    call subcar(intfor(i),trqfor(i))
  enddo
  ! Scale Down internal Forces
  do i=1,ne
    call mulcar(intfor(i),riffac)
  enddo
  ! Add back average Torque Forces and average displacement Forces
  do i=1,ne
    f(i+sr)%x=intfor(i)%x+trqfor(i)%x+avfor%x
    f(i+sr)%y=intfor(i)%y+trqfor(i)%y+avfor%y
    f(i+sr)%z=intfor(i)%z+trqfor(i)%z+avfor%z
  enddo
  ! Save Last Number of Elements
  pne=ne
enddo

if (allocated(intfor)) deallocate(intfor)

end subroutine

