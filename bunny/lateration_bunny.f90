program lateration_bunny
use lateration
implicit none


   type(mesh) :: amesh
   real, dimension(:), allocatable :: dist
   logical :: fast
   
! fast option on
   fast=.true.
! load the Stanford bunny mesh in ascii vtk
   call loadbunny(amesh)
! memory allocation and initialization for distance
   allocate(dist(amesh%Nnodes))
   dist=infinity
! source is on the back of the bunny
   dist(3300)=0.
! computation
   write(0,*) 'starting distance computing'
   call onevsall2d(amesh,dist,fast)
! dumping the result in ascii VTK
   write(*,*) 'dumping distance '
   open(10,file='result.vtk',form = 'formatted')
   call dumpmeshvtk(10,amesh)
   call dumpnodeattributevtk(10,amesh,dist,'distance',.true.)
   close(10)
end program lateration_bunny
!#########################################################################
! load the bunny vtk file
! original version but the nodes not associated to the mesh have been removed.
subroutine loadbunny(amesh)
use lateration
implicit none
  type(mesh) :: amesh

  integer :: i,j,idum

  amesh%Nnodes=34834
  allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
  open(12,file="bunny.vtk")
  do i=1,5
     read(12,*)
  enddo
  do i=1,amesh%Nnodes
     read(12,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  amesh%Ncells=69451
  allocate(amesh%cell(amesh%Ncells,3))
  read(12,*)
  do i=1,amesh%Ncells
     read(12,*) idum,(amesh%cell(i,j),j=1,3)
  enddo
  amesh%cell=amesh%cell+1
  return
end subroutine loadbunny
!###############################################################################
