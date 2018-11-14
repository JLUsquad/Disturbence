subroutine BCMcustom(filename,E,Distance,q0)

!In this function POSCAR is read, and what to do with the BCM q matrix is setted

  implicit none

  integer(kind=4)             :: i, j, nele, natom, Ntype, l
  integer(kind=4),allocatable :: ele_num(:)
  real(kind=8)                :: lat_matrix(3,3),q0(11,1,1),E,Distance,sumdif
  real(kind=8),allocatable    :: q(:,:,:),tau(:,:)
  character(100)              :: c_ele_num, filename
  character(3)                :: c1, c2, c3
  

  open(100,file=trim(filename))
  do i=1, 6
  read(100,*)
  end do
  read(100,'(a100)')c_ele_num
  nele = 0
  do i=1, 100
      if(c_ele_num(i:i) == " ") then
      else
        if(c_ele_num(i-1:i-1) == " ") then
        nele = nele + 1
        else
        end if
      end if
  end do
  rewind(100)
  
  allocate(ele_num(nele))
  allocate(q(11,nele,nele))
  Ntype = nele*(nele-1)/2 + nele

  read(100,*)E
  read(100,*)
  read(100,*)

  do i = 1, 3
      read(100,*) (lat_matrix(i,j),j = 1, 3)
  end do

  read(100,*)ele_num(:)
  read(100,*)
  natom = sum(ele_num)

  allocate(tau(3,natom))

  do i = 1, natom
     read(100,*)(tau(j,i),j=1,3)
  end do
  close(100)
  


  call BCM(nele,natom,tau,lat_matrix,ele_num,q)
  
  sumdif = 0
  do l = 1, 11
     do i = 1, nele
        do j = i, nele
           sumdif = sumdif + (q(l,i,j)-q0(l,1,1))**2
        end do
     end do
  end do

  Distance = sqrt(sumdif/Ntype)

end subroutine BCMcustom
