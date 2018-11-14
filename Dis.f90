
program main

 implicit none

  integer(kind=4)             :: i, j, nele, natom, Ntype, l, model, ls,nn, delete
  integer(kind=4)             :: ele_num(1)
  real(kind=8)   ::lat_matrix(3,3),E,q(11,1,1),tau(3,8)
  real(kind=8)    ::tau2(3,8,101,20)
  character(100)              :: c_ele_num, filename
  

open(100,file=trim('POSCAR1'))

  read(100,*)
  read(100,*)
  do i = 1, 3
      read(100,*) (lat_matrix(i,j),j = 1, 3)
      !Be very careful with lat_matrix, vector a is lat_matrix(1,:)
  end do
print*,lat_matrix
!  read(100,*)(symb(i), i = 1, nele)
  read(100,*)(ele_num(i),i = 1,1)
  read(100,*)
  natom = sum(ele_num)




  do i = 1, natom
     read(100,*)(tau(j,i),j=1,3)
!tau(:,1) is the coordination of atom 1
  end do

  close(100)

print*,tau,natom

call disturbance(tau,tau2,lat_matrix,natom)








end program main
