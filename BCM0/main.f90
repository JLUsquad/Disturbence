program main

  implicit none

  integer(kind=4)             :: k, a, b, c, nele0, natom0, ele_num0(1)
  real(kind=8)                :: E(5755), Distance(5755), tau0(3,2), lat_matrix0(3,3), q0(11,1,1)
  logical                     :: alive
  character(100)              :: filename
  character(3)                :: c1, c2, c3

  nele0=1
  natom0=2
  tau0(:,1)=0.125
  tau0(:,2)=0.875
  lat_matrix0(1,1)=3.8396000862
  lat_matrix0(1,2)=0
  lat_matrix0(1,3)=0
  lat_matrix0(2,1)=1.9198000431
  lat_matrix0(2,2)=3.3251912150
  lat_matrix0(2,3)=0
  lat_matrix0(3,1)=1.9198000431
  lat_matrix0(3,2)=1.1083970717
  lat_matrix0(3,3)=3.1350203425
  ele_num0(1)=2

  call BCM(nele0,natom0,tau0,lat_matrix0,ele_num0,q0)

k = 0
do a = 1, 20
  write(c1,'(i3)')a
  c1 = adjustl(c1)
  do b = 1, 10
    write(c2,'(i3)')b
    c2 = adjustl(c2)
    open(10,file = 'data'//trim(c1)//trim(c2)//'.txt')
    do c = 1, 100
      write(c3,'(i3)')c
      c3 = adjustl(c3)
      filename=trim('Calypso_opt/calypso'//trim(c1)//'/a_'//trim(c2)//'/'//trim(c3)//'/POSCAR')
      inquire(file=filename,exist=alive)
      if(alive)then
        k = k + 1
        call BCMcustom(filename,E(k),Distance(k),q0)
        if(E(k) > -0.1)then
        else
        write(*,*)E(k),Distance(k),k
        write(10,*)E(k),Distance(k),k
        end if
      else
      end if
    end do
    close(10)
  end do
end do


end program main
