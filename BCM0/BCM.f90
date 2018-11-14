subroutine BCM(nele,natom,tau,lat_matrix,ele_num,q)

  implicit none
  integer(kind=4)             :: i, j, k, ix, iy, iz, i1, i2, i3
  integer(kind=4)             :: l, m, n, j1, j2, j3
  integer(kind=4)             :: nele, natom
  integer(kind=4)             :: ele_num(nele), super_ele_num(nele), counter(nele,nele)
  real(kind=8)                :: lat_matrix(3,3), pi
  real(kind=8)                :: eps,eps1,tmp,temp(3)
  real(kind=8)                :: rms,x,y,z,r,theta,phi
  real(kind=8)                :: pos (3,400*natom,nele), mo(nele,nele), tau1(3,natom), tau(3,natom)
  real(kind=8)                :: dis(nele,nele), q(11,nele,nele)
  complex                     :: Yx(nele,nele),sumY(nele,nele)
  logical                     :: cluster

  pi = 3.1415926

  eps = 1e-1
  eps1 = 1e-4
  Yx  = (0.0,0.0)
  sumY = (0.0,0.0)
  mo = 0.0

  if (cluster) then
      ix = 0
      iy = 0
      iz = 0
  else
      ix = 1
      iy = 1
      iz = 1
  end if
  super_ele_num = 0
  k=0

  do i=1,nele
     do j=1,ele_num(i)
        k=k+1
        do i1 = -ix, ix
           do i2 = -iy, iy
              do i3 = -iz, iz
                 super_ele_num(i) = super_ele_num(i) + 1

                 temp(1) = tau(1,k) + i1
                 temp(2) = tau(2,k) + i2
                 temp(3) = tau(3,k) + i3
                 pos(1,super_ele_num(i),i)=temp(1)*lat_matrix(1,1) + temp(2)*lat_matrix(2,1) +temp(3)*lat_matrix(3,1)
                 pos(2,super_ele_num(i),i)=temp(1)*lat_matrix(1,2) + temp(2)*lat_matrix(2,2) +temp(3)*lat_matrix(3,2)
                 pos(3,super_ele_num(i),i)=temp(1)*lat_matrix(1,3) + temp(2)*lat_matrix(2,3) +temp(3)*lat_matrix(3,3)
              end do
           end do
        end do
     end do
  end do

  k=0
    do i=1,nele
       do j=1,ele_num(i)
       k=k+1
        tau1(1,k)=tau(1,k)*lat_matrix(1,1) + tau(2,k)*lat_matrix(2,1) + tau(3,k)*lat_matrix(3,1)
        tau1(2,k)=tau(1,k)*lat_matrix(1,2) + tau(2,k)*lat_matrix(2,2) + tau(3,k)*lat_matrix(3,2)
        tau1(3,k)=tau(1,k)*lat_matrix(1,3) + tau(2,k)*lat_matrix(2,3) + tau(3,k)*lat_matrix(3,3)
       end do
    end do
  dis=100.0
  do l = 0, 10
     mo = 0.0
     do m = -l, l
        sumY = (0.0, 0.0)
        counter = 0
        k = 0
        do i1=1, nele
           do j1 = 1,ele_num(i1)
              k = k + 1
              do i2=1,nele
                 do j2=1,super_ele_num(i2)
                    x = tau1(1,k) - pos(1,j2,i2)
                    y = tau1(2,k) - pos(2,j2,i2)
                    z = tau1(3,k) - pos(3,j2,i2)
                    call car2sphe(x,y,z,r,theta,phi)
                    if (r>eps .and. r<3.0 ) then
                       counter(i1,i2) = counter(i1,i2) + 1
                       call spherical(l,m,theta,phi,Yx(i1,i2))
                       sumY(i1,i2) = sumY(i1,i2) + Yx(i1,i2)
                    end if
                 end do
              end do
           end do
        end do

        do i1=1,nele
           do i2=1,nele
              if (counter(i1,i2)==0) then
                  sumY(i1,i2) = 0.0
              else
                  sumY(i1,i2) = sumY(i1,i2)/real(counter(i1,i2))
              end if
              mo(i1,i2)=mo(i1,i2)+cabs(sumY(i1,i2))**2
           end do
        end do
     end do

     do i1=1,nele
        do i2=1,nele
           q(l+1,i1,i2) = sqrt ( mo(i1,i2) * 4 * pi / ( 2 * l + 1 ) )
        end do
     end do
  end do

  do l=1,11
     if(mod(l,2) == 0) then
     q(l,:,:) = 0
     else
     end if
  end do


  end subroutine BCM
