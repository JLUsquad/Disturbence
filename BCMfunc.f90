
subroutine getpara(filename,nele,natom)
implicit none
integer(kind=4)   ::i,nele,natom
integer(kind=4),allocatable ::ele_num(:)
character(100)    ::filename,c_ele_num


  open(100,file=trim(filename))
!  print* ,trim(filename)
  do i=1, 5
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

  do i = 1, 5
    read(100,*)
  end do

  read(100,*)ele_num(:)
  read(100,*)
  natom = sum(ele_num)
  close(100)

end subroutine

subroutine fingerprint(filename,nele,natom,q,model)
  implicit none  

  integer(kind=4)             :: i, j, nele, natom, Ntype, l, model, ls, nn, delete
  parameter (ls=1,nn=56)
  integer(kind=4)             :: ele_num(nele),atomid(nele+1)
  real(kind=8)                ::lat_matrix(3,3),E,q(11,nele,nele),fp(ls*(nele+1),natom),xx(nn,1,natom),tau(3,natom),cfr(nele*400)
  character(100)              :: c_ele_num, filename
  character(len=2), dimension(nele) :: symb

!  write(*,*)trim(filename) 
  open(100,file=trim(filename))
!  read(100,*)E
!  read(100,*)(symb(i), i = 1, nele)
!  read(100,*)E
  read(100,*)
  read(100,*)
  do i = 1, 3
      read(100,*) (lat_matrix(i,j),j = 1, 3)
      !Be very careful with lat_matrix, vector a is lat_matrix(1,:)
  end do
!  read(100,*)(symb(i), i = 1, nele)
  read(100,*)(ele_num(i),i = 1,nele)
  read(100,*)
  natom = sum(ele_num)
    
  
 
  
  do i = 1, natom
     read(100,*)(tau(j,i),j=1,3)
!tau(:,1) is the coordination of atom 1
  end do

  close(100)

 
  

!  do i= 1, nele+1
!   if(i==1)then
!   atomid(i)=1
!   endif
!   atomid(i)=ele_num(i-1)+1
!  enddo

  q = 0
  fp = 0
  xx = 0
delete=0

 ! print*,tau,lat_matrix,ele_num,q

  call BCM(nele,natom,tau,lat_matrix,ele_num,q)

!if(delete == 1)then
!print* ,filename
!else
!endif

end subroutine fingerprint





subroutine BCM(nele,natom,tau,lat_matrix,ele_num,q)

  implicit none
  integer(kind=4)             :: i, j, k, ix, iy, iz, i1, i2, i3, k1
  integer(kind=4)             :: l, m, n, j1, j2, j3
  integer(kind=4)             :: nele, natom
  integer(kind=4)             :: ele_num(nele), super_ele_num(nele), counter(nele,nele)
  real(kind=8)                :: lat_matrix(3,3), pi
  real(kind=8)                :: eps,eps1,tmp,temp(3)
  real(kind=8)                :: rms,x,y,z,r,theta,phi,maxr,rp(400*natom),rpmin(natom,nele)
  real(kind=8)                :: pos (3,400*natom,nele), mo(nele,nele), tau1(3,natom), tau(3,natom)
  real(kind=8)                :: dis(nele,nele), q(11,nele,nele)
  complex                     :: Yx(nele,nele),sumY(nele,nele)


  pi = 3.1415926

  eps = 1e-1
  eps1 = 1e-4 
  Yx  = (0.0,0.0)
  sumY = (0.0,0.0)
  mo = 0.0

  ix = 1
  iy = 1
  iz = 1
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
        tau1(1,k)=tau(1,k)*lat_matrix(1,1) + tau(2,k)*lat_matrix(2,1) +tau(3,k)*lat_matrix(3,1)
        tau1(2,k)=tau(1,k)*lat_matrix(1,2) + tau(2,k)*lat_matrix(2,2) +tau(3,k)*lat_matrix(3,2)
        tau1(3,k)=tau(1,k)*lat_matrix(1,3) + tau(2,k)*lat_matrix(2,3) +tau(3,k)*lat_matrix(3,3)
       end do
    end do
  dis=100.0

k=0
!-----------------
rpmin(:,:)=10000.0
rp = 0
       do i1=1, nele
           do j1 = 1,ele_num(i1)
              k = k + 1
              k1 = 0
              rp = 0
              do i2=1,nele
                 do j2=1,super_ele_num(i2)
                    k1 = k1 + 1
                    x = tau1(1,k) - pos(1,j2,i2)
                    y = tau1(2,k) - pos(2,j2,i2)
                    z = tau1(3,k) - pos(3,j2,i2)
                    rp(k1)=sqrt(x*x+y*y+z*z)
                    if (rp(k1) .gt. eps)then
                   
                    if (rpmin(k,i2) .gt. rp(k1))then
                    rpmin(k,i2) = rp(k1) 
                    else
                    endif
                    else
                    endif
                 enddo
              enddo
           enddo
       enddo
do i=1,natom
rpmin(i,:)=rpmin(i,:)+1e-2
enddo

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
                    r = sqrt(x*x+y*y+z*z)
                    if (r>eps .and. r<rpmin(k,i2) ) then
                       call car2sphe(x,y,z,r,theta,phi)
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

  subroutine dist_BCM(Distance,q1,q2,nele)
  implicit none
  integer(kind=4) ::l,i,j,nele,Ntype
  real(kind=8)    ::q1(11,nele,nele),q2(11,nele,nele)
  real(kind=8)    ::sumdif,Distance
  
  
  
  sumdif = 0
  do l = 1, 11
     do i = 1, nele
        do j = i, nele
           sumdif = sumdif + (q1(l,i,j)-q2(l,i,j))**2
        end do
     end do
  end do

  Ntype = nele*(nele-1)/2 + nele
  Distance = sqrt(sumdif/Ntype)
end subroutine dist_BCM

