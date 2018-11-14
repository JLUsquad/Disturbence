program main 
  !
  implicit none
  !
  !local variables
  !
  !
  integer(kind=4)             :: i, j, k, loop, loop0, ix, iy, iz, i1, i2, i3,c
  integer(kind=4)             :: l, m, n, j1, j2, j3
  integer(kind=4)             :: nele, natom, Ntype
  integer(kind=4),allocatable :: ele_num(:), super_ele_num(:), counter(:,:)
  real(kind=8)                :: lat_matrix(3,3), pi
  real(kind=8)                :: eps,eps1,tmp,temp(3),sumdif, Distance(101), averdis(20), dev(20)
  real(kind=8)                :: rms,x,y,z,r,theta,phi
  real(kind=8),allocatable    :: pos (:,:,:), mo(:,:) ,tau0(:,:), tau2(:,:,:,:), tau1(:,:), tau(:,:)
  real(kind=8),allocatable    :: dis(:,:), q(:,:,:,:)
  complex,allocatable         :: Yx(:,:),sumY(:,:)
  logical                     :: cluster
  character(100)              :: c_ele_num
  character                   :: l_cluster
 
 
  call getarg(1,l_cluster)
  if (l_cluster == "t") then
       cluster = .true.
       write(*,*)"cluster = true"
  else
       cluster = .false.
       write(*,*)"cluster = false"
  end if

  nele = 0
  open(100,file="POSCAR")
  do i=1, 5
  read(100,*)
  end do
  read(100,'(a100)')c_ele_num
  do i=1, 100
      if(c_ele_num(i:i) == " ") then
      else 
        if(c_ele_num(i-1:i-1) == " ") then
        nele= nele + 1
        else
        end if
      end if
  end do
  rewind(100)
  write(*,*)"nele = ",nele

  pi = 3.1415926
  lat_matrix = 0.0
  Ntype = nele*(nele-1)/2 + nele

  allocate(ele_num(nele))
  allocate(Yx(nele,nele))
  allocate(sumY(nele,nele))
  allocate(mo(nele,nele))
  allocate(q(11,nele,nele,101))
  allocate(super_ele_num(nele))
  allocate(counter(nele,nele))
  allocate(dis(nele,nele))

  read(100,*)
  read(100,*)
  do i = 1, 3
      read(100,*) (lat_matrix(i,j),j = 1, 3)
  end do   
  read(100,*)ele_num(:)  
  read(100,*)
  natom = sum(ele_num)
  allocate(tau0(3,natom))
  allocate(tau(3,natom))
  allocate(tau1(3,natom))
  allocate(tau2(3,natom,101,40))
  allocate(pos(3,400*natom,nele))

  do i = 1, natom
     read(100,*)(tau0(j,i),j=1,3)
  end do
  close(100)
 

  call disturbance(tau0,tau2,lat_matrix,natom)


do c = 30, 40
do loop = 1, 101
  do i = 1, natom
    do j = 1, 3
       tau(j,i) = tau2(j,i,loop,c)
    end do
  end do  


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
      ix = 3
      iy = 3
      iz = 3
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
  
  i=0
    do i=1,nele
       do j=1,ele_num(i)
       k=k+1
        tau1(1,k)=tau(1,k)*lat_matrix(1,1) + tau(2,k)*lat_matrix(2,1) +tau(3,k)*lat_matrix(3,1)
        tau1(2,k)=tau(1,k)*lat_matrix(1,2) + tau(2,k)*lat_matrix(2,2) +tau(3,k)*lat_matrix(3,2)
        tau1(3,k)=tau(1,k)*lat_matrix(1,3) + tau(2,k)*lat_matrix(2,3) +tau(3,k)*lat_matrix(3,3)
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
                    if (r>eps .and. r<10.0 ) then
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
           q(l+1,i1,i2,loop) = sqrt ( mo(i1,i2) * 4 * pi / ( 2 * l + 1 ) )      
        end do
     end do
  end do

  do l=1,10
     if(mod(l,2) == 0) then
     q(l,:,:,loop) = 0
     else
     end if
  end do

  sumdif = 0
  do l = 1, 11
     do i = 1, nele
        do j = i, nele  
           !sumdif = sumdif + q(l,i,j,loop)**2
           sumdif = sumdif + (q(l,i,j,loop) - q(l,i,j,1))**2
        end do
     end do
  end do
  
  Distance(loop) = sqrt(sumdif/Ntype)         !Distance(1) = 0
  
end do

   
  open(10,file = 'data.txt')
  averdis(c) = sum(Distance)/100
  dev(c) = 0
  do i = 2, 101
    dev(c) = dev(c) + (Distance(i) - averdis(c))**2
  end do
  dev(c) = sqrt(dev(c)/100)
  write(*,*)c*0.05,averdis(c),dev(c)
  write(10,*)c*0.05,averdis(c),dev(c)
end do
close(10)


 
!  same_flag=.false.
!  if ( nsimilar == 0 ) then
!     nsimilar = nsimilar + 1
!     StruFaclist1(:,:,:,nsimilar)=q
!     open (1003,file="results/similar.dat")
!     write(1003,*)nsimilar
!     do i=1,nsimilar
!        write(1003,'("struct",1x,I5)')i
!        do i1=1,nele
!           do i2=i1,nele
!              write(1003,"(2x,11f12.6)")(StruFaclist1(j+1,i1,i2,i),j=0,10)
!           end do
!        end do
!     end do
!     close(1003)
!     return
!  else
!     loop1:  do i=1,nsimilar
!        rms=0.0
!        do i1=1,nele
!           do i2=i1,nele
!              do j=0,10,2
!                 rms = rms+(StruFaclist1(j+1,i1,i2,i)-q(j+1,i1,i2))**2               
!              end do
!           end do
!        end do
!        rms = sqrt( rms )
!        if ( rms < eps1 ) same_flag = .true.
!        if ( same_flag ) exit loop1
!     end do loop1
!  end if
!  if (.not.same_flag) then
!     if (nsimilar >= max_num) nsimilar=0
!     nsimilar=nsimilar+1
!     StruFaclist1(:,:,:,nsimilar)=q
!     open (1003,file="results/similar.dat")
!     write(1003,*)nsimilar
!     do i=1,nsimilar
!        write(1003,'("struct",1x,I5)')i
!        do i1=1,nele
!           do i2=i1,nele
!              write(1003,"(2x,11f12.6)")(StruFaclist1(j+1,i1,i2,i),j=0,10)
!           end do
!        end do
!     end do
!     pstruc(nstru)%simf1(:,:,:) = q
!     close(1003)
!  end if
!  !
!  !
!  !
end program main 
!

