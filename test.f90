program main
  implicit none

  integer(kind=4)  ::i,j,k,a,b,c,nele,natom
  real(kind=8)   :: q(11,1,1),q0(11,1,1),D(101,20),eve(20)
  character(100),dimension(2020) ::filename
  character(100)   ::filenameSi
  character(4)     ::c1

q0=0
filenameSi='POSCAR1'
call getpara(filenameSi,nele,natom)
call fingerprint('POSCAR1',nele,natom,q0,1)
k=0
do i=1,20
do j=1,101
k=k+1
print*,k
write(c1,'(i4)')k
filename(k)='POS'//adjustl(trim(c1))
call getpara(filename(k),nele,natom)
q=0
call fingerprint(filename(k),nele,natom,q,1)
call dist_BCM(D(j,i),q,q0,nele)
print*,D(j,i)
!print*,q
!print*,nele,natom,filename(k)
enddo
eve(i)=sum(D(:,i))/100
enddo
print*,eve
!call dist_BCM(D(1),q0,q,1)
!print* ,'test4'
!call dist_BCM(D(2),q(:,:,:,3),q(:,:,:,2),1)
!print* ,'test5'
!call dist_BCM(D(3),q(:,:,:,1),q(:,:,:,2),1)
!open(100,file = 'SiandGa.txt')

!write(100,*)'Si'
!do i=1,2
!do j=1,2
!write(100,*)q(:,i,j)
!enddo
!enddo
!write(100,*)'D = ',D(1)
!write(100,*)
!write(100,*)'Ga'
!write(100,*)q(:,:,:,2)
!write(100,*)'D = ',D(2)
!write(100,*)'D of the two = ',D(3)
!close(100)


!  open(100,file = 'result01.txt')
!  open(50,file = 'disfrom344')
!  do i = 1, 2857

!  do j = 1, 6
!  if(j == 1)then
!  read(100,*)filenameall(i)
!  elseif(j == 2)then
!  read(100,*)(q0(a,1,1,i),a=1,11)
!  elseif(j == 3)then
!  read(100,*)(q0(a,1,2,i),a=1,11)
!  elseif(j == 4)then1
!  read(100,*)(q0(a,2,1,i),a=1,11)
!  elseif(j == 5)then
!  read(100,*)(q0(a,2,2,i),a=1,11)
!  else
!  read(100,*)lastline(i)
!  endif
!  enddo
!  write(50,*)filenameall(i),lastline(i)
!  enddo

!  close(100)
!  close(50)


end program main

