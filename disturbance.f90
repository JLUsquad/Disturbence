subroutine disturbance(tau0,tau2,lat_matrix,natom)

  implicit none
  integer(kind=4)                ::  natom, i, j, a, c, k
  real(kind=8)     ::  tau2(3,natom,101,20), lat_matrix(3,3), p(3)
  real(kind=8)     :: tau0(3,natom), invlat_matrix(3,3), tau1(3,natom), tau3(3,natom)
  real(kind=8)                   ::  x, t
  character(4)  ::c1


  call inv_mat(lat_matrix,invlat_matrix)

  do k=1,natom
    do j=1,3
        tau1(j,k)=tau0(1,k)*lat_matrix(1,j) + tau0(2,k)*lat_matrix(2,j) +tau0(3,k)*lat_matrix(3,j)
        tau3(j,k)=tau1(j,k)
    end do
  end do
  call random_seed ()
k=0

do c= 1,20
  do a= 1,101
k=k+1
    do i = 1, natom
      do j = 1, 3
      call random_number (x)
      t = x-0.5
      t = t*c/40
      tau3(j,i)=tau1(j,i) + t
      tau2(j,i,a,c)=tau3(1,i)*invlat_matrix(1,j) + tau3(2,i)*invlat_matrix(2,j) + tau3(3,i)*invlat_matrix(3,j)
      
      if(tau2(j,i,a,c) < 0) then
        tau2(j,i,a,c)=tau2(j,i,a,c)+1.0
      else
      end if

      if(tau2(j,i,a,c) > 1.0)  then 
        tau2(j,i,a,c)=tau2(j,i,a,c)-1.0
      else
      end if

      if(a == 1)then
        tau2(j,i,a,c) = tau0(j,i)
      else
      end if

      

      end do
    end do
write(c1,'(i4)')k
open(20,file='POS'//adjustl(trim(c1)))
write(20,*)'Si'
write(20,*)'1.0'
do i=1,3
write(20,*)lat_matrix(i,:)
enddo
write(20,*)natom
write(20,*)'Direct'
do i=1,natom
write(20,*)tau2(:,i,a,c)
enddo
close(20)
  end do
end do



end subroutine disturbance

subroutine inv_mat(A,invA)
   implicit none
   integer(kind = 4)     :: i,j,k,l,is(3),js(3)
   real(kind = 8)        :: A(3,3),invA(3,3),A_old(3,3)
   real(kind = 8)        :: t,d
   A_old=A
   l=1
   do k=1,3
      d=1D-3
      do i=k,3
         do j=k,3
            if (abs(a(i,j)).gt.d) then
               d=abs(a(i,j))
               is(k)=i
               js(k)=j
            endif
         enddo
      enddo
      if ((d+1.0).eq.1.0) then
         l=0
      endif
      do j=1,3
         t=a(k,j)
         a(k,j)=a(is(k),j)
         a(is(k),j)=t
      enddo
      do i=1,3
         t=a(i,k)
         a(i,k)=a(i,js(k))
         a(i,js(k))=t
      enddo
      a(k,k)=1/a(k,k)
      do j=1,3
         if(j.ne.k) then
            a(k,j)=a(k,j)*a(k,k)
         endif
      enddo
      do i=1,3
         if (i.ne.k) then
            do j=1,3
               if (j.ne.k) then
                  a(i,j)=a(i,j)-a(i,k)*a(k,j)
               endif
            enddo
         endif
      enddo
      do i=1,3
         if(i.ne.k) then
            a(i,k)=-a(i,k)*a(k,k)
         endif
      enddo
   enddo
   do k=3,1,-1
      do j=1,3
         t=a(k,j)
         a(k,j)=a(js(k),j)
         a(js(k),j)=t
      enddo
      do i=1,3
         t=a(i,k)
         a(i,k)=a(i,is(k))
         a(i,is(k))=t
      enddo
   enddo
   invA=A
   A=A_old
   end subroutine
