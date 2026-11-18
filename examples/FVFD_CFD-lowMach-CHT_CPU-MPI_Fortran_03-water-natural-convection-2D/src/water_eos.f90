module water_eos
   use iso_fortran_env, only: wp=>real64
   implicit none
   private
   public :: load_water,water,face_conductivity,reference_T,reference_p
   real(wp) :: reference_T,reference_p
   real(wp), allocatable :: tab(:,:)
contains
   subroutine load_water(filename)
      character(*), intent(in) :: filename
      integer :: u,n,i
      open(newunit=u,file=filename,status='old',action='read')
      read(u,*) n,reference_p,reference_T
      allocate(tab(6,n))
      do i=1,n
         read(u,*) tab(:,i)
      end do
      close(u)
   end subroutine
   subroutine water(T,rho,h,cp,k,mu,drho,dh)
      real(wp), intent(in) :: T
      real(wp), intent(out) :: rho,h,cp,k,mu
      real(wp), optional, intent(out) :: drho,dh
      integer :: i,n
      real(wp) :: s,d,values(5)
      n=size(tab,2)
      if (T.lt.tab(1,1).or.T.gt.tab(1,n)) then
         print*, 'WATER_EOS_OUT_OF_RANGE_K',T,tab(1,1),tab(1,n)
         error stop 'Water lookup refuses extrapolation or phase change'
      end if
      d=tab(1,2)-tab(1,1)
      i=min(n-1,max(1,1+int((T-tab(1,1))/d)))
      s=(T-tab(1,i))/(tab(1,i+1)-tab(1,i))
      values=(1.0_wp-s)*tab(2:6,i)+s*tab(2:6,i+1)
      rho=values(1); h=values(2); cp=values(3); k=values(4); mu=values(5)
      if (present(drho)) drho=(tab(2,i+1)-tab(2,i))/(tab(1,i+1)-tab(1,i))
      if (present(dh)) dh=(tab(3,i+1)-tab(3,i))/(tab(1,i+1)-tab(1,i))
   end subroutine
   pure real(wp) function face_conductivity(k1,k2)
      real(wp), intent(in) :: k1,k2
      face_conductivity=2.0_wp/(1.0_wp/k1+1.0_wp/k2)
   end function
end module
