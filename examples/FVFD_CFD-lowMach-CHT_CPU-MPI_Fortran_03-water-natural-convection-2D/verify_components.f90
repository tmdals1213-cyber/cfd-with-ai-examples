program verify_components
   use iso_fortran_env, only: wp=>real64
   use water_eos
   implicit none
   integer, parameter :: n=28
   real(wp) :: k(n),aa(n),bb(n),cc(n),rhs(n),T(n),analytical(n),dx,q,g,r,h,cp,kw,mu,dr,dh,x,err
   integer :: i
   call load_water('data/water_p101325.dat')
   call water(298.15_wp,r,h,cp,kw,mu,dr,dh)
   if (abs(r-997.047636760347_wp).gt.1e-8_wp) error stop 'IAPWS reference density failed'
   if (abs(h).gt.1e-5_wp.or.dr.ge.0.0_wp.or.dh.le.0.0_wp) error stop 'EOS derivative or enthalpy reference failed'
   dx=.0005_wp; k(:8)=400.0_wp; k(9:)=kw
   aa=0; bb=0; cc=0; rhs=0
   do i=1,n
      if (i.eq.1) then
         g=2*k(i)/dx; bb(i)=bb(i)+g; rhs(i)=rhs(i)+g*20.0_wp
      else
         g=face_conductivity(k(i-1),k(i))/dx; bb(i)=bb(i)+g; aa(i)=-g
      end if
      if (i.eq.n) then
         g=2*k(i)/dx; bb(i)=bb(i)+g
      else
         g=face_conductivity(k(i),k(i+1))/dx; bb(i)=bb(i)+g; cc(i)=-g
      end if
   end do
   do i=2,n
      g=aa(i)/bb(i-1); bb(i)=bb(i)-g*cc(i-1); rhs(i)=rhs(i)-g*rhs(i-1)
   end do
   T(n)=rhs(n)/bb(n)
   do i=n-1,1,-1
      T(i)=(rhs(i)-cc(i)*T(i+1))/bb(i)
   end do
   q=20.0_wp/(.004_wp/400.0_wp+.010_wp/kw)
   do i=1,n
      x=(i-.5_wp)*dx
      if (i.le.8) then
         analytical(i)=20.0_wp-q*x/400.0_wp
      else
         analytical(i)=q*(.014_wp-x)/kw
      end if
   end do
   err=maxval(abs(T-analytical))
   if (err.gt.1e-9_wp) error stop 'Composite-layer conduction benchmark failed'
   if (abs(face_conductivity(k(8),k(9))*(T(8)-T(9))/dx-q).gt.1e-6_wp) error stop 'Interface flux benchmark failed'
   print*, 'PASS IAPWS_REFERENCE_AND_DERIVATIVE'
   print*, 'PASS COMPOSITE_LAYER_ANALYTIC_MAX_ERROR_K',err
   print*, 'PASS COMPOSITE_LAYER_HEAT_FLUX_W_M2',q
end program
