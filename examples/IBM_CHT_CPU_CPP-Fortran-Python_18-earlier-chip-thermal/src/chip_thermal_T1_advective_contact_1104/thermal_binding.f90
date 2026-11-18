module lica_thermal_component_binding
  use, intrinsic :: iso_c_binding
  implicit none
  ! Interface declaration only. Existing LICA LES/IBM does not call this module yet.
  interface
    integer(c_int) function thermal_step(n,m,nf,edge,g,c,old,source,h,tin,dt,theta,q,boundary) &
      bind(C,name="lica_thermal_step")
      import c_int,c_double
      integer(c_int),value :: n,m,nf
      integer(c_int),intent(in) :: edge(*)
      real(c_double),intent(in) :: g(*),c(*),old(*),source(*)
      real(c_double),value :: h,tin,dt
      real(c_double),intent(out) :: theta(*),q(*),boundary(*)
    end function
  end interface
end module
