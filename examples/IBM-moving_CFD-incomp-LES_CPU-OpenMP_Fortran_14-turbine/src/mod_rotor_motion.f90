! Prescribed rigid rotation of the STL about the x axis. No structural solver.
module mod_rotor_motion
 use iso_fortran_env, only: real64
 use mod_stl_sdf, only: evaluate_sdf
 implicit none
 real(real64),save :: rotor_re,physical_diameter,rpm,kinematic_nu,inflow_ratio,ramp_time
 real(real64),save :: reference_speed,reference_time,omega_max,omega_now=0D0
 real(real64),save :: angle_now=0D0,ctheta=1D0,stheta=0D0,x_half
contains
 subroutine load_rotor(filename)
  character(len=*),intent(in)::filename
  integer::unit
  open(newunit=unit,file=filename,status='old',action='read')
  read(unit,*) physical_diameter,rpm,kinematic_nu,inflow_ratio,ramp_time
  read(unit,*) x_half
  close(unit)
  if(physical_diameter<=0D0.or.abs(rpm)<1D-12.or.kinematic_nu<=0D0.or.ramp_time<=0D0) error stop 'Invalid rotor inputs'
  reference_speed=abs(rpm)*ACOS(-1D0)*physical_diameter/60D0
  reference_time=physical_diameter/reference_speed
  rotor_re=reference_speed*physical_diameter/kinematic_nu
  omega_max=sign(2D0,rpm)
  write(*,*) 'ROTOR D,RPM,NU,UREF,RE=',physical_diameter,rpm,kinematic_nu,reference_speed,rotor_re
 end subroutine
 subroutine update_rotor(t)
  real(real64),intent(in)::t
  real(real64)::pi,s
  pi=ACOS(-1D0)
  if(t<ramp_time)then
   s=max(0D0,t)
   omega_now=omega_max*.5D0*(1D0-cos(pi*s/ramp_time))
   angle_now=omega_max*.5D0*(s-ramp_time/pi*sin(pi*s/ramp_time))
  else
   omega_now=omega_max
   angle_now=omega_max*(t-.5D0*ramp_time)
  endif
  ctheta=cos(angle_now);stheta=sin(angle_now)
 end subroutine
 function rotor_sdf(x,y,z) result(f)
  real(real64),intent(in)::x,y,z
  real(real64)::f,d,a,w,base
  ! Smooth positive broad-phase extension. All possible solid positions lie
  ! inside d<=1, so its zero set is exactly the rotating sampled STL surface.
  d=max(abs(x)/x_half,2D0*abs(y),2D0*abs(z))
  base=.5D0*(d-1D0)
  if(d>=1.3D0)then
   f=base
  else
   f=evaluate_sdf(x,ctheta*y+stheta*z,-stheta*y+ctheta*z)
   if(d>1.1D0)then
    a=(d-1.1D0)/.2D0;w=a*a*(3D0-2D0*a)
    f=(1D0-w)*f+w*base
   endif
  endif
 end function
end module
