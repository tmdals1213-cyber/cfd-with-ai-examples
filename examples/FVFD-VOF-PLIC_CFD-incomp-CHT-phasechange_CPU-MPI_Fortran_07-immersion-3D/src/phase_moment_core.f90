module phase_moment_core
 use iso_fortran_env,only:real64
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 integer,parameter,public::moment_ok=0,moment_bad_input=1,moment_bad_volume=2,moment_range=3
 public::depleting_phase_centroid
contains
 pure subroutine depleting_phase_centroid(advected_moment,available_volume,remaining_volume,centre,status)
  real(real64),intent(in)::advected_moment(3),available_volume,remaining_volume
  real(real64),intent(out)::centre(3)
  integer,intent(out)::status
  integer::a
  centre=0;status=moment_bad_input
  if(.not.all(ieee_is_finite(advected_moment)))return
  if(.not.all(ieee_is_finite([available_volume,remaining_volume])))return
  status=moment_bad_volume
  if(available_volume.le.0.or.remaining_volume.lt.0.or.remaining_volume.gt.available_volume)return
  status=moment_range
  if(available_volume.lt.1)then
   if(maxval(abs(advected_moment)).gt.huge(1.0_real64)*available_volume)return
  endif
  ! Cell-integrated depletion removes the phase-average specific first moment:
  ! M_new = M_adv + S*(M_adv/V_adv) = V_new*(M_adv/V_adv), S<=0.
  ! No centroid/velocity clipping, source-rate change or discarded energy.
  do a=1,3
   centre(a)=advected_moment(a)/available_volume
  enddo
  if(.not.all(ieee_is_finite(centre)))return
  status=moment_ok
 end subroutine
end module
