module phase_birth_core
 use iso_fortran_env,only:real64
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 integer,parameter,public::dp=real64,ok=0,bad_input=1,phase_overdraw=2,unmatched_volume=3,insufficient_heat=4
 type,public::birth_budget
  real(dp)::liquid_volume=0,vapor_volume=0,thermal_energy=0
  real(dp)::mass_residual=0,volume_residual=0,latent_energy=0,expansion_volume=0
 end type
 public::conservative_phase_update,wall_site_activation
contains
 subroutine wall_site_activation(twall,tliquid,tsat,tnuc,site_enabled,vapor_present,ready,activate,status)
  real(dp),intent(in)::twall,tliquid,tsat,tnuc
  logical,intent(in)::site_enabled,vapor_present,ready
  logical,intent(out)::activate
  integer,intent(out)::status
  activate=.false.;status=bad_input
  if(.not.all(ieee_is_finite([twall,tliquid,tsat,tnuc])))return
  ! This is a prescribed-site closure. Its threshold must be configured and validated externally.
  if(tsat.le.0.or.tnuc.le.tsat)return
  status=ok
  activate=site_enabled.and.ready.and..not.vapor_present.and.twall.ge.tnuc.and.tliquid.gt.tsat
 end subroutine
 subroutine conservative_phase_update(cell,vl,vg,flux_l,flux_g,dm,rl,rg,hfg,heat,a,status)
  real(dp),intent(in)::cell,vl,vg,flux_l,flux_g,dm,rl,rg,hfg,heat
  type(birth_budget),intent(out)::a
  integer,intent(out)::status
  real(dp)::nl,ng,latent,mtol,vtol
  a=birth_budget();status=bad_input
  if(.not.all(ieee_is_finite([cell,vl,vg,flux_l,flux_g,dm,rl,rg,hfg,heat])))return
  if(min(cell,rl,rg,hfg).le.0.or.min(vl,vg,heat).lt.0)return
  vtol=128*epsilon(cell)*cell
  if(abs(vl+vg-cell).gt.vtol)return
  nl=vl+flux_l-dm/rl;ng=vg+flux_g+dm/rg;latent=dm*hfg
  a%volume_residual=nl+ng-cell
  a%expansion_volume=dm*(1/rg-1/rl);a%latent_energy=latent
  a%mass_residual=rl*nl+rg*ng-(rl*vl+rg*vg)-(rl*flux_l+rg*flux_g)
  status=phase_overdraw
  if(min(nl,ng).lt.0)return
  status=unmatched_volume
  if(abs(a%volume_residual).gt.vtol)return
  status=insufficient_heat
  if(latent.gt.heat)return
  a%liquid_volume=nl;a%vapor_volume=ng;a%thermal_energy=heat-latent;status=ok
 end subroutine
end module
