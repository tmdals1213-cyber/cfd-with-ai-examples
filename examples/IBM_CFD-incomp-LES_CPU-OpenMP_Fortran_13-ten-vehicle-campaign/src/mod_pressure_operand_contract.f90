! NOT APPLIED / NOT COMPILED. Fixed nonperiodic scope only.
module mod_pressure_operand_contract
 use mod_common, only: M2,N1,N2,N3,N1M,N2M,N3M,NLEV,IPX,IPZ,IOLDV
 use mod_poiss_mgbicg, only: M1MD,M3MD,IIMG,KKMG
 implicit none
 private
 public :: mg_require_fixed_scope,mg_defined_term,mg_defined_value
contains
 subroutine mg_require_fixed_scope()
  if (N1/=161.or.N2/=81.or.N3/=81.or.NLEV/=3.or.IPX/=0.or.IPZ/=0.or.IOLDV/=0) &
    error stop 'Unreviewed pressure operand topology'
 end subroutine
 recursive function mg_defined_term(GI,linked,coefficient,ii,jj,kk) result(value)
  real*8,intent(in)::GI(0:M1MD,0:M2,0:M3MD),coefficient
  logical,intent(in)::linked
  integer*8,intent(in)::ii,jj,kk
  real*8::value
  if (.not.linked) then
    if (coefficient/=0D0) error stop 'Nonzero coefficient on absent pressure link'
    value=0D0 ! LOCAL absent operator term, never initialization of GI/field.
    return
  endif
  value=coefficient*GI(ii,jj,kk)
 end function
 recursive function mg_defined_value(GI,linked,known_zero_weight,ii,jj,kk) result(value)
  real*8,intent(in)::GI(0:M1MD,0:M2,0:M3MD)
  logical,intent(in)::linked,known_zero_weight
  integer*8,intent(in)::ii,jj,kk
  real*8::value
  if (.not.linked) then
    if (.not.known_zero_weight) error stop 'Unsupported active prolongation term'
    value=0D0 ! LOCAL inactive operand only; zero weight is checked explicitly.
    return
  endif
  value=GI(ii,jj,kk)
 end function
end module
