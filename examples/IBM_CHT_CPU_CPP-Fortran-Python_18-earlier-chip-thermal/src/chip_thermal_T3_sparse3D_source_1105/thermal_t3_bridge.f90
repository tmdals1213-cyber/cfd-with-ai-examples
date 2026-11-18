module thermal_t3_bridge
 use, intrinsic :: iso_c_binding, only: c_int,c_double
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: thermal_t3_checked
 interface
  integer(c_int) function t3_native(meta,ml,ctl,cl,kind,cell,edge,geom,adv,mdot, &
   & boundary,bc,tout,qcond,qadv,qbnd,stats) bind(C,name='thermal_t3_step')
   import c_int,c_double
   integer(c_int), intent(in) :: meta(*),kind(*),edge(*),adv(*),boundary(*)
   integer(c_int), value :: ml,cl
   real(c_double), intent(in) :: ctl(*),cell(*),geom(*),mdot(*),bc(*)
   real(c_double), INTENT(INOUT) :: tout(*),qcond(*),qadv(*),qbnd(*),stats(*)
  end function
 end interface
contains
 integer(c_int) function thermal_t3_checked(meta,ctl,kind,cell,edge,geom,adv, &
  & mdot,boundary,bc,tout,qcond,qadv,qbnd,stats) result(code)
  integer(c_int), intent(in) :: meta(:),kind(:),edge(:),adv(:),boundary(:)
  real(c_double), intent(in) :: ctl(:),cell(:),geom(:),mdot(:),bc(:)
  real(c_double), INTENT(INOUT) :: tout(:),qcond(:),qadv(:),qbnd(:),stats(:)
  integer :: n,ne,na,nb
  code=91_c_int
  if(size(meta)/=23 .or. size(ctl)/=4) return
  if(.not.is_contiguous(meta) .or. .not.is_contiguous(ctl)) return
  n=meta(2);ne=meta(3);na=meta(4);nb=meta(5)
  if(n<1 .or. n>4096 .or. ne<0 .or. ne>16384 .or. na<0 .or. na>16384 .or. nb<0 .or. nb>8192) return
  if(size(kind)/=n .or. size(cell)/=6*n .or. size(edge)/=2*ne .or. size(geom)/=4*ne) return
  if(size(adv)/=2*na .or. size(mdot)/=na .or. size(boundary)/=nb .or. size(bc)/=4*nb) return
  if(size(tout)/=n .or. size(qcond)/=ne .or. size(qadv)/=na .or. size(qbnd)/=nb .or. size(stats)/=8) return
  if(.not.is_contiguous(kind) .or. .not.is_contiguous(cell) .or. .not.is_contiguous(edge) .or. &
   & .not.is_contiguous(geom)) return
  if(.not.is_contiguous(adv) .or. .not.is_contiguous(mdot) .or. .not.is_contiguous(boundary) .or. &
   & .not.is_contiguous(bc)) return
  if(.not.is_contiguous(tout) .or. .not.is_contiguous(qcond) .or. .not.is_contiguous(qadv) .or. &
   & .not.is_contiguous(qbnd) .or. .not.is_contiguous(stats)) return
  if(any(.not.ieee_is_finite(ctl)) .or. any(.not.ieee_is_finite(cell))) return
  if(any(.not.ieee_is_finite(geom)) .or. any(.not.ieee_is_finite(mdot)) .or. any(.not.ieee_is_finite(bc))) return
  if(meta(6)<1 .or. any(meta(7:9)/=meta(6))) return
  ! No assignment to any INOUT buffer occurs before the atomic native call.
  ! Zero-size edge/advective/boundary arrays require a future tested compiler ABI.
  if(ne==0 .or. na==0 .or. nb==0) then
   code=92_c_int
   return
  end if
  code=t3_native(meta,23_c_int,ctl,4_c_int,kind,cell,edge,geom,adv,mdot,boundary,bc,tout,qcond,qadv,qbnd,stats)
 end function
end module
