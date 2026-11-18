module case_line_extension
 use iso_c_binding
 use ieee_arithmetic
 use case_dynamic_api,only:rd_event,rd_segment
 implicit none
 private
 public :: extension_from_exposed_events,extend_grid_point
contains
 integer(c_int) function extension_from_exposed_events(direction,num,si,ei,events, &
     mask,velocity,anchor,triangle,component) result(code)
  real(c_double),intent(in) :: direction(3)
  integer(c_int),intent(in) :: num,si,ei
  type(rd_event),intent(in) :: events(:)
  integer(c_int),intent(inout) :: mask,triangle,component
  real(c_double),intent(inout) :: velocity(3),anchor(3)
  integer :: j,state,last
  real(c_double) :: previous,dotn
  code=20
  if(si/=0.or.(ei/=0.and.ei/=1)) return
  if(num<0.or.num>size(events)) return
  if(.not.all(ieee_is_finite(direction))) return
  if(dot_product(direction,direction)<=0) return
  state=0;last=0;previous=0
  code=21
  do j=1,num
   if(.not.all(ieee_is_finite(events(j)%data))) return
   if(events(j)%data(1)<=previous.or.events(j)%data(1)>=1) return
   if(events(j)%info(4)/=merge(-1,1,state==0)) return
   if(events(j)%info(1)<0.or.events(j)%info(2)<0) return
   if(abs(sum(events(j)%data(11:13))-1)>1.e-12_c_double) return
   if(abs(sum(events(j)%data(8:10)**2)-1)>1.e-10_c_double) return
   dotn=dot_product(events(j)%data(8:10),direction)
   if(real(events(j)%info(4),c_double)*dotn<=0) return
   previous=events(j)%data(1)
   state=1-state
   if(state==1) last=j
  end do
  if(state/=ei) return
  if(ei==1) then
   if(last==0) return
   velocity=events(last)%data(5:7)
   anchor=events(last)%data(2:4)
   triangle=events(last)%info(1)
   component=events(last)%info(2)
  else
   ! No wall velocity is consumed for mask=0; retain its defined sentinel.
   triangle=-1;component=-1
  end if
  mask=ei
  code=0
 end function

 integer(c_int) function extend_grid_point(handle,outside,point,time,epoch,owner, &
     mask,velocity,anchor,triangle,component) result(code)
  type(c_ptr),value :: handle
  real(c_double),intent(in) :: outside(3),point(3)
  real(c_double),value :: time
  integer(c_int64_t),value :: epoch
  integer(c_int),value :: owner
  integer(c_int),intent(inout) :: mask,triangle,component
  real(c_double),intent(inout) :: velocity(3),anchor(3)
  integer(c_int) :: raw,num,si,ei
  type(rd_event) :: events(256)
  ! Snapshot update/ownership occurs outside OpenMP in the owning fluid caller.
  ! Only component-OR exposed transitions from the qualified native API enter here.
  raw=-1;num=-1;si=-1;ei=-1
  code=rd_segment(handle,outside,point,time,epoch,owner,256_c_int,raw,num,si,ei,events)
  if(code/=0) return
  code=extension_from_exposed_events(point-outside,num,si,ei,events, &
                                   mask,velocity,anchor,triangle,component)
 end function
end module
