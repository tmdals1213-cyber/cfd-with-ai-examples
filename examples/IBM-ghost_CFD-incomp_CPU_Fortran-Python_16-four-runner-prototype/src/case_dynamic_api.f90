module case_dynamic_api
 use iso_c_binding
 implicit none
 type,bind(C)::rd_event
  integer(c_int)::info(4)
  real(c_double)::data(13)
 end type
 interface
  integer(c_int) function rd_segment(handle,a,b,time,epoch,owner,capacity,raw,num,si,ei,events) bind(C,name='rd_segment')
   import c_int,c_double,c_ptr,c_int64_t,rd_event
   type(c_ptr),value::handle
   real(c_double),intent(in)::a(3),b(3)
   real(c_double),value::time
   integer(c_int64_t),value::epoch
   integer(c_int),value::owner,capacity
   integer(c_int),intent(inout)::raw,num,si,ei
   type(rd_event),intent(inout)::events(*)
  end function
 end interface
end module
