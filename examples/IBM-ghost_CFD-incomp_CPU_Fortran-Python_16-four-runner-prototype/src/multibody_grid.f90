module pacemaker_grid_queries
 use iso_c_binding
 use ieee_arithmetic
 use case_line_extension,only:extend_grid_point
 implicit none
contains
 integer(c_int) function pm_geometry(handle,ndof,points,axis,offsets,people,box,time,epoch, &
    chi,wall,owner,triangle,component,hit,badindex,badcode) bind(C,name='pm_geometry') result(rc)
  type(c_ptr),value::handle
  integer(c_int),value::ndof,people
  real(c_double),intent(in)::points(3,*),offsets(3,*),box(3,2)
  integer(c_int),intent(in)::axis(*)
  real(c_double),value::time
  integer(c_int64_t),value::epoch
  real(c_double),intent(inout)::chi(*),wall(*),hit(3,*)
  integer(c_int),intent(inout)::owner(*),triangle(*),component(*),badindex,badcode
  real(c_double),allocatable::ct(:),wt(:),ht(:,:)
  integer(c_int),allocatable::ot(:),tt(:),pt(:)
  real(c_double)::q(3),outside(3),v(3),anchor(3)
  integer(c_int)::code,mask,tr,part
  integer::j,b
  rc=1;badindex=-1;badcode=-1
  if(.not.c_associated(handle).or.ndof<1.or.ndof>1600000.or.people/=4)return
  if(.not.all(ieee_is_finite(box)).or.any(box(:,1)>box(:,2)))return
  allocate(ct(ndof),wt(ndof),ht(3,ndof),ot(ndof),tt(ndof),pt(ndof))
  ct=0;wt=0;ht=0;ot=0;tt=-1;pt=-1
  do j=1,ndof
   if(axis(j)<0.or.axis(j)>2.or..not.all(ieee_is_finite(points(:,j))))then
    badindex=j-1;badcode=1;return
   end if
   do b=1,people
    q=points(:,j)-offsets(:,b)
    ! AABB exclusion is conservative; near-box points still use the qualified
    ! signed-winding/exposed-event provider. No approximate inside classifier.
    if(any(q<box(:,1)-2.e-10_c_double).or.any(q>box(:,2)+2.e-10_c_double))cycle
    outside=q;outside(1)=box(1,1)-1._c_double
    mask=-31;tr=-41;part=-51;v=17._c_double;anchor=-19._c_double
    code=extend_grid_point(handle,outside,q,time,epoch,1_c_int,mask,v,anchor,tr,part)
    if(code/=0.or.(mask/=0.and.mask/=1))then
     badindex=j-1;badcode=code;rc=2;return
    end if
    if(mask==1)then
     if(ot(j)/=0.or.tr<0.or.tr>=49112.or.part<0.or.part>=95.or..not.all(ieee_is_finite(v)))then
      badindex=j-1;badcode=83;rc=3;return
     end if
     ct(j)=1;wt(j)=v(axis(j)+1);ht(:,j)=anchor+offsets(:,b);ot(j)=b;tt(j)=tr;pt(j)=part
    else
     if(tr/=-1.or.part/=-1.or.any(v/=17._c_double).or.any(anchor/=-19._c_double))then
      badindex=j-1;badcode=84;rc=4;return
     end if
    end if
   end do
  end do
  chi(1:ndof)=ct;wall(1:ndof)=wt;hit(:,1:ndof)=ht
  owner(1:ndof)=ot;triangle(1:ndof)=tt;component(1:ndof)=pt
  badindex=-1;badcode=0;rc=0
 end function
end module
