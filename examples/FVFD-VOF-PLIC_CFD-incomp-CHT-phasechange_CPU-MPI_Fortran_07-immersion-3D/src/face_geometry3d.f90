module face_geometry3d
 use precision,only:WP
 implicit none
contains
 pure function clipped_quad_fraction(vertices,plane,inside_point)result(frac)
  real(WP),intent(in)::vertices(3,4),plane(4),inside_point(3)
  real(WP)::frac,p(3,8),q(3,8),normal(3),off,sa,sb,t,aa(3),bb(3),av(3),full,side
  integer::n,m,l,r
  normal=plane(1:3);off=plane(4);side=dot_product(normal,inside_point)-off
  if(side>0)then
   normal=-normal;off=-off
  endif
  p=0;p(:,1:4)=vertices;n=4;m=0
  do l=1,n
   r=mod(l,n)+1;aa=p(:,l);bb=p(:,r);sa=dot_product(normal,aa)-off;sb=dot_product(normal,bb)-off
   if(sa<=0)then
    m=m+1;q(:,m)=aa
   endif
   if((sa<0.and.sb>0).or.(sa>0.and.sb<0))then
    t=sa/(sa-sb);m=m+1;q(:,m)=aa+t*(bb-aa)
   endif
  enddo
  av=0
  do l=1,4
   r=mod(l,4)+1;av=av+cross(vertices(:,l),vertices(:,r))
  enddo
  full=sqrt(sum(av*av));av=0
  if(m>=3)then
   do l=1,m
    r=mod(l,m)+1;av=av+cross(q(:,l),q(:,r))
   enddo
  endif
  frac=0
  if(full>0)frac=sqrt(sum(av*av))/full
 end function
 pure function cross(a,b)result(c)
  real(WP),intent(in)::a(3),b(3)
  real(WP)::c(3)
  c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
 end function
end module
