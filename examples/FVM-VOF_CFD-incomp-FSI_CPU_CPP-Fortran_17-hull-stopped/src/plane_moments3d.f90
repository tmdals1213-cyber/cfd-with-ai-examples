! Original convex-polyhedron integration, not copied Basilisk source.
! Phase={x in[0,1]^3:dot(normal,x)<=intercept}. No vf clipping/mass fix.
module plane_moments3d
  use geometry3d, only: dp,qp
  implicit none
  private
  public :: cube_phase_moments,cube_interface_polygon
contains
  pure function cross_q(a,b) result(c)
    real(qp),intent(in) :: a(3),b(3)
    real(qp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  pure subroutine cube_vertices(vertices)
    real(qp),intent(out) :: vertices(3,0:7)
    integer :: mask,j
    do mask=0,7
      do j=1,3
        vertices(j,mask)=merge(1._qp,0._qp,btest(mask,j-1))
      enddo
    enddo
  end subroutine

  pure subroutine append_unique(point,points,count)
    real(qp),intent(in) :: point(3)
    real(qp),intent(inout) :: points(3,24)
    integer,intent(inout) :: count
    integer :: j
    do j=1,count
      ! Exact endpoint deduplication avoids erasing legitimately tiny cuts.
      if (all(point==points(:,j))) return
    enddo
    count=count+1
    points(:,count)=point
  end subroutine

  pure subroutine interface_polygon_q(n,c,points,count)
    real(qp),intent(in) :: n(3),c
    real(qp),intent(out) :: points(3,24)
    integer,intent(out) :: count
    real(qp) :: vertices(3,0:7),a(3),b(3),sa,sb,point(3),center(3),normal(3),reference(3),basis1(3),basis2(3),angles(24),angle
    integer :: axis,mask,other,j,k,pivot
    count=0;points=0._qp
    if (sum(abs(n))==0._qp) return
    if (c<=sum(min(n,0._qp)).or.c>=sum(max(n,0._qp))) return
    call cube_vertices(vertices)
    do axis=1,3
      do mask=0,7
        if (btest(mask,axis-1)) cycle
        other=ibset(mask,axis-1)
        a=vertices(:,mask);b=vertices(:,other)
        sa=dot_product(n,a)-c;sb=dot_product(n,b)-c
        if (sa==0._qp) call append_unique(a,points,count)
        if (sb==0._qp) call append_unique(b,points,count)
        if ((sa<0._qp.and.sb>0._qp).or.(sa>0._qp.and.sb<0._qp)) then
          point=a+(sa/(sa-sb))*(b-a)
          call append_unique(point,points,count)
        endif
      enddo
    enddo
    if (count<3) return
    center=sum(points(:,1:count),dim=2)/real(count,qp)
    normal=n/sqrt(sum(n*n))
    pivot=minloc(abs(normal),dim=1)
    reference=0._qp;reference(pivot)=1._qp
    basis1=cross_q(normal,reference);basis1=basis1/sqrt(sum(basis1*basis1))
    basis2=cross_q(normal,basis1)
    do j=1,count
      point=points(:,j)-center
      angles(j)=atan2(dot_product(point,basis2),dot_product(point,basis1))
    enddo
    do j=2,count
      point=points(:,j);angle=angles(j);k=j-1
      do while(k>=1)
        if (angles(k)<=angle) exit
        points(:,k+1)=points(:,k);angles(k+1)=angles(k);k=k-1
      enddo
      points(:,k+1)=point;angles(k+1)=angle
    enddo
    ! Ascending basis angles give outward +normal polygon winding.
  end subroutine

  pure subroutine tetra_triangle(a,b,c,anchor,volume,moment)
    real(qp),intent(in) :: a(3),b(3),c(3),anchor(3)
    real(qp),intent(inout) :: volume,moment(3)
    real(qp) :: signed_volume
    signed_volume=dot_product(a-anchor,cross_q(b-anchor,c-anchor))/6._qp
    volume=volume+signed_volume
    moment=moment+signed_volume*(anchor+a+b+c)/4._qp
  end subroutine

  pure subroutine smaller_phase_q(n,c,volume,moment)
    real(qp),intent(in) :: n(3),c
    real(qp),intent(out) :: volume,moment(3)
    integer,parameter :: faces(4,6)=reshape([0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6],[4,6])
    real(qp) :: vertices(3,0:7),polygon(3,12),cap(3,24),a(3),b(3),anchor(3),sa,sb
    integer :: face,j,previous,count,cap_count
    logical :: a_inside,b_inside
    call cube_vertices(vertices)
    ! The minimizing cube vertex lies inside the phase. Using it as the
    ! tetrahedral anchor avoids subtracting large terms for tiny corner cuts.
    anchor=merge(1._qp,0._qp,n<0._qp)
    volume=0._qp;moment=0._qp
    do face=1,6
      count=0
      do j=1,4
        previous=modulo(j-2,4)+1
        a=vertices(:,faces(previous,face));b=vertices(:,faces(j,face))
        sa=dot_product(n,a)-c;sb=dot_product(n,b)-c
        a_inside=sa<=0._qp;b_inside=sb<=0._qp
        if (a_inside.neqv.b_inside) then
          count=count+1;polygon(:,count)=a+(sa/(sa-sb))*(b-a)
        endif
        if (b_inside) then
          count=count+1;polygon(:,count)=b
        endif
      enddo
      do j=2,count-1
        call tetra_triangle(polygon(:,1),polygon(:,j),polygon(:,j+1),anchor,volume,moment)
      enddo
    enddo
    call interface_polygon_q(n,c,cap,cap_count)
    do j=2,cap_count-1
      call tetra_triangle(cap(:,1),cap(:,j),cap(:,j+1),anchor,volume,moment)
    enddo
  end subroutine

  subroutine cube_phase_moments(normal,intercept,volume,moments)
    real(dp),intent(in) :: normal(3)
    real(qp),intent(in) :: intercept
    real(dp),intent(out) :: volume,moments(3)
    real(qp) :: n(3),c,scale,lower,upper,v,m(3)
    n=real(normal,qp);scale=sum(abs(n))
    if (scale==0._qp) then
      if (intercept>=0._qp) then
        volume=1._dp;moments=.5_dp
      else
        volume=0._dp;moments=0._dp
      endif
      return
    endif
    n=n/scale;c=intercept/scale
    lower=sum(min(n,0._qp));upper=sum(max(n,0._qp))
    if (c<=lower) then
      volume=0._dp;moments=0._dp
    else if (c>=upper) then
      volume=1._dp;moments=.5_dp
    else
      if (c>(lower+upper)/2._qp) then
        call smaller_phase_q(-n,-c,v,m)
        v=1._qp-v;m=.5_qp-m
      else
        call smaller_phase_q(n,c,v,m)
      endif
      volume=real(v,dp);moments=real(m,dp)
    endif
  end subroutine

  subroutine cube_interface_polygon(normal,intercept,has_interface,area,centroid,second_moments)
    real(dp),intent(in) :: normal(3)
    real(qp),intent(in) :: intercept
    logical,intent(out) :: has_interface
    real(dp),intent(out) :: area,centroid(3)
    ! Optional exact raw surface second moments, integral(x_i*x_j dA).
    ! Needed if a later diagnostic squares height; centroid squared is not it.
    real(dp),optional,intent(out) :: second_moments(3,3)
    real(qp) :: n(3),c,scale,points(3,24),a(3),b(3),d(3),cross(3),triangle_area,total,first(3),second(3,3),sums(3)
    integer :: count,j,k,l
    area=0._dp;centroid=0._dp;has_interface=.false.
    if (present(second_moments)) second_moments=0._dp
    n=real(normal,qp);scale=sum(abs(n))
    if (scale==0._qp) return
    n=n/scale;c=intercept/scale
    call interface_polygon_q(n,c,points,count)
    if (count<3) return
    total=0._qp;first=0._qp;second=0._qp
    a=points(:,1)
    do j=2,count-1
      b=points(:,j);d=points(:,j+1)
      cross=cross_q(b-a,d-a);triangle_area=sqrt(sum(cross*cross))/2._qp
      total=total+triangle_area
      sums=a+b+d
      first=first+triangle_area*sums/3._qp
      do k=1,3
        do l=1,3
          second(k,l)=second(k,l)+triangle_area*(sums(k)*sums(l)+a(k)*a(l)+b(k)*b(l)+d(k)*d(l))/12._qp
        enddo
      enddo
    enddo
    if (total>0._qp) then
      has_interface=.true.;area=real(total,dp);centroid=real(first/total,dp)
      if (present(second_moments)) second_moments=real(second,dp)
    endif
  end subroutine
end module plane_moments3d
