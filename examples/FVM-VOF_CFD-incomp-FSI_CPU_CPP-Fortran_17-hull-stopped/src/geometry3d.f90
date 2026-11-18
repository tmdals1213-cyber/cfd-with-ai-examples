! Original standalone implementation, 2026-10-29. No Basilisk source copied.
! Phase convention: {x in [0,1]^3 : dot_product(normal,x) <= intercept}.
! See README.md for derivation, primary references, and numerical limits.
module geometry3d
  use iso_fortran_env, only: real64, real128
  implicit none
  integer, parameter :: dp=real64, qp=real128
contains
  recursive pure function volume_q(normal, intercept) result(v)
    real(qp), intent(in) :: normal(3), intercept
    real(qp) :: v, a(3), b(3), t, scale, total, denom, term, accum, tmp
    integer :: d, j, k, mask, parity
    logical :: complement
    a=abs(normal)
    scale=sum(a)
    if (scale == 0._qp) then
      if (intercept >= 0._qp) then
        v=1._qp
      else
        v=0._qp
      endif
      return
    endif
    t=(intercept-sum(min(normal,0._qp)))/scale
    a=a/scale
    ! Negligible axes are integrated at their midpoint. The normalized
    ! plane offset perturbation is at most 1.5e-18, and is documented.
    d=0
    do j=1,3
      if (a(j) > 1.e-18_qp) then
        d=d+1
        b(d)=a(j)
      else
        t=t-a(j)/2._qp
      endif
    enddo
    total=sum(b(1:d))
    if (t <= 0._qp) then
      v=0._qp
      return
    else if (t >= total) then
      v=1._qp
      return
    endif
    complement=t > total/2._qp
    if (complement) t=total-t
    do j=1,d
      do k=j+1,d
        if (b(k) < b(j)) then
          tmp=b(j); b(j)=b(k); b(k)=tmp
        endif
      enddo
    enddo
    ! Central linear interval when the largest coefficient dominates:
    ! integrating the smaller directions analytically avoids cancellation.
    if (d == 1) then
      v=t/b(1)
    else if (b(d) >= sum(b(1:d-1)) .and. t >= sum(b(1:d-1))) then
      v=(t-sum(b(1:d-1))/2._qp)/b(d)
    else
      ! Inclusion-exclusion on a d-simplex with the upper box limits.
      accum=0._qp
      do mask=0,2**d-1
        term=t
        parity=0
        do j=1,d
          if (btest(mask,j-1)) then
            term=term-b(j)
            parity=parity+1
          endif
        enddo
        if (term > 0._qp) then
          if (mod(parity,2)==0) then
            accum=accum+term**d
          else
            accum=accum-term**d
          endif
        endif
      enddo
      denom=product(b(1:d))
      do j=2,d
        denom=denom*real(j,qp)
      enddo
      v=accum/denom
    endif
    if (complement) v=1._qp-v
    ! Deliberately no numerical clipping.
  end function volume_q

  recursive pure function cube_volume(normal, intercept) result(v)
    real(dp), intent(in) :: normal(3), intercept
    real(dp) :: v
    v=real(volume_q(real(normal,qp),real(intercept,qp)),dp)
  end function cube_volume

  recursive pure function plane_offset(normal, fraction) result(intercept)
    real(dp), intent(in) :: normal(3), fraction
    real(qp) :: intercept, lo, hi, mid, nq(3), target
    integer :: it
    nq=real(normal,qp)
    lo=sum(min(nq,0._qp)); hi=sum(max(nq,0._qp))
    if (fraction<=0._dp) then
      intercept=lo
      return
    else if (fraction>=1._dp) then
      intercept=hi
      return
    endif
    target=real(fraction,qp)
    do it=1,58
      mid=(lo+hi)/2._qp
      if (volume_q(nq,mid) < target) then
        lo=mid
      else
        hi=mid
      endif
    enddo
    intercept=(lo+hi)/2._qp
  end function plane_offset

  recursive pure function box_volume(normal, intercept, lower, upper) result(v)
    real(dp), intent(in) :: normal(3), lower(3), upper(3)
    real(qp), intent(in) :: intercept
    real(dp) :: v
    real(qp) :: width(3), nq(3), lo(3)
    width=real(upper-lower,qp); lo=real(lower,qp); nq=real(normal,qp)
    v=real(product(width)*volume_q(nq*width,intercept-sum(nq*lo)),dp)
  end function box_volume

  recursive pure function slab_volume(normal, intercept, axis, lower, upper) result(v)
    real(dp), intent(in) :: normal(3), lower, upper
    real(qp), intent(in) :: intercept
    integer, intent(in) :: axis
    real(dp) :: v, lo(3), hi(3)
    lo=0._dp; hi=1._dp
    lo(axis)=lower; hi(axis)=upper
    v=box_volume(normal,intercept,lo,hi)
  end function slab_volume

  recursive pure function signed_flux(normal, intercept, axis, courant) result(flux)
    real(dp), intent(in) :: normal(3), courant
    real(qp), intent(in) :: intercept
    integer, intent(in) :: axis
    real(dp) :: flux
    if (courant>=0._dp) then
      flux=slab_volume(normal,intercept,axis,1._dp-courant,1._dp)
    else
      flux=-slab_volume(normal,intercept,axis,0._dp,-courant)
    endif
  end function signed_flux
end module geometry3d
