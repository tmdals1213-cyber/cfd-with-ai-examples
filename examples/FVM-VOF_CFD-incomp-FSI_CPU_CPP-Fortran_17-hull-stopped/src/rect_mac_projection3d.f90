! Standalone axis-constant rectangular MAC pressure operator; not LICA's
! within-axis nonuniform grid, VOF, momentum, wave or body implementation.
module rect_mac_projection3d
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  integer,parameter::dp=real64
  type::rect_grid
    integer::nx=0,ny=0,nz=0
    real(dp)::h(3)=0,origin(3)=0,length(3)=0,cell_volume=0
  end type
  type::boundary_beta
    real(dp),allocatable::xlo(:,:),xhi(:,:),ylo(:,:),yhi(:,:),zlo(:,:),zhi(:,:)
  end type
  type::solve_stats
    integer::flag=0,iterations=0,restarts=0
    real(dp)::rhs_rms=0,rhs_mean=0,compatibility_tolerance=0
    real(dp)::absres=0,original_absres=0,relres=0,stop_ratio=0,pressure_mean=0
  end type
contains
  subroutine configure_grid(nx,ny,nz,spacing,origin,grid)
    integer,intent(in)::nx,ny,nz
    real(dp),intent(in)::spacing(3),origin(3)
    type(rect_grid),intent(out)::grid
    if(min(nx,ny,nz)<2.or.any(spacing<=0).or..not.all(ieee_is_finite(spacing)).or. &
       .not.all(ieee_is_finite(origin)))error stop 'invalid rectangular grid'
    grid%nx=nx;grid%ny=ny;grid%nz=nz;grid%h=spacing;grid%origin=origin
    grid%length=spacing*real([nx,ny,nz],dp);grid%cell_volume=product(spacing)
  end subroutine
  subroutine allocate_boundary(grid,b)
    type(rect_grid),intent(in)::grid
    type(boundary_beta),intent(out)::b
    allocate(b%xlo(grid%ny,grid%nz),b%xhi(grid%ny,grid%nz),b%ylo(grid%nx,grid%nz),b%yhi(grid%nx,grid%nz), &
             b%zlo(grid%nx,grid%ny),b%zhi(grid%nx,grid%ny))
    b%xlo=1;b%xhi=1;b%ylo=1;b%yhi=1;b%zlo=1;b%zhi=1
  end subroutine
  subroutine build_faces(grid,rho,b,dirichlet,periodic,bx,by,bz,diagonal,rho_effective)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::rho(grid%nx,grid%ny,grid%nz)
    type(boundary_beta),intent(in)::b
    logical,intent(in)::dirichlet(6),periodic(3)
    real(dp),intent(out)::bx(grid%nx+1,grid%ny,grid%nz),by(grid%nx,grid%ny+1,grid%nz),bz(grid%nx,grid%ny,grid%nz+1)
    real(dp),intent(out)::diagonal(grid%nx,grid%ny,grid%nz)
    real(dp),intent(in),optional::rho_effective(grid%nx,grid%ny,grid%nz)
    real(dp)::r(grid%nx,grid%ny,grid%nz)
    integer::nx,ny,nz,c
    nx=grid%nx;ny=grid%ny;nz=grid%nz
    if(any(rho<=0).or..not.all(ieee_is_finite(rho)))error stop 'invalid material rho'
    r=rho
    if(present(rho_effective))r=rho_effective
    if(any(r<=0).or..not.all(ieee_is_finite(r)))error stop 'invalid coefficient inertia'
    if(.not.allocated(b%xlo).or..not.allocated(b%xhi).or..not.allocated(b%ylo).or..not.allocated(b%yhi).or. &
       .not.allocated(b%zlo).or..not.allocated(b%zhi))error stop 'allocate boundary planes first'
    if(any(shape(b%xlo)/=[ny,nz]).or.any(shape(b%xhi)/=[ny,nz]).or. &
       any(shape(b%ylo)/=[nx,nz]).or.any(shape(b%yhi)/=[nx,nz]).or. &
       any(shape(b%zlo)/=[nx,ny]).or.any(shape(b%zhi)/=[nx,ny]))error stop 'boundary plane shape'
    if(any(b%xlo<=0).or.any(b%xhi<=0).or.any(b%ylo<=0).or.any(b%yhi<=0).or.any(b%zlo<=0).or.any(b%zhi<=0)) &
      error stop 'nonpositive boundary beta'
    if(.not.all(ieee_is_finite(b%xlo)).or..not.all(ieee_is_finite(b%xhi)).or. &
       .not.all(ieee_is_finite(b%ylo)).or..not.all(ieee_is_finite(b%yhi)).or. &
       .not.all(ieee_is_finite(b%zlo)).or..not.all(ieee_is_finite(b%zhi)))error stop 'nonfinite boundary beta'
    do c=1,3
      if(periodic(c).and.any(dirichlet(2*c-1:2*c)))error stop 'periodic axis has Dirichlet pressure flag'
    enddo
    bx(2:nx,:,:)=2/(r(1:nx-1,:,:)+r(2:nx,:,:))
    by(:,2:ny,:)=2/(r(:,1:ny-1,:)+r(:,2:ny,:))
    bz(:,:,2:nz)=2/(r(:,:,1:nz-1)+r(:,:,2:nz))
    if(periodic(1))then
      bx(1,:,:)=2/(r(nx,:,:)+r(1,:,:));bx(nx+1,:,:)=bx(1,:,:)
    else
      bx(1,:,:)=b%xlo;bx(nx+1,:,:)=b%xhi
    endif
    if(periodic(2))then
      by(:,1,:)=2/(r(:,ny,:)+r(:,1,:));by(:,ny+1,:)=by(:,1,:)
    else
      by(:,1,:)=b%ylo;by(:,ny+1,:)=b%yhi
    endif
    if(periodic(3))then
      bz(:,:,1)=2/(r(:,:,nz)+r(:,:,1));bz(:,:,nz+1)=bz(:,:,1)
    else
      bz(:,:,1)=b%zlo;bz(:,:,nz+1)=b%zhi
    endif
    diagonal=(bx(1:nx,:,:)+bx(2:nx+1,:,:))/grid%h(1)**2 &
            +(by(:,1:ny,:)+by(:,2:ny+1,:))/grid%h(2)**2 &
            +(bz(:,:,1:nz)+bz(:,:,2:nz+1))/grid%h(3)**2
    if(.not.periodic(1))then
      diagonal(1,:,:)=diagonal(1,:,:)+merge(1._dp,-1._dp,dirichlet(1))*bx(1,:,:)/grid%h(1)**2
      diagonal(nx,:,:)=diagonal(nx,:,:)+merge(1._dp,-1._dp,dirichlet(2))*bx(nx+1,:,:)/grid%h(1)**2
    endif
    if(.not.periodic(2))then
      diagonal(:,1,:)=diagonal(:,1,:)+merge(1._dp,-1._dp,dirichlet(3))*by(:,1,:)/grid%h(2)**2
      diagonal(:,ny,:)=diagonal(:,ny,:)+merge(1._dp,-1._dp,dirichlet(4))*by(:,ny+1,:)/grid%h(2)**2
    endif
    if(.not.periodic(3))then
      diagonal(:,:,1)=diagonal(:,:,1)+merge(1._dp,-1._dp,dirichlet(5))*bz(:,:,1)/grid%h(3)**2
      diagonal(:,:,nz)=diagonal(:,:,nz)+merge(1._dp,-1._dp,dirichlet(6))*bz(:,:,nz+1)/grid%h(3)**2
    endif
    ! Mixed periodic + pressure-Dirichlet is supported. The only unanchored
    ! topology tested/supported here is fully periodic, not pure Neumann.
    if(.not.any(dirichlet).and..not.all(periodic))error stop 'unanchored non-fully-periodic topology unsupported'
    if(any(diagonal<=0))error stop 'invalid rectangular Jacobi diagonal'
  end subroutine
  subroutine weighted_gradient(grid,p,bx,by,bz,dirichlet,periodic,gx,gy,gz)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::p(grid%nx,grid%ny,grid%nz),bx(grid%nx+1,grid%ny,grid%nz), &
                        by(grid%nx,grid%ny+1,grid%nz),bz(grid%nx,grid%ny,grid%nz+1)
    logical,intent(in)::dirichlet(6),periodic(3)
    real(dp),intent(out)::gx(grid%nx+1,grid%ny,grid%nz),gy(grid%nx,grid%ny+1,grid%nz),gz(grid%nx,grid%ny,grid%nz+1)
    integer::nx,ny,nz
    nx=grid%nx;ny=grid%ny;nz=grid%nz;gx=0;gy=0;gz=0
    gx(2:nx,:,:)=bx(2:nx,:,:)*(p(2:nx,:,:)-p(1:nx-1,:,:))/grid%h(1)
    gy(:,2:ny,:)=by(:,2:ny,:)*(p(:,2:ny,:)-p(:,1:ny-1,:))/grid%h(2)
    gz(:,:,2:nz)=bz(:,:,2:nz)*(p(:,:,2:nz)-p(:,:,1:nz-1))/grid%h(3)
    if(periodic(1))then
      gx(1,:,:)=bx(1,:,:)*(p(1,:,:)-p(nx,:,:))/grid%h(1);gx(nx+1,:,:)=gx(1,:,:)
    else
      if(dirichlet(1))gx(1,:,:)=bx(1,:,:)*p(1,:,:)/(grid%h(1)*.5_dp)
      if(dirichlet(2))gx(nx+1,:,:)=-bx(nx+1,:,:)*p(nx,:,:)/(grid%h(1)*.5_dp)
    endif
    if(periodic(2))then
      gy(:,1,:)=by(:,1,:)*(p(:,1,:)-p(:,ny,:))/grid%h(2);gy(:,ny+1,:)=gy(:,1,:)
    else
      if(dirichlet(3))gy(:,1,:)=by(:,1,:)*p(:,1,:)/(grid%h(2)*.5_dp)
      if(dirichlet(4))gy(:,ny+1,:)=-by(:,ny+1,:)*p(:,ny,:)/(grid%h(2)*.5_dp)
    endif
    if(periodic(3))then
      gz(:,:,1)=bz(:,:,1)*(p(:,:,1)-p(:,:,nz))/grid%h(3);gz(:,:,nz+1)=gz(:,:,1)
    else
      if(dirichlet(5))gz(:,:,1)=bz(:,:,1)*p(:,:,1)/(grid%h(3)*.5_dp)
      if(dirichlet(6))gz(:,:,nz+1)=-bz(:,:,nz+1)*p(:,:,nz)/(grid%h(3)*.5_dp)
    endif
  end subroutine
  subroutine divergence(grid,u,v,w,d)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::u(grid%nx+1,grid%ny,grid%nz),v(grid%nx,grid%ny+1,grid%nz),w(grid%nx,grid%ny,grid%nz+1)
    real(dp),intent(out)::d(grid%nx,grid%ny,grid%nz)
    d=(u(2:grid%nx+1,:,:)-u(1:grid%nx,:,:))/grid%h(1) &
     +(v(:,2:grid%ny+1,:)-v(:,1:grid%ny,:))/grid%h(2) &
     +(w(:,:,2:grid%nz+1)-w(:,:,1:grid%nz))/grid%h(3)
  end subroutine
  subroutine apply_operator(grid,p,bx,by,bz,dirichlet,periodic,ap)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::p(grid%nx,grid%ny,grid%nz),bx(grid%nx+1,grid%ny,grid%nz), &
                        by(grid%nx,grid%ny+1,grid%nz),bz(grid%nx,grid%ny,grid%nz+1)
    logical,intent(in)::dirichlet(6),periodic(3)
    real(dp),intent(out)::ap(grid%nx,grid%ny,grid%nz)
    real(dp)::gx(grid%nx+1,grid%ny,grid%nz),gy(grid%nx,grid%ny+1,grid%nz),gz(grid%nx,grid%ny,grid%nz+1)
    call weighted_gradient(grid,p,bx,by,bz,dirichlet,periodic,gx,gy,gz)
    call divergence(grid,gx,gy,gz,ap);ap=-ap
  end subroutine
  subroutine remove_mean(a)
    real(dp),intent(inout)::a(:,:,:)
    a=a-sum(a)/real(size(a),dp)
  end subroutine
  subroutine pcg(grid,b,bx,by,bz,diagonal,dirichlet,periodic,p,rtol,atol,maxiter,s)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::b(grid%nx,grid%ny,grid%nz),bx(grid%nx+1,grid%ny,grid%nz),by(grid%nx,grid%ny+1,grid%nz), &
                        bz(grid%nx,grid%ny,grid%nz+1),diagonal(grid%nx,grid%ny,grid%nz),rtol,atol
    logical,intent(in)::dirichlet(6),periodic(3)
    real(dp),intent(inout)::p(grid%nx,grid%ny,grid%nz)
    integer,intent(in)::maxiter
    type(solve_stats),intent(out)::s
    real(dp)::b0(grid%nx,grid%ny,grid%nz),r(grid%nx,grid%ny,grid%nz),z(grid%nx,grid%ny,grid%nz), &
              direction(grid%nx,grid%ny,grid%nz),ad(grid%nx,grid%ny,grid%nz),true_r(grid%nx,grid%ny,grid%nz)
    real(dp)::threshold,rz,rz_new,denom,lambda,recursive_rms,cells
    integer::k
    logical::singular
    cells=real(grid%nx*grid%ny*grid%nz,dp);singular=all(periodic)
    s=solve_stats();s%rhs_rms=sqrt(sum(b*b)/cells);s%rhs_mean=sum(b)/cells
    threshold=max(atol,rtol*s%rhs_rms);b0=b
    if(singular)then
      s%compatibility_tolerance=1e-12_dp*max(1._dp,s%rhs_rms)
      if(abs(s%rhs_mean)>s%compatibility_tolerance)then
        s%flag=4;call apply_operator(grid,p,bx,by,bz,dirichlet,periodic,ad)
        s%absres=sqrt(sum((b-ad)**2)/cells);s%original_absres=s%absres
        s%relres=s%absres/max(s%rhs_rms,atol);s%stop_ratio=s%absres/threshold;s%pressure_mean=sum(p)/cells;return
      endif
      b0=b-s%rhs_mean;call remove_mean(p)
    endif
    call apply_operator(grid,p,bx,by,bz,dirichlet,periodic,ad);r=b0-ad
    s%absres=sqrt(sum(r*r)/cells);s%original_absres=sqrt(sum((b-ad)**2)/cells)
    if(s%absres<=threshold)then
      s%relres=s%absres/max(s%rhs_rms,atol);s%stop_ratio=s%absres/threshold;s%pressure_mean=sum(p)/cells;return
    endif
    if(singular)call remove_mean(r)
    z=r/diagonal;if(singular)call remove_mean(z)
    direction=z;rz=sum(r*z);s%flag=1
    do k=1,maxiter
      s%iterations=k;call apply_operator(grid,direction,bx,by,bz,dirichlet,periodic,ad);denom=sum(direction*ad)
      if(denom<=0.or.rz<=0.or..not.ieee_is_finite(denom).or..not.ieee_is_finite(rz))then
        s%flag=2;exit
      endif
      lambda=rz/denom;p=p+lambda*direction;r=r-lambda*ad
      recursive_rms=sqrt(sum(r*r)/cells)
      if(recursive_rms<=threshold)then
        if(singular)call remove_mean(p)
        call apply_operator(grid,p,bx,by,bz,dirichlet,periodic,ad);true_r=b0-ad;s%absres=sqrt(sum(true_r*true_r)/cells)
        if(s%absres<=threshold)then
          s%flag=0;exit
        endif
        r=true_r;if(singular)call remove_mean(r)
        z=r/diagonal;if(singular)call remove_mean(z)
        direction=z;rz=sum(r*z);s%restarts=s%restarts+1;cycle
      endif
      z=r/diagonal;if(singular)call remove_mean(z)
      rz_new=sum(r*z)
      if(.not.ieee_is_finite(rz_new))then
        s%flag=2;exit
      endif
      direction=z+(rz_new/rz)*direction
      if(singular)call remove_mean(direction)
      rz=rz_new
    enddo
    if(singular)call remove_mean(p)
    call apply_operator(grid,p,bx,by,bz,dirichlet,periodic,ad);true_r=b0-ad
    s%absres=sqrt(sum(true_r*true_r)/cells);s%original_absres=sqrt(sum((b-ad)**2)/cells)
    s%relres=s%absres/max(s%rhs_rms,atol);s%stop_ratio=s%absres/threshold;s%pressure_mean=sum(p)/cells
    if(.not.ieee_is_finite(s%absres))s%flag=3
  end subroutine
  function face_inner(grid,a,b,c,x,y,z,bx,by,bz)result(value)
    type(rect_grid),intent(in)::grid
    real(dp),intent(in)::a(grid%nx+1,grid%ny,grid%nz),b(grid%nx,grid%ny+1,grid%nz),c(grid%nx,grid%ny,grid%nz+1), &
                        x(grid%nx+1,grid%ny,grid%nz),y(grid%nx,grid%ny+1,grid%nz),z(grid%nx,grid%ny,grid%nz+1), &
                        bx(grid%nx+1,grid%ny,grid%nz),by(grid%nx,grid%ny+1,grid%nz),bz(grid%nx,grid%ny,grid%nz+1)
    real(dp)::value
    value=grid%cell_volume*(sum(a*x/bx)-.5_dp*sum(a(1,:,:)*x(1,:,:)/bx(1,:,:)) &
      -.5_dp*sum(a(grid%nx+1,:,:)*x(grid%nx+1,:,:)/bx(grid%nx+1,:,:)) &
      +sum(b*y/by)-.5_dp*sum(b(:,1,:)*y(:,1,:)/by(:,1,:))-.5_dp*sum(b(:,grid%ny+1,:)*y(:,grid%ny+1,:)/by(:,grid%ny+1,:)) &
      +sum(c*z/bz)-.5_dp*sum(c(:,:,1)*z(:,:,1)/bz(:,:,1))-.5_dp*sum(c(:,:,grid%nz+1)*z(:,:,grid%nz+1)/bz(:,:,grid%nz+1)))
  end function
end module rect_mac_projection3d
