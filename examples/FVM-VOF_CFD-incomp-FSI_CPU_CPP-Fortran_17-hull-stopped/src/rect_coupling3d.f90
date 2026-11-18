! Actual phase -> mass -> donor-cell MAC momentum -> updated-rho projection.
module rect_coupling3d
  use rect_transport3d
  use rect_mac_projection3d,only:boundary_beta,solve_stats,allocate_boundary,build_faces,weighted_gradient, &
                                divergence,apply_operator,pcg,face_inner
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type::step_stats
    real(dp)::amin=0,amax=0,cell_closure=0,dual_closure=0,rho_alpha=0,momentum_local=0,momentum_budget=0
    real(dp)::phase_budget=0,mass_budget=0,phase_boundary=0,mass_boundary=0
    real(dp)::source_phase=0,source_mass=0,source_momentum=0,source_identity_phase=0,source_identity_mass=0,source_identity_momentum=0
    real(dp)::source_API=0,sweep_source_max=0,phase_activity(3)=0,volume_activity(3)=0,adv_div=0
    real(dp)::ke_old=0,ke_adv=0,ke_gravity=0,ke_projected=0,gravity_work_new=0,gravity_BE_defect=0,gravity_identity=0
    real(dp)::projection_correction_ke=0,projection_identity=0,projection_orthogonality=0,weighted_adjoint=0
    real(dp)::projection_momentum_identity=0,projection_momentum_change(3)=0
    real(dp)::div_before=0,div_after=0,div_residual_identity=0,beta_density_defect=0,seam_error=0,fullmac_ke_defect=0,liftmax=0
    integer::tiny_calls=0
    type(solve_stats)::solver
  end type
contains
  subroutine pack(g,a,u,v,w)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::a(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(out)::u(g%nx+1,g%ny,g%nz),v(g%nx,g%ny+1,g%nz),w(g%nx,g%ny,g%nz+1)
    u=a(0:g%nx,1:g%ny,1:g%nz,1);v=a(1:g%nx,0:g%ny,1:g%nz,2);w=a(1:g%nx,1:g%ny,0:g%nz,3)
  end subroutine
  subroutine unpack(g,u,v,w,a)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::u(g%nx+1,g%ny,g%nz),v(g%nx,g%ny+1,g%nz),w(g%nx,g%ny,g%nz+1)
    real(dp),intent(out)::a(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    a=0;a(0:g%nx,1:g%ny,1:g%nz,1)=u;a(1:g%nx,0:g%ny,1:g%nz,2)=v;a(1:g%nx,1:g%ny,0:g%nz,3)=w
  end subroutine
  subroutine padded_divergence(g,a,d)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::a(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(out)::d(g%nx,g%ny,g%nz)
    d=(a(1:g%nx,1:g%ny,1:g%nz,1)-a(0:g%nx-1,1:g%ny,1:g%nz,1))/g%h(1) &
     +(a(1:g%nx,1:g%ny,1:g%nz,2)-a(1:g%nx,0:g%ny-1,1:g%nz,2))/g%h(2) &
     +(a(1:g%nx,1:g%ny,1:g%nz,3)-a(1:g%nx,1:g%ny,0:g%nz-1,3))/g%h(3)
  end subroutine
  subroutine project(g,periodic,rho,dt,m,p,s,correction,rhs,true_r)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::rho(g%nx,g%ny,g%nz),dt
    real(dp),intent(inout)::m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),p(g%nx,g%ny,g%nz)
    type(step_stats),intent(inout)::s
    real(dp),intent(out)::correction(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),rhs(g%nx,g%ny,g%nz),true_r(g%nx,g%ny,g%nz)
    type(boundary_beta)::b
    real(dp)::rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),before(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
              after(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),lift(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
              grad(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),beta(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::u(g%nx+1,g%ny,g%nz),v(g%nx,g%ny+1,g%nz),w(g%nx,g%ny,g%nz+1), &
              bx(g%nx+1,g%ny,g%nz),by(g%nx,g%ny+1,g%nz),bz(g%nx,g%ny,g%nz+1), &
              gx(g%nx+1,g%ny,g%nz),gy(g%nx,g%ny+1,g%nz),gz(g%nx,g%ny,g%nz+1)
    real(dp)::diag(g%nx,g%ny,g%nz),d(g%nx,g%ny,g%nz),dlift(g%nx,g%ny,g%nz),ap(g%nx,g%ny,g%nz),mom0(3),mom1(3)
    integer::c,lo(3),hi(3)
    logical::bc(6)
    if(any(periodic).and..not.all(periodic))error stop 'transport stage tests fully periodic or tank only'
    call velocities(g,periodic,rho,m,before,rf);call pack(g,before,u,v,w)
    call integrated_momentum(g,m,mom0);call allocate_boundary(g,b)
    b%xlo=1/rf(0,1:g%ny,1:g%nz,1);b%xhi=1/rf(g%nx,1:g%ny,1:g%nz,1)
    b%ylo=1/rf(1:g%nx,0,1:g%nz,2);b%yhi=1/rf(1:g%nx,g%ny,1:g%nz,2)
    b%zlo=1/rf(1:g%nx,1:g%ny,0,3);b%zhi=1/rf(1:g%nx,1:g%ny,g%nz,3)
    bc=.false.;if(.not.all(periodic))bc(4)=.true.
    call build_faces(g,rho,b,bc,periodic,bx,by,bz,diag)
    call unpack(g,bx,by,bz,beta)
    do c=1,3
      call bounds(g,c,lo,hi)
      s%beta_density_defect=max(s%beta_density_defect,maxval(abs( &
        beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)*rf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)-1)))
    enddo
    lift=0
    if(.not.all(periodic))then
      lift(0,1:g%ny,1:g%nz,1)=before(0,1:g%ny,1:g%nz,1)/dt
      lift(g%nx,1:g%ny,1:g%nz,1)=before(g%nx,1:g%ny,1:g%nz,1)/dt
      lift(1:g%nx,0,1:g%nz,2)=before(1:g%nx,0,1:g%nz,2)/dt
      lift(1:g%nx,1:g%ny,0,3)=before(1:g%nx,1:g%ny,0,3)/dt
      lift(1:g%nx,1:g%ny,g%nz,3)=before(1:g%nx,1:g%ny,g%nz,3)/dt
    endif
    s%liftmax=maxval(abs(lift));call padded_divergence(g,before,d);call padded_divergence(g,lift,dlift)
    s%div_before=maxval(abs(d));rhs=-d/dt+dlift
    s%fullmac_ke_defect=abs(.5_dp*face_inner(g,u,v,w,u,v,w,bx,by,bz)-s%ke_gravity)
    call pcg(g,rhs,bx,by,bz,diag,bc,periodic,p,1e-12_dp,1e-12_dp,4000,s%solver)
    if(s%solver%flag/=0)return
    call weighted_gradient(g,p,bx,by,bz,bc,periodic,gx,gy,gz);call unpack(g,gx,gy,gz,grad)
    correction=grad+lift;after=before-dt*correction
    s%projection_momentum_identity=maxval(abs(rf*after-(m-dt*rf*correction)))
    m=rf*after;call integrated_momentum(g,m,mom1);s%projection_momentum_change=mom1-mom0
    call padded_divergence(g,after,d);s%div_after=maxval(abs(d))
    call apply_operator(g,p,bx,by,bz,bc,periodic,ap);true_r=rhs-ap
    s%div_residual_identity=maxval(abs(d+dt*true_r))
    s%ke_projected=.5_dp*weighted_dot(g,rf,after,after)
    s%projection_correction_ke=.5_dp*weighted_dot(g,rf,after-before,after-before)
    s%projection_identity=s%ke_gravity-s%ke_projected-s%projection_correction_ke+dt*sum(p*d)*g%cell_volume
    s%projection_orthogonality=weighted_dot(g,rf,after,after-before)
    s%weighted_adjoint=weighted_dot(g,rf,after,grad)+sum(p*d)*g%cell_volume
    if(periodic(1))s%seam_error=max(s%seam_error,maxval(abs(after(0,1:g%ny,1:g%nz,1)-after(g%nx,1:g%ny,1:g%nz,1))))
    if(periodic(2))s%seam_error=max(s%seam_error,maxval(abs(after(1:g%nx,0,1:g%nz,2)-after(1:g%nx,g%ny,1:g%nz,2))))
    if(periodic(3))s%seam_error=max(s%seam_error,maxval(abs(after(1:g%nx,1:g%ny,0,3)-after(1:g%nx,1:g%ny,g%nz,3))))
  end subroutine
  subroutine coupled_step(g,periodic,a,rho,m,adv,dt,rhog,drho,gravity,order,p,tr,s,sunit,step,time,mpredict,correction,rhs,true_r)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    integer,intent(in)::order(3),sunit,step
    real(dp),intent(inout)::a(g%nx,g%ny,g%nz),rho(g%nx,g%ny,g%nz),m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),p(g%nx,g%ny,g%nz)
    real(dp),intent(in)::adv(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),dt,rhog,drho,gravity,time
    type(rect_transport),intent(inout)::tr
    type(step_stats),intent(out)::s
    real(dp),intent(out)::mpredict(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),correction(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
                         rhs(g%nx,g%ny,g%nz),true_r(g%nx,g%ny,g%nz)
    type(sweep_stats)::sw
    real(dp)::rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),freeze(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
              uadv(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),ug(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
              grav(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),momsource(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::div(g%nx,g%ny,g%nz),expected(g%nx,g%ny,g%nz),dv,sm,expect_m,source_p(g%nx,g%ny,g%nz),source_m(g%nx,g%ny,g%nz)
    integer::sweep,d,c,i,j,k,pos(3),left(3),lo(3),hi(3)
    s=step_stats();call allocate_transport(g,tr);tr%cc=merge(1._dp,0._dp,a>.5_dp);tr%order=order
    tr%source=0;tr%mass_source=0;tr%phase=0;tr%volume=0;momsource=0
    call velocities(g,periodic,rho,m,freeze,rf);s%ke_old=.5_dp*weighted_dot(g,rf,freeze,freeze)
    s%amin=minval(a);s%amax=maxval(a);call padded_divergence(g,adv,div);s%adv_div=maxval(abs(div))
    do sweep=1,3
      d=order(sweep)
      call coupled_sweep(g,periodic,a,rho,m,adv,freeze,tr%cc,dt,rhog,drho,d, &
        tr%volume(:,:,:,d),tr%phase(:,:,:,d),tr%mass_source(:,:,:,d),sw)
      do k=1,g%nz;do j=1,g%ny;do i=1,g%nx
        pos=[i,j,k];left=pos;left(d)=pos(d)-1
        dv=tr%volume(i,j,k,d)-tr%volume(left(1),left(2),left(3),d)
        tr%source(i,j,k,d)=tr%cc(i,j,k)*dv
        s%source_API=max(s%source_API,abs(tr%mass_source(i,j,k,d)-(rhog*dv+drho*tr%source(i,j,k,d))))
      enddo;enddo;enddo
      do c=1,3
        call bounds(g,c,lo,hi)
        do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
          pos=[i,j,k];sm=dual_source(g,periodic,tr%mass_source(:,:,:,d),c,pos)
          momsource(i,j,k,c)=momsource(i,j,k,c)+sm*freeze(i,j,k,c)
        enddo;enddo;enddo
      enddo
      s%amin=min(s%amin,sw%amin);s%amax=max(s%amax,sw%amax)
      s%cell_closure=max(s%cell_closure,sw%cell_closure);s%dual_closure=max(s%dual_closure,sw%dual_closure)
      s%rho_alpha=max(s%rho_alpha,sw%rho_alpha);s%momentum_local=max(s%momentum_local,sw%momentum_local)
      s%momentum_budget=max(s%momentum_budget,sw%momentum_budget)
      s%phase_budget=max(s%phase_budget,abs(sw%phase_budget));s%mass_budget=max(s%mass_budget,abs(sw%mass_budget))
      s%phase_boundary=s%phase_boundary+sw%phase_boundary;s%mass_boundary=s%mass_boundary+sw%mass_boundary
      s%sweep_source_max=max(s%sweep_source_max,sw%source_density);s%tiny_calls=s%tiny_calls+sw%tiny_calls
      s%phase_activity(d)=sw%phase_activity;s%volume_activity(d)=sw%volume_activity
      write(sunit,'(*(es25.16e3,:,","))')real(step,dp),time,real(sweep,dp),real(d,dp),sw%amin,sw%amax, &
        sw%cell_closure,sw%dual_closure,sw%rho_alpha,sw%momentum_local,sw%phase_budget,sw%mass_budget,sw%momentum_budget, &
        sw%phase_source,sw%mass_source,sw%phase_boundary,sw%mass_boundary,sw%momentum_source,sw%momentum_flux,sw%wall_reaction, &
        sw%source_density,sw%phase_activity,sw%volume_activity,real(sw%tiny_calls,dp)
    enddo
    source_p=sum(tr%source,dim=4);source_m=sum(tr%mass_source,dim=4)
    s%source_phase=maxval(abs(source_p))/g%cell_volume;s%source_mass=maxval(abs(source_m))/g%cell_volume
    expected=tr%cc*dt*div*g%cell_volume
    s%source_identity_phase=maxval(abs(source_p-expected))/g%cell_volume
    expected=(rhog+drho*tr%cc)*dt*div*g%cell_volume
    s%source_identity_mass=maxval(abs(source_m-expected))/g%cell_volume
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        pos=[i,j,k];expect_m=dual_source(g,periodic,expected,c,pos)*freeze(i,j,k,c)
        s%source_momentum=max(s%source_momentum,abs(momsource(i,j,k,c))/(g%cell_volume*dual_scale(g,periodic,c,pos)))
        s%source_identity_momentum=max(s%source_identity_momentum, &
          abs(momsource(i,j,k,c)-expect_m)/(g%cell_volume*dual_scale(g,periodic,c,pos)))
      enddo;enddo;enddo
    enddo
    if(any(rho<=0).or..not.all(ieee_is_finite(m)))error stop 'invalid transported material state'
    call velocities(g,periodic,rho,m,uadv,rf);s%ke_adv=.5_dp*weighted_dot(g,rf,uadv,uadv)
    grav=0;grav(1:g%nx,0:g%ny,1:g%nz,2)=-gravity
    m=m+dt*rf*grav;mpredict=m;call velocities(g,periodic,rho,m,ug,rf)
    s%ke_gravity=.5_dp*weighted_dot(g,rf,ug,ug)
    s%gravity_work_new=dt*weighted_dot(g,rf,ug,grav)
    s%gravity_BE_defect=.5_dp*weighted_dot(g,rf,ug-uadv,ug-uadv)
    s%gravity_identity=s%ke_gravity-s%ke_adv-s%gravity_work_new+s%gravity_BE_defect
    call project(g,periodic,rho,dt,m,p,s,correction,rhs,true_r)
  end subroutine
end module rect_coupling3d
