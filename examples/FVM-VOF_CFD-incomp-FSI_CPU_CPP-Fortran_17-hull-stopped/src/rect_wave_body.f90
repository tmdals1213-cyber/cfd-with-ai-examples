! Rectangular physical PLIC moments/reduced gravity and ONE same-step
! stationary implicit Brinkman / pressure solve. No moving body or FSI.
module rect_wave_body
  use rect_transport3d
  use rect_coupling3d,only:pack,unpack,padded_divergence
  use plane_moments3d,only:cube_phase_moments,cube_interface_polygon
  use rect_mac_projection3d,only:boundary_beta,solve_stats,allocate_boundary,build_faces,pcg,weighted_gradient,apply_operator
  implicit none
  real(dp),parameter::rg=1._dp,dr=999._dp,gravity=9.81_dp,kappa=1.e6_dp
  type::body_stats
    real(dp)::beta_audit=0,diag_audit=0,momentum_residual=0,momentum_scale=1,divergence=0,div_residual=0,wall_normal=0,liftmax=0
    real(dp)::KE_adv=0,KE_new=0,BE_defect=0,red_work=0,p_work=0,penalty_work=0,penalty_dissipation=0,energy_residual=0
    real(dp)::inside_RMS=0,inside_max=0,IB_force(3)=0,IB_Mz=0,IB_max=0
  end type
contains
  subroutine geometric_state(g,a,Z,pe,phi,valid,volume_error)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::a(g%nx,g%ny,g%nz),Z
    real(dp),intent(out)::pe,phi(g%nx,g%ny,g%nz),volume_error
    logical,intent(out)::valid(g%nx,g%ny,g%nz)
    integer::i,j,k
    real(dp)::n(3),mom(3),V,area,centroid(3),y0,flat
    real(qp)::offset
    logical::has_interface,periodic(3)
    periodic=.false.;pe=0;phi=0;valid=.false.;volume_error=0
    do k=1,g%nz;do j=1,g%ny;do i=1,g%nx
      y0=g%origin(2)+(j-1)*g%h(2)
      flat=max(0._dp,min(1._dp,(Z-y0)/g%h(2))) ! ANALYTIC flat reference only.
      if(a(i,j,k)>pure_tol.and.a(i,j,k)<1-pure_tol)then
        call normal_at(g,periodic,a,[i,j,k],n);offset=plane_offset(n,a(i,j,k))
        call cube_phase_moments(n,offset,V,mom);volume_error=max(volume_error,abs(V-a(i,j,k)))
        call cube_interface_polygon(n,offset,has_interface,area,centroid)
        if(has_interface)then
          phi(i,j,k)=dr*gravity*(y0+g%h(2)*centroid(2)-Z);valid(i,j,k)=.true.
        endif
      else
        mom=.5_dp*a(i,j,k) ! Retained tiny tails, not reset/clipping.
      endif
      pe=pe+g%cell_volume*((y0-Z)*(a(i,j,k)-flat)+g%h(2)*(mom(2)-.5_dp*flat))
    enddo;enddo;enddo
    pe=dr*gravity*pe
  end subroutine
  subroutine reduced_force(g,a,phi,valid,f)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::a(g%nx,g%ny,g%nz),phi(g%nx,g%ny,g%nz)
    logical,intent(in)::valid(g%nx,g%ny,g%nz)
    real(dp),intent(out)::f(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    integer::i,j,k,c,p(3),right(3),dims(3)
    real(dp)::potential
    dims=[g%nx,g%ny,g%nz];f=0
    do c=1,3
      do k=1,g%nz;do j=1,g%ny;do i=1,g%nx
        p=[i,j,k];if(p(c)==dims(c))cycle
        right=p;right(c)=right(c)+1;potential=0
        if(valid(i,j,k).and.valid(right(1),right(2),right(3)))then
          potential=.5_dp*(phi(i,j,k)+phi(right(1),right(2),right(3)))
        else if(valid(i,j,k))then
          potential=phi(i,j,k)
        else if(valid(right(1),right(2),right(3)))then
          potential=phi(right(1),right(2),right(3))
        endif
        f(i,j,k,c)=potential*(a(right(1),right(2),right(3))-a(i,j,k))/g%h(c)
      enddo;enddo;enddo
    enddo
  end subroutine
  function face_sum(g,a)result(value)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::a(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::value
    integer::c,i,j,k,lo(3),hi(3)
    value=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        value=value+g%cell_volume*weight(g,c,[i,j,k])*a(i,j,k,c)
      enddo;enddo;enddo
    enddo
  end function
  subroutine body_pressure(g,rho,chi,dt,m,fred,q,cg,solver,s)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::rho(g%nx,g%ny,g%nz),chi(g%nx,g%ny,g%nz),dt,cg(3),fred(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(inout)::m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),q(g%nx,g%ny,g%nz)
    type(solve_stats),intent(out)::solver
    type(body_stats),intent(out)::s
    real(dp)::rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),chif(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      re(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),old(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      R(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),ut(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      new(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),lift(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      bg(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),beta(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      gradq(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),ib(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),tmp(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::rhoe(g%nx,g%ny,g%nz),diag(g%nx,g%ny,g%nz),diagcheck(g%nx,g%ny,g%nz),rhs(g%nx,g%ny,g%nz), &
      div(g%nx,g%ny,g%nz),dl(g%nx,g%ny,g%nz),ap(g%nx,g%ny,g%nz)
    real(dp)::bx(g%nx+1,g%ny,g%nz),by(g%nx,g%ny+1,g%nz),bz(g%nx,g%ny,g%nz+1), &
      gx(g%nx+1,g%ny,g%nz),gy(g%nx,g%ny+1,g%nz),gz(g%nx,g%ny,g%nz+1)
    type(boundary_beta)::b
    integer::c,i,j,k,lo(3),hi(3)
    real(dp)::massmask,x,y
    logical::periodic(3),bc(6)
    periodic=.false.;bc=[.false.,.false.,.false.,.true.,.false.,.false.];s=body_stats()
    call velocities(g,periodic,rho,m,old,rf);call density_faces(g,periodic,chi,chif)
    R=m+dt*fred;re=rf+dt*kappa*chif;rhoe=rho+dt*kappa*chi
    ut=0;call allocate_boundary(g,b)
    b%xlo=1/re(0,1:g%ny,1:g%nz,1);b%xhi=1/re(g%nx,1:g%ny,1:g%nz,1)
    b%ylo=1/re(1:g%nx,0,1:g%nz,2);b%yhi=1/re(1:g%nx,g%ny,1:g%nz,2)
    b%zlo=1/re(1:g%nx,1:g%ny,0,3);b%zhi=1/re(1:g%nx,1:g%ny,g%nz,3)
    call build_faces(g,rho,b,bc,periodic,bx,by,bz,diag,rho_effective=rhoe)
    call unpack(g,bx,by,bz,beta);lift=0
    do c=1,3
      call bounds(g,c,lo,hi)
      ut(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)=R(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)/re(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)
      s%beta_audit=max(s%beta_audit,maxval(abs(beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)* &
        re(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)-1)))
    enddo
    diagcheck=0
    do k=1,g%nz;do j=1,g%ny;do i=1,g%nx
      if(i>1)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i-1,j,k,1)/g%h(1)**2
      if(i<g%nx)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i,j,k,1)/g%h(1)**2
      if(j>1)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i,j-1,k,2)/g%h(2)**2
      if(j<g%ny)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i,j,k,2)/g%h(2)**2
      if(k>1)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i,j,k-1,3)/g%h(3)**2
      if(k<g%nz)diagcheck(i,j,k)=diagcheck(i,j,k)+1/re(i,j,k,3)/g%h(3)**2
      if(j==g%ny)diagcheck(i,j,k)=diagcheck(i,j,k)+2/re(i,j,k,2)/g%h(2)**2
    enddo;enddo;enddo
    s%diag_audit=maxval(abs(diag-diagcheck))/maxval(abs(diagcheck))
    lift(0,1:g%ny,1:g%nz,1)=ut(0,1:g%ny,1:g%nz,1)/dt;lift(g%nx,1:g%ny,1:g%nz,1)=ut(g%nx,1:g%ny,1:g%nz,1)/dt
    lift(1:g%nx,0,1:g%nz,2)=ut(1:g%nx,0,1:g%nz,2)/dt
    lift(1:g%nx,1:g%ny,0,3)=ut(1:g%nx,1:g%ny,0,3)/dt;lift(1:g%nx,1:g%ny,g%nz,3)=ut(1:g%nx,1:g%ny,g%nz,3)/dt
    s%liftmax=maxval(abs(lift));call padded_divergence(g,ut,div);call padded_divergence(g,lift,dl);rhs=-div/dt+dl
    call pcg(g,rhs,bx,by,bz,diag,bc,periodic,q,1e-12_dp,1e-13_dp,4000,solver)
    if(solver%flag/=0)return
    call weighted_gradient(g,q,bx,by,bz,bc,periodic,gx,gy,gz);call unpack(g,gx,gy,gz,bg);bg=bg+lift
    gradq=0
    do c=1,3
      call bounds(g,c,lo,hi)
      gradq(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)=bg(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)/beta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)
    enddo
    new=ut-dt*bg;ib=-kappa*chif*new
    tmp=rf*new-R+dt*gradq-dt*ib;s%momentum_residual=maxval(abs(tmp))
    s%momentum_scale=max(1._dp,maxval(abs(R)),dt*maxval(abs(gradq)),dt*maxval(abs(ib)))
    call padded_divergence(g,new,div);s%divergence=maxval(abs(div))
    call apply_operator(g,q,bx,by,bz,bc,periodic,ap);s%div_residual=maxval(abs(div+dt*(rhs-ap)))
    s%wall_normal=max(maxval(abs(new(0,1:g%ny,1:g%nz,1))),maxval(abs(new(g%nx,1:g%ny,1:g%nz,1))), &
      maxval(abs(new(1:g%nx,0,1:g%nz,2))),maxval(abs(new(1:g%nx,1:g%ny,0,3))),maxval(abs(new(1:g%nx,1:g%ny,g%nz,3))))
    s%KE_adv=.5_dp*weighted_dot(g,rf,old,old);s%KE_new=.5_dp*weighted_dot(g,rf,new,new)
    s%BE_defect=.5_dp*weighted_dot(g,rf,new-old,new-old)
    s%red_work=dt*face_sum(g,fred*new);s%p_work=-dt*face_sum(g,gradq*new)
    s%penalty_work=dt*face_sum(g,ib*new);s%penalty_dissipation=dt*kappa*face_sum(g,chif*new**2)
    s%energy_residual=s%KE_new-s%KE_adv+s%BE_defect-s%red_work-s%p_work-s%penalty_work
    massmask=sum(chi*rho)*g%cell_volume
    if(massmask>0)then
      s%inside_RMS=sqrt(face_sum(g,chif*rf*new**2)/massmask)
      s%inside_max=maxval(abs(new),mask=chif>0)
    endif
    s%IB_max=maxval(abs(ib))
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        s%IB_force(c)=s%IB_force(c)+weight(g,c,[i,j,k])*ib(i,j,k,c)*g%cell_volume
        x=g%origin(1)+(i-.5_dp)*g%h(1);y=g%origin(2)+(j-.5_dp)*g%h(2)
        if(c==1)s%IB_Mz=s%IB_Mz-(y-cg(2))*weight(g,c,[i,j,k])*ib(i,j,k,c)*g%cell_volume
        if(c==2)s%IB_Mz=s%IB_Mz+(x-cg(1))*weight(g,c,[i,j,k])*ib(i,j,k,c)*g%cell_volume
      enddo;enddo;enddo
    enddo
    m=rf*new ! Material momentum, never effective inertia times new velocity.
  end subroutine
end module rect_wave_body
