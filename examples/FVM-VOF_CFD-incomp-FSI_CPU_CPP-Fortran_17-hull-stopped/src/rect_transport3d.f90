! Rectangular constant-per-axis spacing full-MAC phase/mass/momentum transport.
! No clipping, redistribution, velocity reset or production LICA source.
module rect_transport3d
  use geometry3d,only:dp,qp,plane_offset,signed_flux
  use rect_mac_projection3d,only:rect_grid
  implicit none
  real(dp),parameter::pure_tol=1.e-14_dp
  type::rect_transport
    real(dp),allocatable::volume(:,:,:,:),phase(:,:,:,:),source(:,:,:,:),mass_source(:,:,:,:),cc(:,:,:)
    integer::order(3)=0
  end type
  type::sweep_stats
    real(dp)::amin=0,amax=0,cell_closure=0,dual_closure=0,rho_alpha=0,momentum_local=0
    real(dp)::phase_budget=0,mass_budget=0,momentum_budget=0,phase_source=0,mass_source=0
    real(dp)::phase_boundary=0,mass_boundary=0,momentum_source(3)=0,momentum_flux(3)=0,wall_reaction(3)=0
    real(dp)::phase_activity=0,volume_activity=0,source_density=0
    integer::tiny_calls=0
  end type
contains
  subroutine bounds(g,c,lo,hi)
    type(rect_grid),intent(in)::g
    integer,intent(in)::c
    integer,intent(out)::lo(3),hi(3)
    lo=1;hi=[g%nx,g%ny,g%nz];lo(c)=0
  end subroutine
  pure integer function wrap(i,n)result(r)
    integer,intent(in)::i,n
    r=modulo(i-1,n)+1
  end function
  pure integer function reflect(i,n)result(r)
    integer,intent(in)::i,n
    r=i
    if(i<1)r=1-i
    if(i>n)r=2*n+1-i
  end function
  pure real(dp) function weight(g,c,p)result(w)
    type(rect_grid),intent(in)::g
    integer,intent(in)::c,p(3)
    integer::dims(3)
    dims=[g%nx,g%ny,g%nz];w=1
    if(p(c)==0.or.p(c)==dims(c))w=.5_dp
  end function
  pure real(dp) function dual_scale(g,periodic,c,p)result(w)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    integer,intent(in)::c,p(3)
    w=1
    if(.not.periodic(c))w=weight(g,c,p)
  end function
  pure logical function fixed_face(g,periodic,c,p)result(fixed)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    integer,intent(in)::c,p(3)
    integer::dims(3)
    dims=[g%nx,g%ny,g%nz]
    fixed=.not.periodic(c).and.(p(c)==0.or.p(c)==dims(c))
    if(c==2.and.p(2)==g%ny)fixed=.false. ! Top atmospheric pressure outlet.
  end function
  function phase_at(g,periodic,a,p)result(f)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::a(g%nx,g%ny,g%nz)
    integer,intent(in)::p(3)
    real(dp)::f
    integer::r(3),dims(3),d
    dims=[g%nx,g%ny,g%nz];r=p
    if(.not.periodic(2).and.p(2)>g%ny)then
      f=0;return ! Pure-air reservoir, not periodic wrapping.
    endif
    do d=1,3
      if(periodic(d))then
        r(d)=wrap(p(d),dims(d))
      else
        r(d)=reflect(p(d),dims(d))
      endif
    enddo
    f=a(r(1),r(2),r(3))
  end function
  subroutine normal_at(g,periodic,a,p,normal)
    ! Youngs differences are in CELL-INDEX coordinates: these are unit-cube
    ! PLIC coefficients, already proportional to physical gradient*h_axis.
    ! Do NOT divide by h here and pass a physical normal to cube geometry.
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::a(g%nx,g%ny,g%nz)
    integer,intent(in)::p(3)
    real(dp),intent(out)::normal(3)
    real(dp)::w,norm1
    integer::s,t,left(3),right(3),c,d,e
    normal=0
    do c=1,3
      d=mod(c,3)+1;e=mod(c+1,3)+1
      do t=-1,1;do s=-1,1
        left=p;right=p;left(c)=p(c)-1;right(c)=p(c)+1
        left(d)=p(d)+s;right(d)=p(d)+s;left(e)=p(e)+t;right(e)=p(e)+t
        w=real((2-abs(s))*(2-abs(t)),dp)
        normal(c)=normal(c)+w*(phase_at(g,periodic,a,left)-phase_at(g,periodic,a,right))
      enddo;enddo
    enddo
    norm1=sum(abs(normal))
    if(norm1>1.e-30_dp)then
      normal=normal/norm1
    else
      normal=[0._dp,1._dp,0._dp]
    endif
  end subroutine
  subroutine allocate_transport(g,tr)
    type(rect_grid),intent(in)::g
    type(rect_transport),intent(inout)::tr
    if(allocated(tr%volume))return
    allocate(tr%volume(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
      tr%phase(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),tr%source(g%nx,g%ny,g%nz,3), &
      tr%mass_source(g%nx,g%ny,g%nz,3),tr%cc(g%nx,g%ny,g%nz))
  end subroutine
  subroutine density_faces(g,periodic,rho,rf)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::rho(g%nx,g%ny,g%nz)
    real(dp),intent(out)::rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    integer::i,j,k,c,p(3),left(3),right(3),lo(3),hi(3),dims(3)
    dims=[g%nx,g%ny,g%nz];rf=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        p=[i,j,k];left=p;right=p;right(c)=p(c)+1
        if(periodic(c))then
          left(c)=wrap(left(c),dims(c));right(c)=wrap(right(c),dims(c))
        else
          left(c)=max(1,left(c));right(c)=min(dims(c),right(c))
        endif
        rf(i,j,k,c)=.5_dp*(rho(left(1),left(2),left(3))+rho(right(1),right(2),right(3)))
      enddo;enddo;enddo
    enddo
  end subroutine
  subroutine velocities(g,periodic,rho,m,u,rf)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::rho(g%nx,g%ny,g%nz),m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(out)::u(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    integer::c,lo(3),hi(3)
    call density_faces(g,periodic,rho,rf);u=0
    do c=1,3
      call bounds(g,c,lo,hi)
      u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)= &
        m(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)/rf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)
    enddo
  end subroutine
  function weighted_dot(g,rf,a,b)result(value)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::rf(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),a(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
                         b(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::value
    integer::c,i,j,k,lo(3),hi(3)
    value=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        value=value+weight(g,c,[i,j,k])*rf(i,j,k,c)*a(i,j,k,c)*b(i,j,k,c)
      enddo;enddo;enddo
    enddo
    value=value*g%cell_volume
  end function
  subroutine integrated_momentum(g,m,mom)
    type(rect_grid),intent(in)::g
    real(dp),intent(in)::m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(out)::mom(3)
    integer::c,i,j,k,lo(3),hi(3)
    mom=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        mom(c)=mom(c)+weight(g,c,[i,j,k])*m(i,j,k,c)*g%cell_volume
      enddo;enddo;enddo
    enddo
  end subroutine
  function dual_source(g,periodic,s,c,p)result(value)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::s(g%nx,g%ny,g%nz)
    integer,intent(in)::c,p(3)
    real(dp)::value
    integer::left(3),right(3),dims(3)
    dims=[g%nx,g%ny,g%nz];left=p;right=p;right(c)=p(c)+1
    if(periodic(c))then
      left(c)=wrap(left(c),dims(c));right(c)=wrap(right(c),dims(c))
      value=.5_dp*(s(left(1),left(2),left(3))+s(right(1),right(2),right(3)))
    else if(p(c)==0)then
      value=.5_dp*s(right(1),right(2),right(3))
    else if(p(c)==dims(c))then
      value=.5_dp*s(left(1),left(2),left(3))
    else
      value=.5_dp*(s(left(1),left(2),left(3))+s(right(1),right(2),right(3)))
    endif
  end function
  function sample_flux(g,periodic,f,d,p)result(value)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::f(0:g%nx+1,0:g%ny+1,0:g%nz+1)
    integer,intent(in)::d,p(3)
    real(dp)::value
    integer::r(3),dims(3),c
    dims=[g%nx,g%ny,g%nz];r=p
    do c=1,3
      if(periodic(c))then
        if(c==d)then
          r(c)=modulo(p(c),dims(c))
        else
          r(c)=wrap(p(c),dims(c))
        endif
      endif
    enddo
    value=f(r(1),r(2),r(3))
  end function
  function dual_flux(g,periodic,f,c,d,p,positive)result(value)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3),positive
    real(dp),intent(in)::f(0:g%nx+1,0:g%ny+1,0:g%nz+1)
    integer,intent(in)::c,d,p(3)
    real(dp)::value
    integer::face(3),other(3),dims(3),side
    dims=[g%nx,g%ny,g%nz];face=p
    if(c==d)then
      side=merge(1,-1,positive)
      if(.not.periodic(c))then
        if((positive.and.p(c)==dims(c)).or.(.not.positive.and.p(c)==0))then
          value=sample_flux(g,periodic,f,d,p);return
        endif
      endif
      other=p;other(c)=p(c)+side
      value=.5_dp*(sample_flux(g,periodic,f,d,p)+sample_flux(g,periodic,f,d,other))
    else
      if(.not.positive)face(d)=face(d)-1
      other=face;other(c)=p(c)+1
      if(.not.periodic(c).and.p(c)==0)then
        value=.5_dp*sample_flux(g,periodic,f,d,other)
      else if(.not.periodic(c).and.p(c)==dims(c))then
        value=.5_dp*sample_flux(g,periodic,f,d,face)
      else
        value=.5_dp*(sample_flux(g,periodic,f,d,face)+sample_flux(g,periodic,f,d,other))
      endif
    endif
  end function
  function velocity_at(g,periodic,u,c,p)result(value)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::u(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    integer,intent(in)::c,p(3)
    real(dp)::value
    integer::r(3),dims(3),d
    dims=[g%nx,g%ny,g%nz];r=p
    if(.not.periodic(2).and.p(2)>g%ny)then
      if(c==2)then
        r(2)=g%ny;value=u(r(1),r(2),r(3),c)
      else
        value=0 ! Incoming reservoir has zero tangential momentum.
      endif
      return
    endif
    do d=1,3
      if(periodic(d))then
        if(c==d)then
          r(d)=modulo(p(d),dims(d))
        else
          r(d)=wrap(p(d),dims(d))
        endif
      else if(c/=d)then
        r(d)=reflect(p(d),dims(d))
      endif
    enddo
    if(any(r<0).or.any(r>dims))error stop 'velocity ghost index'
    value=u(r(1),r(2),r(3),c)
  end function
  subroutine phase_flux(g,periodic,a,adv,dt,d,fvol,fphase,calls)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    real(dp),intent(in)::a(g%nx,g%ny,g%nz),adv(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),dt
    integer,intent(in)::d
    real(dp),intent(out)::fvol(0:g%nx+1,0:g%ny+1,0:g%nz+1),fphase(0:g%nx+1,0:g%ny+1,0:g%nz+1)
    integer,intent(out)::calls
    integer::lo(3),hi(3),i,j,k,p(3),donor(3),dims(3),c
    real(dp)::courant,q,normal(3)
    real(qp)::offset
    fvol=0;fphase=0;calls=0;dims=[g%nx,g%ny,g%nz]
    call bounds(g,d,lo,hi)
    do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
      p=[i,j,k]
      if(fixed_face(g,periodic,d,p))cycle
      courant=adv(i,j,k,d)*dt/g%h(d)
      if(abs(courant)>.5_dp+1e-12_dp)error stop 'directional Courant exceeds .5'
      donor=p;if(courant<0)donor(d)=donor(d)+1
      fvol(i,j,k)=adv(i,j,k,d)*dt*(g%cell_volume/g%h(d))
      q=phase_at(g,periodic,a,donor)
      if(courant==0)then
        fphase(i,j,k)=0
      else if(q<=pure_tol.or.q>=1-pure_tol)then
        fphase(i,j,k)=courant*q*g%cell_volume
        if(q/=0.and.q/=1)calls=calls+1 ! Retained raw tails, no clamp.
      else
        do c=1,3
          if(periodic(c))donor(c)=wrap(donor(c),dims(c))
        enddo
        call normal_at(g,periodic,a,donor,normal);offset=plane_offset(normal,q)
        fphase(i,j,k)=signed_flux(normal,offset,d,courant)*g%cell_volume
      endif
    enddo;enddo;enddo
  end subroutine
  subroutine coupled_sweep(g,periodic,a,rho,m,adv,freeze,cc,dt,rhog,drho,d,fvol,fphase,source,s)
    type(rect_grid),intent(in)::g
    logical,intent(in)::periodic(3)
    integer,intent(in)::d
    real(dp),intent(inout)::a(g%nx,g%ny,g%nz),rho(g%nx,g%ny,g%nz),m(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp),intent(in)::adv(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),freeze(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
                         cc(g%nx,g%ny,g%nz),dt,rhog,drho
    real(dp),intent(out)::fvol(0:g%nx+1,0:g%ny+1,0:g%nz+1),fphase(0:g%nx+1,0:g%ny+1,0:g%nz+1),source(g%nx,g%ny,g%nz)
    type(sweep_stats),intent(out)::s
    real(dp)::oldrho(g%nx,g%ny,g%nz),fmass(0:g%nx+1,0:g%ny+1,0:g%nz+1),phase0,mass0,mom0(3),mom1(3)
    real(dp)::rf0(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),rf1(0:g%nx+1,0:g%ny+1,0:g%nz+1,3), &
              u(0:g%nx+1,0:g%ny+1,0:g%nz+1,3),updated(0:g%nx+1,0:g%ny+1,0:g%nz+1,3)
    real(dp)::dv,vol,scale,own,fr,fl,vr,vl,sm,moment_flux,moment_source,local_rhs
    integer::i,j,k,c,p(3),left(3),right(3),lo(3),hi(3),top(3),dims(3)
    s=sweep_stats();vol=g%cell_volume;dims=[g%nx,g%ny,g%nz];oldrho=rho;phase0=sum(a)*vol;mass0=sum(rho)*vol
    call integrated_momentum(g,m,mom0);call velocities(g,periodic,rho,m,u,rf0)
    call phase_flux(g,periodic,a,adv,dt,d,fvol,fphase,s%tiny_calls);fmass=rhog*fvol+drho*fphase
    do k=1,g%nz;do j=1,g%ny;do i=1,g%nx
      p=[i,j,k];left=p;left(d)=p(d)-1
      dv=fvol(i,j,k)-fvol(left(1),left(2),left(3))
      source(i,j,k)=(rhog+drho*cc(i,j,k))*dv
      a(i,j,k)=a(i,j,k)+(fphase(left(1),left(2),left(3))-fphase(i,j,k)+cc(i,j,k)*dv)/vol
      rho(i,j,k)=rho(i,j,k)+(fmass(left(1),left(2),left(3))-fmass(i,j,k)+source(i,j,k))/vol
      s%cell_closure=max(s%cell_closure,abs(rho(i,j,k)-oldrho(i,j,k) &
        +(fmass(i,j,k)-fmass(left(1),left(2),left(3))-source(i,j,k))/vol))
      s%phase_source=s%phase_source+cc(i,j,k)*dv
    enddo;enddo;enddo
    s%mass_source=sum(source);s%source_density=maxval(abs(source))/vol
    call density_faces(g,periodic,rho,rf1);updated=m
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        p=[i,j,k];left=p;right=p;left(d)=p(d)-1;right(d)=p(d)+1
        scale=dual_scale(g,periodic,c,p);own=weight(g,c,p)/scale
        fr=dual_flux(g,periodic,fmass,c,d,p,.true.);fl=dual_flux(g,periodic,fmass,c,d,p,.false.)
        sm=dual_source(g,periodic,source,c,p)
        s%dual_closure=max(s%dual_closure,abs(rf1(i,j,k,c)-rf0(i,j,k,c)+(fr-fl-sm)/(vol*scale)))
        vr=u(i,j,k,c);if(fr<0)vr=velocity_at(g,periodic,u,c,right)
        ! At zero mass flux the donor is immaterial. In particular, do not
        ! ask for a nonexistent normal-velocity ghost beyond a fixed wall.
        vl=u(i,j,k,c);if(fl>0)vl=velocity_at(g,periodic,u,c,left)
        moment_flux=fr*vr-fl*vl;moment_source=sm*freeze(i,j,k,c)
        local_rhs=(moment_flux-moment_source)/(vol*scale)
        if(fixed_face(g,periodic,c,p))then
          s%wall_reaction(c)=s%wall_reaction(c)+own*(moment_flux-moment_source)
        else
          updated(i,j,k,c)=m(i,j,k,c)-local_rhs
          s%momentum_local=max(s%momentum_local,abs(updated(i,j,k,c)-m(i,j,k,c)+local_rhs))
        endif
        s%momentum_flux(c)=s%momentum_flux(c)+own*moment_flux
        s%momentum_source(c)=s%momentum_source(c)+own*moment_source
      enddo;enddo;enddo
    enddo
    m=updated;call integrated_momentum(g,m,mom1)
    s%momentum_budget=maxval(abs(mom1-mom0+s%momentum_flux-s%momentum_source-s%wall_reaction))
    call bounds(g,d,lo,hi);hi(d)=0
    do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
      p=[i,j,k];top=p;top(d)=dims(d)
      s%phase_boundary=s%phase_boundary+fphase(top(1),top(2),top(3))-fphase(i,j,k)
      s%mass_boundary=s%mass_boundary+fmass(top(1),top(2),top(3))-fmass(i,j,k)
    enddo;enddo;enddo
    s%phase_budget=sum(a)*vol-phase0+s%phase_boundary-s%phase_source
    s%mass_budget=sum(rho)*vol-mass0+s%mass_boundary-s%mass_source
    s%rho_alpha=maxval(abs(rho-rhog-drho*a));s%amin=minval(a);s%amax=maxval(a)
    s%phase_activity=sum(abs(fphase));s%volume_activity=sum(abs(fvol))
    ! Correct duplicate seam double count for physical activity totals.
    if(periodic(d))then
      call bounds(g,d,lo,hi);hi(d)=0
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        s%phase_activity=s%phase_activity-abs(fphase(i,j,k));s%volume_activity=s%volume_activity-abs(fvol(i,j,k))
      enddo;enddo;enddo
    endif
  end subroutine
end module rect_transport3d
