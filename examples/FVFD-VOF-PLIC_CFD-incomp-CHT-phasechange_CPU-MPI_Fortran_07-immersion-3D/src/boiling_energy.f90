module boiling_energy
 use precision,only:WP
 use face_geometry3d,only:clipped_quad_fraction
 use geometry,only:cfg,physical_span,Lz,block_height,block_width,block_depth,cold_depth,cold_solid,cold_width,cold_thickness,domain_height,pool_height,cold_film_thickness
 use vfs_class,only:vfs,VFlo
 use tpscalar_class,only:tpscalar
 use irl_fortran_interface,only:calculateCentroid,calculateNormal,getPlane
 use hypre_str_class,only:hypre_str,gmres_pfmg,gmres
 use mpi_f08
 use parallel,only:MPI_REAL_WP
 use messager,only:die
 use phase_birth_core,only:wall_site_activation,birth_ok=>ok
 implicit none
 real(WP),allocatable,target::T(:,:,:),E(:,:,:),qphase(:,:,:),theta_l(:,:,:),theta_g(:,:,:)
 real(WP),allocatable::El(:,:,:),Eg(:,:,:),Elold(:,:,:),Egold(:,:,:),HL(:,:,:),HG(:,:,:),inc(:,:,:,:)
 real(WP)::Tsat,hfg,rl,rg,cpl,cpg,kl,kg,power,thermal_out,phase_power,thermal_balance,ks,rhos,cps,wall_superheat,roundoff_energy
 real(WP),allocatable::save_El(:,:,:),save_Eg(:,:,:),save_E(:,:,:),save_T(:,:,:),save_l(:,:,:),save_g(:,:,:),save_s(:,:,:)
 type(hypre_str)::hsL,hsG,hsS
 real(WP),allocatable::theta_s(:,:,:),Sold(:,:,:),s_prev(:,:,:)
 logical::inventory_correction_enabled=.false.
 real(WP)::inventory_resolved_transfer_local=0
 real(WP),allocatable::inventory_scale(:,:,:),inventory_potential_q(:,:,:),inventory_energy_return(:,:,:)
 real(WP)::cold_temperature,cold_heat_out=0
 real(WP)::qvol,coupling_change,thermal_vf_min
 integer::inactive_thermal_interfaces
 real(WP)::inactive_thermal_surface_m2
 integer::coupling_iterations
 integer::wall_birth_enabled,wall_birth_events=0
 logical::wall_birth_done=.false.,wall_birth_trial=.false.
 real(WP)::wall_birth_superheat,wall_birth_z,wall_birth_x,wall_birth_width,wall_birth_height,wall_birth_sigma
 real(WP)::wall_birth_event_mass=0,wall_birth_latent_J=0,wall_birth_surface_J=0
 real(WP),allocatable::wall_birth_m(:,:,:),wall_birth_point(:,:,:,:)

 type(tpscalar)::adv
contains
 subroutine energy_init(vf)
  use param,only:param_read
  type(vfs),intent(in)::vf
  integer::i,j,k
  real(WP)::subcool
  call param_read('Saturation temperature',Tsat);call param_read('Latent heat',hfg)
  call param_read('Liquid density',rl);call param_read('Gas density',rg)
  call param_read('Liquid heat capacity',cpl);call param_read('Gas heat capacity',cpg)
  call param_read('Liquid conductivity',kl);call param_read('Gas conductivity',kg)
  call param_read('Copper conductivity',ks);call param_read('Copper density',rhos);call param_read('Copper heat capacity',cps)
  call param_read('Heat power',power);qvol=power/(block_width*block_height*block_depth)
  call param_read('Initial solid superheat',wall_superheat)
  call param_read('Cold plate temperature',cold_temperature)
  call param_read('Initial liquid subcooling',subcool)
  call param_read('Wall birth enabled',wall_birth_enabled)
  call param_read('Wall birth superheat',wall_birth_superheat)
  call param_read('Wall birth x',wall_birth_x)
  call param_read('Wall birth z',wall_birth_z)
  call param_read('Wall birth width',wall_birth_width)
  call param_read('Wall birth equivalent height',wall_birth_height)
  call param_read('Surface tension coefficient',wall_birth_sigma)
  if(wall_birth_enabled.lt.0.or.wall_birth_enabled.gt.1)call die('Wall birth enabled must be0 or1')
  if(wall_birth_enabled.eq.1)then
   if(wall_birth_superheat.le.0.or.wall_birth_width.le.0.or.wall_birth_height.le.0)call die('Explicit positive wall-site parameters required')
   if(abs(wall_birth_x)+.5_WP*wall_birth_width.gt..5_WP*block_width)call die('Prescribed birth site outside copper top')
   if(abs(wall_birth_z)+.5_WP*wall_birth_width.gt..5_WP*block_depth)call die('Prescribed birth site outside copper z extent')
   if(wall_birth_height.ge..5_WP*cfg%dy(cfg%jmin))call die('Prototype embryo must be thinner than half a first wall cell')
  endif
  call param_read('Thermal interface tolerance',thermal_vf_min)
  if(thermal_vf_min.lt.VFlo.or.thermal_vf_min.ge.0.5_WP)call die('Invalid thermal interface applicability tolerance')
  allocate(T(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
  allocate(inventory_scale,mold=T);allocate(inventory_potential_q,mold=T);allocate(inventory_energy_return,mold=T)
  inventory_scale=1;inventory_potential_q=0;inventory_energy_return=0
  allocate(wall_birth_m,mold=T)
  allocate(wall_birth_point(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
  wall_birth_m=0;wall_birth_point=0
  allocate(E,mold=T);allocate(qphase,mold=T);allocate(theta_l,mold=T);allocate(theta_g,mold=T)
  allocate(El,mold=T);allocate(Eg,mold=T);allocate(Elold,mold=T);allocate(Egold,mold=T)
  allocate(HL,mold=T);allocate(HG,mold=T)
  allocate(inc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,2))
  allocate(theta_s,mold=T);allocate(Sold,mold=T);allocate(s_prev,mold=T)
  allocate(save_El,mold=T);allocate(save_Eg,mold=T);allocate(save_E,mold=T);allocate(save_T,mold=T)
  allocate(save_l,mold=T);allocate(save_g,mold=T);allocate(save_s,mold=T)
  theta_l=-subcool;theta_g=0;theta_s=wall_superheat;qphase=0;roundoff_energy=0
  where(cfg%VF.eq.0)theta_l=wall_superheat
  where(cfg%VF.eq.0)theta_g=wall_superheat
  do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
   if(cold_solid(i,j,k))then
    theta_s(i,j,k)=cold_temperature-Tsat;theta_l(i,j,k)=theta_s(i,j,k);theta_g(i,j,k)=theta_s(i,j,k)
   else if(cfg%VF(i,j,k).gt.0.and.vf%VF(i,j,k).gt.VFlo.and.cfg%ym(j).gt.pool_height.and.abs(cfg%xm(i)).lt.cold_width/2.and.abs(cfg%zm(k)).lt.cold_depth/2)then
    theta_l(i,j,k)=(cold_temperature-Tsat)*max(0._WP,min(1._WP,(vf%Lbary(2,i,j,k)-(domain_height-cold_thickness-cold_film_thickness))/cold_film_thickness))
   endif
  enddo;enddo;enddo
  coupling_iterations=0;coupling_change=0
  call sensible_from_temperature(vf)
  hsL=hypre_str(cfg=cfg,name='liquid heat',method=gmres_pfmg,nst=7)
  hsG=hypre_str(cfg=cfg,name='vapor heat',method=gmres,nst=7)
  hsS=hypre_str(cfg=cfg,name='copper heat',method=gmres_pfmg,nst=7)
  call initialize_solver(hsL);call initialize_solver(hsG);call initialize_solver(hsS)
  call adv%initialize(cfg=cfg,nscalar=2,name='phase sensible energies')
  adv%phase=[0,1];adv%Prho=[rl,rg];adv%diff=0;call adv%setup()
 end subroutine
 subroutine initialize_solver(hs)
  type(hypre_str),intent(inout)::hs
  hs%stc(1,:)=[0,0,0];hs%stc(2,:)=[1,0,0];hs%stc(3,:)=[-1,0,0]
  hs%stc(4,:)=[0,1,0];hs%stc(5,:)=[0,-1,0];hs%stc(6,:)=[0,0,1];hs%stc(7,:)=[0,0,-1]
  hs%rcvg=1e-11_WP;hs%maxit=1000;hs%opr(1,:,:,:)=1;call hs%init()
 end subroutine
 subroutine ghost(a)
  real(WP),intent(inout)::a(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  integer::i,j,k
  if(.not.cfg%xper.and.cfg%iproc.eq.1)then
   do i=cfg%imino_,cfg%imin-1;a(i,:,:)=a(cfg%imin,:,:);enddo
  endif
  if(.not.cfg%xper.and.cfg%iproc.eq.cfg%npx)then
   do i=cfg%imax+1,cfg%imaxo_;a(i,:,:)=a(cfg%imax,:,:);enddo
  endif
  if(cfg%jproc.eq.1)then
   do j=cfg%jmino_,cfg%jmin-1;a(:,j,:)=a(:,cfg%jmin,:);enddo
  endif
  if(cfg%jproc.eq.cfg%npy)then
   do j=cfg%jmax+1,cfg%jmaxo_;a(:,j,:)=a(:,cfg%jmax,:);enddo
  endif
  if(.not.cfg%zper.and.cfg%kproc.eq.1)then
   do k=cfg%kmino_,cfg%kmin-1;a(:,:,k)=a(:,:,cfg%kmin);enddo
  endif
  if(.not.cfg%zper.and.cfg%kproc.eq.cfg%npz)then
   do k=cfg%kmax+1,cfg%kmaxo_;a(:,:,k)=a(:,:,cfg%kmax);enddo
  endif
  call cfg%sync(a)
 end subroutine
 real(WP) function phase_fraction(vf,i,j,k,ph)result(a)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k,ph
  a=vf%VF(i,j,k);if(ph.eq.1)a=1-a
 end function
 subroutine face_vertices(i,j,k,d,v,axis,h,area)
  integer,intent(in)::i,j,k,d
  real(WP),intent(out)::v(3,4),h,area
  integer,intent(out)::axis
  real(WP)::x0,x1,y0,y1,z0,z1,c
  x0=cfg%x(i);x1=cfg%x(i+1);y0=cfg%y(j);y1=cfg%y(j+1);z0=cfg%z(k);z1=cfg%z(k+1)
  axis=(d+1)/2
  select case(axis)
  case(1)
   c=x0;if(mod(d,2)==1)c=x1
   v(:,1)=[c,y0,z0];v(:,2)=[c,y1,z0];v(:,3)=[c,y1,z1];v(:,4)=[c,y0,z1];h=x1-x0;area=(y1-y0)*(z1-z0)
  case(2)
   c=y0;if(mod(d,2)==1)c=y1
   v(:,1)=[x0,c,z0];v(:,2)=[x1,c,z0];v(:,3)=[x1,c,z1];v(:,4)=[x0,c,z1];h=y1-y0;area=(x1-x0)*(z1-z0)
  case(3)
   c=z0;if(mod(d,2)==1)c=z1
   v(:,1)=[x0,y0,c];v(:,2)=[x1,y0,c];v(:,3)=[x1,y1,c];v(:,4)=[x0,y1,c];h=z1-z0;area=(x1-x0)*(y1-y0)
  case default
   call die('Invalid3D face direction')
  end select
 end subroutine
 real(WP) function face_liquid3d(vf,i,j,k,d)result(f)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k,d
  real(WP)::vertices(3,4),plane(4),h,area
  integer::axis
  if(vf%VF(i,j,k).le.VFlo)then
   f=0;return
  else if(vf%VF(i,j,k).ge.1-VFlo)then
   f=1;return
  endif
  call face_vertices(i,j,k,d,vertices,axis,h,area)
  plane=getPlane(vf%liquid_gas_interface(i,j,k),0)
  f=clipped_quad_fraction(vertices,plane,vf%Lbary(:,i,j,k))
  if(f< -1e-12_WP.or.f>1+1e-12_WP)call die('Invalid3D PLIC face fraction')
  f=max(0.0_WP,min(1.0_WP,f))
 end function
 real(WP) function coupling(vf,i,j,k,d,ph)result(c)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k,d,ph
  integer,parameter::di(6)=[1,-1,0,0,0,0],dj(6)=[0,0,1,-1,0,0],dk(6)=[0,0,0,0,1,-1],opp(6)=[2,1,4,3,6,5]
  integer::ii,jj,kk,axis
  real(WP)::vertices(3,4),center(3),f1,f2,a,kap,bi(3),bj(3),distance,h,hn,area
  ii=i+di(d);jj=j+dj(d);kk=k+dk(d);c=0
  if(ii<cfg%imin.or.ii>cfg%imax.or.jj<cfg%jmin.or.jj>cfg%jmax.or.kk<cfg%kmin.or.kk>cfg%kmax)return
  call face_vertices(i,j,k,d,vertices,axis,h,area);center=sum(vertices,dim=2)/4
  if(phase_fraction(vf,i,j,k,ph)<=VFlo)return
  f1=face_liquid3d(vf,i,j,k,d);if(ph==1)f1=1-f1
  a=f1
  if(cfg%VF(ii,jj,kk)>0)then
   if(phase_fraction(vf,ii,jj,kk,ph)<=VFlo)return
   f2=face_liquid3d(vf,ii,jj,kk,opp(d));if(ph==1)f2=1-f2
   a=min(f1,f2)
  endif
  if(a<=0)return
  kap=kl;bi=vf%Lbary(:,i,j,k)
  if(ph==1)then
   kap=kg;bi=vf%Gbary(:,i,j,k)
  endif
  if(cfg%VF(ii,jj,kk)==0)then
   hn=cfg%dx(ii);if(axis==2)hn=cfg%dy(jj);if(axis==3)hn=cfg%dz(kk)
   distance=abs(bi(axis)-center(axis))
   c=a*area/(cfg%vol(i,j,k)*(distance/kap+.5_WP*hn/ks))
  else
   bj=vf%Lbary(:,ii,jj,kk);if(ph==1)bj=vf%Gbary(:,ii,jj,kk)
   distance=abs(bi(axis)-bj(axis))
   if(distance<=0)call die('Nonpositive3D same-phase thermal donor distance')
   c=kap*a*area/(cfg%vol(i,j,k)*distance)
  endif
 end function
 subroutine interface_H(vf)
  type(vfs),intent(in)::vf
  integer::i,j,k,ierr,local_inactive
  real(WP)::p(3),n(3),dl,dg,nn,local_surface
  HL=0;HG=0;local_inactive=0;local_surface=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0.or.vf%VF(i,j,k).le.VFlo.or.vf%VF(i,j,k).ge.1-VFlo.or.vf%SD(i,j,k).le.0)cycle
   if(vf%VF(i,j,k).le.thermal_vf_min.or.vf%VF(i,j,k).ge.1-thermal_vf_min)then
    local_inactive=local_inactive+1;local_surface=local_surface+vf%SD(i,j,k)*cfg%vol(i,j,k)*physical_span/Lz
    cycle
   endif
   p=calculateCentroid(vf%interface_polygon(1,i,j,k));n=calculateNormal(vf%interface_polygon(1,i,j,k));nn=sqrt(sum(n*n))
   if(nn.le.0)call die('Invalid PLIC thermal normal')
   n=n/nn;dl=abs(dot_product(p-vf%Lbary(:,i,j,k),n));dg=abs(dot_product(vf%Gbary(:,i,j,k)-p,n))
   if(min(dl,dg).le.0)call die('Nonpositive phase-to-interface thermal distance')
   HL(i,j,k)=kl*vf%SD(i,j,k)/dl;HG(i,j,k)=kg*vf%SD(i,j,k)/dg
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local_inactive,inactive_thermal_interfaces,1,MPI_INTEGER,MPI_SUM,cfg%comm,ierr)
  call MPI_ALLREDUCE(local_surface,inactive_thermal_surface_m2,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
 end subroutine
 subroutine solve_phase(dt,vf,ph,hs,old,a,H)
  real(WP),intent(in)::dt,old(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),H(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  real(WP),intent(inout)::a(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  type(vfs),intent(in)::vf
  integer,intent(in)::ph
  type(hypre_str),intent(inout)::hs
  integer,parameter::di(6)=[1,-1,0,0,0,0],dj(6)=[0,0,1,-1,0,0],dk(6)=[0,0,0,0,1,-1],slot(6)=[2,3,4,5,6,7]
  integer::i,j,k,d,ii,jj,kk
  real(WP)::Cvol,c,diag
  hs%opr=0;hs%rhs=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)then
    hs%opr(1,i,j,k)=1;hs%rhs(i,j,k)=theta_s(i,j,k);cycle
   endif
   if(phase_fraction(vf,i,j,k,ph).le.VFlo)then
    hs%opr(1,i,j,k)=1;cycle
   endif
   Cvol=phase_fraction(vf,i,j,k,ph)*rl*cpl
   if(ph.eq.1)Cvol=phase_fraction(vf,i,j,k,ph)*rg*cpg
   hs%opr(1,i,j,k)=Cvol/dt+H(i,j,k);hs%rhs(i,j,k)=old(i,j,k)/dt
   do d=1,6
    ii=i+di(d);jj=j+dj(d);kk=k+dk(d);c=coupling(vf,i,j,k,d,ph)
    hs%opr(1,i,j,k)=hs%opr(1,i,j,k)+c
    if(cfg%VF(ii,jj,kk).eq.0)then
     hs%rhs(i,j,k)=hs%rhs(i,j,k)+c*theta_s(ii,jj,kk)
    else
     hs%opr(slot(d),i,j,k)=-c
    endif
   enddo
   diag=hs%opr(1,i,j,k);hs%opr(:,i,j,k)=hs%opr(:,i,j,k)/diag;hs%rhs(i,j,k)=hs%rhs(i,j,k)/diag
  enddo;enddo;enddo
  call hs%setup();hs%sol=a;call hs%solve()
  if(hs%rerr.gt.1.01_WP*hs%rcvg)call die('Two-temperature phase solve residual; case tolerance')
  a=hs%sol;call ghost(a)
 end subroutine
 subroutine sensible_from_temperature(vf)
  type(vfs),intent(in)::vf
  El=vf%VF*rl*cpl*theta_l;Eg=(1-vf%VF)*rg*cpg*theta_g
  where(cfg%VF.eq.0)El=rhos*cps*theta_s
  where(cfg%VF.eq.0)Eg=0
  E=El+Eg;T=Tsat+vf%VF*theta_l+(1-vf%VF)*theta_g
  where(cfg%VF.eq.0)T=Tsat+theta_s
 end subroutine
 subroutine solve_solid(dt,vf)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  integer,parameter::di(6)=[1,-1,0,0,0,0],dj(6)=[0,0,1,-1,0,0],dk(6)=[0,0,0,0,1,-1],slot(6)=[2,3,4,5,6,7],opp(6)=[2,1,4,3,6,5]
  integer::i,j,k,d,ii,jj,kk,ph,axis
  real(WP)::c,h,hn,diag,tau,area,vertices(3,4)
  hsS%opr=0;hsS%rhs=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k)>0)then
    hsS%opr(1,i,j,k)=1;cycle
   endif
   if(cold_solid(i,j,k))then
    hsS%opr(1,i,j,k)=1;hsS%rhs(i,j,k)=cold_temperature-Tsat;cycle
   endif
   hsS%opr(1,i,j,k)=rhos*cps/dt;hsS%rhs(i,j,k)=Sold(i,j,k)/dt+qvol
   do d=1,6
    ii=i+di(d);jj=j+dj(d);kk=k+dk(d)
    if(ii<cfg%imin.or.ii>cfg%imax.or.jj<cfg%jmin.or.jj>cfg%jmax.or.kk<cfg%kmin.or.kk>cfg%kmax)cycle
    if(cfg%VF(ii,jj,kk)==0)then
     call face_vertices(i,j,k,d,vertices,axis,h,area)
     hn=cfg%dx(ii);if(axis==2)hn=cfg%dy(jj);if(axis==3)hn=cfg%dz(kk)
     c=ks*area/(cfg%vol(i,j,k)*(.5_WP*(h+hn)))
     hsS%opr(1,i,j,k)=hsS%opr(1,i,j,k)+c;hsS%opr(slot(d),i,j,k)=-c
    else
     do ph=0,1
      c=coupling(vf,ii,jj,kk,opp(d),ph)*cfg%vol(ii,jj,kk)/cfg%vol(i,j,k)
      tau=theta_l(ii,jj,kk);if(ph==1)tau=theta_g(ii,jj,kk)
      hsS%opr(1,i,j,k)=hsS%opr(1,i,j,k)+c;hsS%rhs(i,j,k)=hsS%rhs(i,j,k)+c*tau
     enddo
    endif
   enddo
   diag=hsS%opr(1,i,j,k);hsS%opr(:,i,j,k)=hsS%opr(:,i,j,k)/diag;hsS%rhs(i,j,k)=hsS%rhs(i,j,k)/diag
  enddo;enddo;enddo
  call hsS%setup();hsS%sol=theta_s;call hsS%solve()
  if(hsS%rerr>1.01_WP*hsS%rcvg)call die('Solid3D heat solve residual')
  theta_s=hsS%sol;call ghost(theta_s)
 end subroutine
 subroutine energy_step(dt,vf,U,V,mdot)
  real(WP),intent(in)::dt,U(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),V(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  type(vfs),intent(in)::vf
  real(WP),intent(out)::mdot(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  integer,parameter::di(6)=[1,-1,0,0,0,0],dj(6)=[0,0,1,-1,0,0],dk(6)=[0,0,0,0,1,-1]
  integer::i,j,k,ierr,it,d,ii,jj,kk
  real(WP)::local(4),global(4),vol,change,c
  Elold=El;Egold=Eg;Sold=rhos*cps*theta_s;call interface_H(vf)
  do it=1,30
   s_prev=theta_s
   call solve_phase(dt,vf,0,hsL,Elold,theta_l,HL);call solve_phase(dt,vf,1,hsG,Egold,theta_g,HG)
   call solve_solid(dt,vf)
   change=0
   do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
    if(cfg%VF(i,j,k).eq.0)change=max(change,abs(theta_s(i,j,k)-s_prev(i,j,k)))
   enddo;enddo;enddo
   call MPI_ALLREDUCE(change,coupling_change,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
   if(coupling_change.lt.1e-10_WP)exit
  enddo
  coupling_iterations=it
  if(it.gt.30)call die('Three-field thermal block iteration failed; case coupling tolerance 1e-10 K')
  call sensible_from_temperature(vf);qphase=HL*theta_l+HG*theta_g;mdot=0
  where(vf%SD.gt.0.and.cfg%VF.gt.0)mdot=qphase/(hfg*vf%SD)
  call cfg%sync(mdot);call cfg%sync(qphase);local=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   vol=cfg%vol(i,j,k)*physical_span/Lz
   local(1)=local(1)+(El(i,j,k)+Eg(i,j,k)-Elold(i,j,k)-Egold(i,j,k))*vol/dt
   if(cfg%VF(i,j,k).gt.0)then
    local(2)=local(2)+qphase(i,j,k)*vol
   else
    if(.not.cold_solid(i,j,k))local(3)=local(3)+qvol*vol
   endif
  enddo;enddo;enddo
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   vol=cfg%vol(i,j,k)*physical_span/Lz
   do d=1,6
    ii=i+di(d);jj=j+dj(d);kk=k+dk(d)
    if(.not.cold_solid(ii,jj,kk))cycle
    c=coupling(vf,i,j,k,d,0);local(4)=local(4)+c*(theta_l(i,j,k)-(cold_temperature-Tsat))*vol
    c=coupling(vf,i,j,k,d,1);local(4)=local(4)+c*(theta_g(i,j,k)-(cold_temperature-Tsat))*vol
   enddo
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,4,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  phase_power=global(2);cold_heat_out=global(4);thermal_out=cold_heat_out;thermal_balance=global(1)+global(2)-global(3)+cold_heat_out
 end subroutine
 subroutine save_thermal_trial()
  save_El=El;save_Eg=Eg;save_E=E;save_T=T;save_l=theta_l;save_g=theta_g;save_s=theta_s
 end subroutine
 subroutine restore_thermal_trial()
  El=save_El;Eg=save_Eg;E=save_E;T=save_T;theta_l=save_l;theta_g=save_g;theta_s=save_s
 end subroutine
 subroutine phase_conversion_cfl(dt,vf,cfl)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  real(WP),intent(out)::cfl
  real(WP)::local,a,da,m,expansion
  integer::i,j,k,ierr
  local=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0.or.qphase(i,j,k).eq.0)cycle
   m=qphase(i,j,k)/hfg;da=dt*m/rl
   a=vf%VF(i,j,k)
   if(m.lt.0)then
    a=1-a;da=dt*m/rg
   endif
   expansion=dt*abs(m*(1/rg-1/rl))
   local=max(local,expansion)
   if(a.le.0)call die('Nonzero resolved conversion source in absent phase before trial commit')
   local=max(local,abs(da)/a)
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,cfl,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
 end subroutine
 subroutine recover_phase_temperature(vf)
  type(vfs),intent(in)::vf
  integer::i,j,k,ierr
  real(WP)::a,b,bound,local_round,global_round,scale_l,scale_g
  scale_l=max(1.0_WP,maxval(abs(theta_l)));scale_g=max(1.0_WP,maxval(abs(theta_g)));local_round=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   a=vf%VF(i,j,k);b=1-a
   if(a.le.VFlo)then
    bound=2*VFlo*rl*cpl*scale_l
    if(abs(El(i,j,k)).gt.bound)call die('Resolved liquid energy in a phase removed by conversion')
    local_round=local_round+abs(El(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    Eg(i,j,k)=Eg(i,j,k)+El(i,j,k);El(i,j,k)=0
   endif
   if(b.le.VFlo)then
    bound=2*VFlo*rg*cpg*scale_g
    if(abs(Eg(i,j,k)).gt.bound)call die('Resolved gas energy in a phase removed by conversion')
    local_round=local_round+abs(Eg(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    El(i,j,k)=El(i,j,k)+Eg(i,j,k);Eg(i,j,k)=0
   endif
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local_round,global_round,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);roundoff_energy=roundoff_energy+global_round
  theta_l=0;theta_g=0
  where(vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_l=El/(vf%VF*rl*cpl)
  where(1-vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_g=Eg/((1-vf%VF)*rg*cpg)
  where(cfg%VF.eq.0)theta_l=theta_s
  where(cfg%VF.eq.0)theta_g=theta_s
  call ghost(theta_l);call ghost(theta_g);E=El+Eg;call cfg%sync(E)
  T(:,:,:)=Tsat+vf%VF*theta_l+(1-vf%VF)*theta_g;where(cfg%VF.eq.0)T=Tsat+theta_s
 end subroutine
 subroutine energy_remap(dt,vf,carrier_u,carrier_v,carrier_w)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  real(WP),intent(in)::carrier_u(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),carrier_v(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),carrier_w(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  real(WP),allocatable::u(:,:,:),v(:,:,:),w(:,:,:)
  integer::i,j,k,ierr
  real(WP)::a,b,local,global,bound,local_round,global_round,scale_l,scale_g
  allocate(u,mold=T);allocate(v,mold=T);allocate(w,mold=T)
  u=carrier_u;v=carrier_v;w=carrier_w
  adv%SC(:,:,:,1)=rl*cpl*theta_l;adv%SC(:,:,:,2)=rg*cpg*theta_g
  where(cfg%VF.eq.0)adv%SC(:,:,:,1)=0
  where(cfg%VF.eq.0)adv%SC(:,:,:,2)=0
  if(cfg%jproc.eq.cfg%npy)then
   j=cfg%jmax
   do k=cfg%kmino_,cfg%kmaxo_;do i=cfg%imino_,cfg%imaxo_
    if(v(i,j+1,k).lt.0)adv%SC(i,j+1:,k,:)=0
   enddo;enddo
  endif
  adv%SCold=adv%SC
  call adv%get_dSCdt(inc,U=u,V=v,W=w,VFold=vf%VFold,VF=vf%VF,detailed_face_flux=vf%detailed_face_flux,dt=dt)
  scale_l=max(1.0_WP,maxval(abs(theta_l)));scale_g=max(1.0_WP,maxval(abs(theta_g)))
  local=0;local_round=0;inventory_resolved_transfer_local=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   local=local+(inc(i,j,k,1)+inc(i,j,k,2))*cfg%vol(i,j,k)*physical_span/Lz
   El(i,j,k)=El(i,j,k)+dt*inc(i,j,k,1);Eg(i,j,k)=Eg(i,j,k)+dt*inc(i,j,k,2)
   a=vf%VF(i,j,k);b=1-a
   if(a.le.VFlo)then
    bound=2*VFlo*rl*cpl*scale_l
    if(abs(El(i,j,k)).gt.bound)then
     if(.not.(inventory_correction_enabled.and.vf%VFold(i,j,k).gt.VFlo.and.(qphase(i,j,k).gt.0.or.inventory_scale(i,j,k).lt.1)))call die('Resolved liquid sensible energy left in absent phase')
     inventory_resolved_transfer_local=inventory_resolved_transfer_local+El(i,j,k)*cfg%vol(i,j,k)*physical_span/Lz
     if(cfg%amRoot)print*,'RESOLVED_PHASE_DISAPPEARANCE_ENTHALPY_TRANSFER',i,j,k,El(i,j,k)
    endif
    if(abs(El(i,j,k)).le.bound)local_round=local_round+abs(El(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    Eg(i,j,k)=Eg(i,j,k)+El(i,j,k);El(i,j,k)=0;theta_l(i,j,k)=0
   else
    theta_l(i,j,k)=El(i,j,k)/(a*rl*cpl)
   endif
   if(b.le.VFlo)then
    bound=2*VFlo*rg*cpg*scale_g
    if(abs(Eg(i,j,k)).gt.bound)then
     if(.not.(inventory_correction_enabled.and.1-vf%VFold(i,j,k).gt.VFlo.and.(qphase(i,j,k).lt.0.or.inventory_scale(i,j,k).lt.1)))call die('Resolved vapor sensible energy left in absent phase')
     inventory_resolved_transfer_local=inventory_resolved_transfer_local+Eg(i,j,k)*cfg%vol(i,j,k)*physical_span/Lz
     if(cfg%amRoot)print*,'RESOLVED_PHASE_DISAPPEARANCE_ENTHALPY_TRANSFER',i,j,k,Eg(i,j,k)
    endif
    if(abs(Eg(i,j,k)).le.bound)local_round=local_round+abs(Eg(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    El(i,j,k)=El(i,j,k)+Eg(i,j,k);Eg(i,j,k)=0;theta_g(i,j,k)=0
   else
    theta_g(i,j,k)=Eg(i,j,k)/(b*rg*cpg)
   endif
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);thermal_out=cold_heat_out-global
  call MPI_ALLREDUCE(local_round,global_round,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);roundoff_energy=roundoff_energy+global_round
  where(vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_l=El/(vf%VF*rl*cpl)
  where(1-vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_g=Eg/((1-vf%VF)*rg*cpg)
  call ghost(theta_l);call ghost(theta_g);E=El+Eg;call cfg%sync(E)
  T(:,:,:)=Tsat+vf%VF*theta_l+(1-vf%VF)*theta_g;where(cfg%VF.eq.0)T=Tsat+theta_s
  deallocate(u,v,w)
 end subroutine

 subroutine wall_birth_thermal_trial(dt,vf)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  real(WP)::dm,qlat,vol,local_mass,global_mass,local_power,global_power
  integer::i,j,k,stat,ierr
  logical::activate
  wall_birth_m=0;wall_birth_point=0;wall_birth_trial=.false.;wall_birth_event_mass=0
  if(wall_birth_enabled.eq.0.or.wall_birth_done)return
  local_mass=0;local_power=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(j.le.cfg%jmin.or.cfg%VF(i,j,k).le.0.or.cfg%VF(i,j-1,k).ne.0)cycle
   if(abs(cfg%ym(j)-block_height-.5_WP*cfg%dy(j)).gt.64*epsilon(1._WP)*cfg%dy(j))cycle
   if((cfg%xm(i)-wall_birth_x)**2+(cfg%zm(k)-wall_birth_z)**2>(.5_WP*wall_birth_width)**2)cycle
   call wall_site_activation(Tsat+theta_s(i,j-1,k),Tsat+theta_l(i,j,k),Tsat,Tsat+wall_birth_superheat,.true.,1-vf%VF(i,j,k).gt.1e-12_WP,.true.,activate,stat)
   if(stat.ne.birth_ok)call die('Invalid prescribed wall birth policy parameters')
   if(.not.activate)cycle
   ! Finite prescribed embryo budget. This is not microscopic nucleation prediction.
   dm=rg*wall_birth_height/cfg%dy(j)*cfg%vol(i,j,k)
   qlat=dm*hfg/cfg%vol(i,j-1,k)
   if(rhos*cps*theta_s(i,j-1,k).lt.qlat)call die('Insufficient local copper heat for embryo latent cost; no clipping')
   theta_s(i,j-1,k)=theta_s(i,j-1,k)-qlat/(rhos*cps)
   wall_birth_m(i,j,k)=dm/(dt*cfg%vol(i,j,k))
   wall_birth_point(:,i,j,k)=[cfg%xm(i),block_height+.5_WP*wall_birth_height,cfg%zm(k)]
   qphase(i,j,k)=qphase(i,j,k)+wall_birth_m(i,j,k)*hfg
   local_mass=local_mass+dm*physical_span/Lz
   local_power=local_power+dm*hfg/dt*physical_span/Lz
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local_mass,global_mass,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  call MPI_ALLREDUCE(local_power,global_power,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  wall_birth_event_mass=global_mass;wall_birth_trial=global_mass.gt.0
  if(wall_birth_trial)then
   phase_power=phase_power+global_power
   call ghost(theta_s);call sensible_from_temperature(vf);call cfg%sync(qphase)
  endif
 end subroutine
 subroutine wall_birth_commit()
  if(.not.wall_birth_trial)return
  wall_birth_done=.true.;wall_birth_events=wall_birth_events+1
  wall_birth_latent_J=wall_birth_latent_J+wall_birth_event_mass*hfg
 end subroutine
 subroutine wall_birth_surface_cost(dt,vf,step,tnow)
  real(WP),intent(in)::dt,tnow
  type(vfs),intent(in)::vf
  integer,intent(in)::step
  real(WP)::local_area,area,work,volf,vols,share,energy
  integer::i,j,k,ierr,u
  if(.not.wall_birth_trial)return
  local_area=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   ! Existing pool/cold-film interfaces are not created by this local birth event.
   if(abs(cfg%xm(i)-wall_birth_x).gt..5_WP*wall_birth_width+cfg%dx(i))cycle
   if(cfg%ym(j).gt.block_height+2*cfg%dy(j).or.cfg%ym(j).lt.block_height)cycle
   if(cfg%VF(i,j,k).gt.0)local_area=local_area+vf%SD(i,j,k)*cfg%vol(i,j,k)*physical_span/Lz
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local_area,area,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  work=wall_birth_sigma*area
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(wall_birth_m(i,j,k).le.0)cycle
   volf=cfg%vol(i,j,k)*physical_span/Lz;vols=cfg%vol(i,j-1,k)*physical_span/Lz
   share=wall_birth_m(i,j,k)*dt*volf/wall_birth_event_mass;energy=work*share/vols
   if(rhos*cps*theta_s(i,j-1,k).lt.energy)call die('Insufficient local copper heat for initial interface surface cost')
   theta_s(i,j-1,k)=theta_s(i,j-1,k)-energy/(rhos*cps)
  enddo;enddo;enddo
  wall_birth_surface_J=wall_birth_surface_J+work
  call ghost(theta_s)
  where(cfg%VF.eq.0)El=rhos*cps*theta_s
  E=El+Eg;where(cfg%VF.eq.0)T=Tsat+theta_s
  call cfg%sync(E);call cfg%sync(T)
  if(cfg%amRoot)then
   open(newunit=u,file='WALL_BIRTH_EVENTS.csv',status='new')
   write(u,'(a)')'step,time_s,event_count,created_mass_kg,latent_cost_J,initial_surface_area_m2,surface_cost_J'
   write(u,'(i0,",",es25.16e3,",",i0,4(",",es25.16e3))')step,tnow,wall_birth_events,wall_birth_event_mass,wall_birth_event_mass*hfg,area,work
   close(u)
   print*,'PRESCRIBED_WALL_BIRTH_COMMITTED',step,tnow,wall_birth_event_mass,work
  endif
 end subroutine


 subroutine inventory_limit_prepare(dt,vf)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  integer::i,j,k,ierr
  real(WP)::local,global
  if(.not.inventory_correction_enabled)return
  inventory_potential_q=qphase;inventory_energy_return=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0.or.wall_birth_m(i,j,k).ne.0)cycle
   if(inventory_scale(i,j,k).lt.0.or.inventory_scale(i,j,k).gt.1)call die('Invalid inventory active-set scale')
   qphase(i,j,k)=inventory_potential_q(i,j,k)*inventory_scale(i,j,k)
   inventory_energy_return(i,j,k)=dt*(inventory_potential_q(i,j,k)-qphase(i,j,k))
  enddo;enddo;enddo
  call cfg%sync(qphase);call cfg%sync(inventory_energy_return)
  local=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).gt.0)local=local+qphase(i,j,k)*cfg%vol(i,j,k)*physical_span/Lz
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);phase_power=global
 end subroutine
 pure function bounded_conversion_scale(available,requested)result(scale)
  real(WP),intent(in)::available,requested
  real(WP)::scale
  scale=1.0_WP
  if(requested.gt.max(0.0_WP,available))scale=max(0.0_WP,available)/requested
 end function
 subroutine inventory_limit_update(vf,dt,reduce_only,changed)
  type(vfs),intent(in)::vf
  real(WP),intent(in)::dt
  logical,intent(in)::reduce_only
  logical,intent(out)::changed
  integer::i,j,k,ierr,ph
  real(WP)::local,global,potential,available,rho,desired,delta_fraction,requested
  integer::local_range,global_range
  real(WP),allocatable::target(:,:,:)
  allocate(target,mold=inventory_scale);target=inventory_scale;local=0;local_range=0
  if(cfg%amRoot)then
   print*,'INVENTORY_UPDATE_ENTER',dt,reduce_only
   flush(6)
  endif
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0.or.wall_birth_m(i,j,k).ne.0)cycle
   potential=inventory_potential_q(i,j,k)/hfg
   if(potential.eq.0)cycle
   if(.not.reduce_only.and.inventory_scale(i,j,k).eq.1)cycle
   ph=1;rho=rl
   if(potential.lt.0)then
    ph=2;rho=rg
   endif
   available=vf%phase_available_before_source(ph,i,j,k)
   if(available.lt.-VFlo)cycle
   requested=dt*abs(potential)/rho
   ! Compare before division. The uncapped quotient may overflow even though
   ! its mathematical bounded value is exactly one. No physical threshold.
   desired=bounded_conversion_scale(available,requested)
   if(requested.eq.0.or.(requested.lt.tiny(1.0_WP).and.available.gt.0))local_range=local_range+1
   if(reduce_only)desired=min(desired,inventory_scale(i,j,k))
   target(i,j,k)=desired
   delta_fraction=requested*abs(desired-inventory_scale(i,j,k))
   local=max(local,delta_fraction)
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
  call MPI_ALLREDUCE(local_range,global_range,1,MPI_INTEGER,MPI_SUM,cfg%comm,ierr)
  if(cfg%amRoot)then
   print*,'INVENTORY_UPDATE_REDUCED',global,global_range
   flush(6)
  endif
  ! Match the already-used native representable phase-volume tolerance, not a
  ! new physical cutoff. Reprojection/flux are recomputed if the active set changes.
  changed=global.gt.VFlo
  if(changed)then
   inventory_scale=target;call cfg%sync(inventory_scale)
  endif
  deallocate(target)
 end subroutine
 subroutine inventory_return_sensible(vf,step,t)
  type(vfs),intent(in)::vf
  integer,intent(in)::step
  real(WP),intent(in)::t
  integer::i,j,k,ierr,u,local_count,global_count
  logical::exists
  real(WP)::delta,local,global,transfer_global,local_remaining,global_remaining
  if(.not.inventory_correction_enabled)return
  local=0;local_count=0;local_remaining=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   delta=inventory_energy_return(i,j,k)
   if(delta.eq.0)cycle
   ! The heat solve removed the potential interface heat. Conversion uses the
   ! available mass only. The difference remains sensible energy, never disappears.
   if(inventory_potential_q(i,j,k).lt.0)then
    if(vf%VF(i,j,k).gt.VFlo)then
     El(i,j,k)=El(i,j,k)+delta
    else
     Eg(i,j,k)=Eg(i,j,k)+delta
    endif
   else
    if(1-vf%VF(i,j,k).gt.VFlo)then
     Eg(i,j,k)=Eg(i,j,k)+delta
    else
     El(i,j,k)=El(i,j,k)+delta
    endif
   endif
   local=local+delta*cfg%vol(i,j,k)*physical_span/Lz;local_count=local_count+1
   if(inventory_potential_q(i,j,k).lt.0)then
    local_remaining=max(local_remaining,1-vf%VF(i,j,k))
   else
    local_remaining=max(local_remaining,vf%VF(i,j,k))
   endif
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  call MPI_ALLREDUCE(local_count,global_count,1,MPI_INTEGER,MPI_SUM,cfg%comm,ierr)
  call MPI_ALLREDUCE(inventory_resolved_transfer_local,transfer_global,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  call MPI_ALLREDUCE(local_remaining,global_remaining,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
  call recover_phase_temperature(vf)
  if(cfg%amRoot)then
   inquire(file='INVENTORY_ENTHALPY_RETURN.csv',exist=exists)
   open(newunit=u,file='INVENTORY_ENTHALPY_RETURN.csv',status='unknown',position='append')
   if(.not.exists)write(u,'(a)')'step,time_s,limited_cells,sensible_return_J,disappearance_sensible_transfer_J,max_remaining_consumed_fraction'
   write(u,'(i0,",",es25.16e3,",",i0,3(",",es25.16e3))')step,t,global_count,global,transfer_global,global_remaining
   close(u)
  endif
 end subroutine

end module
