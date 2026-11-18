module boiling_energy
 use precision,only:WP
 use geometry,only:cfg,physical_span,Lz,block_height
 use vfs_class,only:vfs,VFlo
 use tpscalar_class,only:tpscalar
 use irl_fortran_interface,only:calculateCentroid,calculateNormal,getPlane
 use hypre_str_class,only:hypre_str,gmres_pfmg
 use mpi_f08
 use parallel,only:MPI_REAL_WP
 use messager,only:die
 implicit none
 real(WP),allocatable::T(:,:,:),E(:,:,:),qphase(:,:,:),theta_l(:,:,:),theta_g(:,:,:)
 real(WP),allocatable::El(:,:,:),Eg(:,:,:),Elold(:,:,:),Egold(:,:,:),HL(:,:,:),HG(:,:,:),inc(:,:,:,:)
 real(WP)::Tsat,hfg,rl,rg,cpl,cpg,kl,kg,power,thermal_out,phase_power,thermal_balance,ks,rhos,cps,wall_superheat,roundoff_energy
 type(hypre_str)::hsL,hsG
 type(tpscalar)::adv
contains
 subroutine energy_init(vf)
  use param,only:param_read
  type(vfs),intent(in)::vf
  integer::i,j,k
  real(WP)::lam,delta,ag,tshift,yref
  call param_read('Saturation temperature',Tsat);call param_read('Latent heat',hfg)
  call param_read('Liquid density',rl);call param_read('Gas density',rg)
  call param_read('Liquid heat capacity',cpl);call param_read('Gas heat capacity',cpg)
  call param_read('Liquid conductivity',kl);call param_read('Gas conductivity',kg)
  call param_read('Copper conductivity',ks);call param_read('Copper density',rhos);call param_read('Copper heat capacity',cps)
  call param_read('Heat power',power);if(power.ne.0)call die('This two-temperature branch is fixed-wall reference only')
  call param_read('Initial solid superheat',wall_superheat)
  call param_read('Reference lambda',lam);call param_read('Initial vapor thickness',delta)
  ag=kg/(rg*cpg);tshift=(delta/(2*lam))**2/ag
  allocate(T(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
  allocate(E,mold=T);allocate(qphase,mold=T);allocate(theta_l,mold=T);allocate(theta_g,mold=T)
  allocate(El,mold=T);allocate(Eg,mold=T);allocate(Elold,mold=T);allocate(Egold,mold=T)
  allocate(HL,mold=T);allocate(HG,mold=T)
  allocate(inc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,2))
  theta_l=0;theta_g=0;qphase=0;roundoff_energy=0
  do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
   if(cfg%VF(i,j,k).eq.0)then
    theta_l(i,j,k)=wall_superheat;theta_g(i,j,k)=wall_superheat
   else if(1-vf%VF(i,j,k).gt.VFlo)then
    yref=max(vf%Gbary(2,i,j,k)-block_height,0.0_WP)
    theta_g(i,j,k)=wall_superheat*(1-erf(yref/(2*sqrt(ag*tshift)))/erf(lam))
   endif
  enddo;enddo;enddo
  call sensible_from_temperature(vf)
  hsL=hypre_str(cfg=cfg,name='liquid heat',method=gmres_pfmg,nst=7)
  hsG=hypre_str(cfg=cfg,name='vapor heat',method=gmres_pfmg,nst=7)
  call initialize_solver(hsL);call initialize_solver(hsG)
  call adv%initialize(cfg=cfg,nscalar=2,name='phase sensible energies')
  adv%phase=[0,1];adv%Prho=[rl,rg];adv%diff=0;call adv%setup()
 end subroutine
 subroutine initialize_solver(hs)
  type(hypre_str),intent(inout)::hs
  hs%stc(1,:)=[0,0,0];hs%stc(2,:)=[1,0,0];hs%stc(3,:)=[-1,0,0]
  hs%stc(4,:)=[0,1,0];hs%stc(5,:)=[0,-1,0];hs%stc(6,:)=[0,0,1];hs%stc(7,:)=[0,0,-1]
  hs%rcvg=1e-11_WP;hs%maxit=250;hs%opr(1,:,:,:)=1;call hs%init()
 end subroutine
 subroutine ghost(a)
  real(WP),intent(inout)::a(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  integer::i,j
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
  call cfg%sync(a)
 end subroutine
 real(WP) function phase_fraction(vf,i,j,k,ph)result(a)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k,ph
  a=vf%VF(i,j,k);if(ph.eq.1)a=1-a
 end function
 real(WP) function face_liquid(vf,i,j,k,p1,p2)result(f)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k
  real(WP),intent(in)::p1(3),p2(3)
  real(WP)::plane(4),a,b,side
  if(vf%VF(i,j,k).le.VFlo)then
   f=0;return
  else if(vf%VF(i,j,k).ge.1-VFlo)then
   f=1;return
  endif
  plane=getPlane(vf%liquid_gas_interface(i,j,k),0)
  a=dot_product(plane(1:3),p1)-plane(4);b=dot_product(plane(1:3),p2)-plane(4)
  side=dot_product(plane(1:3),vf%Lbary(:,i,j,k))-plane(4)
  if(side.gt.0)then
   a=-a;b=-b
  endif
  if(a.le.0.and.b.le.0)then
   f=1
  else if(a.ge.0.and.b.ge.0)then
   f=0
  else if(a.lt.0)then
   f=abs(a)/(abs(a)+abs(b))
  else
   f=abs(b)/(abs(a)+abs(b))
  endif
 end function
 real(WP) function coupling(vf,i,j,k,d,ph)result(c)
  type(vfs),intent(in)::vf
  integer,intent(in)::i,j,k,d,ph
  integer,parameter::di(4)=[1,-1,0,0],dj(4)=[0,0,1,-1]
  integer::ii,jj,axis
  real(WP)::p1(3),p2(3),center(3),f1,f2,a,kap,bi(3),bj(3),distance,h
  ii=i+di(d);jj=j+dj(d);c=0
  if(jj.lt.cfg%jmin.or.jj.gt.cfg%jmax)return
  if(.not.cfg%xper.and.(ii.lt.cfg%imin.or.ii.gt.cfg%imax))return
  axis=1;if(d.gt.2)axis=2
  if(axis.eq.1)then
   h=cfg%dx(i);p1=[cfg%x(i+max(di(d),0)),cfg%y(j),cfg%zm(k)];p2=p1;p2(2)=cfg%y(j+1)
  else
   h=cfg%dy(j);p1=[cfg%x(i),cfg%y(j+max(dj(d),0)),cfg%zm(k)];p2=p1;p2(1)=cfg%x(i+1)
  endif
  center=.5_WP*(p1+p2)
  if(phase_fraction(vf,i,j,k,ph).le.VFlo)return
  f1=face_liquid(vf,i,j,k,p1,p2);if(ph.eq.1)f1=1-f1
  if(cfg%VF(ii,jj,k).eq.0)then
   a=f1
  else
   if(phase_fraction(vf,ii,jj,k,ph).le.VFlo)return
   f2=face_liquid(vf,ii,jj,k,p1,p2);if(ph.eq.1)f2=1-f2
   a=min(f1,f2)
  endif
  if(a.le.0)return
  kap=kl;bi=vf%Lbary(:,i,j,k)
  if(ph.eq.1)then
   kap=kg;bi=vf%Gbary(:,i,j,k)
  endif
  if(cfg%VF(ii,jj,k).eq.0)then
   distance=abs(bi(axis)-center(axis))
   c=a/(h*(distance/kap+.5_WP*h/ks))
  else
   bj=vf%Lbary(:,ii,jj,k);if(ph.eq.1)bj=vf%Gbary(:,ii,jj,k)
   distance=abs(bi(axis)-bj(axis))
   if(axis.eq.1.and.cfg%xper.and.distance.gt.cfg%xL/2)distance=abs(distance-cfg%xL)
   if(distance.le.0)call die('Nonpositive same-phase heat donor distance')
   c=kap*a/(h*distance)
  endif
 end function
 subroutine interface_H(vf)
  type(vfs),intent(in)::vf
  integer::i,j,k
  real(WP)::p(3),n(3),dl,dg,nn
  HL=0;HG=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0.or.vf%VF(i,j,k).le.VFlo.or.vf%VF(i,j,k).ge.1-VFlo.or.vf%SD(i,j,k).le.0)cycle
   p=calculateCentroid(vf%interface_polygon(1,i,j,k));n=calculateNormal(vf%interface_polygon(1,i,j,k));nn=sqrt(sum(n*n))
   if(nn.le.0)call die('Invalid PLIC thermal normal')
   n=n/nn;dl=abs(dot_product(p-vf%Lbary(:,i,j,k),n));dg=abs(dot_product(vf%Gbary(:,i,j,k)-p,n))
   if(min(dl,dg).le.0)call die('Nonpositive phase-to-interface thermal distance')
   HL(i,j,k)=kl*vf%SD(i,j,k)/dl;HG(i,j,k)=kg*vf%SD(i,j,k)/dg
  enddo;enddo;enddo
 end subroutine
 subroutine solve_phase(dt,vf,ph,hs,old,a,H)
  real(WP),intent(in)::dt,old(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),H(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  real(WP),intent(inout)::a(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  type(vfs),intent(in)::vf
  integer,intent(in)::ph
  type(hypre_str),intent(inout)::hs
  integer,parameter::di(4)=[1,-1,0,0],dj(4)=[0,0,1,-1],slot(4)=[2,3,4,5]
  integer::i,j,k,d,ii,jj
  real(WP)::Cvol,c,diag
  hs%opr=0;hs%rhs=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)then
    hs%opr(1,i,j,k)=1;hs%rhs(i,j,k)=wall_superheat;cycle
   endif
   if(phase_fraction(vf,i,j,k,ph).le.VFlo)then
    hs%opr(1,i,j,k)=1;cycle
   endif
   Cvol=phase_fraction(vf,i,j,k,ph)*rl*cpl
   if(ph.eq.1)Cvol=phase_fraction(vf,i,j,k,ph)*rg*cpg
   hs%opr(1,i,j,k)=Cvol/dt+H(i,j,k);hs%rhs(i,j,k)=old(i,j,k)/dt
   do d=1,4
    ii=i+di(d);jj=j+dj(d);c=coupling(vf,i,j,k,d,ph)
    hs%opr(1,i,j,k)=hs%opr(1,i,j,k)+c
    if(cfg%VF(ii,jj,k).eq.0)then
     hs%rhs(i,j,k)=hs%rhs(i,j,k)+c*wall_superheat
    else
     hs%opr(slot(d),i,j,k)=-c
    endif
   enddo
   diag=hs%opr(1,i,j,k);hs%opr(:,i,j,k)=hs%opr(:,i,j,k)/diag;hs%rhs(i,j,k)=hs%rhs(i,j,k)/diag
  enddo;enddo;enddo
  call hs%setup();hs%sol=a;call hs%solve()
  if(hs%rerr.gt.1.01_WP*hs%rcvg)call die('Two-temperature linear solve residual; reference-case tolerance')
  a=hs%sol;call ghost(a)
 end subroutine
 subroutine sensible_from_temperature(vf)
  type(vfs),intent(in)::vf
  El=vf%VF*rl*cpl*theta_l;Eg=(1-vf%VF)*rg*cpg*theta_g
  where(cfg%VF.eq.0)El=rhos*cps*theta_l
  where(cfg%VF.eq.0)Eg=0
  E=El+Eg;T=Tsat+vf%VF*theta_l+(1-vf%VF)*theta_g
  where(cfg%VF.eq.0)T=Tsat+theta_l
 end subroutine
 subroutine energy_step(dt,vf,U,V,mdot)
  real(WP),intent(in)::dt,U(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),V(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  type(vfs),intent(in)::vf
  real(WP),intent(out)::mdot(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
  integer,parameter::di(4)=[1,-1,0,0],dj(4)=[0,0,1,-1]
  integer::i,j,k,d,ii,jj,ierr
  real(WP)::local(3),global(3),vol,c
  Elold=El;Egold=Eg;call interface_H(vf)
  call solve_phase(dt,vf,0,hsL,Elold,theta_l,HL);call solve_phase(dt,vf,1,hsG,Egold,theta_g,HG)
  call sensible_from_temperature(vf);qphase=HL*theta_l+HG*theta_g;mdot=0
  where(vf%SD.gt.0.and.cfg%VF.gt.0)mdot=qphase/(hfg*vf%SD)
  call cfg%sync(mdot);call cfg%sync(qphase);local=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   vol=cfg%vol(i,j,k)*physical_span/Lz
   local(1)=local(1)+(El(i,j,k)+Eg(i,j,k)-Elold(i,j,k)-Egold(i,j,k))*vol/dt
   local(2)=local(2)+qphase(i,j,k)*vol
   do d=1,4
    ii=i+di(d);jj=j+dj(d)
    if(cfg%VF(ii,jj,k).ne.0)cycle
    c=coupling(vf,i,j,k,d,0);local(3)=local(3)+c*(wall_superheat-theta_l(i,j,k))*vol
    c=coupling(vf,i,j,k,d,1);local(3)=local(3)+c*(wall_superheat-theta_g(i,j,k))*vol
   enddo
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,3,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
  power=global(3);phase_power=global(2);thermal_out=0;thermal_balance=global(1)+global(2)-global(3)
 end subroutine
 subroutine energy_remap(dt,vf)
  real(WP),intent(in)::dt
  type(vfs),intent(in)::vf
  real(WP),allocatable::u(:,:,:),v(:,:,:),w(:,:,:)
  integer::i,j,k,ierr
  real(WP)::a,b,local,global,bound,local_round,global_round
  allocate(u,mold=T);allocate(v,mold=T);allocate(w,mold=T)
  u=vf%UFl(1,:,:,:)+vf%UFg(1,:,:,:);v=vf%UFl(2,:,:,:)+vf%UFg(2,:,:,:);w=0
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
  local=0;local_round=0
  do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
   if(cfg%VF(i,j,k).eq.0)cycle
   local=local+(inc(i,j,k,1)+inc(i,j,k,2))*cfg%vol(i,j,k)*physical_span/Lz
   El(i,j,k)=El(i,j,k)+dt*inc(i,j,k,1);Eg(i,j,k)=Eg(i,j,k)+dt*inc(i,j,k,2)
   a=vf%VF(i,j,k);b=1-a
   if(a.le.VFlo)then
    bound=2*VFlo*rl*cpl*max(abs(wall_superheat),1.0_WP)
    if(abs(El(i,j,k)).gt.bound)call die('Resolved liquid sensible energy left in absent phase')
    local_round=local_round+abs(El(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    Eg(i,j,k)=Eg(i,j,k)+El(i,j,k);El(i,j,k)=0;theta_l(i,j,k)=0
   else
    theta_l(i,j,k)=El(i,j,k)/(a*rl*cpl)
   endif
   if(b.le.VFlo)then
    bound=2*VFlo*rg*cpg*max(abs(wall_superheat),1.0_WP)
    if(abs(Eg(i,j,k)).gt.bound)call die('Resolved vapor sensible energy left in absent phase')
    local_round=local_round+abs(Eg(i,j,k))*cfg%vol(i,j,k)*physical_span/Lz
    El(i,j,k)=El(i,j,k)+Eg(i,j,k);Eg(i,j,k)=0;theta_g(i,j,k)=0
   else
    theta_g(i,j,k)=Eg(i,j,k)/(b*rg*cpg)
   endif
  enddo;enddo;enddo
  call MPI_ALLREDUCE(local,global,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);thermal_out=-global
  call MPI_ALLREDUCE(local_round,global_round,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr);roundoff_energy=roundoff_energy+global_round
  where(vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_l=El/(vf%VF*rl*cpl)
  where(1-vf%VF.gt.VFlo.and.cfg%VF.gt.0)theta_g=Eg/((1-vf%VF)*rg*cpg)
  call ghost(theta_l);call ghost(theta_g);E=El+Eg;call cfg%sync(E)
  T=Tsat+vf%VF*theta_l+(1-vf%VF)*theta_g;where(cfg%VF.eq.0)T=Tsat+theta_l
  deallocate(u,v,w)
 end subroutine
end module
