program rect_hull_wave_driver
  use rect_wave_body
  use rect_mac_projection3d,only:configure_grid
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(rect_grid)::g
  type(rect_transport)::tr
  type(sweep_stats)::sw
  type(body_stats)::body
  type(solve_stats)::solver
  integer::nx,ny,nz,step,unit,hunit,sunit,snapunit,d,c,sweep,order(3),i,j,k,pidx(3),left(3),lo(3),hi(3),right(3),failures,reg,mx,mz
  real(dp)::h(3),origin(3),cg(3),Z,Htop,Tend,dtmax,dt,t,gravity_dt,speed,CFL,pe,pe0,peold,pe_cell,pcold,phase0,mass0
  real(dp)::maxcell,maxdual,maxmoment,mombudget,rhoalpha,amin,amax,bphase,bmass,bvol,cphase,cmass,csourcep,csourcem
  real(dp)::sphase,smass,smax,sourcepmax,sourcemmax,momsourceMax,geomerr,phasebudget,massbudget,sourcePE,masswork,boundaryPE,boundaryKE,cumBPE,cumBKE,PEidentity
  real(dp)::kold,kadv,kkick,knew,forcework,forceidentity,vol,pi,bxcell,bzcell,amp,meanheight,region_speed(3,3),region_ke(3),globalmax(3)
  real(dp)::fr,vr,own,sm,inside_min,inside_max,maskvol
  real(dp),allocatable::a(:,:,:),rho(:,:,:),chi(:,:,:),q(:,:,:),phi(:,:,:),basis(:,:),height(:,:),div(:,:,:),dv(:,:,:),sourceSum(:,:,:)
  real(dp),allocatable::m(:,:,:,:),u(:,:,:,:),rf(:,:,:,:),adv(:,:,:,:),freeze(:,:,:,:),force(:,:,:,:),kick(:,:,:,:),ukick(:,:,:,:),usweep(:,:,:,:),fmass(:,:,:),momsource(:,:,:,:)
  logical,allocatable::valid(:,:,:)
  logical::periodic(3)
  character(1024)::config,inputfile,prefix
  if(command_argument_count()/=1)error stop 'config required'
  call get_command_argument(1,config)
  open(newunit=unit,file=trim(config),status='old',action='read')
  read(unit,*)nx,ny,nz;read(unit,*)h;read(unit,*)origin;read(unit,*)Z,Htop,mx,mz
  read(unit,*)Tend,dtmax;read(unit,*)cg;read(unit,'(a)')inputfile;read(unit,'(a)')prefix;close(unit)
  if(nx*ny*nz>16384.or.Tend>.5_dp.or.dtmax<=0)error stop 'bounded rectangular wave stage'
  call configure_grid(nx,ny,nz,h,origin,g);periodic=.false.;vol=g%cell_volume;pi=acos(-1._dp)
  allocate(a(nx,ny,nz),rho(nx,ny,nz),chi(nx,ny,nz),q(nx,ny,nz),phi(nx,ny,nz),valid(nx,ny,nz), &
    basis(nx,nz),height(nx,nz),div(nx,ny,nz),dv(nx,ny,nz),sourceSum(nx,ny,nz), &
    m(0:nx+1,0:ny+1,0:nz+1,3),u(0:nx+1,0:ny+1,0:nz+1,3),rf(0:nx+1,0:ny+1,0:nz+1,3), &
    adv(0:nx+1,0:ny+1,0:nz+1,3),freeze(0:nx+1,0:ny+1,0:nz+1,3),force(0:nx+1,0:ny+1,0:nz+1,3), &
    kick(0:nx+1,0:ny+1,0:nz+1,3),ukick(0:nx+1,0:ny+1,0:nz+1,3),usweep(0:nx+1,0:ny+1,0:nz+1,3), &
    fmass(0:nx+1,0:ny+1,0:nz+1),momsource(0:nx+1,0:ny+1,0:nz+1,3))
  open(newunit=unit,file=trim(inputfile),access='stream',form='unformatted',status='old',action='read');read(unit)a,chi;close(unit)
  do k=1,nz;do i=1,nx
    bxcell=(sin(mx*pi*i/nx)-sin(mx*pi*(i-1)/nx))/(mx*pi/nx)
    bzcell=1
    if(mz/=0)bzcell=(sin(mz*pi*k/nz)-sin(mz*pi*(k-1)/nz))/(mz*pi/nz)
    basis(i,k)=bxcell*bzcell
  enddo;enddo
  rho=rg+dr*a;m=0;u=0;q=0;step=0;t=0;dt=0;failures=0
  phase0=sum(a)*vol;mass0=sum(rho)*vol;maskvol=sum(chi)*vol
  cphase=0;cmass=0;csourcep=0;csourcem=0;cumBPE=0;cumBKE=0
  call allocate_transport(g,tr);call geometric_state(g,a,Z,pe,phi,valid,geomerr);pe0=pe
  call velocities(g,periodic,rho,m,u,rf);call cell_PE_value(pe_cell)
  gravity_dt=.2_dp*sqrt(minval(h)/gravity)
  open(newunit=hunit,file=trim(prefix)//'_history.csv',status='replace',action='write')
  write(hunit,'(a)')'step,time,dt,CFL,mode_amplitude,mean_height,phase_volume,mixture_mass,phase_budget,mass_budget,'// &
    'alpha_min,alpha_max,max_divergence,PCG_true_RMS,PCG_original_RMS,PCG_stop_ratio,PCG_iterations,PCG_flag,'// &
    'component_max_u,component_max_v,component_max_w,water_u,water_v,water_w,air_u,air_v,air_w,mixed_u,mixed_v,mixed_w,'// &
    'KE,PE_PLIC,PE_cell,KE_plus_PE_PLIC,cell_mass_closure,dual_mass_closure,rho_alpha_defect,momentum_local_residual,momentum_sweep_budget,'// &
    'fullstep_mass_source_Linf,fullstep_phase_source_Linf,fullstep_momentum_source_Linf,top_volume_step,top_phase_step,top_mass_step,'// &
    'cumulative_top_phase,cumulative_top_mass,cumulative_phase_source,cumulative_mass_source,boundary_KE_step,boundary_PE_step,'// &
    'cumulative_boundary_KE,cumulative_boundary_PE,source_PE,mass_flux_gravity_work,cell_PE_identity,KE_old,KE_advected,KE_kick,'// &
    'redforce_midpoint_predictor_work,force_work_identity,BE_defect,redforce_endpoint_work,q_pressure_work,penalty_work,penalty_dissipation,'// &
    'body_block_energy_identity,body_momentum_residual,body_momentum_scale,effective_beta_audit,effective_diagonal_audit,'// &
    'inside_MAC_vector_RMS,inside_component_max,body_on_fluid_Fx,body_on_fluid_Fy,body_on_fluid_Fz,body_on_fluid_Mz,'// &
    'chi_cell_alpha_min,chi_cell_alpha_max,solid_volume,PLIC_volume_error,water_KE,air_KE,mixed_KE,div_true_residual_identity,five_wall_normal_error'
  open(newunit=sunit,file=trim(prefix)//'_sweeps.csv',status='replace',action='write')
  write(sunit,'(a)')'step,time,sweep,axis,alpha_min,alpha_max,cell_mass_closure,dual_mass_closure,rho_alpha_defect,momentum_local_residual,'// &
    'phase_budget,mass_budget,momentum_budget,phase_source,mass_source,phase_boundary_out,mass_boundary_out,'// &
    'momentum_source_x,momentum_source_y,momentum_source_z,momentum_flux_x,momentum_flux_y,momentum_flux_z,'// &
    'wall_reaction_x,wall_reaction_y,wall_reaction_z,source_density,phase_activity,volume_activity,tiny_calls'
  open(newunit=snapunit,file=trim(prefix)//'_snapshots.f64',access='stream',form='unformatted',status='replace',action='write')
  amin=minval(a);amax=maxval(a);body=body_stats();solver=solve_stats();adv=0
  maxcell=0;maxdual=0;rhoalpha=0;maxmoment=0;mombudget=0;smax=0;sourcepmax=0;momsourceMax=0
  bphase=0;bmass=0;bvol=0;boundaryKE=0;boundaryPE=0;sourcePE=0;masswork=0;PEidentity=0;kold=0;kadv=0;kkick=0;forcework=0;forceidentity=0
  call snapshot();call diagnostics()
  do while(t<Tend-1e-13_dp)
    if(step>=2000)then
      failures=failures+1;exit
    endif
    call velocities(g,periodic,rho,m,u,rf);speed=maxval(abs(u));dt=min(dtmax,gravity_dt,Tend-t)
    do c=1,3
      call bounds(g,c,lo,hi)
      dt=min(dt,.3_dp*h(c)/max(1e-14_dp,maxval(abs(u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)))))
    enddo
    kold=.5_dp*weighted_dot(g,rf,u,u);peold=pe;pcold=pe_cell
    adv=u;freeze=u;tr%cc=merge(1._dp,0._dp,a>.5_dp);order=[1,2,3];if(modulo(step,2)==1)order=[3,2,1]
    tr%order=order;sourceSum=0;momsource=0;sphase=0;smass=0;maxcell=0;maxdual=0;maxmoment=0;mombudget=0;rhoalpha=0
    amin=minval(a);amax=maxval(a);bphase=0;bmass=0;bvol=0;sourcePE=0;masswork=0;boundaryKE=0;boundaryPE=0
    do sweep=1,3
      d=order(sweep);call velocities(g,periodic,rho,m,usweep,rf)
      ! Three actual phase/mass/MAC sweeps, no copied coupled_step/project call.
      call coupled_sweep(g,periodic,a,rho,m,adv,freeze,tr%cc,dt,rg,dr,d,tr%volume(:,:,:,d),tr%phase(:,:,:,d),tr%mass_source(:,:,:,d),sw)
      fmass=rg*tr%volume(:,:,:,d)+dr*tr%phase(:,:,:,d)
      do k=1,nz;do j=1,ny;do i=1,nx
        pidx=[i,j,k];left=pidx;left(d)=pidx(d)-1
        dv(i,j,k)=tr%volume(i,j,k,d)-tr%volume(left(1),left(2),left(3),d)
        tr%source(i,j,k,d)=tr%cc(i,j,k)*dv(i,j,k)
      enddo;enddo;enddo
      sphase=sphase+sum(tr%source(:,:,:,d));smass=smass+sum(tr%mass_source(:,:,:,d));sourceSum=sourceSum+tr%mass_source(:,:,:,d)
      do c=1,3
        call bounds(g,c,lo,hi)
        do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
          pidx=[i,j,k];sm=dual_source(g,periodic,tr%mass_source(:,:,:,d),c,pidx)
          momsource(i,j,k,c)=momsource(i,j,k,c)+sm*freeze(i,j,k,c)
          if(d==2.and.j==ny)then
            fr=dual_flux(g,periodic,fmass,c,d,pidx,.true.);vr=usweep(i,j,k,c)
            right=pidx;right(2)=ny+1
            if(fr<0)vr=velocity_at(g,periodic,usweep,c,right)
            own=weight(g,c,pidx)/dual_scale(g,periodic,c,pidx)
            boundaryKE=boundaryKE+.5_dp*own*fr*vr**2
          endif
        enddo;enddo;enddo
      enddo
      do j=1,ny
        sourcePE=sourcePE+gravity*(origin(2)+(j-.5_dp)*h(2)-Z)*sum(tr%mass_source(:,j,:,d))
      enddo
      if(d==2)then
        bvol=sum(tr%volume(1:nx,ny,1:nz,d));bphase=sum(tr%phase(1:nx,ny,1:nz,d));bmass=sum(fmass(1:nx,ny,1:nz))
        boundaryPE=gravity*(Htop-Z)*bmass
        masswork=gravity*h(2)*(sum(fmass(1:nx,1:ny-1,1:nz))+.5_dp*sum(fmass(1:nx,0,1:nz)+fmass(1:nx,ny,1:nz)))
      endif
      maxcell=max(maxcell,sw%cell_closure);maxdual=max(maxdual,sw%dual_closure);rhoalpha=max(rhoalpha,sw%rho_alpha)
      maxmoment=max(maxmoment,sw%momentum_local);mombudget=max(mombudget,sw%momentum_budget)
      amin=min(amin,sw%amin);amax=max(amax,sw%amax)
      write(sunit,'(*(es25.16e3,:,","))')real(step+1,dp),t,real(sweep,dp),real(d,dp),sw%amin,sw%amax,sw%cell_closure, &
        sw%dual_closure,sw%rho_alpha,sw%momentum_local,sw%phase_budget,sw%mass_budget,sw%momentum_budget,sw%phase_source,sw%mass_source, &
        sw%phase_boundary,sw%mass_boundary,sw%momentum_source,sw%momentum_flux,sw%wall_reaction,sw%source_density,sw%phase_activity,sw%volume_activity,real(sw%tiny_calls,dp)
    enddo
    cphase=cphase+bphase;cmass=cmass+bmass;csourcep=csourcep+sphase;csourcem=csourcem+smass;cumBPE=cumBPE+boundaryPE;cumBKE=cumBKE+boundaryKE
    smax=maxval(abs(sourceSum))/vol;sourcepmax=maxval(abs(sum(tr%source,dim=4)))/vol;momsourceMax=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        momsourceMax=max(momsourceMax,abs(momsource(i,j,k,c))/(vol*dual_scale(g,periodic,c,[i,j,k])))
      enddo;enddo;enddo
    enddo
    call cell_PE_value(pe_cell);PEidentity=pe_cell-pcold+boundaryPE-masswork-sourcePE
    call velocities(g,periodic,rho,m,u,rf);kadv=.5_dp*weighted_dot(g,rf,u,u)
    call geometric_state(g,a,Z,pe,phi,valid,geomerr);call reduced_force(g,a,phi,valid,force)
    kick=m+dt*force;call velocities(g,periodic,rho,kick,ukick,rf);kkick=.5_dp*weighted_dot(g,rf,ukick,ukick)
    forcework=dt*face_sum(g,force*.5_dp*(u+ukick));forceidentity=kkick-kadv-forcework
    ! ONE combined solve: (rho_f+dt*kappa*chi_f)u_new=m_adv+dt*f_red-dt*Gq.
    ! No full rho*g kick and no earlier pressure projection.
    call body_pressure(g,rho,chi,dt,m,force,q,cg,solver,body)
    call velocities(g,periodic,rho,m,u,rf);t=t+dt;step=step+1
    call diagnostics();flush(hunit);flush(sunit)
    if(step==1.or.modulo(step,4)==0.or.t>=Tend-1e-13_dp)call snapshot()
    if(solver%flag/=0.or.amin< -1e-10_dp.or.amax>1+1e-10_dp.or..not.all(ieee_is_finite(u)))then
      failures=failures+1;exit
    endif
    if(maxval(abs(u))>1.or.(pe0>1e-12_dp.and.knew+pe>1.10_dp*pe0).or.maxval(abs(a(:,ny,:)))>1e-10_dp)then
      failures=failures+1;exit
    endif
  enddo
  if(modulo(step,4)/=0)call snapshot()
  close(hunit);close(sunit);close(snapunit)
  open(newunit=unit,file=trim(prefix)//'_checkpoint.f64',access='stream',form='unformatted',status='replace',action='write')
  write(unit)real(step,dp),t,dt,a,rho,q,chi,u,m;close(unit)
  write(*,'(a,i0)')'FAILURES ',failures;write(*,'(a,i0)')'STEPS ',step
  write(*,'(a,es25.16e3)')'TIME_REACHED ',t;write(*,'(a,es25.16e3)')'INITIAL_PLIC_PE ',pe0
  if(failures/=0)stop 3
contains
  subroutine cell_PE_value(value)
    real(dp),intent(out)::value
    integer::jj
    real(dp)::flat,y
    value=0
    do jj=1,ny
      y=origin(2)+(jj-.5_dp)*h(2);flat=rg+dr*merge(1._dp,0._dp,y<Z)
      value=value+gravity*(y-Z)*sum(rho(:,jj,:)-flat)*vol
    enddo
  end subroutine
  subroutine snapshot()
    real(dp)::ux(nx+1,ny,nz),vy(nx,ny+1,nz),wz(nx,ny,nz+1)
    call pack(g,u,ux,vy,wz)
    write(snapunit)real(step,dp),t,dt,a,rho,q,chi,ux,vy,wz;flush(snapunit)
  end subroutine
  subroutine diagnostics()
    real(dp)::fraction
    globalmax=0;region_speed=0;region_ke=0
    do c=1,3
      call bounds(g,c,lo,hi)
      do k=lo(3),hi(3);do j=lo(2),hi(2);do i=lo(1),hi(1)
        fraction=(rf(i,j,k,c)-rg)/dr;reg=3
        if(fraction<1e-8_dp)reg=2
        if(fraction>1-1e-8_dp)reg=1
        globalmax(c)=max(globalmax(c),abs(u(i,j,k,c)));region_speed(reg,c)=max(region_speed(reg,c),abs(u(i,j,k,c)))
        region_ke(reg)=region_ke(reg)+.5_dp*vol*weight(g,c,[i,j,k])*rf(i,j,k,c)*u(i,j,k,c)**2
      enddo;enddo;enddo
    enddo
    height=origin(2)+sum(a,dim=2)*h(2);meanheight=sum(height)/real(nx*nz,dp)
    amp=sum((height-Z)*basis)/sum(basis*basis);knew=sum(region_ke)
    phasebudget=sum(a)*vol-phase0+cphase-csourcep;massbudget=sum(rho)*vol-mass0+cmass-csourcem
    inside_min=0;inside_max=0
    if(any(chi>0))then
      inside_min=minval(a,mask=chi>0);inside_max=maxval(a,mask=chi>0)
    endif
    CFL=0
    do c=1,3
      call bounds(g,c,lo,hi)
      CFL=max(CFL,maxval(abs(adv(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),c)))*dt/h(c))
    enddo
    write(hunit,'(*(es25.16e3,:,","))')real(step,dp),t,dt,CFL,amp,meanheight,sum(a)*vol,sum(rho)*vol,phasebudget,massbudget, &
      amin,amax,body%divergence,solver%absres,solver%original_absres,solver%stop_ratio,real(solver%iterations,dp),real(solver%flag,dp), &
      globalmax,region_speed(1,:),region_speed(2,:),region_speed(3,:),knew,pe,pe_cell,knew+pe,maxcell,maxdual,rhoalpha,maxmoment,mombudget, &
      smax,sourcepmax,momsourceMax,bvol,bphase,bmass,cphase,cmass,csourcep,csourcem,boundaryKE,boundaryPE,cumBKE,cumBPE,sourcePE,masswork, &
      PEidentity,kold,kadv,kkick,forcework,forceidentity,body%BE_defect,body%red_work,body%p_work,body%penalty_work,body%penalty_dissipation, &
      body%energy_residual,body%momentum_residual,body%momentum_scale,body%beta_audit,body%diag_audit,body%inside_RMS,body%inside_max, &
      body%IB_force,body%IB_Mz,inside_min,inside_max,maskvol,geomerr,region_ke,body%div_residual,body%wall_normal
  end subroutine
end program rect_hull_wave_driver
