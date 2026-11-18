
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Lz
   use hypre_str_class,   only: hypre_str
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use evap_class,        only: evap
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use boiling_energy
   use mpi_f08
   use parallel, only: MPI_REAL_WP
   use messager, only: die
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,psL
   type(tpns),        public :: fs,fsL
   type(vfs),         public :: vf
   type(evap),        public :: evp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,mfileL,cflfile,evpfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Ui_L,Vi_L,Wi_L
   
   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius
   real(WP) :: mdotdp
   real(WP), allocatable :: source_local(:,:),source_global(:,:),pressure_keep(:,:,:)
   integer :: history_unit
   real(WP) :: input_integral=0,phase_integral=0,out_integral=0,initial_sensible,initial_vapor,previous_vapor
   real(WP) :: top_volume_out,top_gas_out,source_error
   real(WP) :: rad_drop
   
contains


   !> Function that defines a level set function for a drop problem
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=xyz(2)-center(2)
   end function levelset_drop


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x- boundary for scalar fields
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y- boundary for scalar fields
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc
   
   
   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator


   !> Function that localizes z- boundary for scalar fields
   function zm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin-1) isIn=.true.
   end function zm_locator_sc
   
   
   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,flux_storage,neumann
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,name='VOF')
         vf%cons_correct=.true.
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         ! call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         ! call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize to a droplet
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
         ! Droplet size
         call vf%get_max()
         rad_drop=sqrt(vf%VFint/(Lz*Pi))
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         use tpns_class,      only: dirichlet,neumann,bcond
         type(bcond), pointer :: mybc
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Static contact angle',fs%contact_angle)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Boundary conditions
         if(.not.cfg%xper)call fs%add_bcond(name='xm',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         if(.not.cfg%xper)call fs%add_bcond(name='xp',type=dirichlet,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=neumann,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         ! call fs%add_bcond(name='zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)
         ! call fs%add_bcond(name='zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.false.,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create and initialize an evp object
      create_evp: block
         ! Create the object
         call evp%initialize(cfg=cfg,vf=vf,itp_x=fs%itpr_x,itp_y=fs%itpr_y,itp_z=fs%itpr_z,div_x=fs%divp_x,div_y=fs%divp_y,div_z=fs%divp_z,name='liquid gas pc')
         evp%rho_l=fs%rho_l;evp%rho_g=fs%rho_g
         allocate(pressure_keep,mold=fs%P)
         call energy_init(vf)
         allocate(source_local(cfg%nx,cfg%ny),source_global(cfg%nx,cfg%ny))
         call global_stores(initial_sensible,initial_vapor)
         previous_vapor=initial_vapor
      end block create_evp
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='VaporizingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         call ens_out%add_scalar('mflux',evp%mflux)
         call ens_out%add_scalar('evp_div',evp%div_src)
         call ens_out%add_scalar('mdotdp',evp%mdotdp)
         call ens_out%add_scalar('mfluxL',evp%mfluxL)
         call ens_out%add_scalar('mfluxG',evp%mfluxG)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_vector('vel_itf',evp%U_itf,evp%V_itf,evp%W_itf)
         call ens_out%add_scalar('temperature',T)
         call ens_out%add_scalar('latent_sink',qphase)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      if(cfg%amRoot)then
       open(newunit=history_unit,file='history.csv',status='new')
       write(history_unit,'(a)') 'step,time_s,dt_s,solid_Tmax_K,fluid_Tmax_K,speed_max_m_s,gas_volume_m3,phase_power_W,sensible_heat_J,latent_change_J,input_heat_J,out_heat_J,total_energy_error_J,thermal_balance_W,mass_balance_relative,div_minus_source_s_1,pressure_relative_error,source_integral_error_s_1m3,top_volume_out_m3_s,top_gas_volume_out_m3_s,physical_Tmin_K'
      endif
      call output_fields()
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
    use tpns_class,only:harmonic_visc,static_contact
    real(WP)::q,dt,source_peak,global_peak
    integer::ierr
    do while(.not.time%done())
     call fs%get_cfl(time%dt,time%cfl);call time%adjust_dt();call time%increment();dt=time%dt
     vf%VFold=vf%VF
     call energy_step(dt,vf,fs%U,fs%V,evp%mdotdp)
     call distribute_expansion()
     fs%Uold=fs%U;fs%Vold=fs%V;fs%Wold=fs%W
     call wall_values()
     call fs%get_olddensity(vf=vf)
     call fs%prepare_advection_upwind(dt=dt)
     ! Project the old carrier to the current thermal expansion source for VOF.
     pressure_keep=fs%P;fs%DPjx=0;fs%DPjy=0;fs%DPjz=0
     call fs%get_div(src=evp%div_src);call project_open(dt)
     fs%P=pressure_keep
     call evp%get_vel_pc()
     evp%U_itf=fs%U-evp%vel_pc(:,:,:,1)
     evp%V_itf=fs%V-evp%vel_pc(:,:,:,2)
     evp%W_itf=0
     call open_vof()
     call vf%advance(dt=dt,U=evp%U_itf,V=evp%V_itf,W=evp%W_itf)
     call vf%apply_bcond(time%t,dt)
     call fs%get_viscosity(vf=vf,strat=harmonic_visc)
     ! Native two-phase midpoint momentum uses old/new VOF face densities.
     do while(time%it.le.time%itmax)
      fs%U=.5_WP*(fs%U+fs%Uold);fs%V=.5_WP*(fs%V+fs%Vold);fs%W=0
      call fs%prepare_advection_upwind(dt=dt)
      call fs%get_dmomdt(resU,resV,resW)
      call fs%addsrc_gravity(resU,resV,resW)
      resU=-2*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+dt*resU
      resV=-2*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+dt*resV
      fs%U=2*fs%U-fs%Uold+resU/fs%rho_U
      fs%V=2*fs%V-fs%Vold+resV/fs%rho_V;fs%W=0
      call fs%apply_bcond(time%t,dt);call wall_values()
      call fs%get_div(src=evp%div_src)
      call fs%add_surface_tension_jump(dt=dt,div=fs%div,vf=vf,contact_model=static_contact)
      call project_open(dt)
      time%it=time%it+1
     enddo
     call energy_remap(dt,vf)
     call fs%interp_vel(Ui,Vi,Wi);call fs%get_max();call vf%get_max()
     input_integral=input_integral+power*dt
     phase_integral=phase_integral+phase_power*dt
     out_integral=out_integral+thermal_out*dt+rg*hfg*top_gas_out*dt
     call record_history(dt)
     if(ens_evt%occurs())then
      call vf%update_surfmesh(smesh);call ens_out%write_data(time%t)
     endif
    enddo
    call output_fields()
    if(cfg%amRoot)close(history_unit)
   end subroutine

   subroutine wall_values()
    integer::n,i,j
    where(fs%umask.eq.1)fs%U=0
    where(fs%vmask.eq.1)fs%V=0
    if(.not.cfg%xper.and.cfg%iproc.eq.1)then
     fs%U(cfg%imin,:,:)=0
     do n=1,cfg%no
      i=cfg%imin-n;fs%V(i,:,:)=-fs%V(cfg%imin+n-1,:,:)
     enddo
    endif
    if(.not.cfg%xper.and.cfg%iproc.eq.cfg%npx)then
     fs%U(cfg%imax+1,:,:)=0
     do n=1,cfg%no
      i=cfg%imax+n;fs%V(i,:,:)=-fs%V(cfg%imax-n+1,:,:)
     enddo
    endif
    if(cfg%jproc.eq.1)then
     fs%V(:,cfg%jmin,:)=0
     do n=1,cfg%no
      j=cfg%jmin-n;fs%U(:,j,:)=-fs%U(:,cfg%jmin+n-1,:)
     enddo
    endif
   end subroutine

   subroutine open_vof()
    integer::i,j,k
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmino_,cfg%kmaxo_;do i=cfg%imino_,cfg%imaxo_
      if(fs%V(i,j+1,k).lt.0)vf%VF(i,j+1:,k)=1.0_WP
     enddo;enddo
    endif
    call cfg%sync(vf%VF)
   end subroutine

   subroutine distribute_expansion()
    integer::i,j,k,di,dj,ii,jj,ierr
    real(WP)::weights,weight,m,local_total,global_total,local_src,global_src
    source_local=0;evp%div_src=0;evp%mflux=evp%mdotdp*vf%SD
    local_total=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     m=evp%mflux(i,j,k)
     if(m.eq.0)cycle
     weights=0
     do dj=-3,3;do di=0,0
      ii=i+di;jj=j+dj
      if(ii.lt.cfg%imin.or.ii.gt.cfg%imax.or.jj.lt.cfg%jmin.or.jj.gt.cfg%jmax)cycle
      if(cfg%VF(ii,jj,k).eq.0.or.vf%VF(ii,jj,k).gt.1e-10_WP)cycle
      weight=exp(-0.5_WP*real(di*di+dj*dj,WP));weights=weights+weight
     enddo;enddo
     if(weights.le.0)then
      print*,'NO_NATIVE_GAS_EXPANSION_SUPPORT',i,j,cfg%xm(i),cfg%ym(j),vf%VF(i,j,k),m
      call die('No pure-gas support within3 cells; no solid/old-value fallback')
     endif
     do dj=-3,3;do di=0,0
      ii=i+di;jj=j+dj
      if(ii.lt.cfg%imin.or.ii.gt.cfg%imax.or.jj.lt.cfg%jmin.or.jj.gt.cfg%jmax)cycle
      if(cfg%VF(ii,jj,k).eq.0.or.vf%VF(ii,jj,k).gt.1e-10_WP)cycle
      weight=exp(-0.5_WP*real(di*di+dj*dj,WP))/weights
      source_local(ii,jj)=source_local(ii,jj)+m*weight
     enddo;enddo
     local_total=local_total+m*cfg%vol(i,j,k)
    enddo;enddo;enddo
    call MPI_ALLREDUCE(source_local,source_global,cfg%nx*cfg%ny,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     if(i.lt.cfg%imin.or.i.gt.cfg%imax.or.j.lt.cfg%jmin.or.j.gt.cfg%jmax)cycle
     evp%mfluxG(i,j,k)=source_global(i,j)
     evp%div_src(i,j,k)=source_global(i,j)*(1/rg-1/rl)
    enddo;enddo;enddo
    call cfg%sync(evp%div_src)
    call MPI_ALLREDUCE(local_total,global_total,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    local_src=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     local_src=local_src+evp%div_src(i,j,k)*cfg%vol(i,j,k)
    enddo;enddo;enddo
    call MPI_ALLREDUCE(local_src,global_src,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    source_error=abs(global_src-global_total*(1/rg-1/rl))*physical_span/Lz
   end subroutine

   subroutine project_open(dt)
    real(WP),intent(in)::dt
    integer::i,j,k,ierr
    real(WP)::local(2),global(2),facearea,gasface,local_rhs,global_rhs
    call fs%update_laplacian()
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmin_,cfg%kmax_;do i=cfg%imin_,cfg%imax_
      ps%opr(1,i,j,k)=ps%opr(1,i,j,k)+2*cfg%vol(i,j,k)/(fs%rho_V(i,j+1,k)*cfg%dy(j)**2)
     enddo;enddo
    endif
    call ps%setup();ps%rhs=-cfg%vol*fs%div/dt;ps%sol=0
    local_rhs=maxval(abs(ps%rhs))
    call MPI_ALLREDUCE(local_rhs,global_rhs,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    if(global_rhs.eq.0)then
     ps%rerr=0;ps%it=0
    else
     call ps%solve()
    endif
    if(ps%rerr.gt.1.01_WP*ps%rcvg)call die('Open pressure linear residual; custom input tolerance')
    call fs%get_pgrad(ps%sol,resU,resV,resW)
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmin_,cfg%kmax_;do i=cfg%imin_,cfg%imax_
      resV(i,j+1,k)=-2*ps%sol(i,j,k)/cfg%dy(j)
     enddo;enddo
    endif
    call cfg%sync(resV)
    fs%U=fs%U-dt*resU/fs%rho_U;fs%V=fs%V-dt*resV/fs%rho_V;fs%W=0
    fs%P=fs%P+ps%sol
    ! Do not reapply the Neumann normal-face copy after pressure correction.
    where(fs%umask.eq.1)fs%U=0
    where(fs%vmask.eq.1)fs%V=0
    if(.not.cfg%xper.and.cfg%iproc.eq.1)fs%U(cfg%imin,:,:)=0
    if(.not.cfg%xper.and.cfg%iproc.eq.cfg%npx)fs%U(cfg%imax+1,:,:)=0
    if(cfg%jproc.eq.1)fs%V(:,cfg%jmin,:)=0
    call wall_values()
    call cfg%sync(fs%U);call cfg%sync(fs%V)
    call fs%get_div(src=evp%div_src)
    local=0
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmin_,cfg%kmax_;do i=cfg%imin_,cfg%imax_
      facearea=cfg%dx(i)*physical_span
      local(1)=local(1)+fs%V(i,j+1,k)*facearea
      gasface=0
      if(fs%V(i,j+1,k).gt.0)gasface=1-vf%VF(i,j,k)
      local(2)=local(2)+fs%V(i,j+1,k)*gasface*facearea
     enddo;enddo
    endif
    call MPI_ALLREDUCE(local,global,2,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    top_volume_out=global(1);top_gas_out=global(2)
   end subroutine

   subroutine global_stores(sensible,gasvolume)
    real(WP),intent(out)::sensible,gasvolume
    real(WP)::local(2),global(2),volume
    integer::i,j,k,ierr
    local=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     volume=cfg%vol(i,j,k)*physical_span/Lz
     local(1)=local(1)+E(i,j,k)*volume
     if(cfg%VF(i,j,k).gt.0)local(2)=local(2)+(1-vf%VF(i,j,k))*volume
    enddo;enddo;enddo
    call MPI_ALLREDUCE(local,global,2,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    sensible=global(1);gasvolume=global(2)
   end subroutine

   subroutine record_history(dt)
    real(WP),intent(in)::dt
    real(WP)::sensible,gasvolume,local(3),global(3),masserr,totalerr,speed
    integer::i,j,k,ierr
    call global_stores(sensible,gasvolume);local=-huge(1.0_WP)
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     local(3)=max(local(3),-T(i,j,k))
     if(cfg%VF(i,j,k).eq.0)then
      local(1)=max(local(1),T(i,j,k))
     else
      local(2)=max(local(2),T(i,j,k))
     endif
    enddo;enddo;enddo
    call MPI_ALLREDUCE(local,global,3,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    masserr=(rg*(gasvolume-previous_vapor)/dt+rg*top_gas_out-phase_power/hfg)/max(abs(phase_power/hfg),1e-14_WP)
    totalerr=sensible-initial_sensible+rg*hfg*(gasvolume-initial_vapor)+out_integral-input_integral
    speed=sqrt(fs%Umax**2+fs%Vmax**2)
    if(.not.ieee_is_finite(speed+sum(global)+gasvolume+sensible+phase_power))call die('Nonfinite physical field')
    if(cfg%amRoot)then
     write(history_unit,'(i0,20(",",es25.16e3))')time%n,time%t,dt,global(1:2),speed,gasvolume,phase_power,sensible,rg*hfg*(gasvolume-initial_vapor),input_integral,out_integral,totalerr,thermal_balance,masserr,fs%divmax,ps%rerr,source_error,top_volume_out,top_gas_out,-global(3)
     flush(history_unit)
    endif
    previous_vapor=gasvolume
   end subroutine

   subroutine output_fields()
    integer::u,i,j,k
    character(len=100)::name
    write(name,'("field_step",i8.8,"_rank",i4.4,".csv")')time%n,cfg%rank
    open(newunit=u,file=trim(name),status='new')
    write(u,'(a)')'x_m,y_m,T_K,liquid_fraction,u_m_s,v_m_s,solid,phase_heat_W_m3,liquid_T_K,vapor_T_K'
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     write(u,'(9(es25.16e3,","),es25.16e3)')cfg%xm(i),cfg%ym(j),T(i,j,k),vf%VF(i,j,k),Ui(i,j,k),Vi(i,j,k),1-cfg%VF(i,j,k),qphase(i,j,k),Tsat+theta_l(i,j,k),Tsat+theta_g(i,j,k)
    enddo;enddo;enddo
    close(u)
    if(cfg%amRoot)then
     open(newunit=u,file='PHASE_ENERGY_ROUNDOFF.txt',status='replace')
     write(u,*) roundoff_energy
     close(u)
    endif
   end subroutine

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,Ui_L,Vi_L,Wi_L)

   end subroutine simulation_final
   
   
end module simulation