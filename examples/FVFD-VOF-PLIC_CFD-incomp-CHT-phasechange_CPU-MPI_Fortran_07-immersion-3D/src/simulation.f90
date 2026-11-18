! 3D_DRIVER_EXTENSION: separate unqualified candidate, fresh3D identity required.

module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Lz,domain_depth,block_depth,cold_depth,domain_height,cold_thickness,cold_width,pool_height,cold_film_thickness,cold_solid,top_port,domain_width,port_width
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
   use irl_fortran_interface,only:ByteBuffer_type
   use mpi_f08
   use parallel, only: MPI_REAL_WP
   use messager, only: die
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,psL
   type(tpns), target, public :: fs,fsL
   type(vfs), target, public :: vf
   type(evap), target, public :: evp
   type(timetracker), target, public :: time
   
   !> Ensight postprocessing
   type(surfmesh), target :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt,diagnostic_evt
   integer::last_field_output_step=-1
   integer::trace_unit=-1
   integer::checkpoint_step,stop_after_checkpoint,last_checkpoint_step=-1
   character(len=1024)::restart_prefix
   character(len=64)::checkpoint_identity
   
   !> Simulation monitor file
   type(monitor) :: mfile,mfileL,cflfile,evpfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable, target :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Ui_L,Vi_L,Wi_L
   
   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius
   real(WP) :: mdotdp
   real(WP), allocatable :: pressure_keep(:,:,:)
   integer :: history_unit
   real(WP)::port_velocity_factor=0,port_volume_rate=0
   real(WP)::phase_cfl_limit,phase_cfl_actual
   character(len=32)::phase_budget_mode
   integer::thermal_trial_count,thermal_replay_test,full_step_replay_test
   real(WP) :: input_integral=0,phase_integral=0,out_integral=0,initial_sensible,initial_vapor,previous_vapor
   real(WP) :: top_volume_out,top_gas_out,source_error
   real(WP) :: rad_drop
   
   real(WP),allocatable::accepted_snapshot(:)
   type(ByteBuffer_type),allocatable::accepted_irl
   real(WP)::snapshot_time(7),snapshot_ledger(22)
   integer::snapshot_int(6)
   logical::snapshot_birth_done,snapshot_birth_trial

contains


   !> Function that defines a level set function for a drop problem
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=pool_height-xyz(2)
      G=max(G,min(cold_width/2-abs(xyz(1)),cold_depth/2-abs(xyz(3)),xyz(2)-(domain_height-cold_thickness-cold_film_thickness)))
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
         vf%use_phase_source=.true.
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize to a droplet
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
          block
           real(WP)::lo,hi,dl,fl,moment,aa,bb,dy
           lo=cfg%y(j);hi=cfg%y(j+1);dy=hi-lo
           aa=lo;bb=min(hi,pool_height);dl=max(0.0_WP,bb-aa);moment=dl*.5_WP*(aa+bb)
           if(abs(cfg%xm(i))<cold_width/2.and.abs(cfg%zm(k))<cold_depth/2)then
            aa=max(lo,domain_height-cold_thickness-cold_film_thickness);bb=min(hi,domain_height-cold_thickness)
            fl=max(0.0_WP,bb-aa);dl=dl+fl;moment=moment+fl*.5_WP*(aa+bb)
           endif
           vf%VF(i,j,k)=dl/dy
           if(cfg%VF(i,j,k)==0)vf%VF(i,j,k)=1
           if(vf%VF(i,j,k)<-VFlo.or.vf%VF(i,j,k)>1+VFlo)call die('Invalid analytic3D phase initialization')
           vf%VF(i,j,k)=max(0.0_WP,min(1.0_WP,vf%VF(i,j,k)))
           vf%Lbary(:,i,j,k)=[cfg%xm(i),cfg%ym(j),cfg%zm(k)];vf%Gbary(:,i,j,k)=vf%Lbary(:,i,j,k)
           if(vf%VF(i,j,k)>VFlo.and.vf%VF(i,j,k)<VFhi)then
            vf%Lbary(2,i,j,k)=moment/dl
            vf%Gbary(:,i,j,k)=([cfg%xm(i),cfg%ym(j),cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1-vf%VF(i,j,k))
           endif
          end block
         enddo;enddo;enddo
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
         rad_drop=(3*vf%VFint/(4*Pi))**(1.0_WP/3)
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
         call fs%add_bcond(name='xm',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=dirichlet,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=neumann,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         call fs%add_bcond(name='zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)
         call fs%add_bcond(name='zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         if(cfg%jproc.eq.cfg%npy)then
          block
           integer::i,k
           do k=cfg%kmino_,cfg%kmaxo_;do i=cfg%imino_,cfg%imaxo_
            fs%vmask(i,cfg%jmax+1,k)=1
           enddo;enddo
          end block
         endif
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
         call param_read('Phase conversion CFL',phase_cfl_limit)
         call param_read('Phase budget mode',phase_budget_mode)
         if(trim(phase_budget_mode).ne.'joint_volume_probe'.and.trim(phase_budget_mode).ne.'isolated_inventory'.and.trim(phase_budget_mode).ne.'bounded_inventory_enthalpy')call die('Explicit phase budget mode required')
         inventory_correction_enabled=trim(phase_budget_mode).eq.'bounded_inventory_enthalpy'
         call param_read('Thermal trial replay',thermal_replay_test)
         call param_read('Full step replay test',full_step_replay_test)
         if(full_step_replay_test.ne.0.and.full_step_replay_test.ne.1)call die('Full step replay must be0 or1')
         if(phase_cfl_limit.le.0.or.phase_cfl_limit.ge.1)call die('Invalid phase conversion CFL')
         phase_cfl_actual=0;thermal_trial_count=0
         call global_stores(initial_sensible,initial_vapor)
         previous_vapor=initial_vapor
      end block create_evp
      

      if(vf%two_planes)call die('Checkpoint format qualified only for this one-plane LVIRA case')
      call param_read('Restart prefix',restart_prefix)
      call param_read('Checkpoint step',checkpoint_step)
      call param_read('Stop after checkpoint',stop_after_checkpoint)
      call param_read('Checkpoint identity',checkpoint_identity)
      if(len_trim(checkpoint_identity).ne.64)call die('Checkpoint physics identity must be SHA256')
      ! This new geometry has an explicit vapor headspace. Local wall-site activation guards are retained; global no-vapor restriction belonged to the old fully liquid case.
      if(trim(restart_prefix).ne.'none')call read_checkpoint(trim(restart_prefix))
      if(trim(restart_prefix).eq.'none')then
       block
        real(WP)::initial_zero_dt
        initial_zero_dt=0._WP
        ! Define auxiliary masses/interpolators for fresh-run diagnostics only.
        ! The first actual step recomputes these from the unchanged geometry.
        call fs%get_olddensity(vf=vf)
        call fs%prepare_advection_upwind(dt=initial_zero_dt)
       end block
      endif
      call save_accepted_snapshot();call restore_accepted_snapshot()
      call trace_stage('restored',time%dt,.false.)

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
         diagnostic_evt=event(time=time,name='Double precision flow snapshots')
         call param_read('CSV output period',diagnostic_evt%tper)
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
       write(history_unit,'(a)') 'step,time_s,dt_s,solid_Tmax_K,fluid_Tmax_K,speed_max_m_s,gas_volume_m3,phase_power_W,sensible_heat_J,latent_change_J,input_heat_J,out_heat_J,total_energy_error_J,thermal_balance_W,mass_balance_relative,div_minus_source_s_1,pressure_relative_error,source_integral_error_s_1m3,top_volume_out_m3_s,top_gas_volume_out_m3_s,physical_Tmin_K,thermal_coupling_iterations,thermal_coupling_change_K,inactive_thermal_interface_cells,inactive_thermal_interface_area_m2,phase_conversion_cfl,thermal_trial_count,gas_centroid_x_m,gas_centroid_y_m,heater_gas_contact_area_m2,wall_site_temperature_K,wall_birth_events,wall_birth_surface_J,evaporation_power_W,condensation_power_W,cold_plate_heat_out_W,w_max_m_s,gas_centroid_z_m'
      endif
      call output_fields()
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine save_accepted_snapshot()
    use irl_fortran_interface,only:new,serializeAndPack
    integer::off,n,i,j,k,total
    total=size(wall_birth_m)+size(wall_birth_point)+size(fs%U)+size(fs%V)+size(fs%W)+size(fs%P)+size(fs%rho_U)+size(fs%rho_V)+size(fs%rho_W)+size(fs%DPjx)+size(fs%DPjy)+size(fs%DPjz)+size(fs%div)+size(fs%Pjx)+size(fs%Pjy)+size(fs%Pjz)+size(vf%VF)+size(vf%VFold)+size(vf%Lbary)+size(vf%Gbary)+size(vf%SD)+size(vf%G)+size(vf%curv)+size(vf%Lvol)+size(vf%Gvol)+size(T)+size(E)+size(El)+size(Eg)+size(theta_l)+size(theta_g)+size(theta_s)+size(qphase)+size(evp%div_src)+size(evp%mdotdp)+size(evp%mflux)+size(evp%mfluxL)+size(evp%mfluxG)+size(evp%U_itf)+size(evp%V_itf)+size(evp%W_itf)
    if(.not.allocated(accepted_snapshot))allocate(accepted_snapshot(total))
    if(size(accepted_snapshot).ne.total)call die('Accepted snapshot size changed')
    off=0
    n=size(wall_birth_m);accepted_snapshot(off+1:off+n)=reshape(wall_birth_m,[n]);off=off+n
    n=size(wall_birth_point);accepted_snapshot(off+1:off+n)=reshape(wall_birth_point,[n]);off=off+n
    n=size(fs%U);accepted_snapshot(off+1:off+n)=reshape(fs%U,[n]);off=off+n
    n=size(fs%V);accepted_snapshot(off+1:off+n)=reshape(fs%V,[n]);off=off+n
    n=size(fs%W);accepted_snapshot(off+1:off+n)=reshape(fs%W,[n]);off=off+n
    n=size(fs%P);accepted_snapshot(off+1:off+n)=reshape(fs%P,[n]);off=off+n
    n=size(fs%rho_U);accepted_snapshot(off+1:off+n)=reshape(fs%rho_U,[n]);off=off+n
    n=size(fs%rho_V);accepted_snapshot(off+1:off+n)=reshape(fs%rho_V,[n]);off=off+n
    n=size(fs%rho_W);accepted_snapshot(off+1:off+n)=reshape(fs%rho_W,[n]);off=off+n
    n=size(fs%DPjx);accepted_snapshot(off+1:off+n)=reshape(fs%DPjx,[n]);off=off+n
    n=size(fs%DPjy);accepted_snapshot(off+1:off+n)=reshape(fs%DPjy,[n]);off=off+n
    n=size(fs%DPjz);accepted_snapshot(off+1:off+n)=reshape(fs%DPjz,[n]);off=off+n
    n=size(fs%div);accepted_snapshot(off+1:off+n)=reshape(fs%div,[n]);off=off+n
    n=size(fs%Pjx);accepted_snapshot(off+1:off+n)=reshape(fs%Pjx,[n]);off=off+n
    n=size(fs%Pjy);accepted_snapshot(off+1:off+n)=reshape(fs%Pjy,[n]);off=off+n
    n=size(fs%Pjz);accepted_snapshot(off+1:off+n)=reshape(fs%Pjz,[n]);off=off+n
    n=size(vf%VF);accepted_snapshot(off+1:off+n)=reshape(vf%VF,[n]);off=off+n
    n=size(vf%VFold);accepted_snapshot(off+1:off+n)=reshape(vf%VFold,[n]);off=off+n
    n=size(vf%Lbary);accepted_snapshot(off+1:off+n)=reshape(vf%Lbary,[n]);off=off+n
    n=size(vf%Gbary);accepted_snapshot(off+1:off+n)=reshape(vf%Gbary,[n]);off=off+n
    n=size(vf%SD);accepted_snapshot(off+1:off+n)=reshape(vf%SD,[n]);off=off+n
    n=size(vf%G);accepted_snapshot(off+1:off+n)=reshape(vf%G,[n]);off=off+n
    n=size(vf%curv);accepted_snapshot(off+1:off+n)=reshape(vf%curv,[n]);off=off+n
    n=size(vf%Lvol);accepted_snapshot(off+1:off+n)=reshape(vf%Lvol,[n]);off=off+n
    n=size(vf%Gvol);accepted_snapshot(off+1:off+n)=reshape(vf%Gvol,[n]);off=off+n
    n=size(T);accepted_snapshot(off+1:off+n)=reshape(T,[n]);off=off+n
    n=size(E);accepted_snapshot(off+1:off+n)=reshape(E,[n]);off=off+n
    n=size(El);accepted_snapshot(off+1:off+n)=reshape(El,[n]);off=off+n
    n=size(Eg);accepted_snapshot(off+1:off+n)=reshape(Eg,[n]);off=off+n
    n=size(theta_l);accepted_snapshot(off+1:off+n)=reshape(theta_l,[n]);off=off+n
    n=size(theta_g);accepted_snapshot(off+1:off+n)=reshape(theta_g,[n]);off=off+n
    n=size(theta_s);accepted_snapshot(off+1:off+n)=reshape(theta_s,[n]);off=off+n
    n=size(qphase);accepted_snapshot(off+1:off+n)=reshape(qphase,[n]);off=off+n
    n=size(evp%div_src);accepted_snapshot(off+1:off+n)=reshape(evp%div_src,[n]);off=off+n
    n=size(evp%mdotdp);accepted_snapshot(off+1:off+n)=reshape(evp%mdotdp,[n]);off=off+n
    n=size(evp%mflux);accepted_snapshot(off+1:off+n)=reshape(evp%mflux,[n]);off=off+n
    n=size(evp%mfluxL);accepted_snapshot(off+1:off+n)=reshape(evp%mfluxL,[n]);off=off+n
    n=size(evp%mfluxG);accepted_snapshot(off+1:off+n)=reshape(evp%mfluxG,[n]);off=off+n
    n=size(evp%U_itf);accepted_snapshot(off+1:off+n)=reshape(evp%U_itf,[n]);off=off+n
    n=size(evp%V_itf);accepted_snapshot(off+1:off+n)=reshape(evp%V_itf,[n]);off=off+n
    n=size(evp%W_itf);accepted_snapshot(off+1:off+n)=reshape(evp%W_itf,[n]);off=off+n
    snapshot_time=[time%t,time%dt,time%told,time%dtold,time%tmid,time%dtmid,time%cfl]
    snapshot_ledger=[input_integral,phase_integral,out_integral,initial_sensible,initial_vapor,previous_vapor,top_volume_out,top_gas_out,source_error,phase_cfl_actual,roundoff_energy,phase_power,thermal_out,thermal_balance,coupling_change,inactive_thermal_surface_m2,wall_birth_event_mass,wall_birth_latent_J,wall_birth_surface_J,cold_heat_out,vf%liquid_roundoff_volume,vf%gas_roundoff_volume]
    snapshot_int=[time%n,time%it,thermal_trial_count,coupling_iterations,inactive_thermal_interfaces,wall_birth_events]
    snapshot_birth_done=wall_birth_done;snapshot_birth_trial=wall_birth_trial
    if(allocated(accepted_irl))deallocate(accepted_irl)
    allocate(accepted_irl);call new(accepted_irl)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call serializeAndPack(vf%liquid_gas_interface(i,j,k),accepted_irl)
    enddo;enddo;enddo
   end subroutine
   subroutine restore_accepted_snapshot()
    use irl_fortran_interface,only:unpackAndStore,resetBufferPointer,new,serializeAndPack,getSize,dataPtr
    integer::off,n,i,j,k,ierr
    real(WP)::local_delta,global_delta
    type(ByteBuffer_type)::verify_buffer
    integer(kind=1),pointer::saved_bytes(:),verify_bytes(:)
    off=0
    n=size(wall_birth_m);wall_birth_m=reshape(accepted_snapshot(off+1:off+n),shape(wall_birth_m));off=off+n
    n=size(wall_birth_point);wall_birth_point=reshape(accepted_snapshot(off+1:off+n),shape(wall_birth_point));off=off+n
    n=size(fs%U);fs%U=reshape(accepted_snapshot(off+1:off+n),shape(fs%U));off=off+n
    n=size(fs%V);fs%V=reshape(accepted_snapshot(off+1:off+n),shape(fs%V));off=off+n
    n=size(fs%W);fs%W=reshape(accepted_snapshot(off+1:off+n),shape(fs%W));off=off+n
    n=size(fs%P);fs%P=reshape(accepted_snapshot(off+1:off+n),shape(fs%P));off=off+n
    n=size(fs%rho_U);fs%rho_U=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_U));off=off+n
    n=size(fs%rho_V);fs%rho_V=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_V));off=off+n
    n=size(fs%rho_W);fs%rho_W=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_W));off=off+n
    n=size(fs%DPjx);fs%DPjx=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjx));off=off+n
    n=size(fs%DPjy);fs%DPjy=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjy));off=off+n
    n=size(fs%DPjz);fs%DPjz=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjz));off=off+n
    n=size(fs%div);fs%div=reshape(accepted_snapshot(off+1:off+n),shape(fs%div));off=off+n
    n=size(fs%Pjx);fs%Pjx=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjx));off=off+n
    n=size(fs%Pjy);fs%Pjy=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjy));off=off+n
    n=size(fs%Pjz);fs%Pjz=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjz));off=off+n
    n=size(vf%VF);vf%VF=reshape(accepted_snapshot(off+1:off+n),shape(vf%VF));off=off+n
    n=size(vf%VFold);vf%VFold=reshape(accepted_snapshot(off+1:off+n),shape(vf%VFold));off=off+n
    n=size(vf%Lbary);vf%Lbary=reshape(accepted_snapshot(off+1:off+n),shape(vf%Lbary));off=off+n
    n=size(vf%Gbary);vf%Gbary=reshape(accepted_snapshot(off+1:off+n),shape(vf%Gbary));off=off+n
    n=size(vf%SD);vf%SD=reshape(accepted_snapshot(off+1:off+n),shape(vf%SD));off=off+n
    n=size(vf%G);vf%G=reshape(accepted_snapshot(off+1:off+n),shape(vf%G));off=off+n
    n=size(vf%curv);vf%curv=reshape(accepted_snapshot(off+1:off+n),shape(vf%curv));off=off+n
    n=size(vf%Lvol);vf%Lvol=reshape(accepted_snapshot(off+1:off+n),shape(vf%Lvol));off=off+n
    n=size(vf%Gvol);vf%Gvol=reshape(accepted_snapshot(off+1:off+n),shape(vf%Gvol));off=off+n
    n=size(T);T=reshape(accepted_snapshot(off+1:off+n),shape(T));off=off+n
    n=size(E);E=reshape(accepted_snapshot(off+1:off+n),shape(E));off=off+n
    n=size(El);El=reshape(accepted_snapshot(off+1:off+n),shape(El));off=off+n
    n=size(Eg);Eg=reshape(accepted_snapshot(off+1:off+n),shape(Eg));off=off+n
    n=size(theta_l);theta_l=reshape(accepted_snapshot(off+1:off+n),shape(theta_l));off=off+n
    n=size(theta_g);theta_g=reshape(accepted_snapshot(off+1:off+n),shape(theta_g));off=off+n
    n=size(theta_s);theta_s=reshape(accepted_snapshot(off+1:off+n),shape(theta_s));off=off+n
    n=size(qphase);qphase=reshape(accepted_snapshot(off+1:off+n),shape(qphase));off=off+n
    n=size(evp%div_src);evp%div_src=reshape(accepted_snapshot(off+1:off+n),shape(evp%div_src));off=off+n
    n=size(evp%mdotdp);evp%mdotdp=reshape(accepted_snapshot(off+1:off+n),shape(evp%mdotdp));off=off+n
    n=size(evp%mflux);evp%mflux=reshape(accepted_snapshot(off+1:off+n),shape(evp%mflux));off=off+n
    n=size(evp%mfluxL);evp%mfluxL=reshape(accepted_snapshot(off+1:off+n),shape(evp%mfluxL));off=off+n
    n=size(evp%mfluxG);evp%mfluxG=reshape(accepted_snapshot(off+1:off+n),shape(evp%mfluxG));off=off+n
    n=size(evp%U_itf);evp%U_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%U_itf));off=off+n
    n=size(evp%V_itf);evp%V_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%V_itf));off=off+n
    n=size(evp%W_itf);evp%W_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%W_itf));off=off+n
    time%t=snapshot_time(1)
    time%dt=snapshot_time(2)
    time%told=snapshot_time(3)
    time%dtold=snapshot_time(4)
    time%tmid=snapshot_time(5)
    time%dtmid=snapshot_time(6)
    time%cfl=snapshot_time(7)
    input_integral=snapshot_ledger(1)
    phase_integral=snapshot_ledger(2)
    out_integral=snapshot_ledger(3)
    initial_sensible=snapshot_ledger(4)
    initial_vapor=snapshot_ledger(5)
    previous_vapor=snapshot_ledger(6)
    top_volume_out=snapshot_ledger(7)
    top_gas_out=snapshot_ledger(8)
    source_error=snapshot_ledger(9)
    phase_cfl_actual=snapshot_ledger(10)
    roundoff_energy=snapshot_ledger(11)
    phase_power=snapshot_ledger(12)
    thermal_out=snapshot_ledger(13)
    thermal_balance=snapshot_ledger(14)
    coupling_change=snapshot_ledger(15)
    inactive_thermal_surface_m2=snapshot_ledger(16)
    wall_birth_event_mass=snapshot_ledger(17)
    wall_birth_latent_J=snapshot_ledger(18)
    wall_birth_surface_J=snapshot_ledger(19)
    cold_heat_out=snapshot_ledger(20)
    vf%liquid_roundoff_volume=snapshot_ledger(21)
    vf%gas_roundoff_volume=snapshot_ledger(22)
    time%n=snapshot_int(1)
    time%it=snapshot_int(2)
    thermal_trial_count=snapshot_int(3)
    coupling_iterations=snapshot_int(4)
    inactive_thermal_interfaces=snapshot_int(5)
    wall_birth_events=snapshot_int(6)
    wall_birth_done=snapshot_birth_done;wall_birth_trial=snapshot_birth_trial
    call resetBufferPointer(accepted_irl)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call unpackAndStore(vf%liquid_gas_interface(i,j,k),accepted_irl)
    enddo;enddo;enddo
    call vf%update_band();call vf%polygonalize_interface()
    call resetBufferPointer(accepted_irl)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call unpackAndStore(vf%liquid_gas_interface(i,j,k),accepted_irl)
    enddo;enddo;enddo
    ! Polygon object recreation must not change the exact saved separator/SD or any accepted physical field.
    off=0
    n=size(wall_birth_m);wall_birth_m=reshape(accepted_snapshot(off+1:off+n),shape(wall_birth_m));off=off+n
    n=size(wall_birth_point);wall_birth_point=reshape(accepted_snapshot(off+1:off+n),shape(wall_birth_point));off=off+n
    n=size(fs%U);fs%U=reshape(accepted_snapshot(off+1:off+n),shape(fs%U));off=off+n
    n=size(fs%V);fs%V=reshape(accepted_snapshot(off+1:off+n),shape(fs%V));off=off+n
    n=size(fs%W);fs%W=reshape(accepted_snapshot(off+1:off+n),shape(fs%W));off=off+n
    n=size(fs%P);fs%P=reshape(accepted_snapshot(off+1:off+n),shape(fs%P));off=off+n
    n=size(fs%rho_U);fs%rho_U=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_U));off=off+n
    n=size(fs%rho_V);fs%rho_V=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_V));off=off+n
    n=size(fs%rho_W);fs%rho_W=reshape(accepted_snapshot(off+1:off+n),shape(fs%rho_W));off=off+n
    n=size(fs%DPjx);fs%DPjx=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjx));off=off+n
    n=size(fs%DPjy);fs%DPjy=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjy));off=off+n
    n=size(fs%DPjz);fs%DPjz=reshape(accepted_snapshot(off+1:off+n),shape(fs%DPjz));off=off+n
    n=size(fs%div);fs%div=reshape(accepted_snapshot(off+1:off+n),shape(fs%div));off=off+n
    n=size(fs%Pjx);fs%Pjx=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjx));off=off+n
    n=size(fs%Pjy);fs%Pjy=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjy));off=off+n
    n=size(fs%Pjz);fs%Pjz=reshape(accepted_snapshot(off+1:off+n),shape(fs%Pjz));off=off+n
    n=size(vf%VF);vf%VF=reshape(accepted_snapshot(off+1:off+n),shape(vf%VF));off=off+n
    n=size(vf%VFold);vf%VFold=reshape(accepted_snapshot(off+1:off+n),shape(vf%VFold));off=off+n
    n=size(vf%Lbary);vf%Lbary=reshape(accepted_snapshot(off+1:off+n),shape(vf%Lbary));off=off+n
    n=size(vf%Gbary);vf%Gbary=reshape(accepted_snapshot(off+1:off+n),shape(vf%Gbary));off=off+n
    n=size(vf%SD);vf%SD=reshape(accepted_snapshot(off+1:off+n),shape(vf%SD));off=off+n
    n=size(vf%G);vf%G=reshape(accepted_snapshot(off+1:off+n),shape(vf%G));off=off+n
    n=size(vf%curv);vf%curv=reshape(accepted_snapshot(off+1:off+n),shape(vf%curv));off=off+n
    n=size(vf%Lvol);vf%Lvol=reshape(accepted_snapshot(off+1:off+n),shape(vf%Lvol));off=off+n
    n=size(vf%Gvol);vf%Gvol=reshape(accepted_snapshot(off+1:off+n),shape(vf%Gvol));off=off+n
    n=size(T);T=reshape(accepted_snapshot(off+1:off+n),shape(T));off=off+n
    n=size(E);E=reshape(accepted_snapshot(off+1:off+n),shape(E));off=off+n
    n=size(El);El=reshape(accepted_snapshot(off+1:off+n),shape(El));off=off+n
    n=size(Eg);Eg=reshape(accepted_snapshot(off+1:off+n),shape(Eg));off=off+n
    n=size(theta_l);theta_l=reshape(accepted_snapshot(off+1:off+n),shape(theta_l));off=off+n
    n=size(theta_g);theta_g=reshape(accepted_snapshot(off+1:off+n),shape(theta_g));off=off+n
    n=size(theta_s);theta_s=reshape(accepted_snapshot(off+1:off+n),shape(theta_s));off=off+n
    n=size(qphase);qphase=reshape(accepted_snapshot(off+1:off+n),shape(qphase));off=off+n
    n=size(evp%div_src);evp%div_src=reshape(accepted_snapshot(off+1:off+n),shape(evp%div_src));off=off+n
    n=size(evp%mdotdp);evp%mdotdp=reshape(accepted_snapshot(off+1:off+n),shape(evp%mdotdp));off=off+n
    n=size(evp%mflux);evp%mflux=reshape(accepted_snapshot(off+1:off+n),shape(evp%mflux));off=off+n
    n=size(evp%mfluxL);evp%mfluxL=reshape(accepted_snapshot(off+1:off+n),shape(evp%mfluxL));off=off+n
    n=size(evp%mfluxG);evp%mfluxG=reshape(accepted_snapshot(off+1:off+n),shape(evp%mfluxG));off=off+n
    n=size(evp%U_itf);evp%U_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%U_itf));off=off+n
    n=size(evp%V_itf);evp%V_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%V_itf));off=off+n
    n=size(evp%W_itf);evp%W_itf=reshape(accepted_snapshot(off+1:off+n),shape(evp%W_itf));off=off+n
    call new(verify_buffer)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call serializeAndPack(vf%liquid_gas_interface(i,j,k),verify_buffer)
    enddo;enddo;enddo
    if(getSize(accepted_irl).ne.getSize(verify_buffer))call die('Accepted IRL snapshot size mismatch')
    saved_bytes=>dataPtr(accepted_irl);verify_bytes=>dataPtr(verify_buffer)
    if(any(saved_bytes.ne.verify_bytes))call die('Accepted IRL snapshot byte mismatch')
    nullify(saved_bytes);nullify(verify_bytes)
    local_delta=0;off=0
    n=size(wall_birth_m);local_delta=max(local_delta,maxval(abs(reshape(wall_birth_m,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(wall_birth_point);local_delta=max(local_delta,maxval(abs(reshape(wall_birth_point,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%U);local_delta=max(local_delta,maxval(abs(reshape(fs%U,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%V);local_delta=max(local_delta,maxval(abs(reshape(fs%V,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%W);local_delta=max(local_delta,maxval(abs(reshape(fs%W,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%P);local_delta=max(local_delta,maxval(abs(reshape(fs%P,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%rho_U);local_delta=max(local_delta,maxval(abs(reshape(fs%rho_U,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%rho_V);local_delta=max(local_delta,maxval(abs(reshape(fs%rho_V,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%rho_W);local_delta=max(local_delta,maxval(abs(reshape(fs%rho_W,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%DPjx);local_delta=max(local_delta,maxval(abs(reshape(fs%DPjx,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%DPjy);local_delta=max(local_delta,maxval(abs(reshape(fs%DPjy,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%DPjz);local_delta=max(local_delta,maxval(abs(reshape(fs%DPjz,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%div);local_delta=max(local_delta,maxval(abs(reshape(fs%div,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%Pjx);local_delta=max(local_delta,maxval(abs(reshape(fs%Pjx,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%Pjy);local_delta=max(local_delta,maxval(abs(reshape(fs%Pjy,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(fs%Pjz);local_delta=max(local_delta,maxval(abs(reshape(fs%Pjz,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%VF);local_delta=max(local_delta,maxval(abs(reshape(vf%VF,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%VFold);local_delta=max(local_delta,maxval(abs(reshape(vf%VFold,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%Lbary);local_delta=max(local_delta,maxval(abs(reshape(vf%Lbary,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%Gbary);local_delta=max(local_delta,maxval(abs(reshape(vf%Gbary,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%SD);local_delta=max(local_delta,maxval(abs(reshape(vf%SD,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%G);local_delta=max(local_delta,maxval(abs(reshape(vf%G,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%curv);local_delta=max(local_delta,maxval(abs(reshape(vf%curv,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%Lvol);local_delta=max(local_delta,maxval(abs(reshape(vf%Lvol,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(vf%Gvol);local_delta=max(local_delta,maxval(abs(reshape(vf%Gvol,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(T);local_delta=max(local_delta,maxval(abs(reshape(T,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(E);local_delta=max(local_delta,maxval(abs(reshape(E,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(El);local_delta=max(local_delta,maxval(abs(reshape(El,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(Eg);local_delta=max(local_delta,maxval(abs(reshape(Eg,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(theta_l);local_delta=max(local_delta,maxval(abs(reshape(theta_l,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(theta_g);local_delta=max(local_delta,maxval(abs(reshape(theta_g,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(theta_s);local_delta=max(local_delta,maxval(abs(reshape(theta_s,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(qphase);local_delta=max(local_delta,maxval(abs(reshape(qphase,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%div_src);local_delta=max(local_delta,maxval(abs(reshape(evp%div_src,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%mdotdp);local_delta=max(local_delta,maxval(abs(reshape(evp%mdotdp,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%mflux);local_delta=max(local_delta,maxval(abs(reshape(evp%mflux,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%mfluxL);local_delta=max(local_delta,maxval(abs(reshape(evp%mfluxL,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%mfluxG);local_delta=max(local_delta,maxval(abs(reshape(evp%mfluxG,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%U_itf);local_delta=max(local_delta,maxval(abs(reshape(evp%U_itf,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%V_itf);local_delta=max(local_delta,maxval(abs(reshape(evp%V_itf,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    n=size(evp%W_itf);local_delta=max(local_delta,maxval(abs(reshape(evp%W_itf,[n])-accepted_snapshot(off+1:off+n))));off=off+n
    call MPI_ALLREDUCE(local_delta,global_delta,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    if(global_delta.ne.0)call die('Accepted snapshot physical fields did not restore exactly')
    call fs%interp_vel(Ui,Vi,Wi);call fs%get_max();call vf%get_max()
    if(cfg%amRoot)print*,'ACCEPTED_SNAPSHOT_RESTORED_EXACT',time%n,global_delta
   end subroutine

   subroutine simulation_run
    use tpns_class,only:harmonic_visc,static_contact
    real(WP)::q,dt,source_peak,global_peak,candidate_dt
    integer::ierr,trial,step_attempt,i,j,k
    logical::phase_admitted,forced_replay,inventory_changed
    real(WP)::retry_dt
    main_steps: do while(.not.time%done())
     call fs%get_cfl(time%dt,time%cfl);call time%adjust_dt();candidate_dt=time%dt
     inventory_scale=1
     call save_accepted_snapshot()
     retry_dt=min(candidate_dt,time%tmax-time%t)
     do step_attempt=1,32
      if(step_attempt.gt.1)call restore_accepted_snapshot()
      dt=retry_dt
     vf%VFold=vf%VF
     call trace_stage('step_begin',dt,.false.)
     call save_thermal_trial()
     do trial=1,20
      if(trial.gt.1)call restore_thermal_trial()
      call energy_step(dt,vf,fs%U,fs%V,evp%mdotdp)
      call wall_birth_thermal_trial(dt,vf)
      call inventory_limit_prepare(dt,vf)
     call phase_conversion_cfl(dt,vf,phase_cfl_actual)
      if(trial.eq.1.and.thermal_replay_test.eq.1)cycle
      ! Diagnostic candidate: actual geometric transport+phase-source positivity still gates commit.
      ! Source-only inventory bound is not necessary when advective phase inflow replenishes the cell.
      if(trim(phase_budget_mode).eq.'joint_volume_probe'.or.inventory_correction_enabled)then
       if(cfg%amRoot)print*,'ISOLATED_PHASE_CFL_DIAGNOSTIC',time%n+1,dt,phase_cfl_actual
       exit
      endif
      if(phase_cfl_actual.le.phase_cfl_limit*(1+32*epsilon(1.0_WP)))exit
      if(cfg%amRoot)print*,'THERMAL_PHASE_DT_RETRY',time%n+1,trial,dt,phase_cfl_actual
      dt=dt*min(0.8_WP,0.8_WP*phase_cfl_limit/phase_cfl_actual)
      if(.not.ieee_is_finite(dt).or.time%t+dt.eq.time%t)call die('Phase thermal trial cannot resolve a positive time increment')
     enddo
     if(trial.gt.20)call die('Phase thermal timestep could not be admitted within20 trials')
     call wall_birth_commit()
     thermal_trial_count=trial;time%cfl=time%cfl*dt/candidate_dt;time%dt=dt;call time%increment()
     call distribute_expansion()
     call trace_stage('after_heat_source',dt,.false.)
     fs%Uold=fs%U;fs%Vold=fs%V;fs%Wold=fs%W
     call wall_values()
     call fs%get_olddensity(vf=vf)
     call fs%prepare_advection_upwind(dt=dt)
     call trace_stage('carrier_density',dt,.false.)
     ! Project the old carrier to the current thermal expansion source for VOF.
     pressure_keep=fs%P;fs%DPjx=0;fs%DPjy=0;fs%DPjz=0
     call fs%get_div(src=evp%div_src);call project_open(dt)
     call trace_stage('carrier_projection',dt,.false.)
     fs%P=pressure_keep
     call evp%get_vel_pc()
     evp%U_itf=fs%U
     evp%V_itf=fs%V
     evp%W_itf=fs%W
     call prepare_local_birth_interface(dt)
     call recover_phase_temperature(vf)
     vf%VFold=vf%VF
     call open_vof()
     call vf%advance(dt=dt,U=evp%U_itf,V=evp%V_itf,W=evp%W_itf,admitted=phase_admitted)
     forced_replay=full_step_replay_test.eq.1.and.time%n.eq.1.and.step_attempt.eq.1.and.phase_admitted
     if(forced_replay)phase_admitted=.false.
     if(.not.phase_admitted)then
      if(inventory_correction_enabled.and.vf%phase_source_can_reduce.and..not.forced_replay)then
       if(cfg%amRoot)then
        print*,'INVENTORY_REJECTED_TRIAL_BEFORE_UPDATE',time%n,step_attempt
        flush(6)
       endif
       call inventory_limit_update(vf,dt,.true.,inventory_changed)
       if(inventory_changed)then
        retry_dt=dt
        if(cfg%amRoot)print*,'INVENTORY_ENTHALPY_SOURCE_RETRY',time%n,dt,step_attempt
        cycle
       endif
      endif
      if(inventory_correction_enabled)inventory_scale=1
      retry_dt=.5_WP*dt
      if(forced_replay)then
       retry_dt=dt
       if(cfg%amRoot)print*,'FULL_STEP_REPLAY_REJECTED',time%n,dt
      endif
      if(.not.ieee_is_finite(retry_dt).or.time%told+retry_dt.eq.time%told)call die('Joint phase budget cannot resolve a positive time increment')
      if(cfg%amRoot)print*,'JOINT_PHASE_STEP_RETRY',time%n,time%t,dt,retry_dt,step_attempt
      cycle
     endif
     if(inventory_correction_enabled)then
      call inventory_limit_update(vf,dt,.false.,inventory_changed)
      if(inventory_changed)then
       retry_dt=dt
       if(cfg%amRoot)print*,'INVENTORY_ACTIVE_SET_REPROJECTION',time%n,dt,step_attempt
       cycle
      endif
     endif
     call vf%apply_bcond(time%t,dt)
     call trace_stage('after_vof',dt,.false.)
     call fs%get_viscosity(vf=vf,strat=harmonic_visc)
     ! Native two-phase midpoint momentum uses old/new VOF face densities.
     do while(time%it.le.time%itmax)
      fs%U=.5_WP*(fs%U+fs%Uold);fs%V=.5_WP*(fs%V+fs%Vold);fs%W=.5_WP*(fs%W+fs%Wold)
      call trace_stage('momentum_midpoint',dt,.false.)
      call fs%prepare_advection_upwind(dt=dt)
      call trace_stage('momentum_density',dt,.false.)
      call fs%get_dmomdt(resU,resV,resW)
      call trace_stage('momentum_rhs',dt,.true.)
      call fs%addsrc_gravity(resU,resV,resW)
      call trace_stage('rhs_with_gravity',dt,.true.)
      ! Normal port velocity is prescribed from the global phase-volume balance.
      resU=-2*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+dt*resU
      resV=-2*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+dt*resV
      resW=-2*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+dt*resW
      do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
       fs%U(i,j,k)=2*fs%U(i,j,k)-fs%Uold(i,j,k)+resU(i,j,k)/fs%rho_U(i,j,k)
       fs%V(i,j,k)=2*fs%V(i,j,k)-fs%Vold(i,j,k)+resV(i,j,k)/fs%rho_V(i,j,k)
       fs%W(i,j,k)=2*fs%W(i,j,k)-fs%Wold(i,j,k)+resW(i,j,k)/fs%rho_W(i,j,k)
      enddo;enddo;enddo
      call wall_values();call cfg%sync(fs%U);call cfg%sync(fs%V);call cfg%sync(fs%W)
      call trace_stage('momentum_predictor',dt,.false.)
      call fs%apply_bcond(time%t,dt);call wall_values()
      call wall_values()
      call fs%get_div(src=evp%div_src)
      call fs%add_surface_tension_jump(dt=dt,div=fs%div,vf=vf,contact_model=static_contact)
      call project_open(dt)
      call trace_stage('momentum_projection',dt,.false.)
      time%it=time%it+1
     enddo
     call energy_remap(dt,vf,evp%U_itf,evp%V_itf,evp%W_itf)
     call inventory_return_sensible(vf,time%n,time%t)
     call wall_birth_surface_cost(dt,vf,time%n,time%t)
     call fs%interp_vel(Ui,Vi,Wi);call fs%get_max();call vf%get_max()
     input_integral=input_integral+power*dt
     phase_integral=phase_integral+phase_power*dt
     out_integral=out_integral+thermal_out*dt+rg*hfg*top_gas_out*dt
     call record_history(dt)
     call trace_stage('accepted',dt,.false.)
     if(max(fs%Umax,fs%Vmax,fs%Wmax).gt.20.0_WP)then
      call trace_dump('accepted_speed_over20',dt)
      call die('DIAGNOSTIC_ACCEPTED_SPEED_OVER20: preserved; not a physical result')
     endif
     if(diagnostic_evt%occurs())call output_fields()
     if(ens_evt%occurs())then
      call vf%update_surfmesh(smesh);call ens_out%write_data(time%t)
     endif
     if(time%n.eq.checkpoint_step.or.time%done().or.mod(time%n,1000).eq.0)then
      call write_checkpoint()
      if(stop_after_checkpoint.eq.1)exit main_steps
     endif
      exit
     enddo
     if(step_attempt.gt.32)call die('Joint phase budget not admitted within32 accepted-state retries')
    enddo main_steps
    call output_fields()
    if(cfg%amRoot)close(history_unit)
   end subroutine

   subroutine prepare_local_birth_interface(dt)
    real(WP),intent(in)::dt
    real(WP)::g_old,g_new,added,cellcenter(3)
    integer::i,j,k
    if(.not.wall_birth_trial)return
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(wall_birth_m(i,j,k).le.0)cycle
     g_old=1-vf%VF(i,j,k);added=dt*wall_birth_m(i,j,k)/rg;g_new=g_old+added
     if(g_new.le.0.or.g_new.ge.1)call die('Prescribed local embryo exceeds available cell volume')
     cellcenter=[cfg%xm(i),cfg%ym(j),cfg%zm(k)]
     vf%Gbary(:,i,j,k)=(g_old*vf%Gbary(:,i,j,k)+added*wall_birth_point(:,i,j,k))/g_new
     vf%Lbary(:,i,j,k)=(cellcenter-g_new*vf%Gbary(:,i,j,k))/(1-g_new)
     vf%VF(i,j,k)=1-g_new
    enddo;enddo;enddo
    ! The native transport band is interface-based. Material birth in an initially pure cell must construct an interface before flux assembly.
    call cfg%sync(vf%VF);call vf%sync_and_clean_barycenters();call vf%update_band();call vf%build_interface();call vf%polygonalize_interface()
    if(vf%two_planes)call vf%sense_interface()
    call vf%distance_from_polygon();call vf%subcell_vol();call vf%get_curvature();call vf%reset_volume_moments()
   end subroutine

   real(WP) function port_shape(i,k)result(w)
    integer,intent(in)::i,k
    real(WP)::xi,zi
    xi=(abs(cfg%xm(i))-(domain_width/2-port_width))/port_width
    zi=(abs(cfg%zm(k))-(domain_depth/2-port_width))/port_width
    w=0
    if(xi>0.and.xi<1)w=w+6*xi*(1-xi)
    if(zi>0.and.zi<1)w=w+6*zi*(1-zi)
   end function
   subroutine update_volume_port_rate()
    real(WP)::local(2),global(2)
    integer::i,j,k,ierr
    local=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(cfg%VF(i,j,k)>0)local(1)=local(1)+evp%div_src(i,j,k)*cfg%vol(i,j,k)
    enddo;enddo;enddo
    if(cfg%jproc==cfg%npy)then
     do k=cfg%kmin_,cfg%kmax_;do i=cfg%imin_,cfg%imax_
      if(top_port(i,k))then
       if(vf%VF(i,cfg%jmax,k)>1e-12_WP.or.vf%VF(i,cfg%jmax-1,k)>1e-12_WP)call die('Liquid reached3D vapor reservoir port')
       local(2)=local(2)+port_shape(i,k)*cfg%dx(i)*cfg%dz(k)
      endif
     enddo;enddo
    endif
    call MPI_ALLREDUCE(local,global,2,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    if(global(2)<=0)call die('No3D vapor port area')
    port_volume_rate=global(1);port_velocity_factor=global(1)/global(2)
   end subroutine
   subroutine wall_values()
    integer::n,i,j,k
    real(WP)::bc
    where(fs%umask==1)fs%U=0
    where(fs%vmask==1)fs%V=0
    where(fs%wmask==1)fs%W=0
    if(cfg%iproc==1)then
     fs%U(cfg%imin,:,:)=0
     do n=1,cfg%no
      i=cfg%imin-n;fs%V(i,:,:)=-fs%V(cfg%imin+n-1,:,:);fs%W(i,:,:)=-fs%W(cfg%imin+n-1,:,:)
     enddo
    endif
    if(cfg%iproc==cfg%npx)then
     fs%U(cfg%imax+1,:,:)=0
     do n=1,cfg%no
      i=cfg%imax+n;fs%V(i,:,:)=-fs%V(cfg%imax-n+1,:,:);fs%W(i,:,:)=-fs%W(cfg%imax-n+1,:,:)
     enddo
    endif
    if(cfg%kproc==1)then
     fs%W(:,:,cfg%kmin)=0
     do n=1,cfg%no
      k=cfg%kmin-n;fs%U(:,:,k)=-fs%U(:,:,cfg%kmin+n-1);fs%V(:,:,k)=-fs%V(:,:,cfg%kmin+n-1)
     enddo
    endif
    if(cfg%kproc==cfg%npz)then
     fs%W(:,:,cfg%kmax+1)=0
     do n=1,cfg%no
      k=cfg%kmax+n;fs%U(:,:,k)=-fs%U(:,:,cfg%kmax-n+1);fs%V(:,:,k)=-fs%V(:,:,cfg%kmax-n+1)
     enddo
    endif
    if(cfg%jproc==1)then
     fs%V(:,cfg%jmin,:)=0
     do n=1,cfg%no
      j=cfg%jmin-n;fs%U(:,j,:)=-fs%U(:,cfg%jmin+n-1,:);fs%W(:,j,:)=-fs%W(:,cfg%jmin+n-1,:)
     enddo
    endif
    if(cfg%jproc==cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmino_,cfg%kmaxo_;do i=cfg%imino_,cfg%imaxo_
      bc=0;if(top_port(i,k))bc=port_velocity_factor*port_shape(i,k)
      fs%V(i,j+1,k)=bc
      do n=1,cfg%no
       if(top_port(i,k).and.bc>=0)then
        fs%U(i,j+n,k)=fs%U(i,j-n+1,k);fs%W(i,j+n,k)=fs%W(i,j-n+1,k)
       else
        fs%U(i,j+n,k)=-fs%U(i,j-n+1,k);fs%W(i,j+n,k)=-fs%W(i,j-n+1,k)
       endif
       if(n<cfg%no)fs%V(i,j+1+n,k)=2*bc-fs%V(i,j+1-n,k)
      enddo
     enddo;enddo
    endif
   end subroutine
   subroutine open_vof()
    integer::i,j,k
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmino_,cfg%kmaxo_;do i=cfg%imino_,cfg%imaxo_
      if(fs%V(i,j+1,k).lt.0)vf%VF(i,j+1:,k)=0.0_WP
     enddo;enddo
    endif
    call cfg%sync(vf%VF)
   end subroutine

   subroutine distribute_expansion()
    use irl_fortran_interface,only:calculateCentroid
    integer::i,j,k,ierr
    real(WP)::m,local(2),global(2)
    ! Phase conversion removes the liquid mass BEFORE geometry transport.
    ! Geometric transport then receives the remaining gas volume S*dt*V.
    ! Together: liquid source=-m/rhol; gas source=+m/rhog, without a pure-gas search.
    evp%mflux=qphase/hfg;evp%div_src=0;evp%mfluxG=evp%mflux
    if(inventory_correction_enabled)evp%mdotdp=0
    vf%liquid_volume_source=0;vf%gas_volume_source=0;vf%phase_source_point=0
    local=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(cfg%VF(i,j,k).eq.0.or.evp%mflux(i,j,k).eq.0)cycle
     m=evp%mflux(i,j,k);evp%div_src(i,j,k)=m*(1/rg-1/rl)
     if(inventory_correction_enabled.and.vf%SD(i,j,k).gt.0)evp%mdotdp(i,j,k)=m/vf%SD(i,j,k)
     vf%liquid_volume_source(i,j,k)=-m/rl+wall_birth_m(i,j,k)/rg
     vf%gas_volume_source(i,j,k)=(m-wall_birth_m(i,j,k))/rg
     if(wall_birth_m(i,j,k).gt.0)then
      vf%phase_source_point(:,i,j,k)=wall_birth_point(:,i,j,k)
     else
      vf%phase_source_point(:,i,j,k)=calculateCentroid(vf%interface_polygon(1,i,j,k))
     endif
     local(1)=local(1)+m*cfg%vol(i,j,k);local(2)=local(2)+evp%div_src(i,j,k)*cfg%vol(i,j,k)
    enddo;enddo;enddo
    call cfg%sync(evp%mflux);call cfg%sync(evp%mfluxG);call cfg%sync(evp%div_src)
    if(inventory_correction_enabled)call cfg%sync(evp%mdotdp)
    call cfg%sync(vf%gas_volume_source)
    call MPI_ALLREDUCE(local,global,2,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    source_error=abs(global(2)-global(1)*(1/rg-1/rl))*physical_span/Lz
    call update_volume_port_rate()
   end subroutine

   subroutine project_open(dt)
    real(WP),intent(in)::dt
    integer::i,j,k,ierr,n,ia,ja,ka
    real(WP)::local(2),global(2),facearea,gasface,local_rhs,global_rhs
    call fs%update_laplacian()
    ! Prescribed boundary velocities imply Neumann pressure increments. Pin one fluid DOF only as a gauge and remove its symmetric neighbor columns.
    ia=cfg%imin;ja=cfg%jmax;ka=cfg%kmin
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(i.eq.ia.and.j.eq.ja.and.k.eq.ka)then
      ps%opr(:,i,j,k)=0;ps%opr(1,i,j,k)=1
     else
      do n=2,size(ps%stc,1)
       if(i+ps%stc(n,1).eq.ia.and.j+ps%stc(n,2).eq.ja.and.k+ps%stc(n,3).eq.ka)ps%opr(n,i,j,k)=0
      enddo
     endif
    enddo;enddo;enddo
    call ps%setup();ps%rhs=-cfg%vol*fs%div/dt;ps%sol=0
    if(ia.ge.cfg%imin_.and.ia.le.cfg%imax_.and.ja.ge.cfg%jmin_.and.ja.le.cfg%jmax_.and.ka.ge.cfg%kmin_.and.ka.le.cfg%kmax_)ps%rhs(ia,ja,ka)=0
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
      resV(i,j+1,k)=0
      ! Neumann pressure increment leaves prescribed normal velocity unchanged.
     enddo;enddo
    endif
    call cfg%sync(resV)
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     fs%U(i,j,k)=fs%U(i,j,k)-dt*resU(i,j,k)/fs%rho_U(i,j,k)
     fs%V(i,j,k)=fs%V(i,j,k)-dt*resV(i,j,k)/fs%rho_V(i,j,k)
     fs%W(i,j,k)=fs%W(i,j,k)-dt*resW(i,j,k)/fs%rho_W(i,j,k)
    enddo;enddo;enddo
    fs%P=fs%P+ps%sol
    ! Do not reapply the Neumann normal-face copy after pressure correction.
    where(fs%umask.eq.1)fs%U=0
    where(fs%vmask.eq.1)fs%V=0
    if(cfg%iproc.eq.1)fs%U(cfg%imin,:,:)=0
    if(cfg%iproc.eq.cfg%npx)fs%U(cfg%imax+1,:,:)=0
    if(cfg%jproc.eq.1)fs%V(:,cfg%jmin,:)=0
    call wall_values()
    call cfg%sync(fs%U);call cfg%sync(fs%V);call cfg%sync(fs%W)
    call fs%get_div(src=evp%div_src)
    local=0
    if(cfg%jproc.eq.cfg%npy)then
     j=cfg%jmax
     do k=cfg%kmin_,cfg%kmax_;do i=cfg%imin_,cfg%imax_
      if(.not.top_port(i,k))cycle
      facearea=cfg%dx(i)*cfg%dz(k)
      local(1)=local(1)+fs%V(i,j+1,k)*facearea
      gasface=1
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

   subroutine partition_phase_power(parts)
    real(WP),intent(out)::parts(2)
    real(WP)::local(2),vol
    integer::i,j,k,ierr
    local=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(cfg%VF(i,j,k).eq.0)cycle
     vol=cfg%vol(i,j,k)*physical_span/Lz
     local(1)=local(1)+max(qphase(i,j,k),0._WP)*vol
     local(2)=local(2)+min(qphase(i,j,k),0._WP)*vol
    enddo;enddo;enddo
    call MPI_ALLREDUCE(local,parts,2,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
   end subroutine
   subroutine record_history(dt)
    real(WP),intent(in)::dt
    real(WP)::sensible,gasvolume,local(3),global(3),masserr,totalerr,speed
    integer::i,j,k,ierr,d,ii,jj,kk,axis
    integer,parameter::di(6)=[1,-1,0,0,0,0],dj(6)=[0,0,1,-1,0,0],dk(6)=[0,0,0,0,1,-1]
    real(WP)::gaslocal(3),gasglobal(3),volume,p1(3),p2(3),area,gasface,wall_local,wall_global
    real(WP)::phase_parts(2),vertices(3,4),hh,zlocal,zglobal
    call partition_phase_power(phase_parts)
    call global_stores(sensible,gasvolume);local=-huge(1.0_WP);gaslocal=0
    wall_local=-huge(1.0_WP);zlocal=0
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     if(cfg%VF(i,j,k).ne.0)cycle
     if(abs(cfg%ym(j)-block_height+.5_WP*cfg%dy(j)).gt.64*epsilon(1._WP)*cfg%dy(j))cycle
     if((cfg%xm(i)-wall_birth_x)**2+(cfg%zm(k)-wall_birth_z)**2>(.5_WP*wall_birth_width)**2)cycle
     wall_local=max(wall_local,Tsat+theta_s(i,j,k))
    enddo;enddo;enddo
    call MPI_ALLREDUCE(wall_local,wall_global,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     local(3)=max(local(3),-T(i,j,k))
     if(cfg%VF(i,j,k).gt.0)then
      volume=cfg%vol(i,j,k)*physical_span/Lz
      gaslocal(1)=gaslocal(1)+(1-vf%VF(i,j,k))*vf%Gbary(1,i,j,k)*volume
      gaslocal(2)=gaslocal(2)+(1-vf%VF(i,j,k))*vf%Gbary(2,i,j,k)*volume
      zlocal=zlocal+(1-vf%VF(i,j,k))*vf%Gbary(3,i,j,k)*volume
      do d=1,6
       ii=i+di(d);jj=j+dj(d);kk=k+dk(d)
       if(ii<cfg%imin.or.ii>cfg%imax.or.jj<cfg%jmin.or.jj>cfg%jmax.or.kk<cfg%kmin.or.kk>cfg%kmax)cycle
       if(cfg%VF(ii,jj,kk)>0.or.cold_solid(ii,jj,kk))cycle
       call face_vertices(i,j,k,d,vertices,axis,hh,area)
       gasface=1-face_liquid3d(vf,i,j,k,d);gaslocal(3)=gaslocal(3)+gasface*area
      enddo
     endif

     if(cfg%VF(i,j,k).eq.0)then
      local(1)=max(local(1),T(i,j,k))
     else
      local(2)=max(local(2),T(i,j,k))
     endif
    enddo;enddo;enddo
    call MPI_ALLREDUCE(local,global,3,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    call MPI_ALLREDUCE(gaslocal,gasglobal,3,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    if(gasvolume.gt.0)gasglobal(1:2)=gasglobal(1:2)/gasvolume
    masserr=(rg*(gasvolume-previous_vapor)/dt+rg*top_gas_out-phase_power/hfg)/max(abs(phase_power/hfg),1e-14_WP)
    totalerr=sensible-initial_sensible+rg*hfg*(gasvolume-initial_vapor)+out_integral-input_integral+wall_birth_surface_J
    call MPI_ALLREDUCE(zlocal,zglobal,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
    if(gasvolume>0)zglobal=zglobal/gasvolume
    speed=sqrt(fs%Umax**2+fs%Vmax**2+fs%Wmax**2)
    if(.not.ieee_is_finite(speed+sum(global)+gasvolume+sensible+phase_power))call die('Nonfinite physical field')
    if(cfg%amRoot)then
     write(history_unit,'(i0,37(",",es25.16e3))')time%n,time%t,dt,global(1:2),speed,gasvolume,phase_power,sensible,rg*hfg*(gasvolume-initial_vapor),input_integral,out_integral,totalerr,thermal_balance,masserr,fs%divmax,ps%rerr,source_error,top_volume_out,top_gas_out,-global(3),real(coupling_iterations,WP),coupling_change,real(inactive_thermal_interfaces,WP),inactive_thermal_surface_m2,phase_cfl_actual,real(thermal_trial_count,WP),gasglobal,wall_global,real(wall_birth_events,WP),wall_birth_surface_J,phase_parts(1),phase_parts(2),cold_heat_out,fs%Wmax,zglobal
     flush(history_unit)
    endif
    previous_vapor=gasvolume
   end subroutine

   subroutine write_checkpoint()
    use irl_fortran_interface,only:ByteBuffer_type,new,serializeAndPack,unpackAndStore,getSize,setSize,dataPtr
    integer::u,i,j,k,ierr,ios,magic,version,bits,bounds(12),partition(3),rank_saved
    integer(kind=8)::nbytes
    integer(kind=1),pointer::bytes(:)
    type(ByteBuffer_type)::buffer
    character(len=1200)::name
    character(len=64)::identity
    if(last_checkpoint_step.eq.time%n)return
    last_checkpoint_step=time%n
    bounds=[cfg%imin_,cfg%imax_,cfg%jmin_,cfg%jmax_,cfg%kmin_,cfg%kmax_,cfg%imino_,cfg%imaxo_,cfg%jmino_,cfg%jmaxo_,cfg%kmino_,cfg%kmaxo_];partition=[cfg%npx,cfg%npy,cfg%npz]
    write(name,'("checkpoint_step",i8.8,"_rank",i4.4,".bin")')time%n,cfg%rank
    call new(buffer)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call serializeAndPack(vf%liquid_gas_interface(i,j,k),buffer)
    enddo;enddo;enddo
    nbytes=getSize(buffer);bytes=>dataPtr(buffer)
    open(newunit=u,file=trim(name),form='unformatted',access='stream',status='new',iostat=ios)
    if(ios.ne.0)call die('Cannot create checkpoint without overwriting')
    write(u)20261116,4,storage_size(1.0_WP),bounds,partition,cfg%rank,checkpoint_identity
    write(u)time%n,time%it,time%t,time%dt,time%told,time%dtold,time%tmid,time%dtmid,time%cfl
    write(u)input_integral,phase_integral,out_integral,initial_sensible,initial_vapor,previous_vapor,top_volume_out,top_gas_out,source_error,phase_cfl_actual,thermal_trial_count,roundoff_energy,phase_power,thermal_out,thermal_balance,coupling_change,inactive_thermal_surface_m2,coupling_iterations,inactive_thermal_interfaces
    ! Version4: accepted event state, physical ledgers, and last source geometry.
    write(u)wall_birth_enabled,wall_birth_events,wall_birth_done,wall_birth_trial,wall_birth_event_mass,wall_birth_latent_J,wall_birth_surface_J
    write(u)wall_birth_m,wall_birth_point
    write(u)fs%U,fs%V,fs%W,fs%P,fs%rho_U,fs%rho_V,fs%rho_W,fs%DPjx,fs%DPjy,fs%DPjz,fs%div
    write(u)fs%Pjx,fs%Pjy,fs%Pjz
    write(u)vf%VF,vf%VFold,vf%Lbary,vf%Gbary,vf%SD,vf%G,vf%curv,vf%Lvol,vf%Gvol
    write(u)T,E,El,Eg,theta_l,theta_g,theta_s,qphase
    write(u)evp%div_src,evp%mdotdp,evp%mflux,evp%mfluxL,evp%mfluxG,evp%U_itf,evp%V_itf,evp%W_itf
    write(u)nbytes
    write(u)bytes
    close(u)
    nullify(bytes)
    call MPI_BARRIER(cfg%comm,ierr)
    if(cfg%amRoot)then
     open(newunit=u,file='CHECKPOINT_INDEX.csv',status='unknown',position='append')
     write(u,'(i0,",",es25.16e3)')time%n,time%t
     close(u)
    endif
   end subroutine

   subroutine read_checkpoint(prefix)
    use tpns_class,only:static_contact
    use irl_fortran_interface,only:ByteBuffer_type,new,serializeAndPack,unpackAndStore,getSize,setSize,dataPtr,resetBufferPointer
    integer::u,i,j,k,ierr,ios,magic,version,bits,bounds(12),partition(3),rank_saved,birth_enabled_saved
    integer(kind=8)::nbytes
    integer(kind=1),pointer::bytes(:)
    type(ByteBuffer_type)::buffer,check_buffer
    integer(kind=1),pointer::check_bytes(:)
    real(WP)::zero_dt
    real(WP),allocatable::saved_sd(:,:,:),keep_x(:,:,:),keep_y(:,:,:),keep_z(:,:,:)
    character(len=1200)::name
    character(len=64)::identity
    character(len=*),intent(in)::prefix
    write(name,'(a,"_rank",i4.4,".bin")')trim(prefix),cfg%rank
    open(newunit=u,file=trim(name),form='unformatted',access='stream',status='old',iostat=ios)
    if(ios.ne.0)call die('Missing exact checkpoint rank file')
    read(u,iostat=ios)magic,version,bits,bounds,partition,rank_saved,identity
    if(ios.ne.0)call die('Truncated checkpoint header')
    if(magic.ne.20261116.or.(version.ne.2.and.version.ne.3.and.version.ne.4).or.bits.ne.storage_size(1.0_WP).or.any(bounds.ne.[cfg%imin_,cfg%imax_,cfg%jmin_,cfg%jmax_,cfg%kmin_,cfg%kmax_,cfg%imino_,cfg%imaxo_,cfg%jmino_,cfg%jmaxo_,cfg%kmino_,cfg%kmaxo_]).or.any(partition.ne.[cfg%npx,cfg%npy,cfg%npz]).or.rank_saved.ne.cfg%rank.or.identity.ne.checkpoint_identity)call die('Checkpoint identity, precision or decomposition mismatch')
    read(u,iostat=ios)time%n,time%it,time%t,time%dt,time%told,time%dtold,time%tmid,time%dtmid,time%cfl
    if(ios.ne.0)call die('Truncated checkpoint time state')
    read(u,iostat=ios)input_integral,phase_integral,out_integral,initial_sensible,initial_vapor,previous_vapor,top_volume_out,top_gas_out,source_error,phase_cfl_actual,thermal_trial_count,roundoff_energy,phase_power,thermal_out,thermal_balance,coupling_change,inactive_thermal_surface_m2,coupling_iterations,inactive_thermal_interfaces
    if(ios.ne.0)call die('Truncated checkpoint energy ledgers')
    if(version.eq.4)then
     read(u,iostat=ios)birth_enabled_saved,wall_birth_events,wall_birth_done,wall_birth_trial,wall_birth_event_mass,wall_birth_latent_J,wall_birth_surface_J
     if(ios.ne.0)call die('Truncated checkpoint wall birth state')
     if(birth_enabled_saved.ne.wall_birth_enabled)call die('Checkpoint wall birth enable mismatch')
     if(wall_birth_events.lt.0.or.wall_birth_events.gt.1.or.(wall_birth_done.neqv.(wall_birth_events.eq.1)))call die('Invalid checkpoint wall birth event state')
     read(u,iostat=ios)wall_birth_m,wall_birth_point
     if(ios.ne.0)call die('Truncated checkpoint wall birth source geometry')
     if(cfg%amRoot)print*,'WALL_BIRTH_STATE_RESTORED',wall_birth_events,wall_birth_done,wall_birth_latent_J,wall_birth_surface_J
    else
     if(wall_birth_enabled.eq.1)call die('Older checkpoint lacks wall birth state')
    endif
    read(u,iostat=ios)fs%U,fs%V,fs%W,fs%P,fs%rho_U,fs%rho_V,fs%rho_W,fs%DPjx,fs%DPjy,fs%DPjz,fs%div
    if(ios.ne.0)call die('Truncated checkpoint array group')
    if(version.ge.3)then
     read(u,iostat=ios)fs%Pjx,fs%Pjy,fs%Pjz
     if(ios.ne.0)call die('Truncated persistent surface pressure jump')
    endif
    read(u,iostat=ios)vf%VF,vf%VFold,vf%Lbary,vf%Gbary,vf%SD,vf%G,vf%curv,vf%Lvol,vf%Gvol
    if(ios.ne.0)call die('Truncated checkpoint array group')
    read(u,iostat=ios)T,E,El,Eg,theta_l,theta_g,theta_s,qphase
    if(ios.ne.0)call die('Truncated checkpoint array group')
    read(u,iostat=ios)evp%div_src,evp%mdotdp,evp%mflux,evp%mfluxL,evp%mfluxG,evp%U_itf,evp%V_itf,evp%W_itf
    if(ios.ne.0)call die('Truncated checkpoint array group')
    read(u,iostat=ios)nbytes
    if(ios.ne.0)call die('Missing checkpoint IRL length')
    if(nbytes.le.0.or.nbytes.gt.128_8*size(vf%VF,kind=8))call die('Invalid checkpoint IRL byte length')
    call new(buffer);call setSize(buffer,nbytes);bytes=>dataPtr(buffer)
    read(u,iostat=ios)bytes
    if(ios.ne.0)call die('Truncated checkpoint IRL payload')
    read(u,iostat=ios)magic
    if(ios.ge.0)call die('Unexpected checkpoint trailing bytes')
    close(u);nullify(bytes)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call unpackAndStore(vf%liquid_gas_interface(i,j,k),buffer)
    enddo;enddo;enddo
    ! Same decomposition: exact halo arrays and separators were saved as doubles/bytes.
    ! Recreate only workload maps and polygon objects, not a new interface fit.
    allocate(saved_sd,mold=vf%SD);saved_sd=vf%SD
    call vf%update_band();call vf%polygonalize_interface()
    ! polygonalize_interface can rewrite separators of full/empty neighbor cells.
    ! Recreate polygon objects, then restore the exact saved separator bytes and SD.
    call new(check_buffer)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call serializeAndPack(vf%liquid_gas_interface(i,j,k),check_buffer)
    enddo;enddo;enddo
    bytes=>dataPtr(buffer);check_bytes=>dataPtr(check_buffer)
    print*,'CHECKPOINT_POLYGON_REBUILD_DELTA',cfg%rank,nbytes,getSize(check_buffer),count(bytes(1:min(size(bytes),size(check_bytes))).ne.check_bytes(1:min(size(bytes),size(check_bytes)))),maxval(abs(saved_sd-vf%SD))
    nullify(bytes);nullify(check_bytes)
    call resetBufferPointer(buffer)
    do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
     call unpackAndStore(vf%liquid_gas_interface(i,j,k),buffer)
    enddo;enddo;enddo
    vf%SD=saved_sd;deallocate(saved_sd)
    ! Version2 checkpoint omitted the persistent Pjx/Pjy/Pjz momentum source.
    ! Reconstruct this static-contact source from the saved interface/curvature,
    ! retaining the already-saved incremental pressure-jump arrays.
    if(version.eq.2)then
    allocate(keep_x,mold=fs%DPjx);allocate(keep_y,mold=fs%DPjy);allocate(keep_z,mold=fs%DPjz)
    keep_x=fs%DPjx;keep_y=fs%DPjy;keep_z=fs%DPjz
    zero_dt=0.0_WP
    call fs%add_surface_tension_jump(dt=zero_dt,div=fs%div,vf=vf,contact_model=static_contact)
    fs%DPjx=keep_x;fs%DPjy=keep_y;fs%DPjz=keep_z
    deallocate(keep_x,keep_y,keep_z)
    endif
    call fs%interp_vel(Ui,Vi,Wi);call fs%get_max();call vf%get_max()
    if(cfg%amRoot)print*,'EXACT_CHECKPOINT_READ',time%n,time%t
   end subroutine


   subroutine trace_stage(stage,dt,force_available)
    character(len=*),intent(in)::stage
    real(WP),intent(in)::dt
    logical,intent(in)::force_available
    real(WP)::local(4),global(4)
    integer::ierr
    local(1)=maxval(abs(fs%U(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)))
    local(2)=maxval(abs(fs%V(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)))
    local(3)=maxval(abs(fs%W(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)))
    local(4)=-min(minval(fs%rho_U(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)),minval(fs%rho_V(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)),minval(fs%rho_W(cfg%imin_:cfg%imax_+1,cfg%jmin_:cfg%jmax_+1,cfg%kmin_:cfg%kmax_+1)))
    call MPI_ALLREDUCE(local,global,4,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
    if(.not.all(ieee_is_finite(global)).or.global(4)>=0)then
     call trace_dump(stage,dt);call die('Nonfinite/nonpositive3D staggered field')
    endif
   end subroutine
   subroutine trace_dump(stage,dt)
    character(len=*),intent(in)::stage
    real(WP),intent(in)::dt
    integer::u,bounds(12)
    character(len=150)::name
    bounds=[cfg%imin_,cfg%imax_,cfg%jmin_,cfg%jmax_,cfg%kmin_,cfg%kmax_,cfg%imino_,cfg%imaxo_,cfg%jmino_,cfg%jmaxo_,cfg%kmino_,cfg%kmaxo_]
    write(name,'("DIAGNOSTIC_STAGE_step",i8.8,"_it",i1,"_rank",i4.4,".bin")')time%n,time%it,cfg%rank
    open(newunit=u,file=trim(name),status='new',form='unformatted',access='stream')
    write(u)20261116,99,storage_size(1._WP),bounds,cfg%rank,time%n,time%it,time%t,dt
    write(u)fs%U,fs%V,fs%W,fs%P,fs%Uold,fs%Vold,fs%Wold,fs%rho_Uold,fs%rho_Vold,fs%rho_Wold,fs%rho_U,fs%rho_V,fs%rho_W,vf%VF,vf%SD,vf%curv,evp%div_src,qphase,fs%Pjx,fs%Pjy,T,theta_l,theta_g,cfg%VF
    close(u)
    if(trace_unit.ne.-1)flush(trace_unit)
    print*,'DIAGNOSTIC_STATE_PRESERVED',trim(stage),trim(name)
   end subroutine

   subroutine output_fields()
    integer::u,i,j,k
    character(len=100)::name
    if(time%n.eq.last_field_output_step)return
    last_field_output_step=time%n
    write(name,'("field_step",i8.8,"_rank",i4.4,".csv")')time%n,cfg%rank
    open(newunit=u,file=trim(name),status='new')
    write(u,'(a)')'x_m,y_m,z_m,T_K,liquid_fraction,u_m_s,v_m_s,w_m_s,solid,phase_heat_W_m3,liquid_T_K,vapor_T_K'
    do k=cfg%kmin_,cfg%kmax_;do j=cfg%jmin_,cfg%jmax_;do i=cfg%imin_,cfg%imax_
     write(u,'(11(es25.16e3,","),es25.16e3)')cfg%xm(i),cfg%ym(j),cfg%zm(k),T(i,j,k),vf%VF(i,j,k),Ui(i,j,k),Vi(i,j,k),Wi(i,j,k),1-cfg%VF(i,j,k),qphase(i,j,k),Tsat+theta_l(i,j,k),Tsat+theta_g(i,j,k)
    enddo;enddo;enddo
    close(u)
    if(cfg%amRoot)then
     open(newunit=u,file='PHASE_ENERGY_ROUNDOFF.txt',status='replace')
     write(u,*)roundoff_energy
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