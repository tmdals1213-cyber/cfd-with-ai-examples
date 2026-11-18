module simulation
   use precision, only: WP
   use geometry, only: cfg,physical_span
   use cht_energy
   use water_eos, only: reference_T
   use lowmach_class, only: lowmach,bcond,neumann,dirichlet
   use hypre_str_class, only: hypre_str,pcg_pfmg
   use timetracker_class, only: timetracker
   use ensight_class, only: ensight
   use event_class, only: event
   use mpi_f08
   use parallel, only: MPI_REAL_WP
   use messager, only: die
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final
   type(lowmach) :: fs
   type(hypre_str) :: ps
   type(timetracker) :: time
   type(ensight) :: ens
   type(event) :: out_event
   real(WP), allocatable :: rhsU(:,:,:),rhsV(:,:,:),rhsW(:,:,:),gx(:,:,:),gy(:,:,:),gz(:,:,:)
   real(WP), allocatable :: rho_old(:,:,:),momU_old(:,:,:),momV_old(:,:,:),momW_old(:,:,:),pold(:,:,:)
   real(WP), allocatable :: drhodt(:,:,:),Tguess(:,:,:),Ui(:,:,:),Vi(:,:,:),Wi(:,:,:)
   real(WP) :: gravity(3),Uin,integrated_out=0.0_WP
   integer :: history_unit,mass_history_unit
contains
   function left(pg,i,j,k) result(yes)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: yes
      yes=i.eq.pg%imin
   end function
   function right(pg,i,j,k) result(yes)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: yes
      yes=i.eq.pg%imax+1
   end function
   function bottom(pg,i,j,k) result(yes)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: yes
      yes=j.eq.pg%jmin
   end function
   function top(pg,i,j,k) result(yes)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: yes
      yes=j.eq.pg%jmax+1
   end function

   subroutine simulation_init
      use param, only: param_read
      type(bcond), pointer :: bc
      integer :: n,i,j,k
      call energy_init()
      fs=lowmach(cfg=cfg,name='IAPWS variable-density water')
      fs%rho=density; fs%visc=viscosity
      call param_read('Gravity',gravity); call param_read('Inlet velocity',Uin)
      if (Uin.eq.0.0_WP) then
         call fs%add_bcond(name='left reservoir',type=neumann,locator=left,face='x',dir=-1,canCorrect=.true.)
      else
         call fs%add_bcond(name='left reservoir',type=dirichlet,locator=left,face='x',dir=-1,canCorrect=.false.)
      end if
      call fs%add_bcond(name='right reservoir',type=neumann,locator=right,face='x',dir=1,canCorrect=.true.)
      call fs%add_bcond(name='bottom reservoir',type=neumann,locator=bottom,face='y',dir=-1,canCorrect=.true.)
      call fs%add_bcond(name='top reservoir',type=neumann,locator=top,face='y',dir=1,canCorrect=.true.)
      ps=hypre_str(cfg=cfg,name='Mass continuity projection',method=pcg_pfmg,nst=7)
      call param_read('Pressure iteration',ps%maxit); call param_read('Pressure tolerance',ps%rcvg)
      call fs%setup(pressure_solver=ps)
      ! The reservoir projection has Neumann pressure conditions. Remove its constant
      ! nullspace by fixing one FLUID pressure-correction degree of freedom to zero.
      ! Remove its column too, preserving symmetry for PCG; the omitted continuity
      ! equation follows from the independently enforced global mass balance.
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         if (i.eq.cfg%imin.and.j.eq.cfg%jmin.and.k.eq.cfg%kmin) then
            ps%opr(:,i,j,k)=0.0_WP; ps%opr(1,i,j,k)=1.0_WP
         end if
         if (i.eq.cfg%imin+1.and.j.eq.cfg%jmin.and.k.eq.cfg%kmin) ps%opr(3,i,j,k)=0.0_WP
         if (i.eq.cfg%imin.and.j.eq.cfg%jmin+1.and.k.eq.cfg%kmin) ps%opr(5,i,j,k)=0.0_WP
      end do; end do; end do
      call ps%setup()
      allocate(rhsU,mold=temperature); allocate(rhsV,mold=temperature); allocate(rhsW,mold=temperature)
      allocate(gx,mold=temperature); allocate(gy,mold=temperature); allocate(gz,mold=temperature)
      allocate(rho_old,mold=temperature); allocate(momU_old,mold=temperature); allocate(momV_old,mold=temperature)
      allocate(momW_old,mold=temperature); allocate(pold,mold=temperature); allocate(drhodt,mold=temperature)
      allocate(Tguess,mold=temperature); allocate(Ui,mold=temperature); allocate(Vi,mold=temperature); allocate(Wi,mold=temperature)
      fs%U=Uin; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
      where (fs%umask.eq.1) fs%U=0.0_WP
      call fs%rho_multiply()
      call fs%apply_bcond(0.0_WP,0.0_WP)
      call fs%get_mfr()
      call fs%interp_vel(Ui,Vi,Wi)
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax); time%dt=time%dtmax
      call param_read('Max cfl number',time%cflmax); call param_read('Max time',time%tmax)
      time%print_info=.false.
      ens=ensight(cfg=cfg,name='water_CHT')
      out_event=event(time=time,name='Fields')
      call param_read('Ensight output period',out_event%tper)
      call ens%add_vector('velocity_m_s',Ui,Vi,Wi)
      call ens%add_scalar('temperature_K',temperature)
      call ens%add_scalar('water_density_kg_m3',density)
      call ens%add_scalar('solid_levelset_m',cfg%Gib)
      call ens%add_scalar('dynamic_pressure_Pa',fs%P)
      if (out_event%occurs()) call ens%write_data(time%t)
      if (cfg%amRoot) then
         open(newunit=history_unit,file='history.csv',status='new')
         open(newunit=mass_history_unit,file='mass_compatibility.csv',status='new')
         write(mass_history_unit,'(a)') 'step,time_s,global_error_before_kg_s,global_error_after_kg_s,dMdt_kg_s'
         write(history_unit,'(a)') 'step,time_s,dt_s,coupling_iterations,Tsolid_max_K,Tfluid_max_K,rho_min_kgm3,rho_max_kgm3,Umax_ms,Vmax_ms,mass_residual_kgm3s,pressure_relative_error,applied_CFL,Esolid_J,Hfluid_J,power_W,Qouter_out_W,Qinto_fluid_W,Qinterface_sum_W,heat_budget_relative_error'
      end if
   end subroutine

   subroutine project_flow(dt)
      real(WP), intent(in) :: dt
      integer :: i,j,k,ierr
      real(WP) :: local_rhs_norm,global_rhs_norm
      fs%rho=density; fs%visc=viscosity
      fs%rhoU=momU_old+dt*rhsU
      fs%rhoV=momV_old+dt*rhsV
      fs%rhoW=momW_old+dt*rhsW
      ! Exact EOS density difference; uniform reference gravity is absorbed in hydrostatic pressure.
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         if (fs%umask(i,j,k).eq.0) fs%rhoU(i,j,k)=fs%rhoU(i,j,k)+dt*gravity(1)* &
            sum(fs%itpr_x(:,i,j,k)*(density(i-1:i,j,k)-rho_reference))
         if (fs%vmask(i,j,k).eq.0) fs%rhoV(i,j,k)=fs%rhoV(i,j,k)+dt*gravity(2)* &
            sum(fs%itpr_y(:,i,j,k)*(density(i,j-1:j,k)-rho_reference))
         if (fs%wmask(i,j,k).eq.0) fs%rhoW(i,j,k)=fs%rhoW(i,j,k)+dt*gravity(3)* &
            sum(fs%itpr_z(:,i,j,k)*(density(i,j,k-1:k)-rho_reference))
      end do; end do; end do
      call fs%rho_divide()
      call fs%apply_bcond(time%t,dt)
      drhodt=(density-rho_old)/dt
      where (cfg%VF.eq.0.0_WP) drhodt=0.0_WP
      call correct_reservoir_mass_flux()
      call fs%get_div(drhodt)
      fs%psolv%rhs=-cfg%vol*fs%div/dt
      if (cfg%imin_.le.cfg%imin.and.cfg%jmin_.le.cfg%jmin.and.cfg%kmin_.le.cfg%kmin) &
         fs%psolv%rhs(cfg%imin,cfg%jmin,cfg%kmin)=0.0_WP
      fs%psolv%sol=0.0_WP
      local_rhs_norm=maxval(abs(fs%psolv%rhs))
      call MPI_ALLREDUCE(local_rhs_norm,global_rhs_norm,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
      if (global_rhs_norm.eq.0.0_WP) then
         ! Exact homogeneous system: avoid an undefined relative residual of zero/zero.
         fs%psolv%it=0; fs%psolv%rerr=0.0_WP
      else
         call fs%psolv%solve()
      end if
      if (fs%psolv%rerr.gt.1.01e-9_WP) call die('Mass-projection linear residual failed')
      call fs%shift_p(fs%psolv%sol)
      call fs%get_pgrad(fs%psolv%sol,gx,gy,gz)
      fs%rhoU=fs%rhoU-dt*gx; fs%rhoV=fs%rhoV-dt*gy; fs%rhoW=fs%rhoW-dt*gz
      fs%P=pold+fs%psolv%sol
      call fs%rho_divide()
      call fs%get_div(drhodt)
      call fs%get_max()
   end subroutine

   subroutine correct_reservoir_mass_flux
      ! Native correct_mfr skips errors below 10*epsilon in dimensional kg/s.
      ! For this small SI domain that leaves a visible pointwise compatibility
      ! error. Enforce the same reservoir balance without that dimensional cutoff.
      type(bcond), pointer :: bc
      integer :: n,i,j,k,pass
      real(WP) :: dMdt,err,before,correction
      call cfg%integrate_without_vf(drhodt,dMdt)
      call fs%get_mfr()
      err=sum(fs%mfr)+dMdt; before=err
      if (fs%correctable_area.le.0.0_WP) call die('No open reservoir area for thermal expansion')
      do pass=1,2
         if (err.eq.0.0_WP) exit
         correction=-err/fs%correctable_area
         bc=>fs%first_bc
         do while (associated(bc))
            if (bc%itr%amIn.and.bc%canCorrect) then
               do n=1,bc%itr%n_
                  i=bc%itr%map(1,n);j=bc%itr%map(2,n);k=bc%itr%map(3,n)
                  select case (bc%face)
                  case ('x');fs%rhoU(i,j,k)=fs%rhoU(i,j,k)+bc%rdir*correction
                  case ('y');fs%rhoV(i,j,k)=fs%rhoV(i,j,k)+bc%rdir*correction
                  case ('z');fs%rhoW(i,j,k)=fs%rhoW(i,j,k)+bc%rdir*correction
                  end select
               end do
            end if
            bc=>bc%next
         end do
         call cfg%sync(fs%rhoU);call cfg%sync(fs%rhoV);call cfg%sync(fs%rhoW)
         call fs%get_mfr()
         err=sum(fs%mfr)+dMdt
      end do
      if (cfg%amRoot) write(mass_history_unit,'(i0,4(",",es25.16e3))') time%n,time%t,before,err,dMdt
   end subroutine

   subroutine simulation_run
      integer :: outer,thermal_iterations,ierr
      real(WP) :: cfl,conv_cfl,delta,global_delta,res,diag(7),budget,lims(4),global_lims(4)
      do while (.not.time%done())
         call fs%get_cfl(time%dt,conv_cfl,time%cfl)
         call time%adjust_dt()
         time%dt=min(time%dt,time%tmax-time%t)
         call fs%get_cfl(time%dt,conv_cfl,cfl)
         call time%increment()
         call save_energy()
         rho_old=density; momU_old=fs%rhoU; momV_old=fs%rhoV; momW_old=fs%rhoW; pold=fs%P
         call fs%get_dmomdt(rhsU,rhsV,rhsW)
         do outer=1,20
            Tguess=temperature
            call solve_energy(time%dt,fs%rhoU,fs%rhoV,thermal_iterations)
            call project_flow(time%dt)
            call energy_residual(time%dt,fs%rhoU,fs%rhoV,res)
            delta=maxval(abs(temperature-Tguess))
            call MPI_ALLREDUCE(delta,global_delta,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
            if (res.lt.1e-6_WP.and.global_delta.lt.1e-8_WP.and.fs%divmax.lt.1e-6_WP) exit
         end do
         if (outer.gt.20) then
            if (cfg%amRoot) print*, 'COUPLING_FAILED',time%n,time%t,res,global_delta,fs%divmax
            call die('EOS/enthalpy/mass-flux coupling did not converge')
         end if
         call energy_diagnostics(time%dt,fs%rhoU,fs%rhoV,diag,budget)
         if (abs(budget).gt.1e-6_WP) call die('Global conjugate enthalpy budget failed')
         if (abs(diag(5)+diag(6)).gt.1e-8_WP) call die('Interface heat was not equal and opposite')
         integrated_out=integrated_out+time%dt*diag(4)
         lims(1)=maxval(temperature,mask=cfg%VF.eq.0.0_WP)
         lims(2)=maxval(temperature,mask=cfg%VF.gt.0.0_WP)
         lims(3)=-minval(density,mask=cfg%VF.gt.0.0_WP)
         lims(4)=maxval(density,mask=cfg%VF.gt.0.0_WP)
         call MPI_ALLREDUCE(lims,global_lims,4,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         call fs%interp_vel(Ui,Vi,Wi)
         if (cfg%amRoot) then
            write(history_unit,'(i0,19(",",es25.16e3))') time%n,time%t,time%dt,real(outer,WP),global_lims(1:2), &
               -global_lims(3),global_lims(4),fs%Umax,fs%Vmax,fs%divmax,fs%psolv%rerr,cfl,diag(1:5),diag(5)+diag(6),budget
            flush(history_unit)
            if (mod(time%n,25).eq.0.or.time%done()) print*, 'WATER_CHT',time%n,time%t,'TsC',global_lims(1)-273.15_WP, &
               'Vmax',fs%Vmax,'mass',fs%divmax,'heat_budget',budget
         end if
         if (out_event%occurs().or.time%done()) call ens%write_data(time%t)
      end do
   end subroutine

   subroutine simulation_final
      integer :: rank,ierr,u,i,j,k
      character(len=50) :: name
      call MPI_COMM_RANK(cfg%comm,rank,ierr)
      write(name,'("field_rank",i4.4,".csv")') rank
      open(newunit=u,file=trim(name),status='new')
      write(u,'(a)') 'x_m,y_m,u_m_s,v_m_s,T_K,rho_water_kg_m3,is_solid'
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         write(u,'(6(es25.16e3,","),i1)') cfg%xm(i),cfg%ym(j),Ui(i,j,k),Vi(i,j,k),temperature(i,j,k),density(i,j,k), &
            merge(1,0,cfg%VF(i,j,k).eq.0.0_WP)
      end do; end do; end do
      close(u)
      if (cfg%amRoot) then
         close(history_unit)
         close(mass_history_unit)
         open(newunit=u,file='CUMULATIVE_BUDGET.txt',status='new')
         write(u,*) 'elapsed_s',time%t
         write(u,*) 'heat_generated_J',power*time%t
         write(u,*) 'heat_out_of_domain_J',integrated_out
         close(u)
      end if
   end subroutine
end module
