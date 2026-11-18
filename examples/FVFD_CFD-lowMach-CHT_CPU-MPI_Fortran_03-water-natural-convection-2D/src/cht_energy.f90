module cht_energy
   use precision, only: WP
   use geometry, only: cfg,block_width,block_height,physical_span
   use water_eos, only: water,load_water,face_conductivity,reference_T
   use hypre_str_class, only: hypre_str,gmres_pfmg
   use mpi_f08
   use parallel, only: MPI_REAL_WP
   use messager, only: die
   implicit none
   private
   public :: energy_init,update_properties,save_energy,solve_energy,energy_diagnostics,energy_residual
   public :: temperature,density,conductivity,viscosity,power,rho_reference,stored_old
   real(WP), allocatable :: temperature(:,:,:),density(:,:,:),conductivity(:,:,:),viscosity(:,:,:)
   real(WP), allocatable :: enthalpy(:,:,:),dh(:,:,:),stored(:,:,:),stored_old(:,:,:),capacity(:,:,:)
   type(hypre_str) :: hs
   real(WP) :: power,rho_reference,rhos,cps,ks,qvol,source_scale
contains
   subroutine energy_init
      use param, only: param_read
      integer :: i,j,k,ierr
      real(WP) :: rr,hh,cc,kk,mm,sv,total_sv
      call load_water('data/water_p101325.dat')
      call water(reference_T,rr,hh,cc,kk,mm); rho_reference=rr
      call param_read('Heat power',power)
      call param_read('Copper density',rhos); call param_read('Copper heat capacity',cps)
      call param_read('Copper conductivity',ks)
      if (min(rhos,cps,ks,physical_span).le.0.0_WP.or.power.lt.0.0_WP) call die('Invalid solid heat input')
      qvol=power/(block_width*block_height*physical_span)
      source_scale=max(qvol,1.0_WP)
      allocate(temperature(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(density,mold=temperature); allocate(conductivity,mold=temperature); allocate(viscosity,mold=temperature)
      allocate(enthalpy,mold=temperature); allocate(dh,mold=temperature)
      allocate(stored,mold=temperature); allocate(stored_old,mold=temperature); allocate(capacity,mold=temperature)
      temperature=reference_T
      call update_properties()
      stored_old=stored
      hs=hypre_str(cfg=cfg,name='Conjugate enthalpy',method=gmres_pfmg,nst=7)
      hs%stc(1,:)=[0,0,0]; hs%stc(2,:)=[1,0,0]; hs%stc(3,:)=[-1,0,0]
      hs%stc(4,:)=[0,1,0]; hs%stc(5,:)=[0,-1,0]; hs%stc(6,:)=[0,0,1]; hs%stc(7,:)=[0,0,-1]
      call param_read('Thermal tolerance',hs%rcvg); call param_read('Thermal iteration',hs%maxit)
      hs%opr(1,:,:,:)=1.0_WP
      call hs%init()
      sv=0.0_WP
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         if (cfg%VF(i,j,k).eq.0.0_WP) sv=sv+cfg%vol(i,j,k)*physical_span/cfg%zL
      end do; end do; end do
      call MPI_ALLREDUCE(sv,total_sv,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      if (abs(total_sv-block_width*block_height*physical_span).gt.1e-13_WP) call die('Solid heat-source volume mismatch')
   end subroutine

   subroutine update_properties
      integer :: i,j,k
      real(WP) :: r,h,cp,kappa,mu,dr,dhdT
      ! Cold exterior reservoir; incoming enthalpy is referenced to this temperature.
      if (cfg%iproc.eq.1) temperature(:cfg%imin-1,:,:)=reference_T
      if (cfg%iproc.eq.cfg%npx) temperature(cfg%imax+1:,:,:)=reference_T
      if (cfg%jproc.eq.1) temperature(:,:cfg%jmin-1,:)=reference_T
      if (cfg%jproc.eq.cfg%npy) temperature(:,cfg%jmax+1:,:)=reference_T
      call cfg%sync(temperature)
      do k=cfg%kmino_,cfg%kmaxo_; do j=cfg%jmino_,cfg%jmaxo_; do i=cfg%imino_,cfg%imaxo_
         if (cfg%VF(i,j,k).eq.0.0_WP) then
            density(i,j,k)=rho_reference
            conductivity(i,j,k)=ks; viscosity(i,j,k)=0.0_WP
            enthalpy(i,j,k)=0.0_WP; dh(i,j,k)=0.0_WP
            stored(i,j,k)=rhos*cps*(temperature(i,j,k)-reference_T)
            capacity(i,j,k)=rhos*cps
         else
            call water(temperature(i,j,k),r,h,cp,kappa,mu,dr,dhdT)
            density(i,j,k)=r; conductivity(i,j,k)=kappa; viscosity(i,j,k)=mu
            enthalpy(i,j,k)=h; dh(i,j,k)=dhdT
            stored(i,j,k)=r*h
            capacity(i,j,k)=r*dhdT+dr*h
            if (capacity(i,j,k).le.0.0_WP) call die('Invalid derivative of liquid enthalpy density')
         end if
      end do; end do; end do
   end subroutine

   subroutine save_energy
      stored_old=stored
   end subroutine

   subroutine build_system(dt,mx,my,diag)
      real(WP), intent(in) :: dt,mx(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),my(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
      real(WP), intent(out) :: diag(7)
      integer :: i,j,k,d,ii,jj,slot,ierr
      integer, parameter :: di(4)=[1,-1,0,0],dj(4)=[0,0,1,-1],slots(4)=[2,3,4,5]
      real(WP) :: r,a,m,hup,slope,volume,dist,Tn,c,thermal_into,adv_out,local(7)
      logical :: exterior,isfluid,otherfluid
      hs%opr=0.0_WP; hs%rhs=0.0_WP
      local=0.0_WP
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         isfluid=cfg%VF(i,j,k).gt.0.0_WP
         volume=cfg%vol(i,j,k)*physical_span/cfg%zL
         r=-(stored(i,j,k)-stored_old(i,j,k))/dt
         hs%opr(1,i,j,k)=capacity(i,j,k)/dt
         if (.not.isfluid) then
            r=r+qvol
            local(3)=local(3)+qvol*volume
            local(1)=local(1)+stored(i,j,k)*volume
         else
            local(2)=local(2)+stored(i,j,k)*volume
         end if
         do d=1,4
            ii=i+di(d); jj=j+dj(d); slot=slots(d)
            if (d.le.2) then
               dist=cfg%dx(i)
               if (d.eq.1) m=mx(i+1,j,k)
               if (d.eq.2) m=-mx(i,j,k)
            else
               dist=cfg%dy(j)
               if (d.eq.3) m=my(i,j+1,k)
               if (d.eq.4) m=-my(i,j,k)
            end if
            exterior=ii.lt.cfg%imin.or.ii.gt.cfg%imax.or.jj.lt.cfg%jmin.or.jj.gt.cfg%jmax
            if (exterior) then
               Tn=reference_T
               c=2.0_WP*conductivity(i,j,k)/dist**2
            else
               Tn=temperature(ii,jj,k)
               c=face_conductivity(conductivity(i,j,k),conductivity(ii,jj,k))/dist**2
            end if
            thermal_into=c*(Tn-temperature(i,j,k))
            r=r+thermal_into
            hs%opr(1,i,j,k)=hs%opr(1,i,j,k)+c
            if (.not.exterior) hs%opr(slot,i,j,k)=hs%opr(slot,i,j,k)-c
            if (exterior) local(4)=local(4)-thermal_into*volume
            otherfluid=cfg%VF(ii,jj,k).gt.0.0_WP
            if (.not.exterior.and.(isfluid.neqv.otherfluid)) then
               if (isfluid) local(5)=local(5)+thermal_into*volume
               if (.not.isfluid) local(6)=local(6)+thermal_into*volume
            end if
            if (.not.isfluid.or..not.otherfluid) then
               if (abs(m).gt.1e-12_WP) call die('Advective mass crosses solid/liquid boundary')
               cycle
            end if
            if (m.ge.0.0_WP) then
               hup=enthalpy(i,j,k); slope=dh(i,j,k)
               hs%opr(1,i,j,k)=hs%opr(1,i,j,k)+m*slope/dist
            else
               if (exterior) then
                  hup=0.0_WP; slope=0.0_WP
               else
                  hup=enthalpy(ii,jj,k); slope=dh(ii,jj,k)
                  hs%opr(slot,i,j,k)=hs%opr(slot,i,j,k)+m*slope/dist
               end if
            end if
            adv_out=m*hup/dist
            r=r-adv_out
            if (exterior) local(4)=local(4)+adv_out*volume
         end do
         hs%rhs(i,j,k)=r
         local(7)=max(local(7),abs(r)/source_scale)
      end do; end do; end do
      call MPI_ALLREDUCE(local(1:6),diag(1:6),6,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(local(7),diag(7),1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
   end subroutine

   subroutine solve_energy(dt,mx,my,iterations)
      real(WP), intent(in) :: dt,mx(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),my(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
      integer, intent(out) :: iterations
      real(WP) :: diag(7),delta,global_delta
      integer :: ierr,it
      do it=1,12
         call update_properties()
         call build_system(dt,mx,my,diag)
         if (diag(7).lt.1e-8_WP) then
            iterations=it-1
            return
         end if
         call hs%setup()
         hs%sol=0.0_WP
         call hs%solve()
         if (hs%rerr.gt.1e-9_WP) call die('Thermal linear solve failed')
         temperature=temperature+hs%sol
      end do
      print*, 'THERMAL_NONLINEAR_RESIDUAL',diag(7)
      call die('Thermal nonlinear solve did not converge')
   end subroutine

   subroutine energy_residual(dt,mx,my,res)
      real(WP), intent(in) :: dt,mx(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),my(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
      real(WP), intent(out) :: res
      real(WP) :: diag(7)
      call build_system(dt,mx,my,diag)
      res=diag(7)
   end subroutine

   subroutine energy_diagnostics(dt,mx,my,diag,budget)
      real(WP), intent(in) :: dt,mx(cfg%imino_:,cfg%jmino_:,cfg%kmino_:),my(cfg%imino_:,cfg%jmino_:,cfg%kmino_:)
      real(WP), intent(out) :: diag(7),budget
      real(WP) :: oldE,global_old
      integer :: i,j,k,ierr
      call build_system(dt,mx,my,diag)
      oldE=0.0_WP
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         oldE=oldE+stored_old(i,j,k)*cfg%vol(i,j,k)*physical_span/cfg%zL
      end do; end do; end do
      call MPI_ALLREDUCE(oldE,global_old,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      budget=((diag(1)+diag(2)-global_old)/dt+diag(4)-diag(3))/max(power,1.0_WP)
   end subroutine
end module
