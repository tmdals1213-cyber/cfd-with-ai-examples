module geometry
   use precision, only: WP
   use ibconfig_class, only: ibconfig
   implicit none
   type(ibconfig), public :: cfg
   real(WP), public :: block_width,block_height,physical_span
contains
   subroutine geometry_init
      use sgrid_class, only: sgrid,cartesian
      use ibconfig_class, only: stair
      use parallel, only: group
      use param, only: param_read
      use messager, only: die
      type(sgrid) :: grid
      integer :: nx,ny,i,j,k,part(3)
      real(WP) :: Lx,Ly,dx,qx,qy
      real(WP), allocatable :: x(:),y(:),z(:)
      call param_read('nx',nx); call param_read('ny',ny)
      call param_read('Lx',Lx); call param_read('Ly',Ly)
      call param_read('Block width',block_width); call param_read('Block height',block_height)
      call param_read('Physical span',physical_span)
      dx=Lx/real(nx,WP)
      allocate(x(nx+1),y(ny+1),z(2))
      do i=1,nx+1
         x(i)=real(i-1,WP)*dx-Lx/2
      end do
      do j=1,ny+1
         y(j)=real(j-1,WP)*Ly/real(ny,WP)-0.02_WP
      end do
      z=[-dx/2,dx/2]
      if (minval(abs(x-block_width/2)).gt.1e-12_WP.or.minval(abs(x+block_width/2)).gt.1e-12_WP) &
         call die('Block x walls must coincide with cell faces')
      if (minval(abs(y-block_height/2)).gt.1e-12_WP.or.minval(abs(y+block_height/2)).gt.1e-12_WP) &
         call die('Block y walls must coincide with cell faces')
      grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='waterCHT2D')
      call param_read('Partition',part,short='p')
      cfg=ibconfig(grp=group,decomp=part,grid=grid)
      do k=cfg%kmino_,cfg%kmaxo_; do j=cfg%jmino_,cfg%jmaxo_; do i=cfg%imino_,cfg%imaxo_
         qx=abs(cfg%xm(i))-block_width/2; qy=abs(cfg%ym(j))-block_height/2
         cfg%Gib(i,j,k)=-(sqrt(max(qx,0.0_WP)**2+max(qy,0.0_WP)**2)+min(max(qx,qy),0.0_WP))
      end do; end do; end do
      call cfg%calculate_normal()
      call cfg%calculate_vf(method=stair,allow_zero_vf=.true.)
      call cfg%sync(cfg%VF)
   end subroutine
end module
