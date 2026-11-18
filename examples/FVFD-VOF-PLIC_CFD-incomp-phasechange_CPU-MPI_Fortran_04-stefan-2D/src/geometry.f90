module geometry
 use precision,only:WP
 use config_class,only:config
 implicit none
 type(config),public::cfg
 real(WP),public::block_width,block_height,physical_span,Lz
contains
 subroutine geometry_init
  use sgrid_class,only:sgrid,cartesian
  use parallel,only:group
  use param,only:param_read
  type(sgrid)::g
  integer::nx,ny,i,j,k,part(3)
  real(WP)::Lx,Ly,dx
  real(WP),allocatable::x(:),y(:),z(:)
  call param_read('nx',nx);call param_read('ny',ny)
  call param_read('Lx',Lx);call param_read('Ly',Ly)
  call param_read('Block width',block_width);call param_read('Block height',block_height)
  call param_read('Physical span',physical_span);call param_read('Partition',part,short='p')
  dx=Lx/nx;Lz=dx
  allocate(x(nx+1),y(ny+1),z(2))
  do i=1,nx+1;x(i)=(i-1)*dx-Lx/2;end do
  do j=1,ny+1;y(j)=(j-1)*Ly/ny;end do
  z=[-dx/2,dx/2]
  g=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='heatedBubbleOpen')
  cfg=config(grp=group,decomp=part,grid=g)
  do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
   if(abs(cfg%xm(i)).lt.block_width/2.and.cfg%ym(j).lt.block_height)cfg%VF(i,j,k)=0.0_WP
  end do;end do;end do
  call cfg%sync(cfg%VF)
 end subroutine
end module
