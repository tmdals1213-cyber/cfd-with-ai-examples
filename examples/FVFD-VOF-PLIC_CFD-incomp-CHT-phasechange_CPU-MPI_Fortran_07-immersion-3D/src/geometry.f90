module geometry
 use precision,only:WP
 use config_class,only:config
 implicit none
 type(config),public::cfg
 real(WP),public::block_width,block_height,block_depth,physical_span,Lz
 real(WP),public::domain_width,domain_height,domain_depth,cold_width,cold_depth,cold_thickness,port_width,pool_height,cold_film_thickness
contains
 subroutine geometry_init
  use sgrid_class,only:sgrid,cartesian
  use parallel,only:group
  use param,only:param_read
  type(sgrid)::g
  integer::nx,ny,nz,i,j,k,part(3)
  character(len=64)::boundary_mode
  real(WP)::Lx,Ly
  real(WP),allocatable::x(:),y(:),z(:)
  call param_read('nx',nx);call param_read('ny',ny);call param_read('nz',nz)
  call param_read('Lx',Lx);call param_read('Ly',Ly);call param_read('Lz',Lz)
  call param_read('Block width',block_width);call param_read('Block height',block_height);call param_read('Block depth',block_depth)
  call param_read('Partition',part,short='p');physical_span=Lz
  domain_width=Lx;domain_height=Ly;domain_depth=Lz
  call param_read('Cold plate width',cold_width);call param_read('Cold plate depth',cold_depth);call param_read('Cold plate thickness',cold_thickness)
  call param_read('Boundary mode',boundary_mode)
  if(trim(boundary_mode).ne.'volume_balance_3d_v1')error stop 'Wrong3D boundary mode'
  call param_read('Pressure port width',port_width);call param_read('Pool height',pool_height);call param_read('Initial cold film thickness',cold_film_thickness)
  if(min(nx,ny,nz)<4.or.min(Lx,Ly,Lz)<=0)error stop 'Invalid3D grid'
  if(cold_width+2*port_width>=Lx.or.cold_depth+2*port_width>=Lz.or.pool_height<=block_height.or.pool_height>=Ly-cold_thickness-cold_film_thickness)error stop 'Invalid3D pool/cold plate/port geometry'
  allocate(x(nx+1),y(ny+1),z(nz+1))
  do i=1,nx+1;x(i)=(i-1)*Lx/nx-Lx/2;enddo
  do j=1,ny+1;y(j)=(j-1)*Ly/ny;enddo
  do k=1,nz+1;z(k)=(k-1)*Lz/nz-Lz/2;enddo
  g=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='chip3DOpen')
  cfg=config(grp=group,decomp=part,grid=g)
  do k=cfg%kmino_,cfg%kmaxo_;do j=cfg%jmino_,cfg%jmaxo_;do i=cfg%imino_,cfg%imaxo_
   if(abs(cfg%xm(i))<block_width/2.and.abs(cfg%zm(k))<block_depth/2.and.cfg%ym(j)<block_height)cfg%VF(i,j,k)=0
   if(abs(cfg%xm(i))<cold_width/2.and.abs(cfg%zm(k))<cold_depth/2.and.cfg%ym(j)>domain_height-cold_thickness)cfg%VF(i,j,k)=0
  enddo;enddo;enddo
  call cfg%sync(cfg%VF)
 end subroutine
 logical function cold_solid(i,j,k)
  integer,intent(in)::i,j,k
  cold_solid=cfg%VF(i,j,k)==0.and.cfg%ym(j)>domain_height-cold_thickness
 end function
 logical function top_port(i,k)
  integer,intent(in)::i,k
  top_port=abs(cfg%xm(i))>domain_width/2-port_width.or.abs(cfg%zm(k))>domain_depth/2-port_width
 end function
end module
