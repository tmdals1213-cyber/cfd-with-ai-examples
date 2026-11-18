program checks
 use precision,only:WP
 use face_geometry3d
 implicit none
 real(WP)::v(3,4),p(4),b(3),f
 v(:,1)=[0.0_WP,0.0_WP,0.0_WP];v(:,2)=[1.0_WP,0.0_WP,0.0_WP];v(:,3)=[1.0_WP,1.0_WP,0.0_WP];v(:,4)=[0.0_WP,1.0_WP,0.0_WP]
 p=[1.0_WP,0.0_WP,0.0_WP,.5_WP];b=[.1_WP,.5_WP,0.0_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f-.5_WP)>1e-12_WP)error stop 'Half XY face'
 p=[1.0_WP,1.0_WP,0.0_WP,.5_WP];b=[.1_WP,.1_WP,0.0_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f-.125_WP)>1e-12_WP)error stop 'Oblique triangular face'
 p=[0.0_WP,0.0_WP,1.0_WP,1.0_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f-1)>1e-12_WP)error stop 'Full face'
 p=[1.0_WP,0.0_WP,0.0_WP,-1.0_WP];b=[-2.0_WP,0.0_WP,0.0_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f)>1e-12_WP)error stop 'Empty face'
 v(:,1)=[.5_WP,0.0_WP,0.0_WP];v(:,2)=[.5_WP,1.0_WP,0.0_WP];v(:,3)=[.5_WP,1.0_WP,1.0_WP];v(:,4)=[.5_WP,0.0_WP,1.0_WP]
 p=[0.0_WP,0.0_WP,1.0_WP,.25_WP];b=[.5_WP,.5_WP,.1_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f-.25_WP)>1e-12_WP)error stop 'YZ face'
 v(:,1)=[0.0_WP,.5_WP,0.0_WP];v(:,2)=[1.0_WP,.5_WP,0.0_WP];v(:,3)=[1.0_WP,.5_WP,1.0_WP];v(:,4)=[0.0_WP,.5_WP,1.0_WP]
 p=[0.0_WP,0.0_WP,1.0_WP,.75_WP];b=[.5_WP,.5_WP,.5_WP];f=clipped_quad_fraction(v,p,b)
 if(abs(f-.75_WP)>1e-12_WP)error stop 'XZ face'
 print*,'SIX_3D_PLIC_FACE_KERNEL_CHECKS_PASS'
end program
