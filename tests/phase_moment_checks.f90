program check_phase_moment
 use iso_fortran_env,only:real64
 use phase_moment_core
 implicit none
 real(real64)::a,r,m(3),c(3),p(3),bad(3),point(3),sink
 integer::s,i,n
 n=0;p=[.1_real64,.2_real64,.3_real64]
 do i=0,6
  a=10._real64**(-i*40);r=a*.125_real64;m=a*p
  call depleting_phase_centroid(m,a,r,c,s)
  if(s/=moment_ok.or.maxval(abs(c-p))>1e-14_real64)error stop 'Scale-independent depletion'
  n=n+1
 enddo
 a=scale(1._real64,-1050);p=[.25_real64,.5_real64,.75_real64];m=a*p
 call depleting_phase_centroid(m,a,a/2,c,s)
 if(s/=moment_ok.or.any(c/=p))error stop 'Subnormal support without reciprocal overflow'
 n=n+1
 call depleting_phase_centroid(m,a,0._real64,c,s)
 if(s/=moment_ok.or.any(0._real64*c/=0._real64))error stop 'Complete depletion moment'
 n=n+1
 call depleting_phase_centroid(m,0._real64,0._real64,c,s)
 if(s==moment_ok)error stop 'Absent support must refuse'
 n=n+1
 a=2.0042764952071856e-14_real64+3.0895066299394900e-18_real64
 r=1.6005796797277583e-19_real64;sink=-2.0045694400733823e-14_real64
 point=[-1.9692524249535072e-3_real64,1.8998696288314791e-2_real64,-2.6247600915200327e-3_real64]
 bad=[-1.2859582640596892_real64,.072893787318207656_real64,.028208716839014523_real64]
 ! Reconstruct advected moment from the captured old point-sink operands.
 m=r*bad-sink*point
 call depleting_phase_centroid(m,a,r,c,s)
 if(s/=moment_ok)error stop 'Captured failure support'
 if(c(1)<-.002_real64.or.c(1)>-.00175_real64.or.c(2)<.01875_real64.or.c(2)>.019_real64.or.c(3)<-.00275_real64.or.c(3)>-.0025_real64)error stop 'Captured cell barycenter support'
 if(maxval(abs(bad))<1._real64)error stop 'Captured original pathology missing'
 print*,'CAPTURED_OLD_POINT_SINK_CENTROID',bad
 print*,'REPAIRED_PHASE_MEAN_CENTROID',c
 n=n+1
 print*,'PHASE_MOMENT_KERNEL_CHECKS_PASS',n
end program
