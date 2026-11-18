! STL-derived continuous function for LICA. No Python dependency at runtime.
! Initialize once with load_sdf BEFORE any OpenMP region or geometry query.
module mod_stl_sdf
  use iso_fortran_env, only: int32, int64, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: load_sdf, evaluate_sdf, sdf_loaded
  logical, save :: sdf_loaded = .false.
  integer(int32), save :: dims(3)
  real(real64), save :: origin(3), spacing(3), upper(3)
  real(real64), allocatable, save :: phi(:,:,:)
contains
  subroutine load_sdf(filename)
    character(len=*), intent(in) :: filename
    character(len=8) :: magic
    integer :: unit, ios, stat
    integer(int64) :: file_bytes, count_values
    open(newunit=unit, file=filename, access='stream', form='unformatted', &
         status='old', action='read', convert='little_endian', iostat=ios)
    if (ios /= 0) error stop 'Cannot open SDF grid'
    inquire(unit=unit, size=file_bytes)
    read(unit,iostat=ios) magic, dims, origin, spacing
    if (ios /= 0 .or. magic /= 'LICASDF1') error stop 'Invalid SDF header'
    if (any(dims < 2) .or. any(dims > 100000)) error stop 'Invalid SDF dimensions'
    count_values = product(int(dims,int64))
    if (count_values > 100000000_int64) error stop 'SDF grid exceeds memory limit'
    if (file_bytes /= 68_int64 + 8_int64*count_values) error stop 'SDF size mismatch'
    if (.not.all(ieee_is_finite(origin)) .or. .not.all(ieee_is_finite(spacing))) error stop 'Nonfinite SDF coordinates'
    if (any(spacing <= 0._real64)) error stop 'Invalid SDF spacing'
    if (allocated(phi)) deallocate(phi)
    allocate(phi(dims(1),dims(2),dims(3)), stat=stat)
    if (stat /= 0) error stop 'SDF allocation failed'
    read(unit,iostat=ios) phi
    close(unit)
    if (ios /= 0 .or. .not.all(ieee_is_finite(phi))) error stop 'Invalid SDF payload'
    if (minval(phi(1,:,:)) <= 0._real64 .or. minval(phi(dims(1),:,:)) <= 0._real64 .or. &
        minval(phi(:,1,:)) <= 0._real64 .or. minval(phi(:,dims(2),:)) <= 0._real64 .or. &
        minval(phi(:,:,1)) <= 0._real64 .or. minval(phi(:,:,dims(3))) <= 0._real64) error stop 'SDF boundary not outside solid'
    upper = origin + spacing*real(dims-1,real64)
    sdf_loaded = .true.
  end subroutine

  function evaluate_sdf(x,y,z) result(value)
    real(real64), intent(in) :: x,y,z
    real(real64) :: value, p(3), pc(3), q(3), w(3), weight
    integer :: i(3), a,b,c
    if (.not.sdf_loaded) error stop 'Call load_sdf before FUNCBODY'
    p = [x,y,z]
    if (.not.all(ieee_is_finite(p))) error stop 'Nonfinite SDF query'
    pc = min(max(p,origin),upper)
    q = (pc-origin)/spacing
    i = min(int(floor(q)),int(dims)-2)
    w = q-real(i,real64)
    value = 0._real64
    do c=0,1
      do b=0,1
        do a=0,1
          weight = merge(w(1),1._real64-w(1),a==1)* &
                   merge(w(2),1._real64-w(2),b==1)* &
                   merge(w(3),1._real64-w(3),c==1)
          value = value + weight*phi(i(1)+a+1,i(2)+b+1,i(3)+c+1)
        end do
      end do
    end do
    value = value + sqrt(sum((p-pc)**2))
  end function
end module
