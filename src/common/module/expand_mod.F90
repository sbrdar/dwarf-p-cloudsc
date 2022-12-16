! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module expand_mod
  use atlas_module
  use atlas_fieldset_module
  use atlas_functionspace_blockstructuredcolumns_module

  use parkind1 , only : jpim, jprb
  use yomphyder, only : state_type

  use cloudsc_mpi_mod, only : irank, numproc
  use file_io_mod, only: input_initialize, load_scalar, load_array

  use, intrinsic :: iso_c_binding, only : c_int, c_long, c_float, c_double, c_bool

  implicit none

  interface expand
     procedure expand_l1, expand_i1, expand_r1, expand_r2, expand_r3
  end interface expand

  interface load_and_expand
     procedure load_and_expand_l1, load_and_expand_i1
     procedure load_and_expand_r1, load_and_expand_r2, load_and_expand_r3
  end interface load_and_expand

contains

  subroutine get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    integer(kind=jpim), intent(inout) :: start, end, size
    integer(kind=jpim), intent(in) :: nlon, ndim, nlev, ngptot
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim) :: rankstride
    logical :: use_offset = .false.

    if (present(ngptotg)) use_offset = nlon >= ngptotg
    if (use_offset) then
      rankstride = (ngptotg - 1) / numproc + 1
      start = irank * rankstride + 1
    else
      start = 1
    end if
    size = min(nlon, ngptot)
    end = start + size - 1
  end subroutine get_offsets

  subroutine load_and_expand_i1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    integer(kind=jpim), allocatable, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    integer(kind=jpim), allocatable :: buffer(:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_i1

  subroutine load_and_expand_l1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    logical, allocatable, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    logical, allocatable :: buffer(:), rbuf(:)
    integer(kind=jpim) :: start, end, size
    integer(kind=4), allocatable :: tmp(:)

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_l1

  subroutine load_and_expand_r1(name, field, nlon, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), allocatable, intent(inout) :: field(:,:)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nblocks))
    allocate(buffer(size))
    call load_array(name, start, end, size, nlon, buffer)
    call expand(buffer, field, size, nproma, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r1

  subroutine loadvar_atlas(fset, name, nlon, ngptotg)
    ! Load into the local memory buffer and expand to global field
    type(atlas_fieldset), intent(inout) :: fset
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: nlon
    integer(kind=jpim), intent(in), optional :: ngptotg

    integer(kind=jpim) :: start, end, size, nlev, nproma, ngptot, nblocks, ndim
    type(atlas_field) :: field
    real(kind=jprb), allocatable :: buffer_r1(:), buffer_r2(:,:), buffer_r3(:,:,:)
    integer(kind=jpim), allocatable :: buffer_i1(:)
    logical, allocatable :: buffer_l1(:)
    real(c_double), pointer :: field_r1(:,:), field_r2(:,:,:), field_r3(:,:,:,:)
    integer(c_int), pointer :: field_i1(:,:)
    logical, pointer :: field_l1(:,:)
    type(atlas_functionspace_blockstructuredcolumns) :: fspace
    logical :: lfield, rfield, ifield

    field = fset%field(name)
    ndim = field%rank()
    lfield = (name == "LDCUM")
    ifield = (name == "KTYPE")
    rfield = ((.not. lfield) .and. (.not. ifield))

    fspace = field%functionspace()
    nlev = field%levels()
    nproma = fspace%block_size(1)
    ngptot = fspace%size()
    nblocks = fspace%nblks()

    if (ndim == 2) then
      call get_offsets(start, end, size, nlon, 1, 1, ngptot, ngptotg)
      if (rfield) then
        allocate(buffer_r1(size))
        call field%data(field_r1)
        call load_array(name, start, end, size, nlon, buffer_r1)
        call expand(buffer_r1, field_r1, size, nproma, ngptot, nblocks)
        deallocate(buffer_r1)
      else if (lfield) then
        allocate(buffer_l1(size))
        call field%data(field_l1)
        call load_array(name, start, end, size, nlon, buffer_l1)
        call expand(buffer_l1, field_l1, size, nproma, ngptot, nblocks)
        deallocate(buffer_l1)
      else
        allocate(buffer_i1(size))
        call field%data(field_i1)
        call load_array(name, start, end, size, nlon, buffer_i1)
        call expand(buffer_i1, field_i1, size, nproma, ngptot, nblocks)
        deallocate(buffer_i1)
      endif
    else if (ndim == 3) then
      call get_offsets(start, end, size, nlon, 1, nlev, ngptot, ngptotg)
      if (rfield) then
        call field%data(field_r2)
        allocate(buffer_r2(size, nlev))
        call load_array(name, start, end, size, nlon, nlev, buffer_r2)
        call expand(buffer_r2, field_r2, size, nproma, nlev, ngptot, nblocks)
        deallocate(buffer_r2)
      else
         call exit(1)
      endif
    else if (ndim == 4) then
      call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
      if (rfield) then
        call field%data(field_r3)
        allocate(buffer_r3(size, nlev, ndim))
        call load_array(name, start, end, size, nlon, nlev, ndim, buffer_r3)
        call expand(buffer_r3, field_r3, size, nproma, nlev, ndim, ngptot, nblocks)
        deallocate(buffer_r3)
      else
         call exit(1)
      endif
    endif
  end subroutine loadvar_atlas

  subroutine load_and_expand_r2(name, field, nlon, nlev, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), allocatable, intent(inout) :: field(:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, 1, nlev, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nlev, nblocks))
    allocate(buffer(size, nlev))
    call load_array(name, start, end, size, nlon, nlev, buffer)
    call expand(buffer, field, size, nproma, nlev, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r2

  subroutine load_and_expand_r3(name, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    real(kind=JPRB), allocatable, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start, end, size

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nlev, ndim, nblocks))
    allocate(buffer(size, nlev, ndim))
    call load_array(name, start, end, size, nlon, nlev, ndim, buffer)
    call expand(buffer, field, size, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)
  end subroutine load_and_expand_r3

  subroutine load_and_expand_state_atlas(name, state, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    type(state_type), pointer, intent(inout) :: state(:)
    real(kind=JPRB), allocatable, target, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start, end, size

    integer :: b

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    if (.not. allocated(field))  allocate(field(nproma, nlev, 3+ndim, nblocks))
    allocate(buffer(size, nlev, 3+ndim))

    call load_array(name//'_T', start, end, size, nlon, nlev, buffer(:,:,1))
    call load_array(name//'_A', start, end, size, nlon, nlev, buffer(:,:,2))
    call load_array(name//'_Q', start, end, size, nlon, nlev, buffer(:,:,3))
    call load_array(name//'_CLD', start, end, size, nlon, nlev, ndim, buffer(:,:,4:))

    call expand(buffer(:,:,1), field(:,:,1,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,2), field(:,:,2,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field(:,:,3,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,4:), field(:,:,4:,:), size, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B) schedule(runtime)
    do b=1, nblocks
       state(b)%t => field(:,:,1,b)
       state(b)%a => field(:,:,2,b)
       state(b)%q => field(:,:,3,b)
       state(b)%cld => field(:,:,4:3+ndim,b)
    end do
!$OMP end parallel do
  end subroutine load_and_expand_state_atlas

  subroutine load_and_expand_state(name, state, field, nlon, nlev, ndim, nproma, ngptot, nblocks, ngptotg)
    ! Load into the local memory buffer and expand to global field
    character(len=*) :: name
    type(state_type), allocatable, intent(inout) :: state(:)
    real(kind=JPRB), allocatable, target, intent(inout) :: field(:,:,:,:)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer(kind=jpim), intent(in), optional :: ngptotg
    real(kind=jprb), allocatable :: buffer(:,:,:)
    integer(kind=jpim) :: start, end, size

    integer :: b

    call get_offsets(start, end, size, nlon, ndim, nlev, ngptot, ngptotg)
    if (.not. allocated(state))  allocate(state(nblocks))
    if (.not. allocated(field))  allocate(field(nproma, nlev, 3+ndim, nblocks))
    allocate(buffer(size, nlev, 3+ndim))

    call load_array(name//'_T', start, end, size, nlon, nlev, buffer(:,:,1))
    call load_array(name//'_A', start, end, size, nlon, nlev, buffer(:,:,2))
    call load_array(name//'_Q', start, end, size, nlon, nlev, buffer(:,:,3))
    call load_array(name//'_CLD', start, end, size, nlon, nlev, ndim, buffer(:,:,4:))

    call expand(buffer(:,:,1), field(:,:,1,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,2), field(:,:,2,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,3), field(:,:,3,:), size, nproma, nlev, ngptot, nblocks)
    call expand(buffer(:,:,4:), field(:,:,4:,:), size, nproma, nlev, ndim, ngptot, nblocks)
    deallocate(buffer)

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(B) schedule(runtime)
    do b=1, nblocks
       state(b)%t => field(:,:,1,b)
       state(b)%a => field(:,:,2,b)
       state(b)%q => field(:,:,3,b)
       state(b)%cld => field(:,:,4:3+ndim,b)
    end do
!$OMP end parallel do
  end subroutine load_and_expand_state

  subroutine expand_l1(buffer, field, nlon, nproma, ngptot, nblocks)
    logical, intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx-1,nlon)+1
       bend = min(nlon,bidx+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = .FALSE.
    end do
!$omp end parallel do
  end subroutine expand_l1

  subroutine expand_i1(buffer, field, nlon, nproma, ngptot, nblocks)
    integer(kind=jpim), intent(inout) :: buffer(nlon), field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx-1,nlon)+1
       bend = min(nlon,bidx+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0_JPIM
    end do
!$omp end parallel do
  end subroutine expand_i1

  subroutine expand_r1(buffer, field, nlon, nproma, ngptot, nblocks)
    real(kind=jprb), intent(inout) :: buffer(nlon)
    real(kind=jprb), intent(inout) :: field(nproma, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx-1,nlon)+1
       bend = min(nlon,bidx+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,b) = buffer(bidx:bend)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,b) = buffer(bidx:bend)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,b) = 0.0_JPRB
    end do
!$omp end parallel do
  end subroutine expand_r1

  subroutine expand_r2(buffer, field, nlon, nproma, nlev, ngptot, nblocks)
          use omp_lib
    real(kind=jprb), intent(inout) :: buffer(nlon, nlev)
    real(kind=jprb), intent(inout) :: field(nproma, nlev, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx-1,nlon)+1
       bend = min(nlon,bidx+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,:,b) = buffer(bidx:bend,:)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,:,b) = buffer(bidx:bend,:)
       end do

       field(bsize+1:nproma,:,b) = 0.0_JPRB
    end do
!$omp end parallel do

  end subroutine expand_r2

  subroutine expand_r3(buffer, field, nlon, nproma, nlev, ndim, ngptot, nblocks)
    real(kind=jprb), intent(inout) :: buffer(nlon, nlev, ndim)
    real(kind=jprb), intent(inout) :: field(nproma, nlev, ndim, nblocks)
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim, nproma, ngptot, nblocks
    integer :: b, gidx, bsize, fidx, fend, bidx, bend

!$omp parallel do default(shared) private(b, gidx, bsize, fidx, fend, bidx, bend) schedule(runtime)
    do b=1, nblocks
       gidx = (b-1)*nproma + 1  ! Global starting index of the block in the general domain
       bsize = min(nproma, ngptot - gidx + 1)  ! Size of the field block

       ! First read, might not be aligned
       bidx = mod(gidx-1,nlon)+1
       bend = min(nlon,bidx+bsize-1)
       fidx = 1
       fend = bend - bidx + 1
       field(fidx:fend,:,:,b) = buffer(bidx:bend,:,:)

       ! Fill block by looping over buffer
       do while (fend < bsize)
         fidx = fend + 1
         bidx = 1
         bend = min(bsize - fidx+1, nlon)
         fend = fidx + bend - 1
         field(fidx:fend,:,:,b) = buffer(bidx:bend,:,:)
       end do

       ! Zero out the remainder of last block
       field(bsize+1:nproma,:,:,b) = 0.0_JPRB
    end do
!$omp end parallel do
  end subroutine expand_r3

end module expand_mod
