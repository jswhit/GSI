module readozobs
!$$$  module documentation block
!
! module: readozobs                    read ozone data from diag_sbuv2* files.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: read ozone data from diag_sbuv2* files written out
!  by GSI forward operator code.
!
! Public Subroutines:
!  get_num_ozobs: determine the number of observations to read.
!  get_ozobs_data: read the data and calculate H(x) for ensemble members.
!  write_ozvobs_data: output diag file with spread
!   
! Public Variables: None
!
! program history log:
!   2009-02-23  Initial version.
!   2016-11-29 shlyaeva - updated read routine to calculate linearized H(x)
!                         added write_ozvobs_data to output ensemble spread
!   2017-12-13  shlyaeva - added netcdf diag read/write capability
!
! attributes:
!   language: f95
!
!$$$

use kinds, only: r_single,i_kind,r_kind,r_double
use params, only: nsats_oz,sattypes_oz,npefiles,modelspace_vloc
use constants, only: deg2rad, zero
implicit none

private
public :: get_num_ozobs, get_ozobs_data, write_ozobs_data

contains

! get number of observations from netcdf file
subroutine get_num_ozobs(obspath,datestring,num_obs_tot,num_obs_totdiag,id)
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_dim
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close
  implicit none

    character(len=500), intent(in)  :: obspath
    character(len=10),  intent(in)  :: datestring
    integer(i_kind),    intent(out) :: num_obs_tot, num_obs_totdiag
    character(len=8),   intent(in)  :: id

    character(len=500) obsfile
    character(len=4) pe_name
    real(r_kind) :: errorlimit,errorlimit2
    integer(i_kind) iunit
    integer(i_kind) :: i, nsat, ipe, nobs_curr
    integer(i_kind):: nread,nkeep
    logical :: fexist

    real(r_single),  allocatable, dimension (:) :: Pressure
    integer(i_kind), allocatable, dimension (:) :: Analysis_Use_Flag
    real(r_single),  allocatable, dimension (:) :: Errinv
    real(r_single),  allocatable, dimension (:) :: Observation

    num_obs_tot = 0
    num_obs_totdiag = 0
!   make consistent with screenobs
    errorlimit=1._r_kind/sqrt(1.e9_r_kind)
    errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)
    do nsat=1,nsats_oz
        nread = 0
        nkeep = 0
        peloop: do ipe=0,npefiles
           write(pe_name,'(i4.4)') ipe
           if (npefiles .eq. 0) then
               ! read diag file (concatenated pe* files)
               obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))//'.nc4'
               inquire(file=obsfile,exist=fexist)
               if (.not. fexist .or. datestring .eq. '0000000000') &
               obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))//'.nc4'
           else ! read raw, unconcatenated pe* files.
               obsfile =&
               trim(adjustl(obspath))//'gsitmp_'//trim(adjustl(id))//'/pe'//pe_name//'.'//trim(sattypes_oz(nsat))//'_01.nc4'
           endif
           inquire(file=obsfile,exist=fexist)
           if (.not. fexist) cycle peloop

           call nc_diag_read_init(obsfile, iunit)

           nobs_curr = nc_diag_read_get_dim(iunit,'nobs')

           if (nobs_curr <= 0) then
              call nc_diag_read_close(obsfile)
              cycle peloop
           endif

           allocate(Pressure(nobs_curr), Analysis_Use_Flag(nobs_curr),     &
                   Errinv(nobs_curr), Observation(nobs_curr))

           call nc_diag_read_get_var(iunit, 'Reference_Pressure', Pressure)
           call nc_diag_read_get_var(iunit, 'Analysis_Use_Flag', Analysis_Use_Flag)
           call nc_diag_read_get_var(iunit, 'Inverse_Observation_Error', Errinv)
           call nc_diag_read_get_var(iunit, 'Observation', Observation)

           call nc_diag_read_close(obsfile)

           num_obs_totdiag = num_obs_totdiag + nobs_curr
           nread = nread + nobs_curr
           do i = 1, nobs_curr
             if (Analysis_Use_Flag(i) < 0) cycle
             if (.not. modelspace_vloc .and. (Pressure(i) <= 0.001 .or. &
                 Pressure(i) > 1200._r_kind)) cycle
             if (Errinv(i) <= errorlimit .or.  &
                 Errinv(i) >= errorlimit2 .or.  &
                 abs(Observation(i)) > 1.e9_r_kind) cycle
              nkeep = nkeep + 1
              num_obs_tot = num_obs_tot + 1
           end do
           if (ipe .eq. npefiles) then
              write(6,100) nsat,trim(sattypes_oz(nsat)),nread,nkeep,num_obs_tot
100           format(2x,i3,2x,a20,2x,'nread= ',i9,2x,'nkeep=',i9,2x,'num_obs_tot= ',i9)
           endif
           deallocate(Pressure, Analysis_Use_Flag, Errinv, Observation)
        enddo peloop ! ipe
    enddo ! satellite
    print *,num_obs_tot,' ozone obs'
    print *,num_obs_totdiag, ' total ozone obs in diag file'
end subroutine get_num_ozobs

! read ozone observation data from netcdf file
subroutine get_ozobs_data(obspath, datestring, nobs_max, nobs_maxdiag, hx_mean, hx_mean_nobc, hx, hx_modens, x_obs, x_err, &
           x_lon, x_lat, x_press, x_time, x_code, x_errorig, x_type, x_used, id, nanal, nmem)
  use nc_diag_read_mod, only: nc_diag_read_get_var
  use nc_diag_read_mod, only: nc_diag_read_get_dim, nc_diag_read_get_global_attr
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_close

  use sparsearr,only:sparr, sparr2, readarray, new, delete, assignment(=), init_raggedarr, raggedarr
  use params,only: nanals, lobsdiag_forenkf, neigv, vlocal_evecs
  use statevec, only: state_d
  use mpisetup, only: mpi_wtime, nproc
  use observer_enkf, only: calc_linhx, calc_linhx_modens, setup_linhx
  implicit none

  character*500, intent(in) :: obspath
  character*10, intent(in)  :: datestring

  integer(i_kind), intent(in) :: nobs_max, nobs_maxdiag
  real(r_single), dimension(nobs_max), intent(out)      :: hx_mean, hx_mean_nobc, hx
  ! hx_modens holds modulated ensemble in ob space (zero size and not referenced if neigv=0)
  real(r_single), dimension(neigv,nobs_max), intent(out) :: hx_modens
  real(r_single), dimension(nobs_max), intent(out)      :: x_obs
  real(r_single), dimension(nobs_max), intent(out)      :: x_err, x_errorig
  real(r_single), dimension(nobs_max), intent(out)      :: x_lon, x_lat
  real(r_single), dimension(nobs_max), intent(out)      :: x_press, x_time
  integer(i_kind), dimension(nobs_max), intent(out)     :: x_code
  character(len=20), dimension(nobs_max), intent(out)   :: x_type
  integer(i_kind), dimension(nobs_maxdiag), intent(out) :: x_used

  character(len=8), intent(in) :: id
  integer(i_kind), intent(in)  :: nanal, nmem

  character*500    :: obsfile, obsfile2
  character(len=8) :: id2
  character(len=4) :: pe_name

  integer(i_kind) :: nobs_curr, nob, nobdiag, i, nsat, ipe, nnz, nind, nprof
  integer(i_kind) :: iunit, iunit2

  real(r_double) t1,t2,tsum

  type(sparr2)  :: dhx_dx_read
  type(sparr)   :: dhx_dx
  type(raggedarr) :: hxpert

  real(r_single),  allocatable, dimension (:) :: Latitude, Longitude, Pressure, Time
  integer(i_kind), allocatable, dimension (:) :: Analysis_Use_Flag
  real(r_single),  allocatable, dimension (:) :: Errinv
  real(r_single),  allocatable, dimension (:) :: Observation
  real(r_single),  allocatable, dimension (:) :: Obs_Minus_Forecast_adjusted, Obs_Minus_Forecast_adjusted2
  real(r_single),  allocatable, dimension (:) :: Obs_Minus_Forecast_unadjusted
  integer(i_kind), allocatable, dimension (:,:) :: Observation_Operator_Jacobian_stind
  integer(i_kind), allocatable, dimension (:,:) :: Observation_Operator_Jacobian_endind
  real(r_single), allocatable, dimension (:,:) :: Observation_Operator_Jacobian_val

  logical fexist
  logical twofiles, fexist2
  real(r_kind) :: errorlimit,errorlimit2

  integer(i_kind) :: ix, iy, it, ixp, iyp, itp
  real(r_kind) :: delx, dely, delxp, delyp, delt, deltp
  real(r_single) :: rlat,rlon,rtim,rlat_prev,rlon_prev,rtim_prev,eps

! make consistent with screenobs
  errorlimit=1._r_kind/sqrt(1.e9_r_kind)
  errorlimit2=1._r_kind/sqrt(1.e-6_r_kind)
  eps = 1.e-3

  twofiles = (.not. lobsdiag_forenkf) .and. (nanal <= nanals)
  id2 = 'ensmean'
  if (nanal <= nanals) then
     write(id2,'(a3,(i3.3))') 'mem',nanal
  endif

  tsum = 0
  nob = 0
  rlat_prev = -1.e30; rlon_prev=-1.e30; rtim_prev = -1.e30
  nobdiag = 0
  x_used = 0
  nprof = 0

  hx = zero

  do nsat=1,nsats_oz
      peloop: do ipe=0,npefiles
         write(pe_name,'(i4.4)') ipe
         if (npefiles .eq. 0) then
             ! read diag file (concatenated pe* files)
             obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))//'.nc4'
             inquire(file=obsfile,exist=fexist)
             if (.not. fexist .or. datestring .eq. '0000000000') &
             obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))//'.nc4'
         else ! read raw, unconcatenated pe* files.
             obsfile =&
             trim(adjustl(obspath))//'gsitmp_'//trim(adjustl(id))//'/pe'//pe_name//'.'//trim(sattypes_oz(nsat))//'_01.nc4'
         endif
         inquire(file=obsfile,exist=fexist)
         if (.not. fexist) then
            cycle peloop
         endif

         call nc_diag_read_init(obsfile, iunit)

         nobs_curr = nc_diag_read_get_dim(iunit,'nobs')

         if (nobs_curr <= 0) then
            call nc_diag_read_close(obsfile)
            cycle peloop
         endif
 
         allocate(Latitude(nobs_curr), Longitude(nobs_curr), Time(nobs_curr), Pressure(nobs_curr), &
                  Analysis_Use_Flag(nobs_curr), Errinv(nobs_curr), Observation(nobs_curr),         &
                  Obs_Minus_Forecast_adjusted(nobs_curr), Obs_Minus_Forecast_unadjusted(nobs_curr))
         call nc_diag_read_get_var(iunit, 'Latitude', Latitude)
         call nc_diag_read_get_var(iunit, 'Longitude', Longitude)
         call nc_diag_read_get_var(iunit, 'Time', Time)
         call nc_diag_read_get_var(iunit, 'Reference_Pressure', Pressure)
         call nc_diag_read_get_var(iunit, 'Analysis_Use_Flag', Analysis_Use_Flag)
         call nc_diag_read_get_var(iunit, 'Inverse_Observation_Error', Errinv)
         call nc_diag_read_get_var(iunit, 'Observation', Observation)
         call nc_diag_read_get_var(iunit, 'Obs_Minus_Forecast_adjusted', Obs_Minus_Forecast_adjusted)
         call nc_diag_read_get_var(iunit, 'Obs_Minus_Forecast_unadjusted', Obs_Minus_Forecast_unadjusted)

         if (lobsdiag_forenkf) then
            call nc_diag_read_get_global_attr(iunit, "jac_nnz", nnz)
            call nc_diag_read_get_global_attr(iunit, "jac_nind", nind)
            allocate(Observation_Operator_Jacobian_stind(nind, nobs_curr))
            allocate(Observation_Operator_Jacobian_endind(nind, nobs_curr))
            allocate(Observation_Operator_Jacobian_val(nnz, nobs_curr))
            call nc_diag_read_get_var(iunit,'Observation_Operator_Jacobian_stind', Observation_Operator_Jacobian_stind)
            call nc_diag_read_get_var(iunit,'Observation_Operator_Jacobian_endind', Observation_Operator_Jacobian_endind)
            call nc_diag_read_get_var(iunit,'Observation_Operator_Jacobian_val', Observation_Operator_Jacobian_val)
         endif

         call nc_diag_read_close(obsfile)


         if(twofiles)then
           if (npefiles .eq. 0) then
               ! read diag file (concatenated pe* files)
               obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id2))//'.nc4'
               inquire(file=obsfile2,exist=fexist2)
               if (.not. fexist2 .or. datestring .eq. '0000000000') &
               obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id2))//'.nc4'
           else ! read raw, unconcatenated pe* files.
               obsfile2 =&
               trim(adjustl(obspath))//'gsitmp_'//trim(adjustl(id2))//'/pe'//pe_name//'.'//trim(sattypes_oz(nsat))//'_01.nc4'
           endif

           call nc_diag_read_init(obsfile2, iunit2)

           allocate(Obs_Minus_Forecast_adjusted2(nobs_curr))
           call nc_diag_read_get_var(iunit2, 'Obs_Minus_Forecast_adjusted', Obs_Minus_Forecast_adjusted2)

           call nc_diag_read_close(obsfile2)
 
        end if

        do i = 1, nobs_curr
           nobdiag = nobdiag + 1
           if (Analysis_Use_Flag(i) < 0) cycle
           if (.not. modelspace_vloc .and. (Pressure(i) <= 0.001 .or. &
               Pressure(i) > 1200._r_kind)) cycle

           if (Errinv(i) <= errorlimit .or. Errinv(i) >= errorlimit2 .or.  &
               abs(Observation(i)) > 1.e9_r_kind) cycle
           nob = nob + 1
           x_used(nobdiag) = 1
           x_code(nob) = 700  ! made up code 
           x_lat(nob) = Latitude(i)
           x_lon(nob) = Longitude(i)
           x_press(nob) = Pressure(i)
           x_time(nob) = Time(i)
           x_err(nob) = (1./Errinv(i))**2
           x_errorig(nob) = x_err(nob)
           x_obs(nob) =  Observation(i)
           hx_mean(nob) = Observation(i) - Obs_Minus_Forecast_adjusted(i)
           hx_mean_nobc(nob) = Observation(i) - Obs_Minus_Forecast_unadjusted(i)
           x_type(nob) = ' oz                 '
           if (nanal <= nanals) then
             ! read full Hx from diag file
             if (.not. lobsdiag_forenkf) then
                hx(nob) = Observation(i) - Obs_Minus_Forecast_adjusted2(i)
               ! run linearized Hx
             else
                call new(dhx_dx_read, nnz, nind)
                dhx_dx_read%st_ind = Observation_Operator_Jacobian_stind(:,i)
                dhx_dx_read%end_ind = Observation_Operator_Jacobian_endind(:,i)
                dhx_dx_read%val = Observation_Operator_Jacobian_val(:,i)
                dhx_dx = dhx_dx_read
                t1 = mpi_wtime()
                rlat = x_lat(nob)*deg2rad
                rlon = x_lon(nob)*deg2rad
                rtim = x_time(nob)
                if (nob > 1) then
                   rlat_prev = x_lat(nob-1)*deg2rad
                   rlon_prev = x_lon(nob-1)*deg2rad
                   rtim_prev = x_time(nob-1)
                endif
                if (abs(rlat-rlat_prev) > eps .or. &
                   abs(rlon-rlon_prev) > eps .or. &
                   abs(rtim-rtim_prev) > eps) then
                   call setup_linhx(rlat,rlon,rtim,              &
                                 ix, delx, ixp, delxp, iy, dely,  &
                                 iyp, delyp, it, delt, itp, deltp)
                else
                   nprof = nprof + 1
                endif
                call init_raggedarr(hxpert, dhx_dx%nnz)
                call calc_linhx(hx_mean(nob), state_d(:,:,:,nmem),       &
                                dhx_dx, hxpert, hx(nob),                 &
                                ix, delx, ixp, delxp, iy, dely,   &
                                iyp, delyp, it, delt, itp, deltp)
                ! compute modulated ensemble in obs space
                if (neigv>0) call calc_linhx_modens(hx_mean(nob),dhx_dx,hxpert,hx_modens(:,nob),vlocal_evecs)
                t2 = mpi_wtime()
                tsum = tsum + t2-t1

                call delete(dhx_dx)
                call delete(dhx_dx_read)
             endif
           endif

         end do ! k

         deallocate(Latitude, Longitude, Time, Pressure, Analysis_Use_Flag, Errinv, &
                    Observation, Obs_Minus_Forecast_adjusted,                       &
                    Obs_Minus_Forecast_unadjusted)
         if (twofiles) then
            deallocate(Obs_Minus_Forecast_adjusted2)
         endif
         if (lobsdiag_forenkf) then
            deallocate(Observation_Operator_Jacobian_stind)
            deallocate(Observation_Operator_Jacobian_endind)
            deallocate(Observation_Operator_Jacobian_val)
         endif
      enddo peloop ! ipe
  enddo ! satellite
  if (nanal == nanals) print *,'oz obs profiles, total obs',nprof,nob
  if (nanal == nanals .and. lobsdiag_forenkf) print *, 'time in calc_linhx for oz obs on proc',nproc,' =',tsum

  if (nob /= nobs_max) then
      print *,'number of obs not what expected in get_ozobs_data',nob,nobs_max
      call stop2(93)
  end if
  if (nobdiag /= nobs_maxdiag) then
      print *,'number of total diag obs not what expected in get_ozobs_data',nobdiag,nobs_maxdiag
      call stop2(93)
  end if


end subroutine get_ozobs_data

! writing spread diagnostics to netcdf file
subroutine write_ozobs_data(obspath, datestring, nobs_max, nobs_maxdiag, &
                                 x_fit, x_sprd, x_used, id, gesid)
  use netcdf, only: nf90_inq_dimid, nf90_open, nf90_close, NF90_NETCDF4, &
                    nf90_inquire_dimension, NF90_WRITE, NF90_NOWRITE, nf90_create, nf90_def_dim
  use ncdw_climsg, only: nclayer_check

  use constants, only: r_missing
  implicit none

  character*500,   intent(in) :: obspath
  character*10,    intent(in) :: datestring
  integer(i_kind), intent(in) :: nobs_max, nobs_maxdiag
  real(r_single),  dimension(nobs_max),     intent(in) :: x_fit, x_sprd
  integer(i_kind), dimension(nobs_maxdiag), intent(in) :: x_used
  character(len=8), intent(in) :: id, gesid


  character*500 obsfile, obsfile2
  character(len=4) pe_name

  integer(i_kind) :: iunit, nobsid
  integer(i_kind) :: nob, nobdiag, nobs, ipe, i, nsat
  integer(i_kind), dimension(:), allocatable :: enkf_use_flag
  real(r_single),  dimension(:), allocatable :: enkf_fit, enkf_sprd
  logical :: fexist

  nob  = 0
  nobdiag = 0

  do nsat=1,nsats_oz
      peloop: do ipe=0,npefiles
         write(pe_name,'(i4.4)') ipe
         if (npefiles .eq. 0) then
            ! diag file (concatenated pe* files)
            obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))//'.nc4'
            obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//datestring//'_'//trim(adjustl(id))//'_spread.nc4'
            inquire(file=obsfile,exist=fexist)
            if (.not. fexist .or. datestring .eq. '0000000000') then
               obsfile = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))//'.nc4'
               obsfile2 = trim(adjustl(obspath))//"diag_"//trim(sattypes_oz(nsat))//"_ges."//trim(adjustl(id))//'_spread.nc4'
            endif
         else ! raw, unconcatenated pe* files.
            obsfile = trim(adjustl(obspath))//'gsitmp_'//trim(adjustl(id))//'/pe'//pe_name//'.'//trim(sattypes_oz(nsat))//'_01'//'.nc4'
            obsfile2 = trim(adjustl(obspath))//'gsitmp_'//trim(adjustl(id))//'/pe'//pe_name//'.'//trim(sattypes_oz(nsat))//'_01'//'_spread.nc4'
         endif

         inquire(file=obsfile,exist=fexist)
         if (.not. fexist) cycle peloop

         call nclayer_check(nf90_open(obsfile, NF90_NOWRITE, iunit))
         call nclayer_check(nf90_inq_dimid(iunit, "nobs", nobsid))
         call nclayer_check(nf90_inquire_dimension(iunit, nobsid, len = nobs))
         call nclayer_check(nf90_close(iunit))

         if (nobs <= 0) cycle peloop

         allocate(enkf_use_flag(nobs), enkf_fit(nobs), enkf_sprd(nobs))
         enkf_use_flag = -1
         enkf_fit = r_missing
         enkf_sprd = r_missing

         do i = 1, nobs
           nobdiag = nobdiag + 1
           ! skip if not used in EnKF
           if (x_used(nobdiag) == 1) then
              ! update if it is used in EnKF
              nob = nob + 1
              enkf_use_flag(i) = 1
              enkf_fit(i)  = x_fit(nob)
              enkf_sprd(i) = x_sprd(nob)
           endif
         enddo

         inquire(file=obsfile2,exist=fexist)
         if (.not. fexist) then
            call nclayer_check(nf90_create(trim(obsfile2), NF90_NETCDF4, &
                               iunit))
            call nclayer_check(nf90_def_dim(iunit, "nobs", nobs, nobsid))
         else
            call nclayer_check(nf90_open(obsfile2, NF90_WRITE, iunit))
            call nclayer_check(nf90_inq_dimid(iunit, "nobs", nobsid))
         endif

         call write_ncvar_int(iunit, nobsid, "EnKF_use_flag", enkf_use_flag)
         call write_ncvar_single(iunit, nobsid, "EnKF_fit_"//trim(gesid), enkf_fit)
         call write_ncvar_single(iunit, nobsid, "EnKF_spread_"//trim(gesid), enkf_sprd)
 
         call nclayer_check(nf90_close(iunit))
 
         deallocate(enkf_use_flag, enkf_fit, enkf_sprd)

     enddo peloop ! ipe loop
  enddo

  if (nob .ne. nobs_max) then
      print *,'number of obs not what expected in write_ozobs_data',nob,nobs_max
      call stop2(94)
  end if
  if (nobdiag /= nobs_maxdiag) then
      print *,'number of total obs in diag not what expected in write_ozobs_data',nobdiag, nobs_maxdiag
      call stop2(94)
  endif


  contains
  subroutine write_ncvar_single(iunit, dimid, varname, field)
    use netcdf, only: nf90_def_var, nf90_put_var, nf90_inq_varid,  &
                      nf90_def_var_deflate,NF90_FLOAT, NF90_ENOTVAR
    use ncdw_climsg, only: nclayer_check
    use ncdw_types, only: NLAYER_COMPRESSION
    implicit none
    integer(i_kind), intent(in)  :: iunit, dimid
    character(*), intent(in)     :: varname
    real(r_single), dimension(:), allocatable :: field

    integer :: ierr, varid

    ierr = nf90_inq_varid(iunit, varname, varid)
    if (ierr == NF90_ENOTVAR) then
       call nclayer_check(nf90_def_var(iunit, varname, NF90_FLOAT, dimid, varid))
       call nclayer_check(nf90_def_var_deflate(iunit, varid, 1, 1, int(NLAYER_COMPRESSION)))
    endif
    call nclayer_check(nf90_put_var(iunit, varid, field))
  end subroutine write_ncvar_single

  subroutine write_ncvar_int(iunit, dimid, varname, field)
    use netcdf, only: nf90_def_var, nf90_put_var, nf90_inq_varid,  &
                      nf90_def_var_deflate,NF90_INT, NF90_ENOTVAR
    use ncdw_climsg, only: nclayer_check
    use ncdw_types, only: NLAYER_COMPRESSION
    implicit none
    integer(i_kind), intent(in)  :: iunit, dimid
    character(*), intent(in)     :: varname
    integer(i_kind), dimension(:), allocatable :: field

    integer :: ierr, varid

    ierr = nf90_inq_varid(iunit, varname, varid)
    if (ierr == NF90_ENOTVAR) then
       call nclayer_check(nf90_def_var(iunit, varname, NF90_INT, dimid, varid))
       call nclayer_check(nf90_def_var_deflate(iunit, varid, 1, 1, int(NLAYER_COMPRESSION)))
    endif
    call nclayer_check(nf90_put_var(iunit, varid, field))
  end subroutine write_ncvar_int


end subroutine write_ozobs_data


end module readozobs
