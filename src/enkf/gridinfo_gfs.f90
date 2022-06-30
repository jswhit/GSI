module gridinfo
!$$$  module documentation block
!
! module: gridinfo                     read horizontal (lons, lats) and
!                                      vertical (pressure) information from
!                                      ensemble mean first guess file.
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract: This module reads gfg_YYYYMMDDHH_fhr06_ensmean, and
! extracts information about the analysis grid, including the
! longitudes and latitudes of the analysis grid points and
! the pressure on each grid point/vertical level.
!
! Public Subroutines:
!   getgridinfo: read latitudes, longitudes, pressures and orography for analysis grid,
!    broadcast to each task. Compute spherical cartesian coordinate values
!    for each analysis horizontal grid point.
!   gridinfo_cleanup: deallocate allocated module variables.
!
! Public Variables:
!   npts: number of analysis grid points in the horizontal (from module params).
!   nlevs: number of analysis vertical levels (from module params).
!   ntrac: number of 'tracer' model state variables (3 for GFS,
!    specific humidity, ozone and cloud condensate).
!   ptop: (real scalar) pressure (hPa) at top model layer interface.
!   lonsgrd(npts): real array of analysis grid longitudes (radians).
!   latsgrd(npts): real array of analysis grid latitudes (radians).
!   logp(npts,ndim):  -log(press) for all 2d analysis grids. Assumed invariant
!   in assimilation window, computed fro ensemble mean at middle of window.
!   gridloc(3,npts): spherical cartesian coordinates (x,y,z) for analysis grid.
!
! Modules Used: mpisetup, params, kinds
!
! program history log:
!   2009-02-23  Initial version.
!   2016-05-02: shlyaeva: Modification for reading state vector from table
!   2016-04-20  Modify to handle the updated nemsio sig file (P, DP & DPDT removed)
!   2019-03-13  Add precipitation components
!
! attributes:
!   language: f95
!
!$$$

use mpisetup, only: nproc, mpi_integer, mpi_real4
use mpimod, only: mpi_comm_world
use params, only: datapath,nlevs,nlons,nlats,fgfileprefixes
use kinds, only: r_kind, i_kind, r_double, r_single
use constants, only: one,zero,pi,cp,rd,grav,rearth,max_varname_length
use specmod, only: init_spec_vars, isinitialized, asin_gaulats
use reducedgrid_mod, only: reducedgrid_init, regtoreduced, reducedtoreg,&
                           nptsred, lonsred, latsred
implicit none
private
public :: getgridinfo, gridinfo_cleanup
integer(i_kind),public :: nlevs_pres, idvc
real(r_single),public :: ptop
real(r_single),public, allocatable, dimension(:) :: lonsgrd, latsgrd
! arrays passed to kdtree2 routines must be single
real(r_single),public, allocatable, dimension(:,:) :: gridloc
real(r_single),public, allocatable, dimension(:,:) :: logp
integer,public :: npts
integer,public :: ntrunc
! supported variable names in anavinfo
character(len=max_varname_length),public, dimension(16) :: vars3d_supported = (/'urot', 'vrot', 'u   ', 'v   ', 'tv  ', 'q   ', 'oz  ', 'cw  ', 'tsen', 'prse', &
                                                                                'ql  ', 'qi  ', 'qr  ', 'qs  ', 'qg  ', 'dprs'/) 
character(len=max_varname_length),public, dimension(3)  :: vars2d_supported = (/'ps ', 'pst', 'sst' /)
! supported variable names in anavinfo
contains

subroutine getgridinfo(fileprefix, reducedgrid)
! read latitudes, longitudes and pressures for analysis grid,
! broadcast to each task.
use module_fv3gfs_ncio, only: Dataset, Variable, Dimension, open_dataset,&
                        read_attribute, close_dataset, get_dim, read_vardata 
implicit none

type(Dataset) :: dset
type(Dimension) :: londim,latdim,levdim
character(len=120), intent(in) :: fileprefix
logical, intent(in)            :: reducedgrid

integer(i_kind) nlevsin, ierr, iunit, k, nn, idvc 
character(len=500) filename
integer(i_kind) iret,i,j,nlonsin,nlatsin
real(r_kind), allocatable, dimension(:) :: ak,bk,spressmn
real(r_kind), allocatable, dimension(:,:) :: pressimn,presslmn,values_2d
real(r_kind) kap,kapr,kap1

iunit = 77
kap = rd/cp
kapr = cp/rd
kap1 = kap + one
nlevs_pres=nlevs+1
if (nproc .eq. 0) then
filename = trim(adjustl(datapath))//trim(adjustl(fileprefix))//"ensmean"
dset = open_dataset(filename)
londim = get_dim(dset,'grid_xt'); nlonsin = londim%len
latdim = get_dim(dset,'grid_yt'); nlatsin = latdim%len
levdim = get_dim(dset,'pfull');   nlevsin = levdim%len
idvc = 2; ntrunc = nlatsin-2
if (nlons /= nlonsin .or. nlats /= nlatsin .or. nlevs /= nlevsin) then
  print *,'incorrect dims in netcdf file'
  print *,'expected',nlons,nlats,nlevs
  print *,'got',nlonsin,nlatsin,nlevsin
  call stop2(23)
end if
end if
call mpi_bcast(ntrunc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! initialize spectral module on all tasks.
if (.not. isinitialized) call init_spec_vars(nlons,nlats,ntrunc,4)

if (nproc .eq. 0) then
   ! get pressure, lat/lon information from ensemble mean file.
   allocate(presslmn(nlons*nlats,nlevs))
   allocate(pressimn(nlons*nlats,nlevs+1))
   allocate(spressmn(nlons*nlats))
   call read_vardata(dset, 'pressfc', values_2d,errcode=iret)
   if (iret /= 0) then
      print *,'error reading ps in gridinfo_gfs'
      call stop2(11)
   endif
   ! convert to 1d array, units to millibars, flip so lats go N to S.
   spressmn = 0.01_r_kind*reshape(values_2d,(/nlons*nlats/))
   call read_attribute(dset, 'ak', ak)
   call read_attribute(dset, 'bk', bk)
   call close_dataset(dset)
   ! pressure at interfaces
   do k=1,nlevs+1
      pressimn(:,k) = 0.01_r_kind*ak(nlevs-k+2)+bk(nlevs-k+2)*spressmn(:)
   enddo
   ptop = 0.01_r_kind*ak(1)
   deallocate(ak,bk,values_2d)
   if (reducedgrid) then
      call reducedgrid_init(nlons,nlats,asin_gaulats)
      npts = nptsred
   else
      npts = nlons*nlats
   end if
   allocate(latsgrd(npts),lonsgrd(npts))
   allocate(logp(npts,nlevs_pres)) ! log(ens mean first guess press) on mid-layers
   allocate(gridloc(3,npts))
   !==> pressure at interfaces.
   if (reducedgrid) then
      lonsgrd(:) = lonsred(:)
      latsgrd(:) = latsred(:)
   else
      nn = 0
      do j=1,nlats
         do i=1,nlons
            nn = nn + 1
            lonsgrd(nn) = 2._r_single*pi*float(i-1)/nlons
            latsgrd(nn) = asin_gaulats(j)
         enddo
      enddo
   endif
   do k=1,nlevs
     ! layer pressure from Phillips vertical interpolation.
     presslmn(:,k) = ((pressimn(:,k)**kap1-pressimn(:,k+1)**kap1)/&
                      (kap1*(pressimn(:,k)-pressimn(:,k+1))))**kapr
   end do
   print *,'ensemble mean first guess surface pressure:'
   print *,minval(spressmn),maxval(spressmn)
   !do k=1,nlevs
   !   print *,'min/max ens mean press level',&
   !   k,'=',minval(presslmn(:,k)),maxval(presslmn(:,k))
   !   print *,'min/max ens mean press interface',&
   !   k,'=',minval(pressimn(:,k)),maxval(pressimn(:,k))
   !enddo
   ! logp holds log(pressure) or pseudo-height on grid, for each level/variable.
   do k=1,nlevs
      ! all variables to be updated are on mid-layers, not layer interfaces.
      if (reducedgrid) then
         call regtoreduced(presslmn(:,k),logp(:,k))
         logp(:,k) = -log(logp(:,k))
      else
         logp(:,k) = -log(presslmn(:,k))
      endif
      !print *,'min/max presslmn',k,minval(presslmn(:,k)),maxval(presslmn(:,k)),minval(logp(:,k)),maxval(logp(:,k))
   end do
   if (reducedgrid) then
      call regtoreduced(spressmn,logp(:,nlevs_pres))
      logp(:,nlevs_pres) = -log(logp(:,nlevs_pres))
   else
      logp(:,nlevs_pres) = -log(spressmn(:))
   endif
   deallocate(spressmn,presslmn,pressimn)
end if
call mpi_bcast(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (nproc .ne. 0) then
   ! allocate arrays on other (non-root) tasks
   allocate(latsgrd(npts),lonsgrd(npts))
   allocate(logp(npts,nlevs_pres)) ! log(ens mean first guess press) on mid-layers
   allocate(gridloc(3,npts))
   ! initialize reducedgrid_mod on other tasks.
   if (reducedgrid) then
      call reducedgrid_init(nlons,nlats,asin_gaulats)
   end if
endif
!call mpi_bcast(logp,npts*nlevs_pres,mpi_real4,0,MPI_COMM_WORLD,ierr)
do k=1,nlevs_pres
  call mpi_bcast(logp(1,k),npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
enddo
call mpi_bcast(lonsgrd,npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(latsgrd,npts,mpi_real4,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(ptop,1,mpi_real4,0,MPI_COMM_WORLD,ierr)
!==> precompute cartesian coords of analysis grid points.
do nn=1,npts
   gridloc(1,nn) = cos(latsgrd(nn))*cos(lonsgrd(nn))
   gridloc(2,nn) = cos(latsgrd(nn))*sin(lonsgrd(nn))
   gridloc(3,nn) = sin(latsgrd(nn))
end do

end subroutine getgridinfo

subroutine gridinfo_cleanup()
if (allocated(lonsgrd)) deallocate(lonsgrd)
if (allocated(latsgrd)) deallocate(latsgrd)
if (allocated(logp)) deallocate(logp)
if (allocated(gridloc)) deallocate(gridloc)
end subroutine gridinfo_cleanup

end module gridinfo
