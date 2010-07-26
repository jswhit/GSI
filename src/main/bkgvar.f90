subroutine bkgvar(cvec,sst,slndt,sicet,iflg)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bkgvar      apply background error variances
!   prgmmr: parrish          org: np22                date: 1990-10-06
!
! abstract: apply latitudinal background error variances & manipulate
!           skin temp <--> sst,sfc temp, and ice temp fields
!
! program history log:
!   1990-10-06  parrish
!   2004-08-24  kleist - hoper & htoper replaced
!   2004-11-16  treadon - add longitude dimension to variance array dssv
!   2004-11-22  derber - modify for openMP
!   2005-01-22  parrish - add "use balmod"
!   2005-07-14  wu - add max bound to l2
!   2007-03-13  derber - modify to allow use qvar3d array
!   2007-07-03  kleist - add full 2d error array for surface pressure (global only)
!   2007-11-26  s.liu - correct bug in water point skin temperature variances
!   2010-03-01  zhu   - replace explicit use of each control variable by one array
!                       'cstate' and use nrf* for generalized control variable
!                     - merge global and regional cases
!   2010-05-06  todling - use gsi_bundle
!   2010-06-03  todling - protection for mvars<2
!   2010-07-07  todling - rename cstate to cvec for clarity
!
!   input argument list:
!     t        - t grid values
!     p        - p surface grid values
!     q        - q grid values
!     oz       - ozone grid values
!     skint    - skin temperature grid values
!     cwmr     - cloud water mixing ratio grid values
!     st       - streamfunction grid values
!     vp       - velocity potential grid values
!     sst      - sst grid values
!     slndt    - land surface temperature grid values
!     sicet    - snow/ice covered surface temperature grid values
!     iflg     - flag for skin temperature manipulation
!                0: skint --> sst,slndt,sicet
!                1: sst,slndt,sicet --> skint
!
!   output argument list:
!     t        - t grid values
!     p        - p surface grid values
!     q        - q grid values
!     oz       - ozone grid values
!     skint    - skin temperature grid values
!     cwmr     - cloud water mixing ratio grid values
!     st       - streamfunction grid values
!     vp       - velocity potential grid values
!     sst      - sst grid values
!     slndt    - land surface temperature grid values
!     sicet    - snow/ice covered surface temperature grid values
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: izero,ione,one
  use balmod, only: rllat1,llmax
  use berror, only: dssv,dssvs
  use gridmod, only: nsig,regional,lat2,lon2
  use guess_grids, only: isli2
  use jfunc, only: nval_levs
  use control_vectors, only: mvars, nc2d
  use gsi_bundlemod, only : gsi_bundle
  use gsi_bundlemod, only : gsi_bundlegetpointer
  implicit none

! Declare passed variables
  integer(i_kind),intent(in   ) :: iflg
  type(gsi_bundle),intent(inout) :: cvec
  real(r_kind),dimension(lat2,lon2),intent(inout) :: sst,slndt,sicet

! Declare local variables
  integer(i_kind) i,j,k,n,i_sst,istatus
  real(r_kind) dl1,dl2

! Multipy by variances
!$omp parallel do  schedule(dynamic,1) private(n,k,i,j)
  do n=1,cvec%n3d   ! _RT: must map dssv to this (assumes same order)
     do k=1,nsig
        do i=1,lon2
           do j=1,lat2
              cvec%r3(n)%q(j,i,k)  =cvec%r3(n)%q(j,i,k)*dssv(j,i,k,n)
           end do
        enddo
     enddo
  end do

! Get pointer for SST
  call gsi_bundlegetpointer(cvec,'sst',i_sst,istatus)

! Surface fields
!$omp parallel do  schedule(dynamic,1) private(n,i,j)
  do n=1,cvec%n2d
     if (n/=i_sst) then
        do i=1,lon2
           do j=1,lat2
              cvec%r2(n)%q(j,i)=cvec%r2(n)%q(j,i)*dssvs(j,i,n)
           end do
        end do
     else
        if (mvars>=2) then
           if (iflg == izero) then
!          Break skin temperature into components
               do i=1,lon2
                  do j=1,lat2
                     if(isli2(j,i) == ione) then
                        slndt(j,i)=cvec%r2(n)%q(j,i)*dssvs(j,i,nc2d+1)
                     else if(isli2(j,i) == 2) then
                        sicet(j,i)=cvec%r2(n)%q(j,i)*dssvs(j,i,nc2d+2)
                     else
                        sst(j,i)  =cvec%r2(n)%q(j,i)*dssvs(j,i,n)
                     end if
                  end do
               end do
           else
!          Combine sst,slndt, and sicet into skin temperature field
              do i=1,lon2
                 do j=1,lat2
                    if(isli2(j,i) == ione) then
                       cvec%r2(n)%q(j,i)=slndt(j,i)*dssvs(j,i,nc2d+1)
                    else if(isli2(j,i) == 2) then
                       cvec%r2(n)%q(j,i)=sicet(j,i)*dssvs(j,i,nc2d+2)
                    else
                       cvec%r2(n)%q(j,i)=sst(j,i)*dssvs(j,i,n)
                    end if
                 end do
              end do
           end if
        end if ! mvars
     end if
  end do

  return
end subroutine bkgvar

subroutine bkg_stddev(cvec,svec)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bkg_variances      apply background error variances
!   prgmmr: el akkraoui          org: gmao              date: 2010-06-05
!
! abstract: retrieve background error standard deviations including 
!           flow dependent part
!
! program history log:
!   2010-06-05  el akkraoui
!   2010-07-08  todling - revisit original code
!
!   input argument list:
!     cvec - allocated bundle in control space
!     svec - allocated bundle in state   space
!
!   output argument list:
!     cvec - bundle holding standard deviations in control space
!     svec - bundle holding standard deviations in state   space
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use mpimod, only : mype
  use constants, only: one
  use berror, only: bkgv_flowdep
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_bundlemod, only: assignment(=)
  use gridmod, only: lat2,lon2,nsig

  implicit none
  
! Declare passed variables  
  type(gsi_bundle), intent(inout) :: cvec
  type(gsi_bundle), intent(inout) :: svec

! Declare local variables  	
  integer(i_kind) :: ii,jj,istatus
  real(r_kind),dimension(lat2,lon2) :: sst,slndt,sicet
  real(r_kind),pointer,dimension(:,:)   :: cv_ps
  real(r_kind),pointer,dimension(:,:,:) :: cv_t,cv_sf,cv_vp,cv_rh
  real(r_kind),pointer,dimension(:,:,:) :: sv_tsen,sv_u,sv_v,sv_q,sv_p3d
  logical do_flow_dep,do_getprs_tl,do_normal_rh_to_q,do_tv_to_tsen,do_getuv

! Declare required local control variables
  integer(i_kind), parameter :: ncvars = 5
  integer(i_kind) :: icps(ncvars)
  character(len=3), parameter :: mycvars(ncvars) = (/  &  ! vars from CV needed here
                               'sf ', 'vp ', 'ps ', 't  ', 'q  '/)
  logical lc_sf,lc_vp,lc_t,lc_ps,lc_rh

! Declare required local state variables
  integer(i_kind), parameter :: nsvars = 5
  integer(i_kind)            :: isps(nsvars)
  character(len=4), parameter :: mysvars(nsvars) = (/  &  ! vars from ST needed here
                               'u   ', 'v   ', 'p3d ', 'q   ', 'tsen' /)
  logical ls_u,ls_v,ls_tsen,ls_p3d,ls_q

! Check presence of fields in control bundle
  call gsi_bundlegetpointer (cvec,mycvars,icps,istatus)
  lc_sf =icps(1)>0; lc_vp =icps(2)>0; lc_ps =icps(3)>0
  lc_t  =icps(4)>0; lc_rh =icps(5)>0

! Check presence of fields in state bundle
  call gsi_bundlegetpointer (svec,mysvars,isps,istatus)
  ls_u  =isps(1)>0; ls_v   =isps(2)>0; ls_p3d=isps(3)>0
  ls_q  =isps(4)>0; ls_tsen=isps(5)>0

! Determine what to do given what's available
  do_flow_dep      =lc_sf.and.lc_vp.and.lc_t.and.lc_ps
  do_getprs_tl     =lc_ps.and.lc_t .and.ls_p3d
  do_normal_rh_to_q=lc_rh.and.lc_t .and.ls_p3d.and.ls_q
  do_tv_to_tsen    =lc_t .and.ls_q .and.ls_tsen
  do_getuv         =lc_sf.and.lc_vp.and.ls_u.and.ls_v

  cvec =one

  sst  =one
  slndt=one
  sicet=one

! Get standard deviations (why is this called bkgvar?) in control space
  call bkgvar(cvec,sst,slndt,sicet,1)

  call gsi_bundlegetpointer (cvec,'sf',cv_sf,istatus)
  call gsi_bundlegetpointer (cvec,'vp',cv_vp,istatus)
  call gsi_bundlegetpointer (cvec,'t' ,cv_t ,istatus)
  call gsi_bundlegetpointer (cvec,'ps',cv_ps,istatus)
  call gsi_bundlegetpointer (cvec,'q' ,cv_rh,istatus)

  call gsi_bundlegetpointer (svec,'u'   ,sv_u    ,istatus)
  call gsi_bundlegetpointer (svec,'v'   ,sv_v    ,istatus)
  call gsi_bundlegetpointer (svec,'tsen',sv_tsen ,istatus)
  call gsi_bundlegetpointer (svec,'p3d' ,sv_p3d  ,istatus)
  call gsi_bundlegetpointer (svec,'q'   ,sv_q    ,istatus)

! Add flow dependent part to standard deviations
  if(bkgv_flowdep) then
      if(do_flow_dep) call bkgvar_rewgt(cv_sf,cv_vp,cv_t,cv_ps,mype)
  endif

!  Get 3d pressure
   if(do_getprs_tl) call getprs_tl(cv_ps,cv_t,sv_p3d)

!  Convert input normalized RH to q
   if(do_normal_rh_to_q) call normal_rh_to_q(cv_rh,cv_t,sv_p3d,sv_q)
   
!  Calculate sensible temperature
   if(do_tv_to_tsen) call tv_to_tsen(cv_t,sv_q,sv_tsen)

!  Convert streamfunction and velocity potential to u and v
   if(do_getuv) then
      call getuv(sv_u,sv_v,cv_sf,cv_vp,0)
   end if

!  TO BE DONE: handle the rest of CV and SV fields

end subroutine bkg_stddev
