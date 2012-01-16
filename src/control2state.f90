subroutine control2state_ens4dvar(xhat,sval,bval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    control2state
!   prgmmr: tremolet
!
! abstract:  Converts control variable to physical space
!
! program history log:
!   2007-04-13  tremolet - initial code
!   2008-11-28  todling  - add calc of 3dp; upd rh_to_q (Cucurull 2007-07-26)
!   2009-04-21  derber   - modify call to getuv to getuv(*,0)
!   2009-06-16  parrish  - for l_hyb_ens=.true., add calls to ensemble_forward_model and strong_bk
!   2009-08-14  lueken   - update documentation
!   2009-11-27  parrish  - for uv_hyb_ens=.true., then ensemble perturbations contain u,v instead of st,vp
!                            so introduce extra code to handle this case.
!   2010-02-21  parrish  - introduce changes to allow dual resolution, with ensemble computation on
!                            lower resolution grid compared to analysis grid.
!                            new parameter dual_res=.true. if ensemble grid is different from analysis grid.
!   2010-03-23  zhu      - use cstate for generalizing control variable
!   2010-04-29  todling  - update to use gsi_bundle; some changes toward bypassing standard atmos analysis
!   2010-05-12  todling  - rename cstate to wbundle; state_vector now a bundle
!   2010-05-31  todling  - better consistency checks; add co/co2
!                        - ready to bypass analysis of (any) meteorological fields
!   2010-06-04  parrish  - bug fix: u,v copy to wbundle after getuv for hyb ensemble
!   2010-06-15  todling  - generalized handling of chemistry
!   2011-02-20  zhu      - add gust,vis,pblh
!   2011-05-15  auligne/todling - generalized cloud handling
!   2011-12-20  Ting Lei and Xuguang Wang  - Added Ens4dvar functionality
!                           Contact:  xuguang.wang@ou.edu
!                                     University of Oklahoma      
!
!   input argument list:
!     xhat - Control variable
!     sval - State variable
!     bval - Bias predictors
!
!   output argument list:
!     sval - State variable
!     bval - Bias predictors
!
!$$$ end documentation block
use kinds, only: r_kind,i_kind
use control_vectors, only: control_vector
use control_vectors, only: cvars3d,cvars2d
use bias_predictors, only: predictors
  use mpimod, only: mype
!cltorg use gsi_4dvar, only: nsubwin, nobs_bins, l4dvar, lsqrtb
use gsi_4dvar, only: nsubwin, nobs_bins, l4dvar,lens4dvar,l_hyb_no_bc,ens4dvar_sst_update,ens4dvar_nsubwin,ens4dvar_sst_updt_scheme, k_bk_in_Jo,lsqrtb
use gridmod, only: latlon1n,latlon11
use jfunc, only: nsclen,npclen,nrclen
use hybrid_ensemble_parameters, only: l_hyb_ens,uv_hyb_ens,dual_res,oz_univ_static
use balmod, only: strong_bk
!cltorg use hybrid_ensemble_isotropic, only: ensemble_forward_model,ensemble_forward_model_dual_res
use hybrid_ensemble_isotropic, only: ensemble_forward_model,ensemble_forward_model_dual_res,&
     ensemble_forward_model_ens4dvar,ensemble_forward_model_dual_res_ens4dvar
use gsi_bundlemod, only: gsi_bundlecreate
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundlegetvar
use gsi_bundlemod, only: gsi_bundleputvar
use gsi_bundlemod, only: gsi_bundledestroy
use gsi_bundlemod, only: assignment(=)
use gsi_chemguess_mod, only: gsi_chemguess_get
use gsi_metguess_mod, only: gsi_metguess_get
use mpeu_util, only: getindex
use constants, only : max_varname_length
implicit none
  
! Declare passed variables  
type(control_vector), intent(in   ) :: xhat
!cltorg type(gsi_bundle)    , intent(inout) :: sval(nsubwin)
type(gsi_bundle)    , intent(inout) :: sval(nobs_bins)
type(predictors)    , intent(inout) :: bval

! Declare local variables  	
character(len=*),parameter::myname='control2state'
character(len=max_varname_length),allocatable,dimension(:) :: gases
character(len=max_varname_length),allocatable,dimension(:) :: clouds
!cltorg integer(i_kind) :: i,j,k,ii,jj,im,jm,km,ic,id,ngases,nclouds,istatus
integer(i_kind) :: i,j,k,kk,ii,jj,im,jm,km,ic,id,ngases,nclouds,istatus,istatus4cv_oz
real(r_kind),dimension(:,:,:),allocatable:: u,v
type(gsi_bundle):: wbundle ! work bundle
type(gsi_bundle):: wbundle1 ! work bundle needed for ens4dvar 

! Note: The following does not aim to get all variables in
!       the state and control vectors, but rather the ones
!       this routines knows how to handle.
! Declare required local control variables
integer(i_kind), parameter :: ncvars = 5
integer(i_kind) :: icps(ncvars)
integer(i_kind) :: icpblh,icgust,icvis
character(len=3), parameter :: mycvars(ncvars) = (/  &  ! vars from CV needed here
                               'sf ', 'vp ', 'ps ', 't  ',    &
                               'q  '/)
logical :: lc_sf,lc_vp,lc_ps,lc_t,lc_rh
real(r_kind),pointer,dimension(:,:)   :: cv_ps,cv_vis
real(r_kind),pointer,dimension(:,:)   :: cv_sst
real(r_kind),pointer,dimension(:,:,:) :: cv_sf,cv_vp,cv_t,cv_rh
!clt
real(r_kind),pointer,dimension(:,:,:) :: cv_oz
real(r_kind),pointer,dimension(:,:)   :: cv_ps_s,cv_sst_s
real(r_kind),pointer,dimension(:,:,:) :: cv_sf_s,cv_vp_s,cv_t_s,cv_rh_s
real(r_kind),pointer,dimension(:,:,:) :: cv_vis_s
real(r_kind),pointer,dimension(:,:,:) :: cv_oz_s
!clt real(r_kind),pointer,dimension(:,:,:) :: cv_oz_s,cv_cw_s

! Declare required local state variables
integer(i_kind), parameter :: nsvars = 5
integer(i_kind) :: isps(nsvars)
character(len=4), parameter :: mysvars(nsvars) = (/  &  ! vars from ST needed here
                               'u   ', 'v   ', 'p3d ', 'q   ', 'tsen' /)
logical :: ls_u,ls_v,ls_p3d,ls_q,ls_tsen
real(r_kind),pointer,dimension(:,:)   :: sv_ps,sv_sst
real(r_kind),pointer,dimension(:,:)   :: sv_gust,sv_vis,sv_pblh
real(r_kind),pointer,dimension(:,:,:) :: sv_u,sv_v,sv_p3d,sv_q,sv_tsen,sv_tv,sv_oz
real(r_kind),pointer,dimension(:,:,:) :: sv_rank3
real(r_kind),pointer,dimension(:,:)   :: sv_rank2

logical :: do_strong_bk,do_getprs_tl,do_normal_rh_to_q,do_tv_to_tsen,do_getuv

!******************************************************************************

if (lsqrtb) then
   write(6,*)trim(myname),': not for sqrt(B)'
   call stop2(106)
end if
if (nsubwin/=1 .and. .not.l4dvar) then
   write(6,*)trim(myname),': error 3dvar',nsubwin,l4dvar
   call stop2(107)
end if
if (ens4dvar_nsubwin/=1 .and. .not.lens4dvar) then
   write(6,*)'control2state: error 3dvar',ens4dvar_nsubwin,lens4dvar
   call stop2(107)
end if

im=xhat%step(1)%grid%im
jm=xhat%step(1)%grid%jm
km=xhat%step(1)%grid%km

! Inquire about cloud-vars
call gsi_metguess_get('clouds::3d',nclouds,istatus)
if (nclouds>0) then
    allocate(clouds(nclouds))
    call gsi_metguess_get('clouds::3d',clouds,istatus)
endif

! Inquire about chemistry
call gsi_chemguess_get('dim',ngases,istatus)
if (ngases>0) then
    allocate(gases(ngases))
    call gsi_chemguess_get('gsinames',gases,istatus)
endif

! Since each internal vector of xhat has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (xhat%step(1),mycvars,icps,istatus)
lc_sf =icps(1)>0; lc_vp =icps(2)>0; lc_ps =icps(3)>0
lc_t  =icps(4)>0; lc_rh =icps(5)>0

! Since each internal vector of xhat has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (sval(1),mysvars,isps,istatus)
ls_u  =isps(1)>0; ls_v   =isps(2)>0; ls_p3d=isps(3)>0
ls_q  =isps(4)>0; ls_tsen=isps(5)>0

! Define what to do depending on what's in CV and SV
!cltorg do_strong_bk     =lc_ps.and.lc_sf.and.lc_vp .and.lc_t
do_strong_bk     =lc_ps.and.lc_sf.and.lc_vp .and.lc_t.and.(.not.l_hyb_no_bc)
do_getprs_tl     =lc_ps.and.lc_t .and.ls_p3d
do_normal_rh_to_q=lc_rh.and.lc_t .and.ls_p3d.and.ls_q
do_tv_to_tsen    =lc_t .and.ls_q .and.ls_tsen
do_getuv         =lc_sf.and.lc_vp.and.ls_u.and.ls_v

call gsi_bundlegetpointer (xhat%step(1),'gust',icgust,istatus)
call gsi_bundlegetpointer (xhat%step(1),'vis',icvis,istatus)
call gsi_bundlegetpointer (xhat%step(1),'pblh',icpblh,istatus)

! Loop over control steps
if(.not.lens4dvar) then
do jj=1,nsubwin

!  Create a temporary bundle similar to xhat, and copy contents of xhat into it
   call gsi_bundlecreate ( wbundle, xhat%step(jj), 'control2state work', istatus )
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble creating work bundle'
      call stop2(999)
   endif
   wbundle=xhat%step(jj)

!  Get pointers to required control variables
   call gsi_bundlegetpointer (wbundle,'sf' ,cv_sf ,istatus)
   call gsi_bundlegetpointer (wbundle,'vp' ,cv_vp ,istatus)
   call gsi_bundlegetpointer (wbundle,'ps' ,cv_ps ,istatus)
   call gsi_bundlegetpointer (wbundle,'t'  ,cv_t,  istatus)
   call gsi_bundlegetpointer (wbundle,'q'  ,cv_rh ,istatus)
   if (icvis >0) call gsi_bundlegetpointer (wbundle,'vis',cv_vis,istatus)

!  Get pointers to required state variables
   call gsi_bundlegetpointer (sval(jj),'u'   ,sv_u,   istatus)
   call gsi_bundlegetpointer (sval(jj),'v'   ,sv_v,   istatus)
   call gsi_bundlegetpointer (sval(jj),'ps'  ,sv_ps,  istatus)
   call gsi_bundlegetpointer (sval(jj),'p3d' ,sv_p3d, istatus)
   call gsi_bundlegetpointer (sval(jj),'tv'  ,sv_tv,  istatus)
   call gsi_bundlegetpointer (sval(jj),'tsen',sv_tsen,istatus)
   call gsi_bundlegetpointer (sval(jj),'q'   ,sv_q ,  istatus)
   call gsi_bundlegetpointer (sval(jj),'oz'  ,sv_oz , istatus)
   call gsi_bundlegetpointer (sval(jj),'sst' ,sv_sst, istatus)
   if (icgust>0) call gsi_bundlegetpointer (sval(jj),'gust' ,sv_gust, istatus)
   if (icpblh>0) call gsi_bundlegetpointer (sval(jj),'pblh' ,sv_pblh, istatus)
   if (icvis >0) call gsi_bundlegetpointer (sval(jj),'vis'  ,sv_vis , istatus)

! If this is ensemble run, then add ensemble contribution sum(a_en(k)*xe(k)),  where a_en(k) are the ensemble
!   control variables and xe(k), k=1,n_ens are the ensemble perturbations.
   if(l_hyb_ens) then
      if(uv_hyb_ens) then
!        Convert streamfunction and velocity potential to u,v
         if (do_getuv) then
            allocate(u(im,jm,km))
            allocate(v(im,jm,km))
            call getuv(u,v,cv_sf,cv_vp,0)
            call gsi_bundleputvar ( wbundle, 'sf', u, istatus )
            call gsi_bundleputvar ( wbundle, 'vp', v, istatus )
            deallocate(v)
            deallocate(u)
         endif
      end if
      if(dual_res) then
         call ensemble_forward_model_dual_res(wbundle,xhat%aens(jj,:))
      else
         call ensemble_forward_model(wbundle,xhat%aens(jj,:))
      end if
!     Apply strong constraint to sum of static background and ensemble background combinations to
!     reduce imbalances introduced by ensemble localization in addition to known imbalances from
!     static background
      if(do_strong_bk) call strong_bk(cv_sf,cv_vp,cv_ps,cv_t)
   end if

!  Get 3d pressure
   if(do_getprs_tl) call getprs_tl(cv_ps,cv_t,sv_p3d)

!  Convert input normalized RH to q
   if(do_normal_rh_to_q) call normal_rh_to_q(cv_rh,cv_t,sv_p3d,sv_q)

!  Calculate sensible temperature
   if(do_tv_to_tsen) call tv_to_tsen(cv_t,sv_q,sv_tsen)

!  Convert streamfunction and velocity potential to u,v
   if(do_getuv) then
      if(l_hyb_ens.and.uv_hyb_ens) then
         call gsi_bundlegetvar ( wbundle, 'sf', sv_u, istatus )
         call gsi_bundlegetvar ( wbundle, 'vp', sv_v, istatus )
      else
         call getuv(sv_u,sv_v,cv_sf,cv_vp,0)
      end if
   end if

!  Convert log(vis) to vis
   if (icvis >0)  call logvis_to_vis(cv_vis,sv_vis)

!  Copy other variables
   call gsi_bundlegetvar ( wbundle, 't'  , sv_tv,  istatus )
   call gsi_bundlegetvar ( wbundle, 'oz' , sv_oz,  istatus )
   call gsi_bundlegetvar ( wbundle, 'ps' , sv_ps,  istatus )
   call gsi_bundlegetvar ( wbundle, 'sst', sv_sst, istatus )
   if (icgust>0) call gsi_bundlegetvar ( wbundle, 'gust', sv_gust, istatus )
   if (icpblh>0) call gsi_bundlegetvar ( wbundle, 'pblh', sv_pblh, istatus )

!  Since cloud-vars map one-to-one, take care of them together
   do ic=1,nclouds
      id=getindex(cvars3d,clouds(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),clouds(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle, clouds(ic),sv_rank3,istatus)
      endif
   enddo

!  Same one-to-one map for chemistry-vars; take care of them together 
   do ic=1,ngases
      id=getindex(cvars3d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),gases(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle, gases(ic),sv_rank3,istatus)
      endif
      id=getindex(cvars2d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),gases(ic),sv_rank2,istatus)
          call gsi_bundlegetvar     (wbundle, gases(ic),sv_rank2,istatus)
      endif
   enddo

   call gsi_bundledestroy(wbundle,istatus)
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble destroying work bundle'
      call stop2(999)
   endif

end do
else ! for ens4dvar
!Note, nsubwin = 1 

do jj=1,nsubwin

!  Create a temporary bundle similar to xhat, and copy contents of xhat into it
   call gsi_bundlecreate ( wbundle, xhat%step(jj), 'control2state work', istatus )
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble creating work bundle'
      call stop2(999)
   endif
   wbundle=xhat%step(jj)

!  Get pointers to required control variables
   call gsi_bundlegetpointer (wbundle,'sf' ,cv_sf ,istatus)
   call gsi_bundlegetpointer (wbundle,'vp' ,cv_vp ,istatus)
   call gsi_bundlegetpointer (wbundle,'ps' ,cv_ps ,istatus)
   call gsi_bundlegetpointer (wbundle,'t'  ,cv_t,  istatus)
   call gsi_bundlegetpointer (wbundle,'q'  ,cv_rh ,istatus)
   call gsi_bundlegetpointer (wbundle,'oz'  ,cv_oz ,istatus)
   if (icvis >0) call gsi_bundlegetpointer (wbundle,'vis',cv_vis,istatus)

  call gsi_bundlegetpointer (wbundle,'sst'  ,cv_sst ,istatus) ! clt
!clt 
  do kk=1,nobs_bins
   call gsi_bundlecreate ( wbundle1, xhat%step(jj), 'control2state work1', istatus )
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble creating work bundle'
      call stop2(999)
   endif
   call gsi_bundlegetpointer (wbundle1,'sf' ,cv_sf_s ,istatus)
   call gsi_bundlegetpointer (wbundle1,'vp' ,cv_vp_s ,istatus)
   call gsi_bundlegetpointer (wbundle1,'ps' ,cv_ps_s ,istatus)
   call gsi_bundlegetpointer (wbundle1,'t'  ,cv_t_s,  istatus)
   call gsi_bundlegetpointer (wbundle1,'q'  ,cv_rh_s ,istatus)
   call gsi_bundlegetpointer (wbundle1,'oz'  ,cv_oz_s,  istatus4cv_oz)

   call gsi_bundlegetpointer (wbundle1,'sst'  ,cv_sst_s,  istatus)
   if (icvis>0) call gsi_bundlegetpointer (wbundle1,'vis'  ,cv_vis_s ,istatus)
  cv_sf_s=0
   cv_vp_s=0
   cv_ps_s=0
   cv_t_s=0
   cv_sst_s=0
   if(icvis>0) cv_vis_s=0

!  Get pointers to required state variables
   call gsi_bundlegetpointer (sval(kk),'u'   ,sv_u,   istatus)
   call gsi_bundlegetpointer (sval(kk),'v'   ,sv_v,   istatus)
   call gsi_bundlegetpointer (sval(kk),'ps'  ,sv_ps,  istatus)
   call gsi_bundlegetpointer (sval(kk),'p3d' ,sv_p3d, istatus)
   call gsi_bundlegetpointer (sval(kk),'tv'  ,sv_tv,  istatus)
   call gsi_bundlegetpointer (sval(kk),'tsen',sv_tsen,istatus)
   call gsi_bundlegetpointer (sval(kk),'q'   ,sv_q ,  istatus)
   call gsi_bundlegetpointer (sval(kk),'oz'  ,sv_oz , istatus)
   call gsi_bundlegetpointer (sval(kk),'sst' ,sv_sst, istatus)
   if (icgust>0) call gsi_bundlegetpointer (sval(kk),'gust' ,sv_gust, istatus)
   if (icpblh>0) call gsi_bundlegetpointer (sval(kk),'pblh' ,sv_pblh, istatus)
   if (icvis >0) call gsi_bundlegetpointer (sval(kk),'vis'  ,sv_vis , istatus)
 if(mype.eq.0)  write(6,*)'think-1,in control2state k_bk_in_JO ,kk is ',k_bk_in_Jo,kk
!cltorg   if(k_bk_in_Jo.eq.kk) then
   if(k_bk_in_Jo.eq.2.or.(k_bk_in_Jo-20.eq.kk)) then
 if(mype.eq.0)  write(6,*)'think,in control2state k_bk_in_JO is ',k_bk_in_Jo
      cv_sf_s=cv_sf
      cv_vp_s=cv_vp
      cv_t_s =cv_t
      cv_rh_s=cv_rh
 if(istatus4cv_oz)     cv_oz_s=cv_oz
!clt      cv_cw_s=cv_cw
      cv_ps_s=cv_ps
      cv_sst_s=cv_sst

   else
      cv_sf_s=0
      cv_vp_s=0
      cv_t_s =0
      cv_rh_s=0
if(istatus4cv_oz)      cv_oz_s=0
!clt      cv_cw_s=0
      cv_ps_s=0

 if(ens4dvar_sst_update) then
   if(ens4dvar_sst_updt_scheme.eq.0) then
!clt, now, only sst_updt_scheme =0  is available
   cv_sst_s=cv_sst
   else
      cv_sst_s=0
   endif
  endif

 if(istatus4cv_oz) then
 if(oz_univ_static) then
   cv_oz_s=cv_oz
   else
      cv_oz_s=0
   endif
  endif



   
   endif ! k_bk_in_Jo


! If this is ensemble run, then add ensemble contribution sum(a_en(k)*xe(k)),  where a_en(k) are the ensemble
!   control variables and xe(k), k=1,n_ens are the ensemble perturbations.
   if(l_hyb_ens) then
      if(uv_hyb_ens) then
!        Convert streamfunction and velocity potential to u,v
         if (do_getuv) then
            allocate(u(im,jm,km))
            allocate(v(im,jm,km))
            call getuv(u,v,cv_sf_s,cv_vp_s,0)
            call gsi_bundleputvar ( wbundle1, 'sf', u, istatus )
            call gsi_bundleputvar ( wbundle1, 'vp', v, istatus )
            deallocate(v)
            deallocate(u)
         endif
      end if
      if(dual_res) then
         call ensemble_forward_model_dual_res_ens4dvar(kk,wbundle1,xhat%aens(jj,:))
      else
         call ensemble_forward_model_ens4dvar(kk,wbundle1,xhat%aens(jj,:))
      end if
!     Apply strong constraint to sum of static background and ensemble background combinations to
!     reduce imbalances introduced by ensemble localization in addition to known imbalances from
!     static background
      if(do_strong_bk) call strong_bk(cv_sf_s,cv_vp_s,cv_ps_s,cv_t_s)
   end if

!  Get 3d pressure
   if(do_getprs_tl) call getprs_tl(cv_ps_s,cv_t_s,sv_p3d)

!  Convert input normalized RH to q
   if(do_normal_rh_to_q) call normal_rh_to_q(cv_rh_s,cv_t_s,sv_p3d,sv_q)

!  Calculate sensible temperature
   if(do_tv_to_tsen) call tv_to_tsen(cv_t_s,sv_q,sv_tsen)

!  Convert streamfunction and velocity potential to u,v
   if(do_getuv) then
      if(l_hyb_ens.and.uv_hyb_ens) then
         call gsi_bundlegetvar ( wbundle1, 'sf', sv_u, istatus )
         call gsi_bundlegetvar ( wbundle1, 'vp', sv_v, istatus )
      else
         call getuv(sv_u,sv_v,cv_sf_s,cv_vp_s,0)
      end if
   end if

!  Convert log(vis) to vis
   if (icvis >0)  call logvis_to_vis(cv_vis_s,sv_vis)

!  Copy other variables
   call gsi_bundlegetvar ( wbundle1, 't'  , sv_tv,  istatus )
   call gsi_bundlegetvar ( wbundle1, 'oz' , sv_oz,  istatus )
   call gsi_bundlegetvar ( wbundle1, 'ps' , sv_ps,  istatus )
   call gsi_bundlegetvar ( wbundle1, 'sst', sv_sst, istatus )
   if (icgust>0) call gsi_bundlegetvar ( wbundle1, 'gust', sv_gust, istatus )
   if (icpblh>0) call gsi_bundlegetvar ( wbundle1, 'pblh', sv_pblh, istatus )

!  Since cloud-vars map one-to-one, take care of them together
   do ic=1,nclouds
      id=getindex(cvars3d,clouds(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(kk),clouds(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle1, clouds(ic),sv_rank3,istatus)
      endif
   enddo

!  Same one-to-one map for chemistry-vars; take care of them together 
   do ic=1,ngases
      id=getindex(cvars3d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(kk),gases(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle1, gases(ic),sv_rank3,istatus)
      endif
      id=getindex(cvars2d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(kk),gases(ic),sv_rank2,istatus)
          call gsi_bundlegetvar     (wbundle1, gases(ic),sv_rank2,istatus)
      endif
   enddo

   call gsi_bundledestroy(wbundle1,istatus)
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble destroying work bundle'
      call stop2(999)
   endif

end do ! nobs_bins

   call gsi_bundledestroy(wbundle,istatus)
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble destroying work bundle'
      call stop2(999)
   endif

end do

endif ! for lens4dvar


! Biases
do ii=1,nsclen
   bval%predr(ii)=xhat%predr(ii)
enddo

do ii=1,npclen
   bval%predp(ii)=xhat%predp(ii)
enddo

! Clean up
if (ngases>0) then
    deallocate(gases)
endif

return
end subroutine control2state_ens4dvar
subroutine control2state(xhat,sval,bval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    control2state
!   prgmmr: tremolet
!
! abstract:  Converts control variable to physical space
!
! program history log:
!   2007-04-13  tremolet - initial code
!   2008-11-28  todling  - add calc of 3dp; upd rh_to_q (Cucurull 2007-07-26)
!   2009-04-21  derber   - modify call to getuv to getuv(*,0)
!   2009-06-16  parrish  - for l_hyb_ens=.true., add calls to ensemble_forward_model and strong_bk
!   2009-08-14  lueken   - update documentation
!   2009-11-27  parrish  - for uv_hyb_ens=.true., then ensemble perturbations contain u,v instead of st,vp
!                            so introduce extra code to handle this case.
!   2010-02-21  parrish  - introduce changes to allow dual resolution, with ensemble computation on
!                            lower resolution grid compared to analysis grid.
!                            new parameter dual_res=.true. if ensemble grid is different from analysis grid.
!   2010-03-23  zhu      - use cstate for generalizing control variable
!   2010-04-29  todling  - update to use gsi_bundle; some changes toward bypassing standard atmos analysis
!   2010-05-12  todling  - rename cstate to wbundle; state_vector now a bundle
!   2010-05-31  todling  - better consistency checks; add co/co2
!                        - ready to bypass analysis of (any) meteorological fields
!   2010-06-04  parrish  - bug fix: u,v copy to wbundle after getuv for hyb ensemble
!   2010-06-15  todling  - generalized handling of chemistry
!   2011-02-20  zhu      - add gust,vis,pblh
!   2011-05-15  auligne/todling - generalized cloud handling
!
!   input argument list:
!     xhat - Control variable
!     sval - State variable
!     bval - Bias predictors
!
!   output argument list:
!     sval - State variable
!     bval - Bias predictors
!
!$$$ end documentation block
use kinds, only: r_kind,i_kind
use control_vectors, only: control_vector
use control_vectors, only: cvars3d,cvars2d
use bias_predictors, only: predictors
use gsi_4dvar, only: nsubwin, nobs_bins, l4dvar, lsqrtb
use gridmod, only: latlon1n,latlon11
use jfunc, only: nsclen,npclen,nrclen
use hybrid_ensemble_parameters, only: l_hyb_ens,uv_hyb_ens,dual_res
use balmod, only: strong_bk
use hybrid_ensemble_isotropic, only: ensemble_forward_model,ensemble_forward_model_dual_res
use gsi_bundlemod, only: gsi_bundlecreate
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: gsi_bundlegetvar
use gsi_bundlemod, only: gsi_bundleputvar
use gsi_bundlemod, only: gsi_bundledestroy
use gsi_bundlemod, only: assignment(=)
use gsi_chemguess_mod, only: gsi_chemguess_get
use gsi_metguess_mod, only: gsi_metguess_get
use mpeu_util, only: getindex
use constants, only : max_varname_length
implicit none
  
! Declare passed variables  
type(control_vector), intent(in   ) :: xhat
type(gsi_bundle)    , intent(inout) :: sval(nsubwin)
type(predictors)    , intent(inout) :: bval

! Declare local variables  	
character(len=*),parameter::myname='control2state'
character(len=max_varname_length),allocatable,dimension(:) :: gases
character(len=max_varname_length),allocatable,dimension(:) :: clouds
integer(i_kind) :: i,j,k,ii,jj,im,jm,km,ic,id,ngases,nclouds,istatus
real(r_kind),dimension(:,:,:),allocatable:: u,v
type(gsi_bundle):: wbundle ! work bundle

! Note: The following does not aim to get all variables in
!       the state and control vectors, but rather the ones
!       this routines knows how to handle.
! Declare required local control variables
integer(i_kind), parameter :: ncvars = 5
integer(i_kind) :: icps(ncvars)
integer(i_kind) :: icpblh,icgust,icvis
character(len=3), parameter :: mycvars(ncvars) = (/  &  ! vars from CV needed here
                               'sf ', 'vp ', 'ps ', 't  ',    &
                               'q  '/)
logical :: lc_sf,lc_vp,lc_ps,lc_t,lc_rh
real(r_kind),pointer,dimension(:,:)   :: cv_ps,cv_vis
real(r_kind),pointer,dimension(:,:,:) :: cv_sf,cv_vp,cv_t,cv_rh

! Declare required local state variables
integer(i_kind), parameter :: nsvars = 5
integer(i_kind) :: isps(nsvars)
character(len=4), parameter :: mysvars(nsvars) = (/  &  ! vars from ST needed here
                               'u   ', 'v   ', 'p3d ', 'q   ', 'tsen' /)
logical :: ls_u,ls_v,ls_p3d,ls_q,ls_tsen
real(r_kind),pointer,dimension(:,:)   :: sv_ps,sv_sst
real(r_kind),pointer,dimension(:,:)   :: sv_gust,sv_vis,sv_pblh
real(r_kind),pointer,dimension(:,:,:) :: sv_u,sv_v,sv_p3d,sv_q,sv_tsen,sv_tv,sv_oz
real(r_kind),pointer,dimension(:,:,:) :: sv_rank3
real(r_kind),pointer,dimension(:,:)   :: sv_rank2

logical :: do_strong_bk,do_getprs_tl,do_normal_rh_to_q,do_tv_to_tsen,do_getuv

!******************************************************************************

if (lsqrtb) then
   write(6,*)trim(myname),': not for sqrt(B)'
   call stop2(106)
end if
if (nsubwin/=1 .and. .not.l4dvar) then
   write(6,*)trim(myname),': error 3dvar',nsubwin,l4dvar
   call stop2(107)
end if

im=xhat%step(1)%grid%im
jm=xhat%step(1)%grid%jm
km=xhat%step(1)%grid%km

! Inquire about cloud-vars
call gsi_metguess_get('clouds::3d',nclouds,istatus)
if (nclouds>0) then
    allocate(clouds(nclouds))
    call gsi_metguess_get('clouds::3d',clouds,istatus)
endif

! Inquire about chemistry
call gsi_chemguess_get('dim',ngases,istatus)
if (ngases>0) then
    allocate(gases(ngases))
    call gsi_chemguess_get('gsinames',gases,istatus)
endif

! Since each internal vector of xhat has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (xhat%step(1),mycvars,icps,istatus)
lc_sf =icps(1)>0; lc_vp =icps(2)>0; lc_ps =icps(3)>0
lc_t  =icps(4)>0; lc_rh =icps(5)>0

! Since each internal vector of xhat has the same structure, pointers are
! the same independent of the subwindow jj
call gsi_bundlegetpointer (sval(1),mysvars,isps,istatus)
ls_u  =isps(1)>0; ls_v   =isps(2)>0; ls_p3d=isps(3)>0
ls_q  =isps(4)>0; ls_tsen=isps(5)>0

! Define what to do depending on what's in CV and SV
do_strong_bk     =lc_ps.and.lc_sf.and.lc_vp .and.lc_t
do_getprs_tl     =lc_ps.and.lc_t .and.ls_p3d
do_normal_rh_to_q=lc_rh.and.lc_t .and.ls_p3d.and.ls_q
do_tv_to_tsen    =lc_t .and.ls_q .and.ls_tsen
do_getuv         =lc_sf.and.lc_vp.and.ls_u.and.ls_v

call gsi_bundlegetpointer (xhat%step(1),'gust',icgust,istatus)
call gsi_bundlegetpointer (xhat%step(1),'vis',icvis,istatus)
call gsi_bundlegetpointer (xhat%step(1),'pblh',icpblh,istatus)

! Loop over control steps
do jj=1,nsubwin

!  Create a temporary bundle similar to xhat, and copy contents of xhat into it
   call gsi_bundlecreate ( wbundle, xhat%step(jj), 'control2state work', istatus )
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble creating work bundle'
      call stop2(999)
   endif
   wbundle=xhat%step(jj)

!  Get pointers to required control variables
   call gsi_bundlegetpointer (wbundle,'sf' ,cv_sf ,istatus)
   call gsi_bundlegetpointer (wbundle,'vp' ,cv_vp ,istatus)
   call gsi_bundlegetpointer (wbundle,'ps' ,cv_ps ,istatus)
   call gsi_bundlegetpointer (wbundle,'t'  ,cv_t,  istatus)
   call gsi_bundlegetpointer (wbundle,'q'  ,cv_rh ,istatus)
   if (icvis >0) call gsi_bundlegetpointer (wbundle,'vis',cv_vis,istatus)

!  Get pointers to required state variables
   call gsi_bundlegetpointer (sval(jj),'u'   ,sv_u,   istatus)
   call gsi_bundlegetpointer (sval(jj),'v'   ,sv_v,   istatus)
   call gsi_bundlegetpointer (sval(jj),'ps'  ,sv_ps,  istatus)
   call gsi_bundlegetpointer (sval(jj),'p3d' ,sv_p3d, istatus)
   call gsi_bundlegetpointer (sval(jj),'tv'  ,sv_tv,  istatus)
   call gsi_bundlegetpointer (sval(jj),'tsen',sv_tsen,istatus)
   call gsi_bundlegetpointer (sval(jj),'q'   ,sv_q ,  istatus)
   call gsi_bundlegetpointer (sval(jj),'oz'  ,sv_oz , istatus)
   call gsi_bundlegetpointer (sval(jj),'sst' ,sv_sst, istatus)
   if (icgust>0) call gsi_bundlegetpointer (sval(jj),'gust' ,sv_gust, istatus)
   if (icpblh>0) call gsi_bundlegetpointer (sval(jj),'pblh' ,sv_pblh, istatus)
   if (icvis >0) call gsi_bundlegetpointer (sval(jj),'vis'  ,sv_vis , istatus)

! If this is ensemble run, then add ensemble contribution sum(a_en(k)*xe(k)),  where a_en(k) are the ensemble
!   control variables and xe(k), k=1,n_ens are the ensemble perturbations.
   if(l_hyb_ens) then
      if(uv_hyb_ens) then
!        Convert streamfunction and velocity potential to u,v
         if (do_getuv) then
            allocate(u(im,jm,km))
            allocate(v(im,jm,km))
            call getuv(u,v,cv_sf,cv_vp,0)
            call gsi_bundleputvar ( wbundle, 'sf', u, istatus )
            call gsi_bundleputvar ( wbundle, 'vp', v, istatus )
            deallocate(v)
            deallocate(u)
         endif
      end if
      if(dual_res) then
         call ensemble_forward_model_dual_res(wbundle,xhat%aens(jj,:))
      else
         call ensemble_forward_model(wbundle,xhat%aens(jj,:))
      end if
!     Apply strong constraint to sum of static background and ensemble background combinations to
!     reduce imbalances introduced by ensemble localization in addition to known imbalances from
!     static background
      if(do_strong_bk) call strong_bk(cv_sf,cv_vp,cv_ps,cv_t)
   end if

!  Get 3d pressure
   if(do_getprs_tl) call getprs_tl(cv_ps,cv_t,sv_p3d)

!  Convert input normalized RH to q
   if(do_normal_rh_to_q) call normal_rh_to_q(cv_rh,cv_t,sv_p3d,sv_q)

!  Calculate sensible temperature
   if(do_tv_to_tsen) call tv_to_tsen(cv_t,sv_q,sv_tsen)

!  Convert streamfunction and velocity potential to u,v
   if(do_getuv) then
      if(l_hyb_ens.and.uv_hyb_ens) then
         call gsi_bundlegetvar ( wbundle, 'sf', sv_u, istatus )
         call gsi_bundlegetvar ( wbundle, 'vp', sv_v, istatus )
      else
         call getuv(sv_u,sv_v,cv_sf,cv_vp,0)
      end if
   end if

!  Convert log(vis) to vis
   if (icvis >0)  call logvis_to_vis(cv_vis,sv_vis)

!  Copy other variables
   call gsi_bundlegetvar ( wbundle, 't'  , sv_tv,  istatus )
   call gsi_bundlegetvar ( wbundle, 'oz' , sv_oz,  istatus )
   call gsi_bundlegetvar ( wbundle, 'ps' , sv_ps,  istatus )
   call gsi_bundlegetvar ( wbundle, 'sst', sv_sst, istatus )
   if (icgust>0) call gsi_bundlegetvar ( wbundle, 'gust', sv_gust, istatus )
   if (icpblh>0) call gsi_bundlegetvar ( wbundle, 'pblh', sv_pblh, istatus )

!  Since cloud-vars map one-to-one, take care of them together
   do ic=1,nclouds
      id=getindex(cvars3d,clouds(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),clouds(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle, clouds(ic),sv_rank3,istatus)
      endif
   enddo

!  Same one-to-one map for chemistry-vars; take care of them together 
   do ic=1,ngases
      id=getindex(cvars3d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),gases(ic),sv_rank3,istatus)
          call gsi_bundlegetvar     (wbundle, gases(ic),sv_rank3,istatus)
      endif
      id=getindex(cvars2d,gases(ic))
      if (id>0) then
          call gsi_bundlegetpointer (sval(jj),gases(ic),sv_rank2,istatus)
          call gsi_bundlegetvar     (wbundle, gases(ic),sv_rank2,istatus)
      endif
   enddo

   call gsi_bundledestroy(wbundle,istatus)
   if(istatus/=0) then
      write(6,*) trim(myname), ': trouble destroying work bundle'
      call stop2(999)
   endif

end do

! Biases
do ii=1,nsclen
   bval%predr(ii)=xhat%predr(ii)
enddo

do ii=1,npclen
   bval%predp(ii)=xhat%predp(ii)
enddo

! Clean up
if (ngases>0) then
    deallocate(gases)
endif

return
end subroutine control2state
