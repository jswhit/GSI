module enkf_nlbalmod
 use constants, only: zero,one,cp,fv,rd,omega,tiny_r_kind,max_varname_length,t0c,r0_05,&
                      init_constants,init_constants_derived
 use kinds, only: i_kind,r_double,r_kind,r_single
 use specmod, only: sptezv_s, sptez_s, init_spec_vars, ndimspec => nc, getgrad,&
                    isinitialized, lap, invlap, ncd2, gaulats, areameanwts, imax, jmax
 implicit none
 private
 public :: getbaltps_fromvrt

 contains

 subroutine getbaltps_fromvrt(ug, vg, upsig, vpsig, ak, bk, bk_top, bk_bot, tvg_ref, psg_ref, tvg, psg,&
                              tvgbal, psgbal, nlons, nlats, nlevs, nmem)
 real(r_kind), intent(in), dimension(nlevs+1) :: ak,bk
 real(r_kind), intent(in) :: bk_top,bk_bot
 real(r_kind), intent(in) :: psg_ref(nlons*nlats)
 real(r_kind), intent(in) :: tvg_ref(nlons*nlats,nlevs)
 real(r_kind), intent(in) :: psg(nlons*nlats)
 integer(i_kind), intent(in) :: nlons,nlats,nlevs,nmem
 real(r_kind), intent(out) :: psgbal(nlons*nlats)
 real(r_kind), intent(inout) :: ug(nlons*nlats,nlevs),vg(nlons*nlats,nlevs),tvg(nlons*nlats,nlevs)
 real(r_kind), intent(out)   :: tvgbal(nlons*nlats,nlevs)
 real(r_kind), intent(out)   :: upsig(nlons*nlats,nlevs),vpsig(nlons*nlats,nlevs)

 real(r_kind) taper_pbl(nlevs+1),ptop,psmean
 real(r_kind), dimension(ndimspec)  :: vrtspec,divspec,lnpsspec
 real(r_kind), dimension(nlevs) :: tvmean
 real(r_kind), dimension(nlons*nlats) :: &
    ug2,vg2,lnpsx,lnpsy,lnpsrx,lnpsry
 real(r_kind), dimension(nlons*nlats,nlevs) :: vrtg,phigbal
 real(r_kind), dimension(nlons*nlats,nlevs+1) :: pressi,phigbali
 real(r_kind), dimension(ndimspec,nlevs)  :: Gspec

 !real(r_kind),allocatable,dimension(:,:) :: pressi,phigbali,Gspec,vrtg,phigbal

 integer(i_kind) k,nn,j,i,niter,nlevp
 logical ::  nobal=.false. ! for debugging

 !nobal = .true.

 call init_constants(.false.) ! .false. means global
 call init_constants_derived

 if (.not. isinitialized) then
    print *,'specmod not initialized!'
    stop
 endif

 ! compute taper function for PBL (all flow is unbalanced when taper_pbl -> 0)
 nlevp = 0 ! top of pressure domain
 do k=1,nlevs+1
     if (bk(k) .ge. bk_top) then
        if (bk(k) .le. bk_bot) then
            taper_pbl(k) = 1.-(bk(k)-bk_top)/(bk_bot-bk_top)
        else
            taper_pbl(k) = 0.
        endif
     else
        taper_pbl(k) = 1.
     endif
     if (bk(k) > 1.e-10 .and. nlevp == 0) then
          nlevp = k-1
          ptop = ak(k-1)
     endif
     if (nmem .eq. 1) print *,'k,bk,taper_bpl',k,bk(k),taper_pbl(k)
 enddo
 !if (nmem .eq. 1) print *,'ptop,nlevp = ',ptop,nlevp,ak(nlevp)

 ! compute global mean tv, ps
 do k=1,nlevs
    tvmean(k) = sum(tvg(:,k)*areameanwts)
    !if (nmem .eq. 1) print*,k,tvmean(k)
 enddo
 psmean = sum(psg*areameanwts)

 ! compute abs vorticity on grid, rotational wind
 !     idir     - integer transform flag
 !                (idir>0 for wave to grid, idir<0 for grid to wave)
 do k=1,nlevs
    call sptezv_s(divspec,vrtspec,ug(:,k),vg(:,k),-1)
    call sptez_s(vrtspec,vrtg(:,k),1)
    divspec(:) = 0
    call sptezv_s(divspec,vrtspec,upsig(:,k),vpsig(:,k),1)
    nn = 0
    do j=1,nlats
    do i=1,nlons
       nn = nn + 1
       vrtg(nn,k) = vrtg(nn,k) + 2.*omega*gaulats(j)
    enddo
    enddo
    ug2 = upsig(:,k)*vrtg(:,k); vg2 = vpsig(:,k)*vrtg(:,k)
    call sptezv_s(divspec,vrtspec,ug2,vg2,-1)
    ug2 = 0.5*(upsig(:,k)**2 + vpsig(:,k)**2)
    call sptez_s(divspec,ug2,-1)
    Gspec(:,k) = vrtspec - lap*divspec ! RHS of balance eqn
    call sptez_s(Gspec(:,k),vg2,1)
    !if (nmem .eq. 1) print *,k,minval(vrtg(:,k)),maxval(vrtg(:,k)),minval(vg2),maxval(vg2)
 enddo 

! for testing, return zero tvbal, psbal (plus rotational wind)
 if (nobal) then
    psgbal= 0.
    tvgbal = 0.
    return
 endif
    

 ! initial guess of balanced state
 psgbal = psg_ref
 tvgbal = tvg_ref
 ug2 = log(psg_ref)
 call sptez_s(lnpsspec,ug2,-1)
 call getgrad(lnpsspec, lnpsrx, lnpsry)
 !if (nmem .eq. 1) print *,'lnprsx,lnprsy',minval(lnpsrx),maxval(lnpsry)

 ! determine balanced ps iteratively by solving NBE at lowest model level
 !if (nmem .eq. 1) print*,'min/max psgbal',minval(psgbal),maxval(psgbal)
 !if (nmem .eq. 1) print *,'psmean,tvmean',psmean,tvmean(nlevs)
 if (nmem .eq. 1) print *,'psbal iteration:'
 do niter=1,6
     ug2 = log(psgbal)
     call sptez_s(lnpsspec,ug2,-1)
     call getgrad(lnpsspec, lnpsx, lnpsy)
     ug2 = rd*(tvg(:,nlevs)-tvmean(nlevs))*lnpsx - rd*tvg_ref(:,nlevs)*lnpsrx
     vg2 = rd*(tvg(:,nlevs)-tvmean(nlevs))*lnpsy - rd*tvg_ref(:,nlevs)*lnpsry
     call sptezv_s(divspec,vrtspec,ug2,vg2,-1)
     lnpsspec = invlap*(Gspec(:,nlevs) - divspec)/(rd*tvmean(nlevs)) ! new estimate
     ug2 = psgbal
     call sptez_s(lnpsspec,psgbal,1)
     psgbal = exp(psgbal + log(psmean))
     if (nmem .eq. 1) print *,niter,minval(psgbal),maxval(psgbal),minval(psgbal-ug2),maxval(psgbal-ug2)
 enddo

 ! now determine balanced geopot
 do k=1,nlevs+1
    pressi(:,k) =  ak(k) + bk(k)*psg
    !pressi(:,k) =  ak(k) + bk(k)*psgbal
    phigbali(:,k) =  ak(k) + bk(k)*psg_ref
 enddo
 if (nmem .eq. 1) print *,'pressi',minval(pressi),maxval(pressi)
 if (nmem .eq. 1) print *,'pressir',minval(phigbali),maxval(phigbali)
 if (nmem .eq. 1) print *,'phibal (min/max):'
 do k=1,nlevs
    vg2 = 0.5*(log(pressi(:,k+1))+log(pressi(:,k))) ! press on model levels
    call sptez_s(lnpsspec,vg2,-1)
    call getgrad(lnpsspec, lnpsx, lnpsy)
    !print *,minval(lnpsx),maxval(lnpsx),minval(lnpsy),maxval(lnpsy)
    vg2 = 0.5*(log(phigbali(:,k+1))+log(phigbali(:,k))) ! ref press on model levels
    call sptez_s(lnpsspec,vg2,-1)
    call getgrad(lnpsspec, lnpsrx, lnpsry)
    ug2 = rd*tvg(:,k)*lnpsx - rd*tvg_ref(:,k)*lnpsrx
    vg2 = rd*tvg(:,k)*lnpsy - rd*tvg_ref(:,k)*lnpsry
    call sptezv_s(divspec,vrtspec,ug2,vg2,-1)
    divspec = invlap*(Gspec(:,k) - divspec)
    call sptez_s(divspec,phigbal(:,k),1)
    if (nmem .eq. 1) print *,k,minval(phigbal(:,k)),maxval(phigbal(:,k))
 enddo

 ! infer balanced tv from balanced geopot
 do k=1,nlevs+1
    pressi(:,k) =  ak(k) + bk(k)*psg
 enddo
 phigbali(:,1) = phigbal(:,1)
 phigbali(:,nlevs+1) = phigbal(:,nlevs)
 do k=2,nlevs
    phigbali(:,k) = 0.5*(phigbal(:,k)+phigbal(:,k-1)) ! pressi now holds geopot
 enddo
 if (nmem .eq. 1) print *,'tv (min/max), tvbal (min/max):'
 do k=1,nlevs
    ug2 = log(pressi(:,k+1)/pressi(:,k)) ! dlnp
    vg2 = (phigbali(:,k)-phigbali(:,k+1))/(rd*ug2)
    vg2 = vg2 - sum(vg2*areameanwts)
    tvgbal(:,k) = tvg_ref(:,k)+taper_pbl(k)*vg2
    if (nmem .eq. 1) print *,k,minval(tvg(:,k)),maxval(tvg(:,k)),minval(tvgbal(:,k)),maxval(tvgbal(:,k))
 enddo

 !deallocate(pressi,phigbali,Gspec,vrtg,phigbal)

 end subroutine getbaltps_fromvrt

end module enkf_nlbalmod
