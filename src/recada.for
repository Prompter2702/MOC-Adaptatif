      SUBROUTINE SWEE_TOT_3D_ADAPTIVE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,nb_niv,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir,
       ! inputs: cross section  
     &                          sigt,sigs,
       ! inputs: source and flux moments at previous iteration 
     &                          srcm,asrc,
       ! input/output: boundary flux 
     &                          bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                          mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                          sgnc,sgni,sgne,sgnt,
     &                          ccofglb,icofglb,ecofglb,tcofglb,
       ! input: auxiliar memory for boundary angular fluxes, angular 
       ! source moments and self-scattering xs 
     &                          sigg, tmom,
       ! input: auxiliar memory angular flux & source 
     &                          dsrc, aflx,
       ! input: auxiliar memory to drive angular mirror-reflection 
       ! or rotation/translation b.c.
     &                          rdir,
       ! input: pixel-to-medium array 
     &                          zreg, 
       ! input: pixel-to-medium array 
     &                          dira,dirf,
       ! input: logical to add angular source in case of
       ! time-dependent calculation 
     &                          lgki,
       ! output: flux moments 
     &                          flxm,
       ! Delta on each direction of the region
     &                          delt3,
      ! max inner iterations and tolerance for inner iterations
     &                          maxinner, tolinner, tolcor,
     &                          aflxmean,
     &                          pdslu4,
     &                          aflx0, aflx1,
     &                          asrcm0, asrcm1,
     &                          finc1, fout1,
     &                          finc0, fout0, lastadd)
   
      USE FLGOCT
      USE FLGBCD
      USE FLGINC
      USE SWEEP8ONE
      USE SRCCORR
      
      IMPLICIT NONE
   
      INTEGER, PARAMETER :: noct = 8, n8=8, ns=4
      INTEGER, INTENT(IN) :: nn,ng,nr,nh,nc,nx,ny,nz,nmat,nb_niv
      INTEGER, INTENT(IN) :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nani,nhrm,nd,ndir
      REAL, INTENT(INOUT) :: flxm(ng,nr,nh,nc,2)
      REAL, INTENT(INOUT) :: bflx(nn,nb,nbfx,2,noct),
     &                       bfly(nn,nb,nbfy,2,noct),
     &                       bflz(nn,nb,nbfz,2,noct)
      REAL, INTENT(IN)    :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),
     &                            w(ndir),pisn
      REAL, INTENT(IN)    :: sphr(nhrm,nd)
      INTEGER,INTENT(IN) :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER,INTENT(IN) :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(IN)    :: srcm(ng,nr,nh,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc,noct)
      REAL, INTENT(INOUT) :: aflxmean(nn,nr,noct)
   
      REAL, INTENT(INOUT)    :: asrc(nn,nr,nc,noct)
      REAL, INTENT(INOUT)    :: dsrc(nn,nr,nc,*)
      REAL, INTENT(INOUT)    :: tmom(ng,nr,nh,nc)
      REAL, INTENT(INOUT)    :: sigg(ng,nr)
      INTEGER, INTENT(INOUT) :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL, INTENT(INOUT) :: lgki
      REAL, INTENT(IN)       :: delt3(3)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir) ! pds to compute src 4 \int(L^4)
   
      INTEGER :: xout,yout,zout,fst
      INTEGER :: oc,oct,d,da,c,h
   
      REAL, INTENT(INOUT)  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT)  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3), fout0(nn,nb,3)
   
      REAL, INTENT(IN)  :: ccofglb(nn,nc, nc, nmat,nb_niv),
     &                        icofglb(nn,nc, nbd,nmat,nb_niv)
      REAL, INTENT(IN)  :: ecofglb(nn,nbd,nc,nmat,nb_niv),
     &                        tcofglb(nn,nbd,nbd,nmat,nb_niv)

      REAL  :: errbnd(3)
      REAL(KIND=8)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL  :: ok, okinner
      INTEGER  :: ii,jj,kk,rr 
      INTEGER , INTENT(IN) :: maxinner
      REAL, INTENT(IN) :: tolinner, tolcor
      REAL :: errinner
      INTEGER :: cnt,g,r,nb_cell,a,b,inew,iold, addrflx(2)
      INTEGER, PARAMETER :: olda=1, newa=2
      INTEGER, INTENT(INOUT) :: lastadd
      REAL(KIND=8) :: max_err_inner, tol_tmp
      REAL(KIND=8) ::  errtot
      LOGICAL :: oksrc
   
      REAL :: xstlvl1(ng,n8)

      cnt = 0
      addrflx = (/olda, newa/)
      xstlvl1 = 0.0
      
      okinner = .FALSE.
      errbnd = 0.0
      ! Initialization of the coefficients to compute the 
   
      ! Internal iterations loop
      DO WHILE (.NOT. okinner .AND. cnt<maxinner)
      cnt = cnt + 1
      inew = addrflx(newa)
      iold = addrflx(olda)

   !     Adds the contribution of selfscattering sources "sigs*flxp"
   !     to source moments "tmom"
      CALL GMOM3D(ng,nani,nh,nc,nr,sigs,flxm(1,1,1,1,iold),
     &            srcm,zreg,sigg,tmom)
    
      flxm(:,:,:,:,inew) = 0.0
   
      DO oc=1,8
        ! Loop over octants of angular space in order defined by
        ! octant list "olst3D",
        oct=olst3D(oc)
   
   !       Sweep for all directions in octant "oct".
   !       initial direction of the octant 
        da = (oct-1)*ndir+1
   
   !       Define the directional sources in asrc from tmom and sphr
   
        CALL GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr(1,da),
     &               tmom(1,1,1,1),asrc(1,1,1,oct))
        IF(lgki)CALL SPXPY(nr*nc,dsrc(1,1,1,oct),asrc(1,1,1,oct))
        
   !        Index of outgoing side.
         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
   
         bflx(:,:,:,xout,oct) = bflx(:,:,:,xinc(oct),oct)
         bfly(:,:,:,yout,oct) = bfly(:,:,:,yinc(oct),oct)
         bflz(:,:,:,zout,oct) = bflz(:,:,:,zinc(oct),oct)
         
   !        Sweep for all directions in octant "oct".
   !        initial direction of the octant 
         da = (oct-1)*ndir+1
   
        ! The boundary entries are computed on the exiting part
        ! in order to work only on the exiting part of the bfl
   
   !     Level 0 needs to already be computed before entering recursive 
   !     Function
   
        CALL SRCHOMO0(nn,nr,nx,ny,
     &                 1,nx,1,ny,1,nz,
     &                 zreg,asrc(:,:,:,oct),asrcm0)

        CALL MERGEBOUND0(nn,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx(:,:,:,xinc(oct),oct),
     &                   bfly(:,:,:,yinc(oct),oct),
     &                   bflz(:,:,:,zinc(oct),oct),
     &                   finc0)

         
   
        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                sgnc(:,:,oct),sgni(:,:,oct),
     &                sgne(:,:,oct),sgnt(:,:,oct),   
     &                ccofglb(:,:,:,1,1),icofglb(:,:,:,1,1),
     &                ecofglb(:,:,:,1,1),tcofglb(:,:,:,1,1) )


        CALL ERRSURF(nn,nx,ny,nz,
     &       1,nx,1,ny,1,nz,
     &       finc0, bflx,bfly,bflz,tcofglb(:,:,:,1,1) ,errbnd)

     
       ! RECURIVE SWEEP
        nb_cell = 0
        
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                         nx,ny, nz,nmat, nb_niv,
     &                         xinc(oct), yinc(oct), zinc(oct),oct,
     &                         1,nx,1,ny,1,nz,0,
     &                         sigt,xstlvl1, 
     &                         mu,eta,ksi,w,pdslu4,
     &                         sgnc(:,:,oct),sgni(:,:,oct),
     &                         sgne(:,:,oct),sgnt(:,:,oct),
     &                         zreg, asrc(:,:,:,oct),
     &                         aflx(:,:,:,oct),
     &                         bflx(1,1,1,xout,oct),
     &                         bfly(1,1,1,yout,oct),
     &                         bflz(1,1,1,zout,oct),
     &                         aflx0, aflx1,
     &                         asrcm0, asrcm1,finc0,finc1,
     &                         fout0,fout1,
     &                         ccofglb,icofglb,ecofglb,tcofglb,
     &                         ccof8, icof8,ecof8, tcof8,
     &                         tolcor,errcor,errmul,ok,
     &                         delt3,ii,jj,kk,rr,
     &                         nb_cell,a,b,oksrc,errtot,
     &                         aflxmean(1,1,oct), errbnd )



      print *, xinc(oct), yinc(oct), zinc(oct)
      print *,"nb_cell", nb_cell
   
      
    !   CALL BCD2R2(oct,xout,yout,bflx,bfly,mu,eta,w,pisn,rdir)
   
   !  Compute the flux moments
      DO c=1,nc
      DO h=1,nh
      fst = 0
      DO d = 1,ndir
        flxm(:ng,:nr,h,c,inew) = 
     &    flxm(:ng,:nr,h,c,inew) +  
     &    (sphr(h,da+d-1)*w(d))*aflx(fst+1:fst+ng,:nr,c,oct)
        fst = fst + ng
      ENDDO
      ENDDO
      ENDDO
   
   
      ENDDO ! octant loop
   
   !  Compute the inner error
      okinner = .TRUE.
   
      max_err_inner = 0.0
      tol_tmp = 0.0 
   
      spat :   DO c=1,nc
      harm :   DO h=1,nh
      region : DO r = 1,nr
      groupe : DO g = 1,ng
        errinner = ABS(flxm(g,r,h,c,iold) - flxm(g,r,h,c,inew))
        IF ( errinner + epsilon(1.0) 
     &   > tolinner * ABS( flxm(g,r,h,c,iold) )   )THEN
          okinner = .FALSE.
          EXIT spat
        ENDIF
      ENDDO groupe
      ENDDO region
      ENDDO harm
      ENDDO spat
    
      addrflx = CSHIFT(addrflx, 1)
      print *,"Error inner", errinner
   
   ! End inner iterations
      ENDDO

      lastadd = addrflx(olda)
   
      print *,"Number iteration", cnt
    !   print *,"Error inner", errinner
      
      END


!----------------------------------------------------------------------
!----------------------------------------------------------------------

      
      RECURSIVE SUBROUTINE RECADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &                             imin,imax,jmin,jmax,kmin,kmax,niv,
     &                             sigt, xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,ecofglb,tcofglb,
     &                             ccof8, icof8,ecof8, tcof8,
     &                             tolcor, errcor,errmul,
     &                             ok,delt3,i,j,k,r,
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      USE SWEEP8ONE
      USE SRCCORR 

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh,
     &           imin,imax,jmin,jmax,kmin,kmax,
     &           nx,ny,nz,niv,oct,nmat,nb_niv
      INTEGER, INTENT(IN) :: xinc,yinc,zinc
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)

      REAL, INTENT(IN)    ::  tolcor
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),bfly(nn,nb,nx*nz),
     &           bflz(nn,nb,nx*ny)
     
      REAL, INTENT(INOUT)  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT)  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3),    fout0(nn,nb,3)
      REAL, INTENT(INOUT)  :: aflxmean(nn,nr)

      REAL, INTENT(INOUT) :: xstlvl1(ng,n8)

      REAL, INTENT(INOUT) :: ccof8(nn,nc ,nc), icof8(nn,nc ,nbd),
     &                       ecof8(nn,nbd,nc), tcof8(nn,nbd,nbd)
      
      REAL, INTENT(IN) ::    ccofglb(nn,nc ,nc ,nmat,nb_niv),
     &                       icofglb(nn,nc ,nbd,nmat,nb_niv),
     &                       ecofglb(nn,nbd,nc ,nmat,nb_niv),
     &                       tcofglb(nn,nbd,nbd,nmat,nb_niv)

      REAL(KIND=8),INTENT(INOUT)  :: errcor(ng, ndir),errmul(ng, ndir)
      LOGICAL, INTENT(INOUT)  :: ok, oksrc
      REAL(KIND=8), INTENT(INOUT)  :: errtot
      REAL, INTENT(IN)        :: delt3(3)
      INTEGER, INTENT(INOUT)  :: i,j,k,r, nb_cell,g,cnt
      REAL, INTENT(INOUT) :: errbnd(3)

      !Variables locales 

      INTEGER :: order_ijk_cell(6,n8)
      INTEGER :: order_cell(n8)
      
      REAL :: asrcmtps(ng,ndir,nc,n8)
      REAL :: finc0tps(nn,nb,nb), finc1tps(nn,nb,nb,ns)

      fout1 = 0.0
      order_ijk_cell = 0
      order_cell = 0

      CALL SRCCOR(ng,ndir,delt3/(2**niv),mu,eta,ksi,asrcm0,
     &       sigt(ng,zreg(((imin-1)*ny + (jmin-1))*nx + kmin)),
     &       errcor,pdslu4)
      
      CALL FILLXSLVL1(ng,nx,ny,
     &                imin,imax,jmin,jmax,kmin,kmax,
     &                zreg,xstlvl1,sigt)
      
      errbnd = errbnd*(delt3(3)**2)/(4**niv)

      
      oksrc = .TRUE.
      errtot = 0.0
      
      cnt = 0
      r   = ((imin-1)*ny + (jmin-1))*nx + kmin
      ido  : DO i= kmin,kmax
      DO j= jmin,jmax
      DO k= imin,imax
         cnt = ((i-1)*ny + (j-1))*nx + k 
         if ( zreg(cnt) .NE. zreg(r)) THEN
             oksrc = .FALSE.
             EXIT ido
         ENDIF
      ENDDO
      ENDDO
      ENDDO ido

      errbnd = 0.0

      cnt = 0
      group : DO i = 1, ng
      direc : DO j = 1, ndir
        cnt = cnt + 1
        errtot = ABS( errcor(i,j)) + SUM(errbnd)/3
    !     IF ( errtot  > 3*epsilon(1.0) +
    !  &      tolcor*SUM( ABS(fout0(cnt,1,:)), DIM=1)/3  )THEN
    !         oksrc = .FALSE.
    !         EXIT direc
        IF ( errtot  > 3*epsilon(1.0) + tolcor*1.0  )THEN
            oksrc = .FALSE.
            EXIT group
        ENDIF
      ENDDO direc
      ENDDO group


      IF(imax-imin>0 .AND. jmax-jmin>0 
     & .AND. kmax-kmin>0 .AND.  .NOT. oksrc) THEN
        ok = .FALSE.
      ELSE
        ok = .TRUE.
      ENDIF

      IF (.NOT. ok) THEN 
! If the criterion is not satisfied, compute on lvl 1

        CALL FILL_ORDER_CELL(imin, imax, jmin, jmax, kmin, kmax, nx,ny,
     &                       xinc, yinc, zinc,order_cell,order_ijk_cell) 

        CALL SRCHOMO1(nn,nr,nx,ny,
     &                 imin, imax,jmin,jmax,kmin,kmax,
     &                 zreg, asrc,asrcm1)

        CALL MERGEBOUND1(nn,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, finc1)

        CALL ASSIGN_EIGHT_COEF3D(nn,nx,ny,
     &                    order_ijk_cell, nmat,zreg,
     &                    sgnc,sgni,sgne,sgnt,
     &                    ccofglb(:,:,:,:,niv+2),icofglb(:,:,:,:,niv+2),
     &                    ecofglb(:,:,:,:,niv+2),tcofglb(:,:,:,:,niv+2),
     &                    ccof8,icof8,ecof8,tcof8)   
     
        CALL SWEEP_8REGIONS(nn,2,asrcm1, finc1, aflx1, fout1,
     &                     ccof8,icof8,ecof8,tcof8, xinc, yinc, zinc)

        ! read(*,*)
        ! CALL SRC2LVL(ng, ndir, asrcm1, sigt, delt3/(2**niv), errmul)
        ! errtot = 0.0
        ! DO i = 1, ng
        !     DO j = 1, ndir
        !         IF ( errtot < ABS(errmul(i,j)) ) THEN
        !             errtot = ABS(errmul(i,j))
        !         ENDIF
        !     ENDDO
        ! ENDDO
        ! errtot = SQRT(errtot/norme)
        ! IF(errtot > tolcor/(2**niv)) THEN
        ! !    print *, "Finalement pas ok"
        !    ok = .FALSE.
        !  ELSE
        !    print *,"Finalement ok"
        !    ok = .TRUE.
        ! ENDIF

      ENDIF
      
      IF (ok) THEN

! If the criterion is satisfied, compute on lvl 0 and 
        IF (niv > 0) THEN
            CALL MERGEBOUND0(nn,nb,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx,
     &                   bfly,
     &                   bflz,
     &                   finc0)

        r = zreg(((imin-1)*ny + (jmin-1))*nx + kmin)
        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                   sgnc,sgni,sgne,sgnt,    
     &                   ccofglb(:,:,:,r,niv+1),icofglb(:,:,:,r,niv+1),
     &                   ecofglb(:,:,:,r,niv+1),tcofglb(:,:,:,r,niv+1))
        END IF 
        
        ! - write the angular flux of node-0 on the fine-mesh flux 
        CALL SPLITAFLX0(nn,nr,nc,nx,ny,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   aflx,aflx0)
!      - write the outgoing angular flux of node-0 
!            on the fine-mesh boundary flux
         CALL SPLITBOUND0(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)
        
         !Aflxmean
         DO i=kmin,kmax
         DO j= jmin,jmax
         DO k= imin,imax
            r = ((i-1)*ny + (j-1))*nx + k 
            aflxmean(:,r) = aflx0(:,1)
         ENDDO
         ENDDO
         ENDDO

         nb_cell = nb_cell + 1
         RETURN

      ELSE IF((imax-imin)>1) THEN
!         15/ si ko : 
!        - call recursively the subroutine for the 8 nodes of level-1

        asrcmtps(:,:,:,1) = asrcm1(:,:,:, xinc+2*(yinc-1) + 4*(zinc-1))  
        asrcmtps(:,:,:,2) = asrcm1(:,:,:,3-xinc+2*(yinc-1) + 4*(zinc-1))  
        asrcmtps(:,:,:,3)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(zinc-1))  
        asrcmtps(:,:,:,4)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))  
        asrcmtps(:,:,:,5)= asrcm1(:,:,:,xinc + 2*(yinc-1) + 4*(2-zinc))  
        asrcmtps(:,:,:,6)=asrcm1(:,:,:,3-xinc + 2*(yinc-1) + 4*(2-zinc))  
        asrcmtps(:,:,:,7)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))  
        asrcmtps(:,:,:,8)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(2-zinc))  
        
      finc1tps = finc1

      finc0tps(:,:,1) = finc1tps(:,:,1, yinc + 2*(zinc-1))
      finc0tps(:,:,2) = finc1tps(:,:,2, xinc + 2*(zinc-1))
      finc0tps(:,:,3) = finc1tps(:,:,3, xinc + 2*(yinc-1))
      
      CALL ERRSURF(nn,nx,ny,nz,
     &   order_ijk_cell(1,1),order_ijk_cell(2,1),order_ijk_cell(3,1),
     &   order_ijk_cell(4,1),order_ijk_cell(5,1),order_ijk_cell(6,1),
     &   finc0, bflx,bfly,bflz,
     &   tcofglb(:,:,:,zreg(order_cell(1)),niv+1),errbnd)


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &   order_ijk_cell(1,1),order_ijk_cell(2,1),order_ijk_cell(3,1),
     &   order_ijk_cell(4,1),order_ijk_cell(5,1),order_ijk_cell(6,1),
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,1),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)

      finc1tps(:,:,1, yinc + 2*(zinc-1)) = fout0(:,:,1)
      finc1tps(:,:,2, xinc + 2*(zinc-1)) = fout0(:,:,2)
      finc1tps(:,:,3, xinc + 2*(yinc-1)) = fout0(:,:,3)


      finc0tps(:,:,1) = finc1tps(:,:,1, yinc + 2*(zinc-1))
      finc0tps(:,:,2) = finc1tps(:,:,2, 3-xinc+ 2*(zinc-1))
      finc0tps(:,:,3) = finc1tps(:,:,3, 3-xinc+ 2*(yinc-1))
      CALL ERRSURF(nn,nx,ny,nz,
     &   order_ijk_cell(1,2),order_ijk_cell(2,2),order_ijk_cell(3,2),
     &   order_ijk_cell(4,2),order_ijk_cell(5,2),order_ijk_cell(6,2), 
     &   finc0, bflx,bfly,bflz,
     &   tcofglb(:,:,:,zreg(order_cell(2)),niv+1),errbnd)     

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &   order_ijk_cell(1,2),order_ijk_cell(2,2),order_ijk_cell(3,2),
     &   order_ijk_cell(4,2),order_ijk_cell(5,2),order_ijk_cell(6,2), 
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,2),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      finc1tps(:,:,1, yinc + 2*(zinc-1))  = fout0(:,:,1)
      finc1tps(:,:,2, 3-xinc+ 2*(zinc-1)) = fout0(:,:,2)
      finc1tps(:,:,3, 3-xinc+ 2*(yinc-1)) = fout0(:,:,3)

      finc0(:,:,1) = finc1(:,:,1, 3-yinc + 2*(zinc-1))
      finc0(:,:,2) = finc1(:,:,2, xinc + 2*(zinc-1))
      finc0(:,:,3) = finc1(:,:,3, xinc + 2*(2-yinc))
      CALL ERRSURF(nn,nx,ny,nz,
     &    order_ijk_cell(1,3),order_ijk_cell(2,3),order_ijk_cell(3,3),
     &    order_ijk_cell(4,3),order_ijk_cell(5,3),order_ijk_cell(6,3), 
     &    finc0, bflx,bfly,bflz,
     &    tcofglb(:,:,:,zreg(order_cell(3)),niv+1),errbnd)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &    order_ijk_cell(1,3),order_ijk_cell(2,3),order_ijk_cell(3,3),
     &    order_ijk_cell(4,3),order_ijk_cell(5,3),order_ijk_cell(6,3), 
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,3),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      finc1(:,:,1, 3-yinc + 2*(zinc-1)) = fout0(:,:,1)
      finc1(:,:,2, xinc + 2*(zinc-1)) = fout0(:,:,2)
      finc1(:,:,3, xinc + 2*(2-yinc)) = fout0(:,:,3)

      finc0(:,:,1) = finc1(:,:,1, 3-yinc + 2*(zinc-1))
      finc0(:,:,2) = finc1(:,:,2, 3-xinc + 2*(zinc-1))
      finc0(:,:,3) = finc1(:,:,3, 3-xinc + 2*(2-yinc))
      CALL ERRSURF(nn,nx,ny,nz,
     &    order_ijk_cell(1,4),order_ijk_cell(2,4),order_ijk_cell(3,4),
     &    order_ijk_cell(4,4),order_ijk_cell(5,4),order_ijk_cell(6,4), 
     &    finc0, bflx,bfly,bflz,
     &    tcofglb(:,:,:,zreg(order_cell(4)),niv+1),errbnd)             
      CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &    order_ijk_cell(1,4),order_ijk_cell(2,4),order_ijk_cell(3,4),
     &    order_ijk_cell(4,4),order_ijk_cell(5,4),order_ijk_cell(6,4), 
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,4),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      finc1(:,:,1, 3-yinc + 2*(zinc-1)) = fout0(:,:,1)
      finc1(:,:,2, 3-xinc + 2*(zinc-1)) = fout0(:,:,2)
      finc1(:,:,3, 3-xinc + 2*(2-yinc)) = fout0(:,:,3)
!-----------------------------------------------------------------------
      finc0(:,:,1) = finc1(:,:,1, yinc + 2*(2-zinc))
      finc0(:,:,2) = finc1(:,:,2, xinc + 2*(2-zinc))
      finc0(:,:,3) = finc1(:,:,3, xinc + 2*(yinc-1))
      CALL ERRSURF(nn,nx,ny,nz,
     &    order_ijk_cell(1,5),order_ijk_cell(2,5),order_ijk_cell(3,5),
     &    order_ijk_cell(4,5),order_ijk_cell(5,5),order_ijk_cell(6,5), 
     &    finc0, bflx,bfly,bflz,
     &    tcofglb(:,:,:,zreg(order_cell(5)),niv+1),errbnd) 
         CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &    order_ijk_cell(1,5),order_ijk_cell(2,5),order_ijk_cell(3,5),
     &    order_ijk_cell(4,5),order_ijk_cell(5,5),order_ijk_cell(6,5), 
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,5),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
       finc1(:,:,1, yinc + 2*(2-zinc)) = fout0(:,:,1)
       finc1(:,:,2, xinc + 2*(2-zinc)) = fout0(:,:,2)
       finc1(:,:,3, xinc + 2*(yinc-1)) = fout0(:,:,3)
    
       finc0(:,:,1) = finc1(:,:,1, yinc + 2*(2-zinc))
       finc0(:,:,2) = finc1(:,:,2, 3-xinc + 2*(2-zinc))
       finc0(:,:,3) = finc1(:,:,3, 3-xinc + 2*(yinc-1))
       CALL ERRSURF(nn,nx,ny,nz,
     &   order_ijk_cell(1,6),order_ijk_cell(2,6),order_ijk_cell(3,6),
     &   order_ijk_cell(4,6),order_ijk_cell(5,6),order_ijk_cell(6,6),  
     &   finc0, bflx,bfly,bflz,
     &   tcofglb(:,:,:,zreg(order_cell(6)),niv+1),errbnd) 
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &   order_ijk_cell(1,6),order_ijk_cell(2,6),order_ijk_cell(3,6),
     &   order_ijk_cell(4,6),order_ijk_cell(5,6),order_ijk_cell(6,6),  
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,6),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,   
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      finc1(:,:,1, yinc + 2*(2-zinc)) = fout0(:,:,1)
      finc1(:,:,2, 3-xinc + 2*(2-zinc)) = fout0(:,:,2)
      finc1(:,:,3, 3-xinc + 2*(yinc-1)) = fout0(:,:,3)

      finc0(:,:,1) = finc1(:,:,1, 3-yinc + 2*(2-zinc))
      finc0(:,:,2) = finc1(:,:,2, xinc + 2*(2-zinc))
      finc0(:,:,3) = finc1(:,:,3, xinc + 2*(2-yinc))
      CALL ERRSURF(nn,nx,ny,nz,
     &   order_ijk_cell(1,7),order_ijk_cell(2,7),order_ijk_cell(3,7),
     &   order_ijk_cell(4,7),order_ijk_cell(5,7),order_ijk_cell(6,7),   
     &   finc0, bflx,bfly,bflz,
     &   tcofglb(:,:,:,zreg(order_cell(7)),niv+1),errbnd) 

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &   order_ijk_cell(1,7),order_ijk_cell(2,7),order_ijk_cell(3,7),
     &   order_ijk_cell(4,7),order_ijk_cell(5,7),order_ijk_cell(6,7),
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,7),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,        
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)
      finc1(:,:,1, 3-yinc + 2*(2-zinc)) = fout0(:,:,1)
      finc1(:,:,2, xinc + 2*(2-zinc)) = fout0(:,:,2)
      finc1(:,:,3, xinc + 2*(2-yinc)) = fout0(:,:,3)

      finc0(:,:,1) = finc1(:,:,1, 3-yinc + 2*(2-zinc))
      finc0(:,:,2) = finc1(:,:,2, 3-xinc + 2*(2-zinc))
      finc0(:,:,3) = finc1(:,:,3, 3-xinc + 2*(2-yinc))

      CALL ERRSURF(nn,nx,ny,nz,
     &    order_ijk_cell(1,8),order_ijk_cell(2,8),order_ijk_cell(3,8),
     &    order_ijk_cell(4,8),order_ijk_cell(5,8),order_ijk_cell(6,8),
     &    finc0, bflx,bfly,bflz,
     &    tcofglb(:,:,:,zreg(order_cell(8)),niv+1),errbnd) 
     
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,nmat,nb_niv,
     &                             xinc,yinc,zinc, oct,
     &    order_ijk_cell(1,8),order_ijk_cell(2,8),order_ijk_cell(3,8),
     &    order_ijk_cell(4,8),order_ijk_cell(5,8),order_ijk_cell(6,8),
     &                             niv+1,
     &                             sigt,xstlvl1,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             asrcmtps(:,:,:,8),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccofglb,icofglb,
     &                             ecofglb,tcofglb,
     &                             ccof8,icof8,ecof8,tcof8,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,        
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean,
     &                             errbnd)

    !- write the outgoing angular flux of node-0 on the fine-mesh
    ! boundary flux
         RETURN
      
        ELSE 
!   Case where the error is not guaranted but the pixel-lvl is reached
            CALL SPLITAFLX1(nn,nr,nc,nx,ny,nz,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)

            CALL PROJMEAN1(nn,nr,nc,nx,ny,
     &                     imin,imax,jmin,jmax,kmin,kmax,
     &                     aflxmean, aflx1)

            CALL SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)

            CALL MERGEBOUND0(nn,nb,imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)
          

      ENDIF
      nb_cell = nb_cell + 8

      END SUBROUTINE RECADA_ONE_OCTANT


      SUBROUTINE FILL_ORDER_CELL(imin, imax, jmin, jmax, kmin, kmax,
     &                           nx,ny,xinc, yinc, zinc,
     &                           order_cell,order_ijk_cell) 


      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: imin, imax, jmin, jmax, kmin, kmax,
     &                        xinc, yinc, zinc, nx,ny
      INTEGER, INTENT(OUT) :: order_ijk_cell(6,8), order_cell(8)

      INTEGER :: l

      order_ijk_cell(1,1) = (2*imin + (xinc-1)*(imax-imin+2))/2 
      order_ijk_cell(2,1) = (imin + imax + (xinc-1)*(imax-imin))/2
      order_ijk_cell(3,1) = (2*jmin + (yinc-1)*(jmax-jmin+2))/2
      order_ijk_cell(4,1) = (jmin + jmax + (yinc-1)*(jmax-jmin))/2
      order_ijk_cell(5,1) = (2*kmin + (zinc-1)*(kmax-kmin+2))/2
      order_ijk_cell(6,1) = (kmin + kmax + (zinc-1)*(jmax-jmin))/2

      order_ijk_cell(1,2) = (2*imin + (2-xinc)*(imax-imin+2))/2
      order_ijk_cell(2,2) = (imin + imax + (2-xinc)*(imax-imin))/2
      order_ijk_cell(3,2) = (2*jmin + (yinc-1)*(jmax-jmin+2))/2
      order_ijk_cell(4,2) = (jmin + jmax + (yinc-1)*(jmax-jmin))/2
      order_ijk_cell(5,2) = (2*kmin + (zinc-1)*(kmax-kmin+2))/2
      order_ijk_cell(6,2) = (kmin + kmax + (zinc-1)*(jmax-jmin))/2

      order_ijk_cell(1,3) = (2*imin + (xinc-1)*(imax-imin+2))/2 
      order_ijk_cell(2,3) = (imin + imax + (xinc-1)*(imax-imin))/2
      order_ijk_cell(3,3) = (2*jmin + (2-yinc)*(jmax-jmin+2))/2
      order_ijk_cell(4,3) = (jmin + jmax + (2-yinc)*(jmax-jmin))/2
      order_ijk_cell(5,3) = (2*kmin + (zinc-1)*(kmax-kmin+2))/2
      order_ijk_cell(6,3) = (kmin + kmax + (zinc-1)*(jmax-jmin))/2

      order_ijk_cell(1,4) = (2*imin + (2-xinc)*(imax-imin+2))/2
      order_ijk_cell(2,4) = (imin + imax + (2-xinc)*(imax-imin))/2
      order_ijk_cell(3,4) = (2*jmin + (2-yinc)*(jmax-jmin+2))/2
      order_ijk_cell(4,4) = (jmin + jmax + (2-yinc)*(jmax-jmin))/2
      order_ijk_cell(5,4) = (2*kmin + (zinc-1)*(kmax-kmin+2))/2
      order_ijk_cell(6,4) = (kmin + kmax + (zinc-1)*(jmax-jmin))/2

      order_ijk_cell(1,5) = (2*imin + (xinc-1)*(imax-imin+2))/2 
      order_ijk_cell(2,5) = (imin + imax + (xinc-1)*(imax-imin))/2
      order_ijk_cell(3,5) = (2*jmin + (yinc-1)*(jmax-jmin+2))/2
      order_ijk_cell(4,5) = (jmin + jmax + (yinc-1)*(jmax-jmin))/2
      order_ijk_cell(5,5) = (2*kmin + (2-zinc)*(kmax-kmin+2))/2
      order_ijk_cell(6,5) = (kmin + kmax + (2-zinc)*(kmax-kmin))/2

      order_ijk_cell(1,6) = (2*imin + (2-xinc)*(imax-imin+2))/2
      order_ijk_cell(2,6) = (imin + imax + (2-xinc)*(imax-imin))/2
      order_ijk_cell(3,6) = (2*jmin + (yinc-1)*(jmax-jmin+2))/2
      order_ijk_cell(4,6) = (jmin + jmax + (yinc-1)*(jmax-jmin))/2
      order_ijk_cell(5,6) = (2*kmin + (2-zinc)*(kmax-kmin+2))/2
      order_ijk_cell(6,6) = (kmin + kmax + (2-zinc)*(kmax-kmin))/2

      order_ijk_cell(1,7) = (2*imin + (xinc-1)*(imax-imin+2))/2 
      order_ijk_cell(2,7) = (imin + imax + (xinc-1)*(imax-imin))/2
      order_ijk_cell(3,7) = (2*jmin + (2-yinc)*(jmax-jmin+2))/2
      order_ijk_cell(4,7) = (jmin + jmax + (2-yinc)*(jmax-jmin))/2
      order_ijk_cell(5,7) = (2*kmin + (2-zinc)*(kmax-kmin+2))/2
      order_ijk_cell(6,7) = (kmin + kmax + (2-zinc)*(kmax-kmin))/2

      order_ijk_cell(1,8) = (2*imin + (2-xinc)*(imax-imin+2))/2
      order_ijk_cell(2,8) = (imin + imax + (2-xinc)*(imax-imin))/2
      order_ijk_cell(3,8) = (2*jmin + (2-yinc)*(jmax-jmin+2))/2
      order_ijk_cell(4,8) = (jmin + jmax + (2-yinc)*(jmax-jmin))/2
      order_ijk_cell(5,8) = (2*kmin + (2-zinc)*(kmax-kmin+2))/2
      order_ijk_cell(6,8) = (kmin + kmax + (2-zinc)*(kmax-kmin))/2


      DO l=1,8
         order_cell(l) = order_ijk_cell(1,l) 
     &    + nx*(order_ijk_cell(3,l)-1) + nx*ny*(order_ijk_cell(5,l)-1)
      ENDDO
    

      END SUBROUTINE FILL_ORDER_CELL

 