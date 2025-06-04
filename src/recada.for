       SUBROUTINE SWEE3D_ADAPTIVE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir,
       ! inputs: cross section  
     &                          sigt,sigs,
       ! inputs: source and flux moments at previous iteration 
    !  &                          flxp,srcm,
     &                          srcm,
       ! input/output: boundary flux 
     &                          bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                          mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                          sgnc,sgni,sgne,sgnt,
       ! input: auxiliar memory for boundary angular fluxes, angular 
       ! source moments and self-scattering xs 
     &                          sigg, !tmom,
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
     &                          aflxmean)

      USE FLGOCT
      USE FLGBCD
      USE FLGINC
      USE SWEEP8ONE
      USE SRCCORR
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: noct = 8, n8=8, ns=4
      
      INTEGER, INTENT(IN) :: nn,ng,nr,nh,nc,nx,ny,nz,nmat
      INTEGER, INTENT(IN) :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nani,nhrm,nd,ndir
      REAL, INTENT(INOUT) :: flxm(ng,nr,nh,nc)
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

      REAL    :: asrc(nn,nr,nc,noct)
      REAL    :: flxp(ng,nr,nh,nc)
      REAL    :: dsrc(nn,nr,nc,*)
      REAL    :: tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      INTEGER :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL :: lgki
      REAL, INTENT(IN)     :: delt3(3)
      REAL(KIND=8) :: pdslu4(ndir) ! pds to compute src 4 \int(L^4)

      INTEGER :: xout,yout,zout,fst
      INTEGER :: oc,oct,d,dd,da,c,h,i

      REAL  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL  :: xshom0(ng), xshom1(ng,n8)
      REAL  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL  :: ok, okinner
      INTEGER  :: ii,jj,kk,rr 
      INTEGER , INTENT(IN) :: maxinner
      REAL, INTENT(IN) :: tolinner, tolcor
      REAL :: errinner, tmp
      INTEGER :: cnt,g,r,j, nb_cell, a,b
      REAL :: max_err_inner, tol_tmp, errtot
      LOGICAL :: oksrc

      INTEGER :: posmax(4)

      cnt = 0
      flxp = flxm
      okinner = .FALSE.

      ! Initialization of the coefficients to compute the 
      ! srccor error
      CALL GAUSS_LU4(100,ndir,(/1.0, 1.0, 1.0/), mu,eta,ksi, pdslu4)


      ! Internal iterations loop
      DO WHILE (.NOT. okinner .AND. cnt<maxinner)
      cnt = cnt + 1

!     Adds the contribution of selfscattering sources "sigs*flxp"
!     to source moments "tmom"
      CALL GMOM3D(ng,nani,nh,nc,nr,sigs,flxp,srcm,zreg,sigg,tmom)
    
      flxm = 0.0

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

        CALL XSSRCHOMO0(nn,ng,nr,nx,ny, ndir,
     &                 1,nx,1,ny,1,nz,
     &                 sigt, zreg, aflx(:,:,:,oct), w,
     &                 asrc(:,:,:,oct), xshom0,
     &                 asrcm0)

        CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                delt3,xshom0,
     &                sgnc(:,:,oct),sgni(:,:,oct),
     &                sgne(:,:,oct),sgnt(:,:,oct),
     &                ccof,icof,ecof,tcof)

        CALL MERGEBOUND0(nn,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx(:,:,:,xinc(oct),oct),
     &                   bfly(:,:,:,yinc(oct),oct),
     &                   bflz(:,:,:,zinc(oct),oct),
     &                   finc0)

        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)

       ! RECURIVE SWEEP
        nb_cell = 0
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                         nx,ny, nz,
     &                         xinc(oct), yinc(oct), zinc(oct),oct,
     &                         1,nx,1,ny,1,nz,
     &                         0,sigt, 
     &                         mu,eta,ksi,w,pdslu4,
     &                         sgnc(:,:,oct),sgni(:,:,oct),
     &                         sgne(:,:,oct),sgnt(:,:,oct),
     &                         zreg, asrc(:,:,:,oct),
     &                         aflx(:,:,:,oct),
     &                         bflx(:,:,:,xout,oct),
     &                         bfly(:,:,:,yout,oct),
     &                         bflz(:,:,:,zout,oct),
     &                         aflx0, aflx1,xshom0, xshom1,
     &                         asrcm0, asrcm1,finc0,finc1,
     &                         fout0,fout1,ccof,icof,ecof,tcof,
     &                         ccof8,icof8,ecof8,tcof8,
     &                         tolinner,tolcor,
     &                         errcor,errmul,ok,delt3,ii,jj,kk,rr,
     &                         nb_cell,a,b,oksrc,errtot,aflxmean)

      print *, xinc(oct), yinc(oct), zinc(oct)
      print *,"nb_cell", nb_cell

      
      CALL BCD2R2(oct,xout,yout,bflx,bfly,mu,eta,w,pisn,rdir)

!  Compute the flux moments
      DO c=1,nc
      DO h=1,nh
      fst = 0
      DO d = 1,ndir
        flxm(:ng,:nr,h,c) = flxm(:ng,:nr,h,c) + 
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
        errinner = ABS(flxp(g,r,h,c) - flxm(g,r,h,c)) 
        IF ( errinner  > tolinner * ABS( flxp(g,r,h,c) )
     &                   + epsilon(1.0) )THEN
          okinner = .FALSE.
          EXIT spat
        ENDIF
      ENDDO groupe
      ENDDO region
      ENDDO harm
      ENDDO spat
      
      flxp = flxm
    !   read(*,*)

! End inner iterations
      ENDDO

      print *,"Number iteration", cnt
    !   print *,"Error inner", errinner
      
      END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      RECURSIVE SUBROUTINE RECADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,
     &                             xinc,yinc,zinc, oct,
     &                             imin,imax,jmin,jmax,kmin,kmax,niv,
     &                             sigt, 
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1, xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)
     
      USE SWEEP8ONE
      USE SRCCORR 

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh,
     &           imin,imax,jmin,jmax,kmin,kmax,
     &           nx,ny,nz,niv,oct
      INTEGER, INTENT(IN) :: xinc,yinc,zinc
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(IN)    :: tol, tolcor
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),bfly(nn,nb,nx*nz),
     &           bflz(nn,nb,nx*ny)
     
      REAL, INTENT(INOUT)  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT)  :: xshom0(ng), xshom1(ng,n8)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT)  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3), fout0(nn,nb,3)
      REAL, INTENT(INOUT)  :: aflxmean(nn,nr,n8)


      REAL, INTENT(INOUT)  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL, INTENT(INOUT)  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8),INTENT(INOUT)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL, INTENT(IN)     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL, INTENT(IN)     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL, INTENT(INOUT)  :: ok, oksrc
      REAL, INTENT(INOUT)  :: errtot
      REAL, INTENT(IN)        :: delt3(3)
      INTEGER, INTENT(INOUT)  :: i,j,k,r, nb_cell,g,cnt


      !Variables locales 

      REAL :: asrcmtps(ng,ndir,nc,n8),xshomtps(ng,n8)
      REAL :: ccoftps(nn,nc,nc,n8),icoftps(nn,nc,nbd,n8)
      REAL :: ecoftps(nn,nbd,nc,n8),tcoftps(nn,nbd,nbd,n8)
      
      
      fout1 = 0.0

      CALL SRCCOR(ng,ndir,delt3/(2**niv),mu,eta,ksi,asrcm0,
     &            xshom0,errcor,pdslu4)

      oksrc = .TRUE.
      errtot = 0.0
      cnt = 0

      group1 : DO g=1,ng
      DO i=kmin,kmax
      DO j= jmin,jmax
      DO k= imin,imax
            r = ((i-1)*ny + (j-1))*nx + k 
            if (xshom0(g) .NE. sigt(g,zreg(r))) THEN
                oksrc = .FALSE.
                EXIT group1
            ENDIF
        ENDDO
        ENDDO
        ENDDO
      ENDDO group1


      group : DO i = 1, ng
      direc : DO j = 1, ndir
        cnt = cnt + 1
        errtot = ABS( errcor(i,j)) 
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
        CALL XSSRCHOMO1(nn,ng,nr,nx,ny, ndir,
     &                 imin, imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom1,
     &                 asrcm1)

        CALL MERGEBOUND1(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, finc1)

        CALL EIGHT_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                    delt3/(2**(niv+1)),xshom1,
     &                    sgnc,sgni,sgne,sgnt,
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

        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)
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
        
         !!!Aflxmean
         DO i=kmin,kmax
         DO j= jmin,jmax
         DO k= imin,imax
            r = ((i-1)*ny + (j-1))*nx + k 
            aflxmean(:,r,oct) = aflx0(:,1)
         ENDDO
         ENDDO
         ENDDO


          nb_cell = nb_cell + 1
          RETURN

      ELSE IF((imax-imin)>1) THEN
!         15/ si ko : 
!        - call recursively the subroutine for the 8 nodes of level-1

        asrcmtps(:,:,:,1) = asrcm1(:,:,:, xinc+2*(yinc-1) + 4*(zinc-1))  
        xshomtps(:,1) = xshom1(:,       xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps(:,:,:,2) = asrcm1(:,:,:,3-xinc+2*(yinc-1) + 4*(zinc-1))  
        xshomtps(:,2) = xshom1(:,      3-xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps(:,:,:,3)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(zinc-1))  
        xshomtps(:,3)= xshom1(:,      xinc+2*(2-yinc) + 4*(zinc-1))

        asrcmtps(:,:,:,4)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))  
        xshomtps(:,4)= xshom1(:,      3-xinc + 2*(2-yinc) + 4*(zinc-1))

        asrcmtps(:,:,:,5)= asrcm1(:,:,:,xinc + 2*(yinc-1) + 4*(2-zinc))  
        xshomtps(:,5)= xshom1(:,      xinc + 2*(yinc-1) + 4*(2-zinc))

        asrcmtps(:,:,:,6)=asrcm1(:,:,:,3-xinc + 2*(yinc-1) + 4*(2-zinc))  
        xshomtps(:,6) = xshom1(:,     3-xinc + 2*(yinc-1) + 4*(2-zinc))  

        asrcmtps(:,:,:,7)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))  
        xshomtps(:,7)= xshom1(:,      xinc+2*(2-yinc) + 4*(2-zinc))

        asrcmtps(:,:,:,8)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(2-zinc))  
        xshomtps(:,8)= xshom1(:,      3-xinc + 2*(2-yinc) + 4*(2-zinc))

        
      ccoftps(:,:,:,1) = ccof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      ecoftps(:,:,:,1) = ecof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      tcoftps(:,:,:,1) = tcof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      icoftps(:,:,:,1) = icof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))

      ccoftps(:,:,:,2) = ccof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      ecoftps(:,:,:,2) = ecof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      tcoftps(:,:,:,2) = tcof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      icoftps(:,:,:,2) = icof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))

      ccoftps(:,:,:,3) = ccof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      ecoftps(:,:,:,3) = ecof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      tcoftps(:,:,:,3) = tcof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      icoftps(:,:,:,3) = icof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))

      ccoftps(:,:,:,4) = ccof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      ecoftps(:,:,:,4) = ecof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      tcoftps(:,:,:,4) = tcof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      icoftps(:,:,:,4) = icof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))

      ccoftps(:,:,:,5) = ccof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      ecoftps(:,:,:,5) = ecof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      tcoftps(:,:,:,5) = tcof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      icoftps(:,:,:,5) = icof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))

      ccoftps(:,:,:,6) = ccof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      ecoftps(:,:,:,6) = ecof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      tcoftps(:,:,:,6) = tcof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      icoftps(:,:,:,6) = icof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))

      ccoftps(:,:,:,7) = ccof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      ecoftps(:,:,:,7) = ecof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      tcoftps(:,:,:,7) = tcof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      icoftps(:,:,:,7) = icof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))

      ccoftps(:,:,:,8) = ccof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      ecoftps(:,:,:,8) = ecof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      tcoftps(:,:,:,8) = tcof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      icoftps(:,:,:,8) = icof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,1),xshom1,
     &                             asrcmtps(:,:,:,1),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,1),icoftps(:,:,:,1),
     &                             ecoftps(:,:,:,1),tcoftps(:,:,:,1),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)
      
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,2),xshom1,
     &                             asrcmtps(:,:,:,2),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,2),icoftps(:,:,:,2),
     &                             ecoftps(:,:,:,2),tcoftps(:,:,:,2),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,3),xshom1,
     &                             asrcmtps(:,:,:,3),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,3),icoftps(:,:,:,3),
     &                             ecoftps(:,:,:,3),tcoftps(:,:,:,3),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)
        
      CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,4),xshom1,
     &                             asrcmtps(:,:,:,4),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,4),icoftps(:,:,:,4),
     &                             ecoftps(:,:,:,4),tcoftps(:,:,:,4),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)

!-----------------------------------------------------------------------

         CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,5),xshom1,
     &                             asrcmtps(:,:,:,5),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,5),icoftps(:,:,:,5),
     &                             ecoftps(:,:,:,5),tcoftps(:,:,:,5),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r, 
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,6),xshom1,
     &                             asrcmtps(:,:,:,6),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,6),icoftps(:,:,:,6),
     &                             ecoftps(:,:,:,6),tcoftps(:,:,:,6),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,   
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,7),xshom1,
     &                             asrcmtps(:,:,:,7),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,7),icoftps(:,:,:,7),
     &                             ecoftps(:,:,:,7),tcoftps(:,:,:,7),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,        
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0,aflx1,
     &                             xshomtps(:,8),xshom1,
     &                             asrcmtps(:,:,:,8),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,8),icoftps(:,:,:,8),
     &                             ecoftps(:,:,:,8),tcoftps(:,:,:,8),
     &                             ccof8,icof8,ecof8,tcof8,tol,tolcor,
     &                             errcor,errmul,ok,delt3,i,j,k,r,        
     &                             nb_cell,g,cnt,oksrc,errtot,aflxmean)
    !- write the outgoing angular flux of node-0 on the fine-mesh
    ! boundary flux
         RETURN
      
        ELSE 
!   Case where the error is not guaranted but the pixel-lvl is reached
            CALL SPLITAFLX1(nn,nr,nc,nx,ny,nz,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)

            CALL PROJMEAN1(nn,nr,nc,nx,ny,nz,
     &                     imin,imax,jmin,jmax,kmin,kmax,
     &                     aflxmean, aflx1)

            CALL SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)

      ENDIF
      nb_cell = nb_cell + 8

      END SUBROUTINE RECADA_ONE_OCTANT