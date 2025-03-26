       SUBROUTINE SWEE3D_ADAPTIVE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir,
       ! inputs: cross section  
     &                          sigt,sigs,
       ! inputs: source and flux moments at previous iteration 
     &                          flxp,srcm,
       ! input/output: boundary flux 
     &                          bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                          mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                          sgnc,sgni,sgne,sgnt,
       ! input: auxiliar memory for boundary angular fluxes, angular source moments and self-scattering xs 
     &                          finc,fout,tmom,sigg,
       ! input: auxiliar memory angular flux & source 
     &                          aflx,asrc,dsrc,
       ! input: auxiliar memory for the interface angular flux 
     &                          xaux,yaux,zaux,
       ! input: auxiliar memory to drive angular mirror-reflection or rotation/translation b.c.
     &                          rdir,
       ! input: pixel-to-medium array 
     &                          zreg, 
       ! input: pixel-to-medium array 
     &                          dira,dirf,
       ! input: logical to add angular source in case of time-dependent calculation 
     &                          lgki,
       ! output: flux moments 
     &                          flxm)

      USE FLGOCT
      USE FLGBCD
      USE FLGINC

      IMPLICIT NONE

      INTEGER, PARAMETER :: noct = 8
      
      INTEGER :: nn,ng,nr,nh,nc,nx,ny,nz,nmat
      INTEGER :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER :: nani,nhrm,nd,ndir
      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct),bfly(nn,nb,nbfy,2,noct),
     &           bflz(nn,nb,nbfz,2,noct)
      REAL    :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(nr)
      REAL    :: aflx(nn,nr,nc,noct),asrc(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,*)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      REAL    :: xaux(nn,nb),yaux(nn,nx,nb),zaux(nn,nx,ny,nb)
      INTEGER :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL :: lgki

      INTEGER :: xout,yout,zout,fst
      INTEGER :: oc,oct,d,dd,da,c,h,i


     
     
      CALL GMOM3D(ng,nani,nh,nc,nr,sigs,flxp,srcm,zreg,sigg,tmom)

!     Loop over octants of angular space in order defined by
!     octant list "olst3D".

      DO oc=1,8
         oct=olst3D(oc)
!        Index of outgoing side.
         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
!        Sweep for all directions in octant "oct".
!        initial direction of the octant 

         da = (oct-1)*ndir+1

!        Directional source.

         CALL GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr(1,da),
     &               tmom(1,1,1,1),asrc(1,1,1,oct))
        
         IF(lgki)CALL SPXPY(nr*nc,dsrc(1,1,1,oct),asrc(1,1,1,oct))


!        Sweep.
        ! The boundary entries are computed on the exiting part
        ! in order to work only on the exiting part of the bdfl
        DO i=1,nbfx
            bflx(:,:,i,xout,oct) = bflx(:,:,i,xinc(oct),oct)
        ENDDO
        DO i=1,nbfy
            bfly(:,:,i,yout,oct) = bflx(:,:,i,yinc(oct),oct)
        ENDDO
        DO i=1,nbfz
            bflz(:,:,i,zout,oct) = bflx(:,:,i,zinc(oct),oct)
        ENDDO

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc(oct), yinc(oct), zinc(oct),oct,
     &                             1,nx,1,ny,1,nz,
     &                             0,2.0,sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc(:,:,oct),sgni(:,:,oct),
     &                             sgne(:,:,oct),sgnt(:,:,oct),
     &                             zreg, asrc(:,:,:,oct),
     &                             aflx(:,:,:,oct),
     &                             bflx(:,:,:,xout,oct),
     &                             bfly(:,:,:,yout,oct),
     &                             bflz(:,:,:,zout,oct))      
     
         ! CALL SWEEP
!        Angular moments.
     
         DO c=1,nc
         DO h=1,nh
         fst = 0
         DO d = 1,ndir
            flxm(:ng,:nr,h,c) = flxm(:ng,:nr,h,c) + 
     &        (sphr(h,da)*w(d))*aflx(fst+1:fst+ng,:nr,c,oct)
            fst = fst + ng
         ENDDO
         ENDDO
         ENDDO
 
      ENDDO
     

!     Boundary conditions: specular reflection, albedo
    
      DO oc=1,8
         oct=olst3D(oc)

!        Index of outgoing side.

         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
         CALL BCD2R3(oct,xout,yout,zout,
     &               bflx,bfly,bflz,mu,eta,
     &               ksi,w,pisn,rdir)
     
      ENDDO 
      
      END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      RECURSIVE SUBROUTINE RECADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,
     &                             xinc,yinc,zinc, oct,
     &                             imin,imax,jmin,jmax,kmin,kmax,niv,
     &                             delt, sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)
     
      USE SWEEP8ONE
      USE SRCCORR 

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER :: nn,ng,ndir,nr,nh,
     &           imin,imax,jmin,jmax,kmin,kmax,
     &           nx,ny,nz,niv,oct
      INTEGER, INTENT(IN) :: xinc,yinc,zinc
      REAL, INTENT(IN)    :: delt
    !   REAL, INTENT(IN)    :: flxm(ng,nr,nh,nc)
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),bfly(nn,nb,nx*nz),
     &           bflz(nn,nb,nx*ny)
     
      REAL  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL  :: xshom0(ng)
      REAL  :: xshom1(ng,n8)
      REAL  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL          :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL          :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8)  :: errcor(ng, ndir)
      REAL(KIND=8)  :: errmul(ng, ndir)
      
      REAL          :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL          :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL       ::  ok
      REAL          :: delt3(3)
      INTEGER       :: i,j,k,r 

      delt3 = (/ delt, delt, delt /)/(2**niv)

      fout0 = 0.0
      fout1 = 0.0

      CALL XSSRCHOMO(nn,ng,nr,nx,ny,ndir,
     &               imin,imax,jmin,jmax,kmin,kmax,
     &               sigt, zreg, aflx, w,
     &               asrc, xshom1, xshom0,
     &               asrcm1, asrcm0)
     
! Compute coefficients 
    !   print *,"xshom0", xshom0
    !   print *,"asrcm0", asrcm0
    !   print*, "oct", oct

      CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                delt3,xshom0,
     &                sgnc,sgni,sgne,sgnt,
     &                ccof,icof,ecof,tcof)

      CALL MERGEBOUND(xinc, yinc, zinc, nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, finc1, finc0)

! Compute level-0
      CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)

! Error check by source-correction estimation
      CALL SRCCOR(ng,ndir,delt,mu,eta,ksi,asrcm0, xshom0,errcor)

!     7/ query ok/ko? 
!     Définir un critère 
      ok = .FALSE.
    !   print *,"aflx0", aflx0

      IF (.NOT. ok) THEN 
!         8/ if ko: compute coefficiets for level-1
        CALL EIGHT_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                    delt3/2,xshom1,
     &                    sgnc,sgni,sgne,sgnt,
     &                    ccof8,icof8,ecof8,tcof8)
!         9/ upload the boundary source !         Done 
!         10/ upload the source on level-1!         Done
!         11/ solve node-1

        CALL SWEEP_8REGIONS(nn,2,asrcm1, finc1, aflx1, fout1,
     &                     ccof8,icof8,ecof8,tcof8, xinc, yinc, zinc,
     &                     nx,ny)

    !   CALL SRC2LVL(ng, ndir, srcm1, sigt, errmul)
!         12/ error check by 2nd -order estimation 
!         13/ query ok/ko ? 
!         14/ si ok :
      ENDIF

      IF (ok) THEN 
        print *,"ok"
        ! - write the angular flux of node-0 on the fine-mesh flux 
        DO i=imin,imax
            DO j=jmin,jmax
                DO k=kmin,kmax
                    r = ((k-1)*ny + (j-1))*nx + i
                    aflx(:,r,:) = aflx0(:,:)
                ENDDO
            ENDDO
        ENDDO
!      - write the outgoing angular flux of node-0 
!            on the fine-mesh boundary flux
         CALL SPLITBOUND0(xinc, yinc, zinc, nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)
        
          !         - return 
          RETURN

      ELSE IF((imax-imin)>2 .AND. .NOT. ok ) THEN
!         15/ si ko : 
!        - call recursively the subroutine for the 8 nodes of level-1

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   imin + (xinc-1)*(imax-imin) ,(imin + imax)/2,
     &   jmin + (yinc-1)*(jmax-jmin), (jmin + jmax)/2,
     &   kmin + (zinc-1)*(kmax-kmin),(kmin + kmax)/2,
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (imax+imin)/2 ,imin + (2-xinc)*(imax-imin),
     &   jmin + (yinc-1)*(jmax-jmin),(jmin + jmax)/2,
     &   kmin + (zinc-1)*(kmax-kmin),(kmin + kmax)/2,
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   imin + (xinc-1)*(imax-imin) ,(imin + imax)/2,
     &   (jmin + jmax)/2, jmin + (2-yinc)*(jmax-jmin),
     &   kmin + (zinc-1)*(kmax-kmin),(kmin + kmax)/2,
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &  (imax+imin)/2,imin + (2-xinc)*(imax-imin),
     &  (jmin + jmax)/2,jmin + (2-yinc)*(jmax-jmin),
     &  kmin + (zinc-1)*(kmax-kmin),(kmin + kmax)/2,
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)


!!!!!!!!!!!!!!!!!!!!!!!

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &  imin + (xinc-1)*(imax-imin) ,(imin + imax)/2,
     &  jmin + (yinc-1)*(jmax-jmin),(jmin + jmax)/2,
     &  (kmin + kmax)/2,kmin + (2-zinc)*(kmax-kmin),
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)



        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (imax+imin)/2, imin + (2-xinc)*(imax-imin),
     &   jmin + (yinc-1)*(jmax-jmin),(jmin + jmax)/2,
     &   (kmin + kmax)/2,kmin + (2-zinc)*(kmax-kmin),
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)        

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   imin + (xinc-1)*(imax-imin),(imin + imax)/2,
     &   (jmin + jmax)/2,jmin + (2-yinc)*(jmax-jmin),
     &   (kmin + kmax)/2,kmin + (2-zinc)*(kmax-kmin),
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)        

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (imax+imin)/2,imin + (2-xinc)*(imax-imin),
     &   (jmin + jmax)/2,jmin + (2-yinc)*(jmax-jmin),
     &   (kmin + kmax)/2,kmin + (2-zinc)*(kmax-kmin),
     &                             niv+1,
     &                             delt, sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz)        

    !- write the outgoing angular flux of node-0 on the fine-mesh
    ! boundary flux
         RETURN
      
        ELSE 
!   Case where the error is not guaranted but the pixel-lvl is reached
            CALL SPLITAFLX1(nn,nr,nc,imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)

            CALL SPLITBOUND1(xinc, yinc, zinc,nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)
      ENDIF
      
      END SUBROUTINE RECADA_ONE_OCTANT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE SPLITAFLX1(nn,nr,nc,imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)

        INTEGER, PARAMETER :: n8=8
        INTEGER, INTENT(IN) :: nn,nr,nc,imin,imax,jmin,jmax,kmin,kmax
        REAL, INTENT(IN) :: aflx1(nn,nc,n8)
        REAL, INTENT(OUT) :: aflx(nn,nr,nc)

        INTEGER :: i_half,j_half,k_half,cnt,icnt,jcnt,kcnt
        INTEGER  :: itot,jtot,ktot,ii,jj,kk

        i_half = (imin+imax)/2
        itot   = (imax-imin)/2
        j_half = (jmin+jmax)/2
        jtot   = (jmax-jmin)/2
        k_half = (kmin+kmax)/2
        ktot   = (kmax-kmin)/2

        cnt = 1
        kcnt = kmin
        DO kk =0,1
            jcnt =jmin
            DO jj =0,1 
                icnt = imin
                DO ii =0,1 
                   DO z= kcnt, kcnt+ktot
                   DO y= jcnt, jcnt+jtot
                   DO x= icnt, icnt+itot
                      r = ((z-1)*ny + (y-1))*nx + x 
                      aflx(:,r,:) = aflx1(:,:,cnt)
                   ENDDO
                   ENDDO
                   ENDDO
                   cnt = cnt + 1 
                   icnt = i_half + 1
                ENDDO
               jcnt = j_half + 1
            ENDDO
           kcnt = k_half + 1
        ENDDO  

      END SUBROUTINE SPLITAFLX1


      SUBROUTINE XSSRCHOMO(nn,ng,nr,nx,ny, ndir,
     &                  imin,imax,jmin,jmax,kmin,kmax,
     &                  sigt, zreg, aflx, w, asrc, xshom1, xshom0,
     &                  srchomo1, srchomo0)
      ! inputs
      IMPLICIT NONE

      INTEGER, PARAMETER :: n8 = 8, nc = 4
      INTEGER, INTENT(IN) :: nn,ng,nr,ndir,
     &           imin,imax,jmin,jmax,kmin,kmax,nx,ny
      
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL(KIND=8), INTENT(IN)    :: w(ndir)
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL, INTENT(IN)    :: aflx(ng,ndir,nr,nc)
      REAL, INTENT(IN)    :: asrc(nn,nr,nc)
      ! output 
      REAL, INTENT(OUT)    :: xshom1(ng,n8), srchomo1(nn,nc,n8)
      REAL, INTENT(OUT)    :: xshom0(ng), srchomo0(nn,nc)
      
      ! locals 
      REAL(KIND=8)    :: flxoct(ng,nr)
      INTEGER :: kk,jj,ii,cnt,icnt,jcnt,kcnt,x,y,z,r,d,
     &           i_half,j_half,k_half,itot,ktot,jtot,m
      REAL(KIND=8) :: aux(ng,n8), auy(ng,n8)

      i_half = (imin+imax)/2
      itot   = (imax-imin)/2
      j_half = (jmin+jmax)/2
      jtot   = (jmax-jmin)/2
      k_half = (kmin+kmax)/2
      ktot   = (kmax-kmin)/2

      flxoct = 0.0
      DO z=kmin,kmax
        DO y=jmin,jmax
            DO x=imin,imax
                DO d=1,ndir
                r = ((z-1)*ny + (y-1))*nx + x 
                flxoct(:,r) = flxoct(:,r) + 8.0*w(d)*aflx(:,d,r,1)
                ENDDO
            ENDDO 
        ENDDO
      ENDDO

      aux = 0.0D0
      auy = 0.0D0
      srchomo1 = 0.0D0

      cnt = 1
      kcnt = kmin
      DO kk =0,1
          jcnt =jmin
          DO jj =0,1 
              icnt = imin
              DO ii =0,1 
                 DO z= kcnt, kcnt+ktot
                 DO y= jcnt, jcnt+jtot
                 DO x= icnt, icnt+itot
                    r = ((z-1)*ny + (y-1))*nx + x 
                    m = zreg(r)
                    aux(:,cnt) = aux(:,cnt) + sigt(:,m) * flxoct(:,r)
                    auy(:,cnt) = auy(:,cnt) + flxoct(:,r)
                    srchomo1(:,:,cnt) = srchomo1(:,:,cnt) 
     &                                + asrc(:,r,:)
                 ENDDO
                 ENDDO
                 ENDDO
                 cnt = cnt + 1 
                 icnt = i_half + 1
              ENDDO
             jcnt = j_half + 1
          ENDDO
         kcnt = k_half + 1
      ENDDO
!   &                      + asrc(:,r,4)*4/(kmax-kmin)**2
!   &                      + asrc(:,r,1)*(2*z - 2*kcnt + ktot)

      xshom1(:,:) = aux(:,:)/auy(:,:)
      xshom0 = sum(aux, DIM=2)/sum(auy, DIM=2)

      !!faux
      srchomo1 = srchomo1/((itot+1)*(jtot+1)*(ktot+1))  
      srchomo0(:,1) = sum(srchomo1(:,1,:), DIM=2)

      cnt = 0
      DO ii =0,1
        DO jj = 0,1
          DO kk = 0,1
            cnt = cnt + 1
            srchomo0(:,2) = srchomo0(:,2) + srchomo1(:,1,cnt)*(ii-0.5)
     &                      + srchomo1(:,2,cnt)/2
            srchomo0(:,3) = srchomo0(:,3) + srchomo1(:,1,cnt)*(jj-0.5)
     &                      + srchomo1(:,3,cnt)/2
            srchomo0(:,4) = srchomo0(:,4) + srchomo1(:,1,cnt)*(kk-0.5)
     &                      + srchomo1(:,4,cnt)/2
          ENDDO
        ENDDO
      ENDDO

      srchomo0  = srchomo0/8

      END SUBROUTINE XSSRCHOMO




      SUBROUTINE MERGEBOUND(xinc, yinc, zinc, nn, ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, finc1, finc0)

      IMPLICIT NONE

! Define the projection of the angular flux on the boundary
! on level 1 and level 0
! finc0(nn, nb, nb): on level 0 for each group for each direction
! on each 3 spatial component on each 3 incoming face
! finc1(nn, nb, nb, ns) same but on level 1 and for each 4 subfaces
       
      INTEGER ,PARAMETER :: ns = 4

      INTEGER, INTENT(IN) :: nn, ng,nb,xinc,yinc,zinc  
      INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax, nx,ny,nz
      REAL, INTENT(IN)    :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
    
      REAL, INTENT(INOUT)  :: finc1(nn,nb,nb,ns), finc0(nn,nb,nb)

      INTEGER :: cnt, icnt, jcnt, kcnt, ii,jj, kk,i,j,k
      INTEGER :: i_half, j_half, k_half

      finc1 = 0.0
      finc0 = 0.0

      i_half = (imin+imax)/2
      j_half = (jmin+jmax)/2
      k_half = (kmin+kmax)/2
      
      cnt  = 1
      jcnt = 0
      kcnt = 0
      DO kk=0,1
        jcnt = 0
        DO jj=0,1
          DO k=kcnt, kcnt+k_half
            DO j=jcnt, jcnt+j_half
              finc1(:,:,1,cnt) = finc1(:,:,1,cnt) + bflx(:,:,j,k)
            ENDDO
          ENDDO
          cnt = cnt + 1
          jcnt = jcnt + j_half
        ENDDO
        kcnt = kcnt + k_half
      ENDDO
      
      cnt  = 1
      icnt = 0
      kcnt = 0
      DO kk=0,1
        icnt = 0
        DO ii=0,1
          DO k=kcnt, kcnt+k_half
            DO i=icnt, icnt+i_half
              finc1(:,:,2,cnt) = finc1(:,:,2,cnt) + bfly(:,:,i,k)
            ENDDO
          ENDDO
        icnt = icnt +i_half
        cnt = cnt + 1
        ENDDO
        kcnt = kcnt + k_half
      ENDDO
      
      cnt = 1
      icnt = 0
      jcnt = 0
      DO jj=0,1
        icnt=0
        DO ii=0,1
          DO j=jcnt, jcnt+j_half
            DO i=icnt, icnt+i_half
               finc1(:,:,3,cnt) = finc1(:,:,3,cnt) + bflz(:,:,i,j)
            ENDDO
          ENDDO
          icnt = icnt + i_half
          cnt = cnt + 1
        ENDDO
        jcnt = jcnt+ j_half
      ENDDO

      finc1 = finc1/4
      finc0 = sum(finc1, DIM=4)/4

      END SUBROUTINE MERGEBOUND

!----------------------------------------------------------------------
      
      SUBROUTINE SPLITBOUND0(xinc, yinc, zinc, nn, ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout0)

      IMPLICIT NONE
      INTEGER ,PARAMETER :: ns = 4

      INTEGER, INTENT(IN) :: nn, ng,nb,xinc, yinc, zinc 
      INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax, nx,ny,nz
      REAL, INTENT(IN)    :: fout0(nn,nb,nb)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)

      INTEGER :: i,j,k, xout, yout, zout
      xout= 3-xinc
      yout= 3-yinc
      zout= 3-zinc

      DO k=kmin, kmax
        DO j=jmin, jmax
          bflx(:,:,j,k) = fout0(:,:,1)
        ENDDO
      ENDDO
      DO k=kmin, kmax
        DO i=imin, imax
          bfly(:,:,i,k) = fout0(:,:,2)
        ENDDO
      ENDDO
      DO j=kmin, kmax
        DO i=imin, imax
          bflz(:,:,i,j) = fout0(:,:,3)
        ENDDO
      ENDDO
      
      END SUBROUTINE SPLITBOUND0


      SUBROUTINE SPLITBOUND1(xinc, yinc, zinc,nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)
   
         IMPLICIT NONE
         INTEGER ,PARAMETER :: ns = 4
   
         INTEGER, INTENT(IN) :: nn, ng,nb,xinc,yinc, zinc 
         INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax,nx,ny,nz
         REAL, INTENT(IN)    :: fout1(nn,nb,nb,ns)
         REAL, INTENT(INOUT) :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
   
         INTEGER :: i,j,k, xout, yout, zout
         INTEGER :: cnt, icnt, jcnt, kcnt, ii,jj, kk
         INTEGER :: i_half, j_half, k_half

         xout= 3-xinc
         yout= 3-yinc
         zout= 3-zinc

         i_half = (imin+imax)/2
         j_half = (jmin+jmax)/2
         k_half = (kmin+kmax)/2
         
         cnt  = 1
         jcnt = 0
         kcnt = 0
         DO kk=0,1
           jcnt = 0
           DO jj=0,1
             DO k=kcnt, kcnt+k_half
               DO j=jcnt, jcnt+j_half
                bflx(:,:,j,k) =  fout1(:,:,1,cnt)
               ENDDO
             ENDDO
             cnt = cnt + 1
             jcnt = jcnt + j_half
           ENDDO
           kcnt = kcnt + k_half
         ENDDO
         
         cnt  = 1
         icnt = 0
         kcnt = 0
         DO kk=0,1
           icnt = 0
           DO ii=0,1
             DO k=kcnt, kcnt+k_half
               DO i=icnt, icnt+i_half
                 bfly(:,:,i,k) = fout1(:,:,2,cnt)
               ENDDO
             ENDDO
           icnt = icnt +i_half
           cnt = cnt + 1
           ENDDO
           kcnt = kcnt + k_half
         ENDDO
         
         cnt = 1
         icnt = 0
         jcnt = 0
         DO jj=0,1
           icnt=0
           DO ii=0,1
             DO j=jcnt, jcnt+j_half
               DO i=icnt, icnt+i_half
                  bflz(:,:,i,j) = fout1(:,:,3,cnt)
               ENDDO
             ENDDO
             icnt = icnt + i_half
             cnt = cnt + 1
           ENDDO
           jcnt = jcnt+ j_half
         ENDDO
         
      END SUBROUTINE SPLITBOUND1