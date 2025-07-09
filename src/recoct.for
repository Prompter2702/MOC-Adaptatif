      SUBROUTINE SWEEP_3D_OCTREE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,nb_niv,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir, octree_root,
     &                          tab_morton,
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
     &                          maxinner, tolinner, tolcor, pdslu4,
     &                          aflx0,asrcm0, finc0, fout0, lastadd)
   
      USE FLGOCT
      USE FLGBCD
      USE FLGINC
      USE SWEEP8ONE
      USE SRCCORR
      USE OCTREE
      
      IMPLICIT NONE
   
      INTEGER, PARAMETER :: noct = 8, n8=8, ns=4
      INTEGER, INTENT(IN) :: nn,ng,nr,nh,nc,nmat
      INTEGER(KIND=2), INTENT(IN) :: nx,ny,nz
      INTEGER(KIND=1), INTENT(IN) :: nb_niv
      INTEGER, INTENT(IN) :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER, INTENT(IN) :: nani,nhrm,nd,ndir
      REAL, INTENT(INOUT) :: flxm(ng,nr,nh,nc,2)
      REAL, INTENT(INOUT) :: bflx(nn,nb,nbfx,2,noct),
     &                       bfly(nn,nb,nbfy,2,noct),
     &                       bflz(nn,nb,nbfz,2,noct)
      REAL, INTENT(IN)    :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),
     &                            w(ndir),pisn
      
      TYPE(CellOctree), INTENT(INOUT) :: octree_root
      INTEGER(KIND=4), INTENT(INOUT)  :: tab_morton

      REAL, INTENT(IN)    :: sphr(nhrm,nd)
      INTEGER,INTENT(IN) :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER,INTENT(IN) :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(IN)    :: srcm(ng,nr,nh,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc,noct)
   
      REAL, INTENT(INOUT)    :: asrc(nn,nr,nc,noct)
      REAL, INTENT(INOUT)    :: dsrc(nn,nr,nc,*)
      REAL, INTENT(INOUT)    :: tmom(ng,nr,nh,nc)
      REAL, INTENT(INOUT)    :: sigg(ng,nr)
      INTEGER, INTENT(INOUT) :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL, INTENT(INOUT) :: lgki
      REAL, INTENT(IN)       :: delt3(3)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir) ! pds to compute src 4 \int(L^4)
   
      INTEGER(KIND=1) :: xout,yout,zout, xin,yin,zin
      INTEGER :: fst
      INTEGER :: oc,oct,d,da,c,h
   
      REAL, INTENT(INOUT)  :: aflx0(nn,nc)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3), fout0(nn,nb,3)
   
      REAL, INTENT(IN)  :: ccofglb(nn,nc, nc, nmat,nb_niv),
     &                        icofglb(nn,nc, nbd,nmat,nb_niv)
      REAL, INTENT(IN)  :: ecofglb(nn,nbd,nc,nmat,nb_niv),
     &                        tcofglb(nn,nbd,nbd,nmat,nb_niv)

      REAL  :: errbnd(3)
      REAL(KIND=8)  :: errcor(ng, ndir),errmul(ng, ndir)
      INTEGER :: ordre_local(n8)


      LOGICAL  :: okinner
      INTEGER , INTENT(IN) :: maxinner
      REAL, INTENT(IN) :: tolinner, tolcor
      REAL :: errinner
      INTEGER :: cnt_iter,g,r,nb_cell,inew,iold, addrflx(2),cnt
      INTEGER, PARAMETER :: olda=1, newa=2
      INTEGER, INTENT(INOUT) :: lastadd
      REAL(KIND=8) :: max_err_inner, tol_tmp
      REAL(KIND=8) ::  errtot
   
      cnt_iter = 0
      addrflx = (/olda, newa/)
      
      okinner = .FALSE.
      errbnd = 0.0
      ordre_local = 0
      tmom = 0.0
      flxm = 0.0
      ! Initialization of the coefficients to compute the 
   
      ! Internal iterations loop
      DO WHILE (.NOT. okinner .AND. cnt_iter<maxinner)
      cnt_iter = cnt_iter + 1
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
   
      !Define the directional sources in asrc from tmom and sphr
   
        CALL GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr(1,da),
     &              tmom(1,1,1,1),asrc(1,1,1,oct))
        IF(lgki)CALL SPXPY(nr*nc,dsrc(1,1,1,oct),asrc(1,1,1,oct))
        
   !        Index of outgoing side.
         xin = xinc(oct)
         yin = yinc(oct)
         zin = zinc(oct)
         xout=3_1-xin
         yout=3_1-yin
         zout=3_1-zin
   
         bflx(:,:,:,xout,oct) = bflx(:,:,:,xinc(oct),oct)
         bfly(:,:,:,yout,oct) = bfly(:,:,:,yinc(oct),oct)
         bflz(:,:,:,zout,oct) = bflz(:,:,:,zinc(oct),oct)

   !       Initial direction of the octant
         da = (oct-1)*ndir+1
         nb_cell = 0

         ordre_local(1) = xinc(oct) + 2*(yinc(oct)-1) + 4*(zinc(oct)-1)
         ordre_local(2) = xout      + 2*(yinc(oct)-1) + 4*(zinc(oct)-1)
         ordre_local(3) = xinc(oct) + 2*(yout-1)      + 4*(zinc(oct)-1)
         ordre_local(4) = xout      + 2*(yout-1)      + 4*(zinc(oct)-1)
         ordre_local(5) = xinc(oct) + 2*(yinc(oct)-1) + 4*(zout-1)
         ordre_local(6) = xout      + 2*(yinc(oct)-1) + 4*(zout-1)
         ordre_local(7) = xinc(oct) + 2*(yout-1)      + 4*(zout-1)
         ordre_local(8) = xout      + 2*(yout-1)      + 4*(zout-1)
         
         tab_morton = 0


        CALL RECADA_OCTREE(nn,ng,ndir,nr,nh,
     &                     nx,ny,nz,nmat,nb_niv,
     &                     xin, yin, zin,oct,
     &                     octree_root,
     &                     tab_morton,ordre_local,0_1,
     &                     sigt,mu,eta,ksi,w,pdslu4,
     &                     sgnc(:,:,oct),sgni(:,:,oct),
     &                     sgne(:,:,oct),sgnt(:,:,oct),
     &                     zreg, asrc(:,:,:,oct),
     &                     aflx(:,:,:,oct),flxm(1,1,1,1,iold),
     &                     bflx(1,1,1,xout,oct),
     &                     bfly(1,1,1,yout,oct),
     &                     bflz(1,1,1,zout,oct),
     &                     aflx0,asrcm0,finc0, fout0,
     &                     ccofglb,icofglb,ecofglb,tcofglb,
     &                     tolcor, errcor,
     &                     delt3,r,nb_cell,cnt,
     &                     errtot,errbnd)

     
        print *, xinc(oct), yinc(oct), zinc(oct)
        print *,"nb_cell", nb_cell
        ! read(*,*)
   
      
    !   CALL BCD2R2(oct,xout,yout,bflx,bfly,mu,eta,w,pisn,rdir)
   
    !  Compute the flux moments
        DO c=1,nc
        DO h=1,nh
        fst = 0
        DO d = 1,ndir
            flxm(:ng,:nr,h,c,inew) = flxm(:ng,:nr,h,c,inew) +
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
        IF ( errinner
     &   > tolinner * ABS( flxm(g,r,h,c,iold) ) + epsilon(1.0)    )THEN
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
      print *,"Number iteration", cnt_iter
    !   print *,"Error inner", errinner
      
      END SUBROUTINE SWEEP_3D_OCTREE
      
      
      !-----------------------------------------------------------------
      
      RECURSIVE SUBROUTINE RECADA_OCTREE(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,nmat,nb_niv,
     &                             xinc,yinc,zinc,oct,
     &                             cell,tab_morton, ordre_local, niv,
     &                             sigt,mu,eta,ksi,w,pdslu4,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,flx,
     &                             bflx, bfly, bflz,
     &                             aflx0,asrcm0,finc0,fout0,
     &                             ccofglb,icofglb,ecofglb,tcofglb,
     &                             tolcor, errcor,
     &                             delt3,r,
     &                             nb_cell,cnt,errtot, errbnd)

      USE SWEEP8ONE
      USE SRCCORR
      USE OCTREE

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh,oct,nmat
      INTEGER(KIND=1), INTENT(IN) :: niv, nb_niv
      INTEGER(KIND=2), INTENT(IN) :: nx,ny,nz
      INTEGER(KIND=1), INTENT(IN) :: xinc,yinc,zinc
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)


      Type(CellOctree), INTENT(INOUT) :: cell
      INTEGER, INTENT(IN) :: ordre_local(n8)
      INTEGER(KIND=4), INTENT(INOUT)  :: tab_morton
      INTEGER(KIND=4) :: num
      INTEGER(KIND=2) :: imin,imax,jmin,jmax,kmin,kmax


      REAL, INTENT(IN)    ::  tolcor
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc), flx(ng,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),
     &                       bfly(nn,nb,nx*nz),
     &                       bflz(nn,nb,nx*ny)
     
      REAL, INTENT(INOUT)  :: aflx0(nn,nc)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3),    fout0(nn,nb,3)

      REAL, INTENT(IN) ::    ccofglb(nn,nc ,nc ,nmat,nb_niv),
     &                       icofglb(nn,nc ,nbd,nmat,nb_niv),
     &                       ecofglb(nn,nbd,nc ,nmat,nb_niv),
     &                       tcofglb(nn,nbd,nbd,nmat,nb_niv)

      REAL(KIND=8),INTENT(INOUT)  :: errcor(ng, ndir)
      REAL(KIND=8), INTENT(INOUT)  :: errtot
      REAL, INTENT(IN)        :: delt3(3)
      INTEGER, INTENT(INOUT)  :: r,nb_cell,cnt
      REAL, INTENT(INOUT) :: errbnd(3)

      !Variables locales 
      LOGICAL :: allchilddone, ok, isleaf
      REAL :: norm1
      INTEGER :: i,j

      !--------IF NOT LEAF THEN CALL RECURSIVE SUBROUTINE IMMEDIATELY---

      IF (tab_morton < 0) THEN
         print *, "Error in tab_morton, it should be positive"
         RETURN
      ENDIF 


      IF( ibits(tab_morton, 30,1) == 1) THEN
         print *, "Octree completed"
         RETURN
      END IF
      
      CALL IS_LEAF(cell,isleaf)

      ok = .FALSE.

      IF (isleaf) THEN
         ! Compute, imin, imax, jmin, jmax, kmin, kmax

         CALL BOUND_REG(nx,ny,nz,niv,tab_morton, xinc,yinc,zinc,
     &                  imin,imax,jmin,jmax,kmin,kmax)

         CALL IS_HOMOEGENEOUS(zreg,imin,imax,jmin,jmax,kmin,kmax,
     &                        nx,ny,ok)

         IF (ok) THEN
            !Compute this one
            r = zreg(((imin-1)*ny + (jmin-1))*nx + kmin)

            CALL SRCHOMO0(nn,nr,nx,ny,
     &                    imin,imax,jmin,jmax,kmin,kmax,
     &                    asrc,asrcm0)
            CALL MERGEBOUND0(nn,nb,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx,bfly,bflz,finc0)
            CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                   sgnc,sgni,sgne,sgnt,   
     &                   ccofglb(:,:,:,r,niv+1),icofglb(:,:,:,r,niv+1),
     &                   ecofglb(:,:,:,r,niv+1),tcofglb(:,:,:,r,niv+1))
            CALL ERRSURF(nn,nx,ny,nz,
     &          imin,imax,jmin,jmax,kmin,kmax,
     &          finc0, bflx,bfly,bflz,tcofglb(:,:,:,r,niv+1) ,errbnd)

            CALL SRCCOR(ng,ndir,delt3/(2**niv),mu,eta,ksi,asrcm0,
     &       sigt(ng,r),errcor,pdslu4)

     
            errbnd = errbnd*(delt3(3)**2)/(4**niv)

            !Erreur en norme sup
    !         cnt = 0
    !         group : DO i = 1, ng
    !            DO j = 1, ndir
    !               cnt = cnt + 1
    !               errtot = ABS( errcor(i,j)) + SUM(errbnd)/3
    !               IF ( errtot  > 3*epsilon(1.0) + 
    !  &              tolcor*MAXVAL( ABS( flx(:,r,:) ) )    )THEN
    !                 ok = .FALSE.
    !                 EXIT group
    !               ENDIF
    !            ENDDO
    !         ENDDO group

            !Erreur en norme 1
            cnt = 0
            errtot = 0.0
            norm1  = 0.0
            group : DO i = 1, ng
               DO j = 1, ndir
                  cnt = cnt + 1
                  errtot = errtot + ABS( errcor(i,j))
                  norm1  = norm1 +  ABS( flx(i,r,1) ) + epsilon(1.0)
               ENDDO
            ENDDO group            

            print *, errtot
            print *, tolcor*norm1
            read(*,*)

            IF ( errtot  > tolcor*norm1 ) THEN
                ok = .FALSE.
            ENDIF

         END IF


         IF(ok .OR. (niv.GE.nb_niv-1_1) ) THEN
            tab_morton = tab_morton + 2**(30-3*niv)
            CALL SPLITAFLX0(nn,nr,nc,nx,ny,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   aflx,aflx0)
!      - write the outgoing angular flux of node-0 
!            on the fine-mesh boundary flux
            CALL SPLITBOUND0(nn,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)
            nb_cell = nb_cell + 1
            RETURN 
         ELSE
            CALL ADD_CHILDREN(cell)
         ENDIF

      ENDIF ! END IF LEAF


         
      IF( (.NOT. isleaf) .OR. (niv < nb_niv-1)) THEN
            num = ibits(tab_morton, 30-3*(niv+1_1), 3) 
            allchilddone = .FALSE.

            DO WHILE(.NOT. allchilddone)

               IF (num==7) THEN
                  allchilddone = .TRUE. 
               ENDIF
               CALL RECADA_OCTREE(nn,ng,ndir,nr,nh,
     &                      nx,ny,nz,nmat,nb_niv,
     &                      xinc,yinc,zinc, oct,
     &                      cell%children(ordre_local(num+1))%ptrChild,
     &                      tab_morton,ordre_local,
     &                      niv + 1_1,
     &                      sigt,
     &                      mu,eta,ksi,w,pdslu4,
     &                      sgnc,sgni,sgne,sgnt,
     &                      zreg,asrc,aflx,flx,
     &                      bflx, bfly, bflz,
     &                      aflx0,asrcm0,finc0,fout0,
     &                      ccofglb,icofglb,ecofglb,tcofglb,
     &                      tolcor, errcor,
     &                      delt3,r,
     &                      nb_cell,cnt,errtot,errbnd)

               num = ibits(tab_morton, 30-3*(niv+1_1), 3)

            ENDDO


      END IF

      
      END SUBROUTINE RECADA_OCTREE
      
      !------------------------------------------------------------------


      SUBROUTINE IS_HOMOEGENEOUS(zreg,imin,imax,jmin,jmax,kmin,kmax,
     &                           nx,ny,ok)

      IMPLICIT NONE

      INTEGER(KIND=2), INTENT(IN) :: imin,imax,jmin,jmax,kmin,kmax, 
     &                               nx,ny
      INTEGER, INTENT(IN) :: zreg(*)

      LOGICAL, INTENT(OUT) :: ok

      INTEGER(KIND=2) :: i,j,k
      INTEGER :: cnt,r

      ok = .TRUE.
      cnt = 0
      r   = ((imin-1)*ny + (jmin-1))*nx + kmin
      ido  : DO i= kmin,kmax
        DO j= jmin,jmax
            DO k= imin,imax
                cnt = ((i-1)*ny + (j-1))*nx + k 
                if ( zreg(cnt) .NE. zreg(r)) THEN
                    ok = .FALSE.
                    EXIT ido
                ENDIF
            ENDDO
        ENDDO
      ENDDO ido

      END SUBROUTINE IS_HOMOEGENEOUS

      

      SUBROUTINE BOUND_REG(nx,ny,nz,niv,tab_morton,xi,yi,zi,
     &                     imin,imax,jmin,jmax,kmin,kmax)

      IMPLICIT NONE

      INTEGER(KIND=2), INTENT(IN) :: nx,ny,nz
      INTEGER(KIND=1), INTENT(IN) :: niv
      INTEGER(KIND=4), INTENT(IN) :: tab_morton
      INTEGER(KIND=1), INTENT(IN)  :: xi,yi,zi
      INTEGER(KIND=2), INTENT(INOUT) :: imin,imax,jmin,jmax,kmin,kmax
      
      INTEGER(KIND=1) :: cnt
      INTEGER(KIND=2) :: a,b,c
      
      imin = 1_2
      jmin = 1_2
      kmin = 1_2
      imax = nx
      jmax = ny
      kmax = nz
         
      DO cnt=1,niv
         a = ibits(tab_morton, 30-3*cnt  ,1)
         b = ibits(tab_morton, 30-3*cnt+1,1)
         c = ibits(tab_morton, 30-3*cnt+2,1)
         IF (xi==2)  THEN
            a = 1_2-a
         END IF 
         IF (yi==2)  THEN
            b = 1_2-b
         END IF 
         IF (zi==2)  THEN
            c = 1_2-c
         END IF 
         imin = imin + a*nx/(2_2**cnt)
         jmin = jmin + b*ny/(2_2**cnt)
         kmin = kmin + c*nz/(2_2**cnt)
         imax = imax - (1_2-a)*nx/(2_2**cnt)
         jmax = jmax - (1_2-b)*ny/(2_2**cnt) 
         kmax = kmax - (1_2-c)*nz/(2_2**cnt)
      ENDDO

      END SUBROUTINE BOUND_REG



      subroutine print_binaire(n)
         implicit none
         integer, intent(in) :: n
         integer :: i
         character(len=32) :: bits
         bits = ""

         do i = 31, 0, -1
            if (iand(n, ishft(1, i)) /= 0) then
                  bits(32 - i:32 - i) = '1'
            else
                  bits(32 - i:32 - i) = '0'
            end if
         end do

         print *, trim(bits)
      end subroutine print_binaire

      !-----------------------------------------------------------------

      RECURSIVE SUBROUTINE MEAN_FLUX_CELL_OCTREE(nn,ng,ndir,nr,nh,nhrm,
     &                              nx,ny,nz,nb_niv,
     &                              xinc,yinc,zinc,sphr,w,
     &                              cell,tab_morton,niv,
     &                              delt3,r, aflx,flxmean,depth)

        USE OCTREE

        IMPLICIT NONE

        INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                        ns = 4

        Type(CellOctree), INTENT(INOUT) :: cell
        REAL, INTENT(IN)  :: sphr(nhrm,ndir*8)
        REAL, INTENT(OUT) :: flxmean(ng,nr)
        INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh, nhrm
        INTEGER(KIND=1), INTENT(IN) :: niv, nb_niv
        INTEGER(KIND=2), INTENT(IN) :: nx,ny,nz
        INTEGER(KIND=1), INTENT(IN) :: xinc,yinc,zinc

        REAL(KIND=8), INTENT(IN) :: w(ndir)
        INTEGER(KIND=4), INTENT(INOUT)  :: tab_morton
        REAL, INTENT(INOUT) :: aflx(nn,nr,nc,n8)
        REAL, INTENT(IN)        :: delt3(3)
        LOGICAL, INTENT(IN) :: depth !Put depth to True to have the 
        ! depth of the cell instead of the mean

        INTEGER(KIND=4) :: num,r,da,h,d,fst
        INTEGER(KIND=2) :: imin,imax,jmin,jmax,kmin,kmax,oct,i,j,k
        LOGICAL :: isleaf,allchilddone
        REAL :: flx(ng)

      !--------IF NOT LEAF THEN CALL RECURSIVE SUBROUTINE IMMEDIATELY---

        IF (tab_morton < 0) THEN
            print *, "Error in tab_morton, it should be positive"
            RETURN
        ENDIF 
   
   
        IF( ibits(tab_morton, 30,1) == 1) THEN
           print *, "Octree completed"
           RETURN
        END IF
         
        CALL IS_LEAF(cell,isleaf)
   
        IF (isleaf) THEN
            ! Compute, imin, imax, jmin, jmax, kmin, kmax
   
            CALL BOUND_REG(nx,ny,nz,niv,tab_morton, xinc,yinc,zinc,
     &                     imin,imax,jmin,jmax,kmin,kmax)
               
            r = (imin-1)*ny + (jmin-1)*nx + kmin
            tab_morton = tab_morton + 2**(30-3*niv)
            
            flx = 0.0
            Do oct=1,8
                da = (oct-1)*ndir+1
                DO h=1,nh
                    fst = 0
                    DO d = 1,ndir
                      flx(:ng) = flx(:ng) + 
     &                (sphr(h,da+d-1)*w(d))*aflx(fst+1:fst+ng,r,1,oct)
                       fst = fst + ng
                    ENDDO
                ENDDO
            ENDDO
            

            IF (depth) THEN
                DO i=imin,imax
                    DO j=jmin,jmax
                        DO k=kmin,kmax
                            r = ((k-1)*ny + (j-1))*nx + i
                            flxmean(:,r) = niv
                        ENDDO
                    ENDDO
                ENDDO
            ELSE
                DO i=imin,imax
                    DO j=jmin,jmax
                        DO k=kmin,kmax
                            r = ((k-1)*ny + (j-1))*nx + i
                            flxmean(:,r) = flx(:)
                        ENDDO
                    ENDDO
                ENDDO
            END IF                

            RETURN 

        ELSE ! IF  NOT LEAF
            num = ibits(tab_morton, 30-3*(niv+1_1), 3) 
            allchilddone = .FALSE.
   
            DO WHILE(.NOT. allchilddone)
                IF (num==7) THEN
                   allchilddone = .TRUE. 
                ENDIF
                CALL MEAN_FLUX_CELL_OCTREE(nn,ng,ndir,nr,nh,nhrm,
     &                         nx,ny,nz,nb_niv,
     &                         xinc,yinc,zinc,sphr,w,
     &                         cell%children(num+1)%ptrChild,
     &                         tab_morton,niv+1_1,
     &                         delt3,r, aflx,flxmean,depth)
   
                num = ibits(tab_morton, 30-3*(niv+1_1), 3)
            ENDDO
   
        END IF

      END SUBROUTINE MEAN_FLUX_CELL_OCTREE

