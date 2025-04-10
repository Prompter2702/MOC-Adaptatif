        SUBROUTINE XSSRCHOMO1(nn,ng,nr,nx,ny, ndir,
     &                 imin,imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom1,
     &                 srchomo1)

! Compute the homogenized angular source and cross-section at lvl 1

      IMPLICIT NONE

      INTEGER, PARAMETER :: n8 = 8, nc = 4
      INTEGER, INTENT(IN) :: nn,ng,nr,ndir,
     &                imin,imax,jmin,jmax,kmin,kmax,nx,ny
     
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL(KIND=8), INTENT(IN)    :: w(ndir)
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL, INTENT(IN)    :: aflx(ng,ndir,nr,nc)
      REAL, INTENT(IN)    :: asrc(nn,nr,nc)
! output 
      REAL, INTENT(INOUT) :: xshom1(ng,n8)
      REAL, INTENT(INOUT) :: srchomo1(nn,nc,n8)
! locals 
      REAL(KIND=8)    :: flxoct(ng,nr)
      INTEGER :: kk,jj,ii,cnt,icnt,jcnt,kcnt,x,y,z,r,d,
     &           i_half,j_half,k_half,itot,ktot,jtot,m
      REAL(KIND=8) :: aux(ng,n8), auy(ng,n8), xsmoy(ng,n8)

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
              xsmoy(:,cnt) = xsmoy(:,cnt) + sigt(:,m)
              srchomo1(:,1,cnt) = srchomo1(:,1,cnt) + asrc(:,r,1)
                
              srchomo1(:,2,cnt) = srchomo1(:,2,cnt) + (asrc(:,r,2)
     &        + asrc(:,r,1)*(2*x-2*icnt-itot) )/(itot+1)
     
              srchomo1(:,3,cnt) = srchomo1(:,3,cnt) + ( asrc(:,r,3)
     &        + asrc(:,r,1)*(2*y-2*jcnt-jtot) )/(jtot+1)

              srchomo1(:,4,cnt) = srchomo1(:,4,cnt) + ( asrc(:,r,4)
     &        + asrc(:,r,1)*(2*z-2*kcnt-ktot) )/(ktot+1)
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
     
      xshom1(:,:) = aux(:,:)/auy(:,:)
      DO m=1,ng
        DO kk=1,n8
          IF (auy(m,kk) == 0.0D0) THEN
            xshom1(m,kk) = xsmoy(m,kk)/(itot+1)*(jtot+1)*(ktot+1)
          ENDIF
        ENDDO
      ENDDO
      srchomo1 = srchomo1/((itot+1)*(jtot+1)*(ktot+1))  

      END SUBROUTINE XSSRCHOMO1

!----------------------------------------

      SUBROUTINE XSSRCHOMO0(nn,ng,nr,nx,ny, ndir,
     &                 imin,imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom0,
     &                 srchomo0)

! Compute the homogenized angular source and cross-section at lvl 0

        IMPLICIT NONE
   
        INTEGER, PARAMETER :: nc = 4
        INTEGER, INTENT(IN) :: nn,ng,nr,ndir,
     &                     imin,imax,jmin,jmax,kmin,kmax,nx,ny
        
        INTEGER, INTENT(IN) :: zreg(nr)
        REAL(KIND=8), INTENT(IN)    :: w(ndir)
        REAL, INTENT(IN)    :: sigt(ng,*)
        REAL, INTENT(IN)    :: aflx(ng,ndir,nr,nc)
        REAL, INTENT(IN)    :: asrc(nn,nr,nc)
        ! output 
        REAL, INTENT(OUT)    :: xshom0(ng), srchomo0(nn,nc)
        
        ! locals 
        REAL(KIND=8)    :: flxoct(ng,nr)
        INTEGER :: x,y,z,r,d,m
        REAL(KIND=8) :: aux(ng), auy(ng), xsmoy(ng)
   
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
        srchomo0 = 0.0D0
        xsmoy = 0.0D0
   
        DO z=kmin,kmax
          DO y= jmin,jmax
            DO x= imin,imax
              r = ((z-1)*ny + (y-1))*nx + x 
              m = zreg(r)
              aux(:) = aux(:) + sigt(:,m) * flxoct(:,r)
              xsmoy(:) = xsmoy(:) + sigt(:,m)
              auy(:) = auy(:) + flxoct(:,r)
              srchomo0(:,1) = srchomo0(:,1) + asrc(:,r,1)
              srchomo0(:,2) = srchomo0(:,2) + ( asrc(:,r,2)
     &        + asrc(:,r,1)*(2*x-imax-imin) )/(imax-imin+1)
     
              srchomo0(:,3) = srchomo0(:,3) + ( asrc(:,r,3)
     &        + asrc(:,r,1)*(2*y-jmax-jmin) )/(jmax-jmin+1)

              srchomo0(:,4) = srchomo0(:,4) + ( asrc(:,r,4)
     &        + asrc(:,r,1)*(2*z-kmax-kmin) )/(kmax-kmin+1)

            ENDDO
          ENDDO
        ENDDO

      xshom0= aux/auy
      DO m=1,ng
        IF (auy(m) == 0.0D0) THEN
          xshom0(m) = xsmoy(m)/(kmax-kmin+1)*(jmax-jmin+1)*(imax-imin+1)
        ENDIF
      ENDDO

      srchomo0 = srchomo0/((imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1))  

      END SUBROUTINE XSSRCHOMO0

      SUBROUTINE MERGEBOUND0(nn,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, finc0)
    
      IMPLICIT NONE
   
! Define the projection of the angular flux on the boundary
! on level 1 and level 0
! finc0(nn, nb, nb): on level 0 for each group for each direction
! on each 3 spatial component on each 3 incoming face
! finc1(nn, nb, nb, ns) same but on level 1 and for each 4 subfaces
         
      INTEGER, INTENT(IN) :: nn,nb 
      INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax, nx,ny,nz
      REAL, INTENT(IN)    :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
       
      REAL, INTENT(INOUT)  :: finc0(nn,nb,nb)
      INTEGER :: i,j,k
    
      finc0 = 0.0
    
      DO k=kmin,kmax
        DO j=jmin, jmax
            finc0(:,1,1) = finc0(:,1,1) + bflx(:,1,j,k)
            finc0(:,2,1) = finc0(:,2,1) + (bflx(:,2,j,k) 
     &       + bflx(:,1,j,k)*(2*j-jmax-jmin))/(jmax-jmin+1)
            finc0(:,3,1) = finc0(:,3,1) + (bflx(:,3,j,k)
     &       + bflx(:,1,j,k)*(2*k-kmax-kmin))/(kmax-kmin+1)
        ENDDO
      ENDDO
         
      DO k=kmin,kmax
        DO i=imin, imax
          finc0(:,1,2) = finc0(:,1,2) + bfly(:,1,i,k)
          finc0(:,2,2) = finc0(:,2,2) + (bfly(:,2,i,k)
     &       + bfly(:,1,i,k)*(2*i-imax-imin))/(imax-imin+1)
          finc0(:,3,2) = finc0(:,3,2) + (bfly(:,3,i,k)
     &       + bfly(:,1,i,k)*(2*k-kmax-kmin))/(kmax-kmin+1)
        ENDDO
      ENDDO
         
      DO j=jmin, jmax
        DO i=imin, imax
           finc0(:,1,3) = finc0(:,1,3) + bflz(:,1,i,j)
           finc0(:,2,3) = finc0(:,2,3) + (bflz(:,2,i,j)
     &       + bflz(:,1,i,j)*(2*i-imax-imin))/(imax-imin+1)
           finc0(:,3,3) = finc0(:,3,3) + (bflz(:,3,i,j)
     &       + bflz(:,1,i,j)*(2*j-jmax-jmin))/(jmax-jmin+1)
        ENDDO
      ENDDO
             
      finc0(:,:,1) = finc0(:,:,1)/((kmax-kmin+1)*(jmax-jmin+1))
      finc0(:,:,2) = finc0(:,:,2)/((imax-imin+1)*(kmax-kmin+1))
      finc0(:,:,3) = finc0(:,:,3)/((imax-imin+1)*(jmax-jmin+1))
    
      END SUBROUTINE MERGEBOUND0
   


      SUBROUTINE MERGEBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, finc1)
 
      IMPLICIT NONE

! Define the projection of the angular flux on the boundary
! on level 1 and level 0
! finc0(nn, nb, nb): on level 0 for each group for each direction
! on each 3 spatial component on each 3 incoming face
! finc1(nn, nb, nb, ns) same but on level 1 and for each 4 subfaces
      
      INTEGER ,PARAMETER :: ns = 4
 
      INTEGER, INTENT(IN) :: nn, ng,nb 
      INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax, nx,ny,nz
      REAL, INTENT(IN)    :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
    
      REAL, INTENT(INOUT)  :: finc1(nn,nb,nb,ns)
 
      INTEGER :: cnt, icnt, jcnt, kcnt, ii,jj, kk,i,j,k
      INTEGER :: i_half, j_half, k_half, itot,jtot,ktot
 
      finc1 = 0.0
 
      i_half = (imin+imax)/2
      itot   = (imax-imin)/2
      j_half = (jmin+jmax)/2
      jtot   = (jmax-jmin)/2
      k_half = (kmin+kmax)/2
      ktot   = (kmax-kmin)/2
      
      cnt  = 1
      kcnt = kmin
      DO kk=0,1
        jcnt = jmin
        DO jj=0,1
          DO k=kcnt, kcnt+ktot
            DO j=jcnt, jcnt+jtot
              finc1(:,1,1,cnt) = finc1(:,1,1,cnt) + bflx(:,1,j,k)
              finc1(:,2,1,cnt) = finc1(:,2,1,cnt) + (bflx(:,2,j,k)
     &         + bflx(:,1,j,k)*(2*j-jtot-2*jcnt))/(jtot+1)
              finc1(:,3,1,cnt) = finc1(:,3,1,cnt) + (bflx(:,3,j,k)
     &         + bflx(:,1,j,k)*(2*k-ktot-2*kcnt))/(ktot+1)
            ENDDO
          ENDDO
          cnt = cnt + 1
          jcnt = j_half + 1
        ENDDO
        kcnt = k_half + 1
      ENDDO
      
      cnt  = 1
      kcnt = kmin
      DO kk=0,1
        icnt = imin
        DO ii=0,1
          DO k=kcnt, kcnt+ktot
            DO i=icnt, icnt+itot
              finc1(:,1,2,cnt) = finc1(:,1,2,cnt) + bfly(:,1,i,k)
              finc1(:,2,2,cnt) = finc1(:,2,2,cnt) + (bfly(:,2,i,k)
     &         + bfly(:,1,i,k)*(2*i-itot-2*icnt))/(itot+1)
              finc1(:,3,2,cnt) = finc1(:,3,2,cnt) + (bfly(:,3,i,k)
     &         + bfly(:,1,i,k)*(2*k-ktot-2*kcnt))/(ktot+1)
            ENDDO
          ENDDO
        icnt = i_half + 1 
        cnt = cnt + 1
        ENDDO
        kcnt = k_half + 1
      ENDDO
      
      cnt = 1
      jcnt = jmin
      DO jj=0,1
        icnt=imin
        DO ii=0,1
          DO j=jcnt, jcnt+jtot
            DO i=icnt, icnt+itot
               finc1(:,1,3,cnt) = finc1(:,1,3,cnt) + bflz(:,1,i,j)
               finc1(:,2,3,cnt) = finc1(:,2,3,cnt) + (bflz(:,2,i,j)
     &         + bflz(:,1,i,j)*(2*i-itot-2*icnt))/(itot+1)
               finc1(:,3,3,cnt) = finc1(:,3,3,cnt) + (bflz(:,3,i,j)
     &         + bflz(:,1,i,j)*(2*j-jtot-2*jcnt))/(jtot+1)
            ENDDO
          ENDDO
          icnt = i_half + 1 
          cnt = cnt + 1
        ENDDO
        jcnt = j_half + 1 
      ENDDO

      finc1(:,:,1,:) = finc1(:,:,1,:)/((ktot+1)*(jtot+1))
      finc1(:,:,2,:) = finc1(:,:,2,:)/((itot+1)*(ktot+1))
      finc1(:,:,3,:) = finc1(:,:,3,:)/((itot+1)*(jtot+1))
 
      END SUBROUTINE MERGEBOUND1

!----------------------------------------------------------------------
     
      SUBROUTINE SPLITBOUND0(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout0)


! Split the lvl 0 boundary flux on the outgoing faces 

      IMPLICIT NONE
      INTEGER ,PARAMETER :: ns = 4
 
      INTEGER, INTENT(IN) :: nn, ng,nb
      INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax, nx,ny,nz
      REAL, INTENT(IN)    :: fout0(nn,nb,nb)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
 
      INTEGER :: i,j,k
 
      DO k=kmin, kmax
        DO j=jmin, jmax
          bflx(:,1,j,k) = fout0(:,1,1)
     &            + 3.0*fout0(:,2,1)*(2*j-jmin-jmax)/(jmax-jmin+1.0)
     &            + 3.0*fout0(:,3,1)*(2*k-kmin-kmax)/(kmax-kmin+1.0)
          bflx(:,2,j,k) = fout0(:,2,1)/(jmax-jmin+1)
          bflx(:,3,j,k) = fout0(:,3,1)/(kmax-kmin+1)  
        ENDDO
      ENDDO

      DO k=kmin, kmax
        DO i=imin, imax
          bfly(:,1,i,k) = fout0(:,1,2)
     &           + 3.0*fout0(:,2,2)*(2.0*i-imax-imin)/(imax-imin+1.0)
     &           + 3.0*fout0(:,3,2)*(2.0*k-kmax-kmin)/(kmax-kmin+1.0)
          bfly(:,2,i,k) = fout0(:,2,2)/(imax-imin+1)
          bfly(:,3,i,k) = fout0(:,3,2)/(kmax-kmin+1)
        ENDDO
      ENDDO
      DO j=jmin, jmax
        DO i=imin, imax
          bflz(:,1,i,j) = fout0(:,1,3)
     &           + 3.0*fout0(:,2,3)*(2*i-imax-imin)/(imax-imin+1.0)
     &           + 3.0*fout0(:,3,3)*(2*j-jmax-jmin)/(jmax-jmin+1.0)
          bflz(:,2,i,j) = fout0(:,2,3)/(imax-imin+1)
          bflz(:,3,i,j) = fout0(:,3,3)/(jmax-jmin+1)
        ENDDO
      ENDDO
      
      END SUBROUTINE SPLITBOUND0
 
 
      SUBROUTINE SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)

! Split the lvl 1 boundary flux on the outgoing faces 

         IMPLICIT NONE
         INTEGER ,PARAMETER :: ns = 4
   
         INTEGER, INTENT(IN) :: nn, ng,nb
         INTEGER, INTENT(IN) :: imin, imax,jmin,jmax,kmin,kmax,nx,ny,nz
         REAL, INTENT(IN)    :: fout1(nn,nb,nb,ns)
         REAL, INTENT(INOUT) :: bflx(nn,nb,ny,nz),
     &                       bfly(nn,nb,nx,nz),
     &                       bflz(nn,nb,nx,ny)
   
         INTEGER :: i,j,k
         INTEGER :: cnt, icnt, jcnt, kcnt, ii,jj, kk
         INTEGER :: i_half, j_half, k_half, itot, jtot, ktot
 
         i_half = (imin+imax)/2
         itot   = (imax-imin)/2
         j_half = (jmin+jmax)/2
         jtot   = (jmax-jmin)/2
         k_half = (kmin+kmax)/2
         ktot   = (kmax-kmin)/2
         
         cnt  = 1
         kcnt = kmin
         DO kk=0,1
           jcnt = jmin
           DO jj=0,1
             DO k=kcnt, kcnt+ktot
               DO j=jcnt, jcnt+jtot
                bflx(:,1,j,k) =  fout1(:,1,1,cnt)
     &           + 3.0*fout1(:,2,1,cnt)*(2*j-2*jcnt-jtot)/(jtot+1)
     &           + 3.0*fout1(:,3,1,cnt)*(2*k-2*kcnt-ktot)/(ktot+1)
                bflx(:,2,j,k) =  fout1(:,2,1,cnt)/(jtot+1)
                bflx(:,3,j,k) =  fout1(:,3,1,cnt)/(ktot+1)
               ENDDO
             ENDDO
             cnt = cnt + 1
             jcnt = j_half + 1
           ENDDO
           kcnt = k_half + 1
         ENDDO
         
         cnt  = 1
         kcnt = kmin
         DO kk=0,1
           icnt = imin
           DO ii=0,1
             DO k=kcnt, kcnt+ktot
               DO i=icnt, icnt+itot
                 bfly(:,1,i,k) = fout1(:,1,2,cnt)
     &           + 3.0*fout1(:,2,2,cnt)*(2*i-2*icnt-itot)/(itot+1)
     &           + 3.0*fout1(:,3,2,cnt)*(2*k-2*kcnt-ktot)/(ktot+1)
                 bfly(:,2,i,k) = fout1(:,2,2,cnt)/(itot+1)
                 bfly(:,3,i,k) = fout1(:,3,2,cnt)/(ktot+1)
               ENDDO
             ENDDO
           icnt = i_half + 1
           cnt = cnt + 1
           ENDDO
           kcnt = k_half + 1
         ENDDO
         
         cnt = 1
         jcnt = jmin
         DO jj=0,1
           icnt=imin
           DO ii=0,1
             DO j=jcnt, jcnt+jtot
               DO i=icnt, icnt+itot
                  bflz(:,1,i,j) = fout1(:,1,3,cnt)
     &        + 3.0*fout1(:,2,3,cnt)*(2*i-2*icnt-itot)/(itot+1)
     &        + 3.0*fout1(:,3,3,cnt)*(2*j-2*jcnt-jtot)/(jtot+1)
                  bflz(:,2,i,j) = fout1(:,2,3,cnt)/(itot+1)
                  bflz(:,3,i,j) = fout1(:,3,3,cnt)/(jtot+1)
               ENDDO
             ENDDO
             icnt = i_half + 1
             cnt = cnt + 1
           ENDDO
           jcnt = j_half + 1
         ENDDO
         
      END SUBROUTINE SPLITBOUND1

      SUBROUTINE SPLITAFLX1(nn,nr,nc,nx,ny,nz,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)
 
! Split the lvl 1 anglular flux on the global anglular flux array

          IMPLICIT NONE
   
          INTEGER, PARAMETER :: n8=8
          INTEGER, INTENT(IN) :: nn,nr,nc,imin,imax,jmin,jmax,kmin,kmax
          INTEGER, INTENT(IN) :: nx,ny,nz
          REAL, INTENT(IN) :: aflx1(nn,nc,n8)
          REAL, INTENT(OUT) :: aflx(nn,nr,nc)
   
          INTEGER :: i_half,j_half,k_half,cnt,icnt,jcnt,kcnt
          INTEGER  :: itot,jtot,ktot,ii,jj,kk,r,x,y,z
   
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
                        aflx(:,r,1) = aflx1(:,1,cnt)
     &               + 3.0*aflx1(:,2,cnt)*(2*x-2*icnt-itot)/(itot+1)
     &               + 3.0*aflx1(:,3,cnt)*(2*y-2*jcnt-jtot)/(jtot+1)
     &               + 3.0*aflx1(:,4,cnt)*(2*z-2*kcnt-ktot)/(ktot+1)
                        aflx(:,r,2) = aflx1(:,2,cnt)/(itot+1)
                        aflx(:,r,3) = aflx1(:,3,cnt)/(jtot+1)
                        aflx(:,r,4) = aflx1(:,4,cnt)/(ktot+1)
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


      SUBROUTINE SPLITAFLX0(nn,nr,nc,nx,ny,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx0)
 
! Split the lvl 0 anglular flux on the global anglular flux array

        IMPLICIT NONE
   
        INTEGER, INTENT(IN) :: nx,ny
        INTEGER, INTENT(IN) :: nn,nr,nc,imin,imax,jmin,jmax,kmin,kmax
        REAL, INTENT(IN) :: aflx0(nn,nc)
        REAL, INTENT(OUT) :: aflx(nn,nr,nc)
   
        INTEGER :: i,j,k,r

        DO i=imin,imax
          DO j=jmin,jmax
            DO k=kmin,kmax
                r = ((k-1)*ny + (j-1))*nx + i
                aflx(:,r,1) = aflx0(:,1) 
     &                + 3.0*aflx0(:,2)*(2*i-imax-imin)/(imax-imin+1)
     &                + 3.0*aflx0(:,3)*(2*j-jmax-jmin)/(jmax-jmin+1)
     &                + 3.0*aflx0(:,4)*(2*k-kmax-kmin)/(kmax-kmin+1)
                aflx(:,r,2) = aflx0(:,2)/(imax-imin+1)
                aflx(:,r,3) = aflx0(:,3)/(jmax-jmin+1)
                aflx(:,r,4) = aflx0(:,4)/(kmax-kmin+1)

            ENDDO
          ENDDO
        ENDDO

        END SUBROUTINE SPLITAFLX0