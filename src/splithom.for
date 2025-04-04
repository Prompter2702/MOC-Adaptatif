      SUBROUTINE XSSRCHOMO1(nn,ng,nr,nx,ny, ndir,
     &                 imin,imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom1,
     &                 srchomo1)
! inputs
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
      srchomo1 = srchomo1/((itot+1)*(jtot+1)*(ktot+1))  

    !   cnt = 0
    !   DO ii =0,1
    !     DO jj = 0,1
    !       DO kk = 0,1
    !         cnt = cnt + 1
    !         srchomo0(:,2) = srchomo0(:,2) + srchomo1(:,1,cnt)*(ii-0.5)
    !  &                      + srchomo1(:,2,cnt)/2
    !         srchomo0(:,3) = srchomo0(:,3) + srchomo1(:,1,cnt)*(jj-0.5)
    !  &                      + srchomo1(:,3,cnt)/2
    !         srchomo0(:,4) = srchomo0(:,4) + srchomo1(:,1,cnt)*(kk-0.5)
    !  &                      + srchomo1(:,4,cnt)/2
    !       ENDDO
    !     ENDDO
    !   ENDDO
    !   srchomo0  = srchomo0/8

      END SUBROUTINE XSSRCHOMO1

!----------------------------------------

      SUBROUTINE XSSRCHOMO0(nn,ng,nr,nx,ny, ndir,
     &                 imin,imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom0,
     &                 srchomo0)
   ! inputs
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
        REAL(KIND=8) :: aux(ng), auy(ng)
   
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
   
        DO z=kmin,kmax
            DO y= jmin,jmax
                DO x= imin,imax
                    r = ((z-1)*ny + (y-1))*nx + x 
                    m = zreg(r)
                    aux(:) = aux(:) + sigt(:,m) * flxoct(:,r)
                    auy(:) = auy(:) + flxoct(:,r)
                    srchomo0(:,1) = srchomo0(:,1) + asrc(:,r,1)
                    ! srchomo0(:,2:4) = srchomo0(:,2:4) + asrc(:,r,2:4)
     &                         + 
                ENDDO
            ENDDO
        ENDDO
        srchomo0(:,2:4) = asrc(:,imax+(jmax-1)*nx + (kmax-1)*nx*ny,2:4)
     &   - asrc(:,imin+(jmin-1)*nx + (kmin-1)*nx*ny,2:4)
     
   !   & + asrc(:,r,4)*4/(kmax-kmin)**2
   !   & + asrc(:,r,1)*(2*z - 2*kcnt + ktot)
        
        xshom0= aux/auy
        srchomo0 = srchomo0/((imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1))  
   
            !    srchomo0(:,2) = srchomo0(:,2) + srchomo1(:,1,cnt)*(ii-0.5)
        ! &                      + srchomo1(:,2,cnt)/2
        !        srchomo0(:,3) = srchomo0(:,3) + srchomo1(:,1,cnt)*(jj-0.5)
        ! &                      + srchomo1(:,3,cnt)/2
        !        srchomo0(:,4) = srchomo0(:,4) + srchomo1(:,1,cnt)*(kk-0.5)
        ! &                      + srchomo1(:,4,cnt)/2
   
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
            finc0(:,:,1) = finc0(:,:,1) + bflx(:,:,j,k)
        ENDDO
      ENDDO
         
      DO k=kmin,kmax
        DO i=imin, imax
          finc0(:,:,2) = finc0(:,:,2) + bfly(:,:,i,k)
        ENDDO
      ENDDO
         
      DO j=jmin, jmax
        DO i=imin, imax
           finc0(:,:,3) = finc0(:,:,3) + bflz(:,:,i,j)
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
              finc1(:,:,1,cnt) = finc1(:,:,1,cnt) + bflx(:,:,j,k)
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
              finc1(:,:,2,cnt) = finc1(:,:,2,cnt) + bfly(:,:,i,k)
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
               finc1(:,:,3,cnt) = finc1(:,:,3,cnt) + bflz(:,:,i,j)
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
     
      SUBROUTINE SPLITBOUND0(nn, ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout0)
 
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
 
 
      SUBROUTINE SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)
   
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
                bflx(:,:,j,k) =  fout1(:,:,1,cnt)
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
                 bfly(:,:,i,k) = fout1(:,:,2,cnt)
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
                  bflz(:,:,i,j) = fout1(:,:,3,cnt)
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