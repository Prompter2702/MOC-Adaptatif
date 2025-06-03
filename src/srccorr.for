      MODULE SRCCORR

      IMPLICIT NONE
      PUBLIC :: SRCCOR
      PUBLIC :: ERRSURF
    !   PRIVATE :: GAUSS_LU4
      CONTAINS

      SUBROUTINE SRCCOR(ng, ndir, delt3, mu,eta, ksi,
     &                       srcm, sigt, err,pdslu4)
        
      IMPLICIT NONE

      INTEGER, PARAMETER :: nc = 4

      INTEGER, INTENT(IN) :: ng, ndir
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir), ksi(ndir)
      REAL, INTENT(IN) :: srcm(ng, ndir, nc), sigt(ng)
      REAL, INTENT(IN) :: delt3(3)
      REAL(KIND=8), INTENT(IN) :: pdslu4(ndir)
      REAL(KIND=8), INTENT(INOUT):: err(ng, ndir)

      REAL(KIND=8) :: l4(ndir)
      INTEGER :: d, nquad
      REAL(KIND=8) :: coeff(ng)

      nquad = 10
      err = 0.0D0
      
    !   CALL GAUSS_LU4(nquad,ndir, delt3, mu,eta, ksi, l4)

      coeff(:) = (sigt(:)**2)/72

      DO d=1,ndir
        err(:,d) = ABS( srcm(:,d,2)*mu(d) )
     &           + ABS(srcm(:,d,3)*eta(d) )
     &           + ABS(srcm(:,d,4)*ksi(d) )

        err(:,d) = coeff(:)*err(:,d) * pdslu4(d)*(delt3(1)**5 )
        ! err(:,d) = coeff(:)*err(:,d) * pdslu4(d)*(delt3(1)**3 )
      END DO

      END SUBROUTINE SRCCOR



      SUBROUTINE GAUSS_LU4(nquad,ndir, delt3, mu,eta,ksi, val)
    
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nquad, ndir
      REAL, INTENT(IN)    :: delt3(3)
      REAL(KIND=8), INTENT(IN)    :: mu(ndir), eta(ndir), ksi(ndir)
      REAL(KIND=8), INTENT(INOUT)   :: val(ndir)
      INTEGER :: i,j
      REAL(KIND=8) :: dx,dy,dz ! size of cube
      REAL(KIND=8) :: pds1, pds2, pds3
      
      
      dx = delt3(1)/nquad
      dy = delt3(2)/nquad
      dz = delt3(3)/nquad

      pds1 = dx*dy
      pds2 = dx*dz
      pds3 = dy*dz

      val(:) = 0.0D0

      ! x plan
      DO i=1,nquad
        DO j= 1,nquad
          val(:) = val(:) + pds1*(min(i*dx/ABS(mu(:)),
     &                     j*dy/ABS(eta(:)), delt3(3)/ABS(ksi(:))))**4
          val(:) = val(:) + pds2*(min(i*dx/ABS(mu(:)),
     &                     j*dy/ABS(eta(:)), delt3(3)/ABS(ksi(:))))**4
          val(:) = val(:) + pds3*(min(delt3(1)/ABS(mu(:)),
     &                     i*dy/ABS(eta(:)), j*dz/ABS(ksi(:))))**4
        ENDDO
      ENDDO
  
      END SUBROUTINE GAUSS_LU4


      SUBROUTINE SRC2LVL(ng, ndir, srcm1, sigt, delt3, err)
        
      IMPLICIT NONE

      INTEGER, PARAMETER :: n8 = 8, nc = 4
      INTEGER, INTENT(IN) :: ng, ndir
      REAL, INTENT(IN) :: srcm1(ng,ndir,nc,n8), sigt(ng,n8), delt3(3)
      REAL(KIND=8), INTENT(INOUT) :: err(ng,ndir)
      INTEGER :: dir

      err = 0.0


      DO dir=1,ndir 


      err(:,dir) = 
     &  (ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,2)) )*(delt3(1)**1)*
     &                                    ( sigt(:,1) + sigt(:,2))/2
     & +(ABS( srcm1(:,dir,1,3) - srcm1(:,dir,1,4)) )*delt3(1)*
     &                                    ( sigt(:,3) + sigt(:,4))/2
     & +(ABS( srcm1(:,dir,1,5) - srcm1(:,dir,1,6)) )*(delt3(1)**1)*
     &                                    ( sigt(:,5) + sigt(:,6))/2
     & +(ABS( srcm1(:,dir,1,7) - srcm1(:,dir,1,8)) )*(delt3(1)**1)*
     &                                    ( sigt(:,7) + sigt(:,8))/2
     & +(ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,3)) )*(delt3(2)**2)*
     &                                    ( sigt(:,1) + sigt(:,3))/2
     & +(ABS( srcm1(:,dir,1,2) - srcm1(:,dir,1,4)) )*(delt3(2)**1)*
     &                                    ( sigt(:,2) + sigt(:,4))/2
     & +(ABS( srcm1(:,dir,1,5) - srcm1(:,dir,1,7)) )*(delt3(2)**1)*
     &                                    ( sigt(:,5) + sigt(:,7))/2
     & +(ABS( srcm1(:,dir,1,6) - srcm1(:,dir,1,8)) )*(delt3(2)**1)*
     &                                    ( sigt(:,6) + sigt(:,8))/2
     & +(ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,5)) )*(delt3(3)**1)*
     &                                    ( sigt(:,1) + sigt(:,5))/2
     & +(ABS( srcm1(:,dir,1,2) - srcm1(:,dir,1,6)) )*(delt3(3)**1)*
     &                                    ( sigt(:,2) + sigt(:,6))/2
     & +(ABS( srcm1(:,dir,1,3) - srcm1(:,dir,1,7)) )*(delt3(3)**1)*
     &                                    ( sigt(:,3) + sigt(:,7))/2
     & +(ABS( srcm1(:,dir,1,4) - srcm1(:,dir,1,8)) )*(delt3(3)**1)*
     &                                    ( sigt(:,4) + sigt(:,8))/2

      END DO

      err(:,:)= err(:,:)/1080 ! 1080 = 12*90

      END SUBROUTINE SRC2LVL

      SUBROUTINE ERRSURF(nn, nx,ny,nz,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   flx0, bflx,bfly,bflz,tcof,errbnd)

        IMPLICIT NONE

        INTEGER,PARAMETER :: nbd = 9, nb = 3

        INTEGER, INTENT(IN) :: nn,nx,ny,nz
        INTEGER, INTENT(IN) :: imin, imax, jmin, jmax, kmin, kmax
        REAL, INTENT(IN) :: flx0(nn,nb,nb)
        REAL, INTENT(IN) :: bflx(nn,nb,ny,nz)
        REAL, INTENT(IN) :: bfly(nn,nb,nx,nz)
        REAL, INTENT(IN) :: bflz(nn,nb,nx,ny)

        REAL, INTENT(IN) :: tcof(nn,nbd,nbd)   
        REAL, INTENT(INOUT) :: errbnd(3)
        REAL :: err(nn,nb), errfin(nn,nb)

        INTEGER :: i,j,k,Di,Dj,Dk
        Di = imax - imin + 1
        Dj = jmax - jmin + 1
        Dk = kmax - kmin + 1      
   
        err = 0.0
          DO k=kmin, kmax
          DO j=jmin, jmax
            err(:,1) = err(:,1) + ABS(bflx(:,1,j,k) - flx0(:,1,1)
     &           - 3.0*flx0(:,2,1)*(2*j-jmin-jmax)/Dj
     &           - 3.0*flx0(:,3,1)*(2*k-kmin-kmax)/Dk)
          ENDDO
          ENDDO
  
          DO k=kmin, kmax
          DO i=imin, imax
            err(:,2) = err(:,2) + ABS(bfly(:,1,i,k) - flx0(:,1,2)
     &           - 3.0*flx0(:,2,2)*(2.0*i-imax-imin)/Di
     &           - 3.0*flx0(:,3,2)*(2.0*k-kmax-kmin)/Dk)
          ENDDO
          ENDDO

          DO j=jmin, jmax
          DO i=imin, imax
            err(:,3) = err(:,3) + ABS(bflz(:,1,i,j) - flx0(:,1,3)
     &           - 3.0*flx0(:,2,3)*(2*i-imax-imin)/Di
     &           - 3.0*flx0(:,3,3)*(2*j-jmax-jmin)/Dj)
          ENDDO
          ENDDO

        err(:,1) = err(:,1)/(Dj*Dk)
        err(:,2) = err(:,2)/(Di*Dk)
        err(:,3) = err(:,3)/(Di*Dj)

        errfin = 0.0D0
        
        errfin(:,1) = tcof(:,1,1)*err(:,1)
     &           +    tcof(:,1,4)*err(:,2)
     &           +    tcof(:,1,7)*err(:,3)

        errfin(:,2) = tcof(:,2,1)*err(:,1)
     &           +    tcof(:,2,4)*err(:,2)
     &           +    tcof(:,2,7)*err(:,3)

        errfin(:,3) = tcof(:,3,1)*err(:,1)
     &           +    tcof(:,3,4)*err(:,2)
     &           +    tcof(:,3,7)*err(:,3)

        ! errbnd = SUM(ABS(errfin))/(nn*3)
        errbnd(1) = SUM( ABS(err(:,1)) )/nn
        errbnd(2) = SUM( ABS(err(:,2)) )/nn
        errbnd(3) = SUM( ABS(err(:,3)) )/nn

      END SUBROUTINE ERRSURF

      END MODULE SRCCORR