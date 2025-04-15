      MODULE SRCCORR

      IMPLICIT NONE
      PUBLIC :: SRCCOR
    !   PRIVATE :: GAUSS_LU4
      CONTAINS

      SUBROUTINE SRCCOR(ng, ndir, delt3, mu,eta, ksi, srcm, sigt, err)
        
      IMPLICIT NONE

      INTEGER, PARAMETER :: nc = 4

      INTEGER, INTENT(IN) :: ng, ndir
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir), ksi(ndir)
      REAL, INTENT(IN) :: srcm(ng, ndir, nc), sigt(ng)
      REAL, INTENT(IN) :: delt3(3)
      REAL(KIND=8), INTENT(INOUT):: err(ng, ndir)

      REAL(KIND=8) :: l4(ndir)
      INTEGER :: d, nquad
      REAL(KIND=8) :: coeff(ng)

      nquad = 10
      err = 0.0D0
      
      CALL GAUSS_LU4(nquad,ndir, delt3, mu,eta, ksi, l4)

      coeff(:) = (sigt(:)**2)/24.0

      DO d=1,ndir
        err(:,d) = srcm(:,d,2)*mu(d) 
     &           + srcm(:,d,3)*eta(d) 
     &           + srcm(:,d,4)*ksi(d)

        err(:,d) = coeff(:)*err(:,d) * l4(d)
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
     &  (ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,2)) )*delt3(1)*0.5*
     &                                    ( sigt(:,1) + sigt(:,2))
     & +(ABS( srcm1(:,dir,1,3) - srcm1(:,dir,1,4)) )*delt3(1)*
     &                                    ( sigt(:,3) + sigt(:,4))
     & +(ABS( srcm1(:,dir,1,5) - srcm1(:,dir,1,6)) )*delt3(1)*0.5*
     &                                    ( sigt(:,5) + sigt(:,6))
     & +(ABS( srcm1(:,dir,1,7) - srcm1(:,dir,1,8)) )*delt3(1)*0.5*
     &                                    ( sigt(:,7) + sigt(:,8))
     & +(ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,3)) )*delt3(2)*0.5*
     &                                    ( sigt(:,1) + sigt(:,3))
     & +(ABS( srcm1(:,dir,1,2) - srcm1(:,dir,1,4)) )*delt3(2)*0.5*
     &                                    ( sigt(:,2) + sigt(:,4))
     & +(ABS( srcm1(:,dir,1,5) - srcm1(:,dir,1,7)) )*delt3(2)*0.5*
     &                                    ( sigt(:,5) + sigt(:,7))
     & +(ABS( srcm1(:,dir,1,6) - srcm1(:,dir,1,8)) )*delt3(2)*0.5*
     &                                    ( sigt(:,6) + sigt(:,8))
     & +(ABS( srcm1(:,dir,1,1) - srcm1(:,dir,1,5)) )*delt3(3)*0.5*
     &                                    ( sigt(:,1) + sigt(:,5))
     & +(ABS( srcm1(:,dir,1,2) - srcm1(:,dir,1,6)) )*delt3(3)*0.5*
     &                                    ( sigt(:,2) + sigt(:,6))
     & +(ABS( srcm1(:,dir,1,3) - srcm1(:,dir,1,7)) )*delt3(3)*0.5*
     &                                    ( sigt(:,3) + sigt(:,7))
     & +(ABS( srcm1(:,dir,1,4) - srcm1(:,dir,1,8)) )*delt3(3)*0.5*
     &                                    ( sigt(:,4) + sigt(:,8) )

      END DO

      err(:,:)= err(:,:)/(12.0*90.0)


      END SUBROUTINE SRC2LVL

      END MODULE SRCCORR