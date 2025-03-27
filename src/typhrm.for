      MODULE SNQHRM

      IMPLICIT NONE
!     Coefficients for spherical harmonics expansion

      PUBLIC :: SNQDLFT

      CONTAINS

      SUBROUTINE SNQDLFT(ndir,nd,nani,ndim,
     &                 nhrm,mu,eta,ksi,w,sphr)
!       'LevelSymmetric (Default)'
        INTEGER, INTENT(IN) :: ndir,ndim, nani,nd, nhrm
        !INTEGER, INTENT(INOUT) :: nhrm
        REAL(KIND=8) , INTENT(OUT) :: mu(ndir),eta(ndir),ksi(ndir)
        REAL(KIND=8) , INTENT(OUT) :: w(ndir)
        REAL, INTENT(INOUT) :: sphr(nhrm,nd)

        INTEGER , DIMENSION(8) :: sx=(/ 1,-1,-1, 1, 1,-1,-1, 1/)
        INTEGER , DIMENSION(8) :: sy=(/ 1, 1,-1,-1, 1, 1,-1,-1/)
        INTEGER , DIMENSION(8) :: sz=(/ 1, 1, 1, 1,-1,-1,-1,-1/)
        REAL (KIND=8) , PARAMETER :: pi=3.1415926535897932D0
        REAL (KIND=8) :: phi
        INTEGER :: dd,o,d,noct

        IF(ndir/=3)STOP 'SNQDFLT: ndir must be 3' 
        IF(ndim/=3)STOP 'SNQDFLT: ndim must be 3'

        mu  = (/.30163872, .30163872, .90444905/)
        eta = (/.90444905, .30163872, .30163872/)
        ksi = (/.30163872, .90444905, .30163872/)
        w=1.0D+0/24.

        noct=2**ndim

        SELECT CASE (ndim)
        CASE (1:2)
           STOP "Dimension 1 or not implemented, SNQDLFT"
        CASE (3)
           IF(nhrm/=(1+nani)**2)STOP"ERROR Dimension nhrm, SNQDLFT"
        END SELECT

        SELECT CASE (ndim)
        CASE (2:3)
          dd=0
          DO o=1,noct
            DO d=1,ndir
              dd=dd+1
              phi=ATAN(sy(o)*ksi(d)/eta(d))
              if(sy(o).lt.0)phi=phi+pi
              if(sz(o).lt.0)phi=-phi
              CALL SPHARF(ndim,nani,DBLE(sx(o))*mu(d),
     &                    phi,sphr(1,dd))
              ENDDO
           ENDDO

        CASE(1)
           STOP "Dimension 1 not implemented, typhrm" 

        END SELECT
        

      END SUBROUTINE SNQDLFT


C     Last change:  I    22 Mar 2009   11:24 pm
      SUBROUTINE SPHARF(ndim,lmax,x,phi,coef)
    !     Given the anisotropy order "lmax", direction cosine "x",
    !     and azimuthal angle "phi", computes the coefficients "coef"
    !     for the spherical harmonics expansion of the source terms.
    !     The "coef" array must be dimensioned to:
    !        (lmax+1)**2         for 3D problem,
    !        (lmax+1)(lmax+2)/2  for 2D problem.
    !     Terms are arranged as
    !     3D:
    !     ((C{lm}*P{lm}(x),
    !       C{lm}*P{lm}(x)*Cos(m phi),
    !       C{lm}*P{lm}(x)*Sin(m phi),  m = 0 , l ), l = 0 , lmax )
    !     2D:
    !     ((C{lm}*P{lm}(x),
    !       C{lm}*P{lm}(x)*Cos(m phi),  m = 0 , l ), l = 0 , lmax )
    !     where
    !     P{lm}(x) = (-1)**m (1-x**2)**(m/2) Derivative{m}LegendreP{l}
    !     C{lm}    = (-1)**m Sqrt((2-KronDelta(m,0)(l-m)!/(l+m)!)
    !     Computation of P{lm} uses the following recursion:
    !     P{m+1 m}(x) = (2m+1)*x*P{mm}(x)
    !     with the starting value
    !     P{mm}(x) = (-1)**m (2m+1)!!(1-x*x)**(m/2)

          IMPLICIT NONE  

          INTEGER, INTENT(IN) :: ndim, lmax
          REAL(KIND=8), INTENT(IN) :: x, phi
          REAL , DIMENSION(*) :: coef   
          INTEGER :: h,mm,m,l
          REAL(KIND=8) :: somx2,fact,faci,two,cosa,sina,pmm,facp,facm,
     &                    cfac,plm,pm1,pp1
    
          IF(lmax.GT.0)somx2=SQRT((1.-x)*(1.+x))
          fact=-1
          faci= 1
          two=1
          cosa=1
          sina=1
          mm=0
          pmm=1
          DO m=0,lmax
             IF(m.GT.0)THEN
                faci=faci*m*(4*m-2)
                pmm=pmm*fact*somx2     !pmm=-pmm*fact*somx2
                cosa=COS(m*phi)
                IF(ndim.EQ.3) sina=SIN(m*phi)
             ENDIF
             fact=fact+2
             facp=faci
             cfac=SQRT(two/facp)
             IF(ndim.EQ.2)THEN  
                mm=mm+m+1
                coef(mm)=cfac*pmm*cosa
             ELSE
                mm=mm+2*m+1
                IF(m.GT.0)coef(mm-1)=cfac*pmm*cosa
                          coef(mm  )=cfac*pmm*sina
             ENDIF
             IF(m.LT.lmax)THEN
                facp=facp*(2*m+1)
                facm=1
                cfac=SQRT(two/facp)
                pm1=pmm
                plm=x*(2*m+1)*pmm
                IF(ndim.EQ.2)THEN
                   h=mm+m+1  
                   coef(h)=cfac*plm*cosa
                ELSE
                   h=mm+2*m+1
                   IF(m.GT.0)coef(h-1)=cfac*plm*cosa
                             coef(h  )=cfac*plm*sina
                ENDIF
                DO l=m+2,lmax
                   facp=facp*(l+m)
                   facm=facm*(l-m)
                   cfac=SQRT(two*facm/facp)
                   pp1=(x*(2*l-1)*plm-(l+m-1)*pm1)/(l-m)
                   pm1=plm
                   plm=pp1
                   IF(ndim.EQ.2)THEN
                      h=h+l
                      coef(h)=cfac*plm*cosa
                   ELSE
                      h=h+2*l-1
                      IF(m.GT.0)coef(h-1)=cfac*plm*cosa
                                coef(h  )=cfac*plm*sina
                   ENDIF
                ENDDO
             ENDIF
             two=2
          ENDDO
    
          END
    
      END MODULE SNQHRM