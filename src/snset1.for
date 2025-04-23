      SUBROUTINE SNSET1(x,w,lgdble,n,lgrod)

!     Generates Gauss-Legendre or Double Gauss-Legendre 
!     quadrature set, where

!     "n"      is quadrature order (number of positive directions)
!     "x(1:n)" is direction cosine,
!     "w(1:n)" is wieght,
!     "lgdble" is .TRUE. if Double Gauss-Legendre formula,
!     "lgrod"  is .TRUE. if rod problem.

      IMPLICIT REAL (KIND=8) (A-H,O-Z)
      REAL (KIND=8) :: x(n),w(n)
      LOGICAL       :: lgdble,lgrod
      REAL (KIND=8) , PARAMETER :: eps = 3.D-14
      INTEGER       , PARAMETER :: nit = 30

      IF(lgrod)THEN
         m=1
         x(1)=1.
         w(1)=.5
      ELSE
         IF(lgdble)THEN
            m=(n+1)/2
            nn=n
         ELSE
            m=n
            nn=2*n
         ENDIF

!        Gauss Legendre quadrature set
!        (taken from W.H. Press et all, "Numerical Recipes",
!         Cambridge University Press 1989., Chapter 4.)

         DO i=1,m
            z=COS(3.14159265358979D0*(i-0.25)/(nn+0.5))

            DO it=1,nit
               p1=1.
               p2=0.
               DO j=1,nn
                  p3=p2
                  p2=p1
                  p1=((2*j-1)*z*p2-(j-1)*p3)/j
               ENDDO
               pp=nn*(z*p1-p2)/(z*z-1)
               z1=z
               z=z1-p1/pp
               IF(ABS(z-z1).LE.eps)EXIT
            ENDDO

            x(i)=z
            ww=1/((1-z*z)*pp*pp)
            w(i)=ww
         ENDDO

!        Double Gauss-Legendre quadrature

         IF(lgdble)THEN
            DO i=1,m
               aux=x(i)
               x(i)    =.5*(1.D0+aux)
               x(n-i+1)=.5*(1.D0-aux)
               w(i)    =.5*w(i)
               w(n-i+1)= w(i)
            ENDDO
         ENDIF
      ENDIF

      END
