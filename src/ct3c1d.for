      SUBROUTINE CT3C1D(t,nb,x,y,e,i,j)

!     Computes a part of transmission matrix "t(nb,nb)" corresponding
!     to a diagonal block z->z in simplified linear
!     3D characteristics scheme.
!     Input variables are:
!     x    = ox*lmax (1st)
!     y    = oy*lmax (2nd)
!     e    = EXP(-tmax)
!     i,j  = address shift in "t" matrix

      REAL , DIMENSION(nb,nb):: t
      REAL(KIND=8)           :: e
      REAL                   :: x,y
      INTEGER                :: nb,i,j


      t(i+1,j+1)=0.25*(-2+x)*(-2+y)*e

      t(i+2,j+1)=0.125*x*(-2+x)*(-2+y)*e

      t(i+3,j+1)=0.125*y*(-2+x)*(-2+y)*e

      t(i+1,j+2)=-0.375*x*(-2+x)*(-2+y)*e

      t(i+2,j+2)=-0.125*(4-6*x+x**3)*(-2+y)*e

      t(i+3,j+2)=-0.1875*x*y*(-2+x)*(-2+y)*e

      t(i+1,j+3)=-0.375*y*(-2+x)*(-2+y)*e

      t(i+2,j+3)=-0.1875*x*y*(-2+x)*(-2+y)*e

      t(i+3,j+3)=-0.125*(-2+x)*(4-6*y+y**3)*e

      END SUBROUTINE
