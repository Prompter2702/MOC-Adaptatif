      SUBROUTINE CE3C1Z(e,nb,z,x,y,lmax,f,i,ii,jj)

!     Computes a part of escape matrix "e(nb,nc)" corresponding
!     to one exiting surface in simplified linear
!     3D characteristics scheme.
!     Input variables are:
!     z    = oz*lmax (exiting)
!     x    = ox*lmax (1st)
!     y    = oy*lmax (2nd)
!     f(n) = (1 - Sum[tmax^i/i!,{i,0,n}] Exp[-tmax]) / tmax^(n+1)
!     i,j  = address shift in "e" matrix
!     ii,jj = moments ordering (side dependent)

      REAL , DIMENSION(nb,*):: e
      REAL        :: x,y,z
      INTEGER     :: nb,i
      INTEGER , DIMENSION(*) :: ii,jj

      REAL(KIND=8)                   :: lmax
      REAL(KIND=8) , DIMENSION(0:*)  :: f


      e(i+ii(1),jj(1))=0.5*lmax*(2*f(0)-(x+y)*f(1)+x*y*f(2))

      e(i+ii(1),jj(2))=0.25*3*lmax*x*
     &        (-2*f(1)+2*(x+y)*f(2)-3*x*y*f(3))

      e(i+ii(1),jj(3))=0.25*3*lmax*y*
     &        (-2*f(1)+2*(x+y)*f(2)-3*x*y*f(3))

      e(i+ii(1),jj(4))=0.5*3*lmax*
     &        (2*f(0)+(-x-y-2*z)*f(1)+(2*y*z+x*(y+2*z))*f(2)-
     &        3*x*y*z*f(3))

      e(i+ii(2),jj(1))=0.25*lmax*x*(2*f(1)-2*(x+y)*f(2)+3*x*y*f(3))

      e(i+ii(2),jj(2))=0.5*lmax*
     &        (2*f(0)+(-3*x-y)*f(1)+3*x*y*f(2)+3*x**3*f(3)-
     &        6*x**3*y*f(4))

      e(i+ii(2),jj(3))=0.25*3*lmax*x*y*
     &        (-2*f(2)+3*(x+y)*f(3)-6*x*y*f(4))

      e(i+ii(2),jj(4))=0.25*3*lmax*x*
     &        (2*f(1)-2*(x+y+2*z)*f(2)+(6*y*z+3*x*(y+2*z))*f(3)-
     &        12*x*y*z*f(4))

      e(i+ii(3),jj(1))=0.25*lmax*y*(2*f(1)-2*(x+y)*f(2)+3*x*y*f(3))

      e(i+ii(3),jj(2))=0.25*3*lmax*x*y*(-2*f(2)+3*(x+y)*f(3)-
     &        6*x*y*f(4))

      e(i+ii(3),jj(3))=0.5*lmax*
     &        (2*f(0)+(-x-3*y)*f(1)+3*x*y*f(2)+3*y**3*f(3)-
     &        6*x*y**3*f(4))

      e(i+ii(3),jj(4))=0.25*3*lmax*y*
     &        (2*f(1)-2*(x+y+2*z)*f(2)+(6*y*z+3*x*(y+2*z))*f(3)-
     &        12*x*y*z*f(4))
      
      END SUBROUTINE
