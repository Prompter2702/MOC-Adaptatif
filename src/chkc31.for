      SUBROUTINE CHKC31(lmax,x,y,z,t,ccof,icof,ecof,tcof)

!     Checks the conservation 3D bilinear schemes.

      REAL , DIMENSION(4,4) :: ccof
      REAL , DIMENSION(4,9) :: icof
      REAL , DIMENSION(9,4) :: ecof
      REAL , DIMENSION(9,9) :: tcof
      REAL(KIND=8) , DIMENSION(4,4) :: a,rv
      REAL(KIND=8) , DIMENSION(4,9) :: alpha,beta,rs

      REAL(KIND=8) :: lmax,x,y,z,t


      alpha=0

      alpha(1,1)=x
      alpha(2,1)=x
      alpha(3,2)=x
      alpha(4,3)=x
      alpha(1,4)= y
      alpha(2,5)= y
      alpha(3,4)= y
      alpha(4,6)= y
      alpha(1,7)=  z
      alpha(2,8)=  z
      alpha(3,9)=  z
      alpha(4,7)=  z

      alpha=(0.5/lmax)*alpha

      beta=0      

      beta(1,1)= x
      beta(2,1)=-x
      beta(3,2)= x
      beta(4,3)= x
      beta(1,4)=  y
      beta(2,5)=  y
      beta(3,4)= -y
      beta(4,6)=  y
      beta(1,7)=  z
      beta(2,8)=  z
      beta(3,9)=  z
      beta(4,7)= -z

      beta=(0.5/lmax)*beta

      a=0
      a(1,1)=t
      a(2,2)=t
      a(3,3)=t
      a(4,4)=t
      a(2,1)=-x
      a(3,1)=-y
      a(4,1)=-z
      a=(1/lmax)*a

      rv=MATMUL(alpha,ecof)+MATMUL(a,ccof)
      rs=MATMUL(alpha,tcof)+MATMUL(a,icof)-beta

      CALL CHKCMX(4,9,rv,rs)

      END
