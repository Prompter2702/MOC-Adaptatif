      SUBROUTINE INCEXP(e,f,x,nmax)

!     Computes the functions
!     f_{n}(x) = (1 - Sum[x^i/i!,{i,0,n}] Exp[-x]) / x^(n+1)
!     (n=0,nmax) using recursion or series expansion.

      IMPLICIT REAL(KIND=8)(A-Z)

      REAL(KIND=8) , INTENT(OUT)                  :: e
      REAL(KIND=8) , INTENT(OUT) , DIMENSION(0:*) :: f
      REAL(KIND=8) , INTENT(IN) :: x
      INTEGER      , INTENT(IN) :: nmax
      INTEGER :: n,i

      INTEGER , PARAMETER :: degr=6
      INTEGER , PARAMETER :: nfun=9

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co00 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co00 =
     &   (/
     &    1.000000000000000D+00,   -5.000000000000000D-01,
     &    1.666666666666666D-01,   -4.166666666666666D-02,
     &    8.333333333333330D-03,   -1.388888888888888D-03,
     &    1.984126984126984D-04
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co01 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co01 =
     &   (/
     &    5.000000000000000D-01,   -3.333333333333333D-01,
     &    1.250000000000000D-01,   -3.333333333333333D-02,
     &    6.944444444444444D-03,   -1.190476190476190D-03,
     &    1.736111111111111D-04
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co02 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co02 =
     &   (/
     &    1.666666666666666D-01,   -1.250000000000000D-01,
     &    5.000000000000000D-02,   -1.388888888888888D-02,
     &    2.976190476190476D-03,   -5.208333333333332D-04,
     &    7.716049382716049D-05
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co03 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co03 =
     &   (/
     &    4.166666666666666D-02,   -3.333333333333333D-02,
     &    1.388888888888888D-02,   -3.968253968253968D-03,
     &    8.680555555555560D-04,   -1.543209876543209D-04,
     &    2.314814814814814D-05
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co04 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co04 =
     &   (/
     &    8.333333333333330D-03,   -6.944444444444444D-03,
     &    2.976190476190476D-03,   -8.680555555555560D-04,
     &    1.929012345679012D-04,   -3.472222222222222D-05,
     &    5.260942760942761D-06
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co05 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co05 =
     &   (/
     &    1.388888888888888D-03,   -1.190476190476190D-03,
     &    5.208333333333332D-04,   -1.543209876543209D-04,
     &    3.472222222222222D-05,   -6.313131313131313D-06,
     &    9.645061728395060D-07
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co06 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co06 =
     &   (/
     &    1.984126984126984D-04,   -1.736111111111111D-04,
     &    7.716049382716049D-05,   -2.314814814814814D-05,
     &    5.260942760942761D-06,   -9.645061728395060D-07,
     &    1.483855650522317D-07
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co07 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co07 =
     &   (/
     &    2.480158730158729D-05,   -2.204585537918871D-05,
     &    9.920634920634921D-06,   -3.006253006253006D-06,
     &    6.889329805996473D-07,   -1.271876271876271D-07,
     &    1.968379944570420D-08
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co08 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co08 =
     &   (/
     &    2.755731922398589D-06,   -2.480158730158730D-06,
     &    1.127344877344877D-06,   -3.444664902998236D-07,
     &    7.949226699226700D-08,   -1.476284958427815D-08,
     &    2.296443268665490D-09
     &    /)

!      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr) :: co09 =
      REAL(KIND=8) , PARAMETER , DIMENSION(degr+1) :: co09 =
     &   (/
     &    2.755731922398588D-07,   -2.505210838544171D-07,
     &    1.148221634332745D-07,   -3.532989644100754D-08,
     &    8.201583102376750D-09,   -1.530962179110327D-09,
     &    2.392128404859886D-10
     &    /)

      REAL(KIND=8) , PARAMETER , DIMENSION(0:degr,0:nfun) :: coef =
     &    RESHAPE ((/
     &    co00(:),co01(:),co02(:),co03(:),co04(:),co05(:),
     &    co06(:),co07(:),co08(:),co09(:)

     &    /), (/1+degr,1+nfun/))

      REAL(KIND=8) , PARAMETER :: one = 1.0D+0
      REAL(KIND=8) , PARAMETER :: eps = 1.0D-1

      e=EXP(-x)

      IF (x.GT.eps) THEN

         f(0)=(one-e)/x
         fact=one
         DO i=1,nmax
            fact=fact*i
            f(i)=(f(i-1)-e/fact)/x
         ENDDO

      ELSE

         DO n=0,nmax
            f(n)=coef(degr,n)
            DO i=degr-1,0,-1
               f(n)=f(n)*x+coef(i,n)
            ENDDO
         ENDDO

      ENDIF

      END SUBROUTINE
