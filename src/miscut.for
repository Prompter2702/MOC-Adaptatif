      SUBROUTINE ARRPOS(p,n,d,s)

!     Given the position "p" of an array element in a sequence
!     representing the array element order, returns the subscripts
!     for all array dimensions in the vector "s".
!     The sequence is that of common fortran convention (the leftmost
!     dimension changes the most rapidly).
!     Scalar "n" and vector "d" are input respectively as rank of
!     the array and size of each dimension.
!     Lower bound of each dimension is asumed to be 1.

      INTEGER , INTENT(IN) :: p,n
      INTEGER , INTENT(IN)  , DIMENSION(*) :: d
      INTEGER , INTENT(OUT) , DIMENSION(*) :: s

      INTEGER :: k,i

      k=p
      DO i=1,n
         s(i)=MOD(k,d(i))
         IF(s(i).EQ.0)s(i)=d(i)
         k=(k-s(i))/d(i)+1
      ENDDO

      END SUBROUTINE ARRPOS
