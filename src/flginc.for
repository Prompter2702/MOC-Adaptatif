      MODULE FLGINC

!     Index of incoming boundary, perpendiclular to respectively
!     x-, y-, and z-axis, per octant.
!     1 means lower, 2 means upper boundary coordinate.

      INTEGER , DIMENSION(8) , PARAMETER :: xinc=(/1,2,2,1,1,2,2,1/)
      INTEGER , DIMENSION(8) , PARAMETER :: yinc=(/1,1,2,2,1,1,2,2/)
      INTEGER , DIMENSION(8) , PARAMETER :: zinc=(/1,1,1,1,2,2,2,2/)

      INTEGER , PARAMETER , DIMENSION(8,3) :: oinc = RESHAPE(
     &   (/1,2,2,1,1,2,2,1,  1,1,2,2,1,1,2,2,  1,1,1,1,2,2,2,2/),
     &   (/8,3/))

!     Sweeping order for the Octant index  
      INTEGER , PARAMETER , DIMENSION(8) :: olst3D=(/7,8,6,3,5,4,2,1/)
      INTEGER , PARAMETER , DIMENSION(4) :: olst2D=(/3,4,2,1/)
      INTEGER , PARAMETER , DIMENSION(2) :: olst1D=(/2,1/)

      END MODULE
