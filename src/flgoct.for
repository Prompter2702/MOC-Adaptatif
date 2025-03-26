      MODULE FLGOCT

!     Arrival octant after reflection
!     on respectively x-, y-, z-plane per departure octant.

      INTEGER , DIMENSION(8) , PARAMETER :: xoct=(/2,1,4,3,6,5,8,7/)
      INTEGER , DIMENSION(8) , PARAMETER :: yoct=(/4,3,2,1,8,7,6,5/)
      INTEGER , DIMENSION(8) , PARAMETER :: zoct=(/5,6,7,8,1,2,3,4/)

      INTEGER , PARAMETER , DIMENSION(8,3) :: roct = RESHAPE(
     &   (/2,1,4,3,6,5,8,7, 4,3,2,1,8,7,6,5 , 5,6,7,8,1,2,3,4/),(/8,3/))


!     Arrival octant after exiting a domain having Pi/2 rotational
!     symmetry (arround z-axis only),
!     on respectively x-, y-, z-plane per departure octant.
!     "zoc4" is dummy.

      INTEGER , DIMENSION(8) , PARAMETER :: xoc4=(/4,1,2,3,8,5,6,7/)
      INTEGER , DIMENSION(8) , PARAMETER :: yoc4=(/2,3,4,1,6,7,8,5/)
      INTEGER , DIMENSION(8) , PARAMETER :: zoc4=(/1,1,1,1,1,1,1,1/)
      
      INTEGER , DIMENSION(3) , PARAMETER :: nino=(/1,2,4/) ! number of inc. onctant for 1D, 2D, 3D

      INTEGER , PARAMETER , DIMENSION(8,3) :: roc4 = RESHAPE(
     &   (/4,1,2,3,8,5,6,7, 2,3,4,1,6,7,8,5 , 1,1,1,1,1,1,1,1/),(/8,3/))
     
      INTEGER , DIMENSION(4,2) , PARAMETER ::  octx = RESHAPE(
     &   (/1,4,5,8,3,2,7,6/),(/4,2/))
      INTEGER , DIMENSION(4,2) , PARAMETER ::  octy = RESHAPE(
     &   (/1,2,5,6,3,4,7,8/),(/4,2/))
      INTEGER , DIMENSION(4,2) , PARAMETER ::  octz = RESHAPE(
     &   (/1,2,3,4,5,6,7,8/),(/4,2/))
      INTEGER , DIMENSION(4,2,3), PARAMETER ::  octd = RESHAPE(
     &  (/1,4,5,8, 3,2,7,6,   
     &    1,2,5,6, 3,4,7,8,
     &    1,2,3,4, 5,6,7,8/),(/4,2,3/))

      END MODULE
