      MODULE FLGBCD

!     Boundary conditions flags.

      INTEGER , PARAMETER :: Vacuum              =  1
      INTEGER , PARAMETER :: Translation         =  2
      INTEGER , PARAMETER :: HalfPiRotation      =  3
      INTEGER , PARAMETER :: PiRotation          =  4
      INTEGER , PARAMETER :: SpecularReflection  =  5
      INTEGER , PARAMETER :: IsotropicReflection =  6
      INTEGER , PARAMETER :: SingleAlbedo        =  7
      INTEGER , PARAMETER :: DiagonalAlbedo      =  8
      INTEGER , PARAMETER :: MultigroupAlbedo    =  9
      INTEGER , PARAMETER :: TranslationOfHalfPi = 10
      INTEGER , PARAMETER :: AxialSymmetry       = 11
      INTEGER , PARAMETER :: DiagonalSymmetry    = 12
      INTEGER , PARAMETER :: HalfPiRotationN     = 13
      INTEGER , PARAMETER :: TranslationOfHalfPiN= 14
      INTEGER , PARAMETER :: AxialHalfPiRotation = 15
      
      
      INTEGER , PARAMETER :: nbcdw = 15

      CHARACTER(LEN=32) , PARAMETER , DIMENSION(nbcdw) :: bcdw = (/
     &'Vacuum              ','Translation         ',
     &'HalfPiRotation      ','PiRotation          ',
     &'SpecularReflection  ','IsotropicReflection ',
     &'SingleAlbedo        ','DiagonalAlbedo      ',
     &'MultigroupAlbedo    ','TranslationOfHalfPi ',
     &'AxialSymmetry       ','DiagonalSymmetry    ',
     &'HalfPiRotationN     ','TranslationOfHalfPiN',
     &'AxialHalfPiRotation '
     & /)

      END MODULE
