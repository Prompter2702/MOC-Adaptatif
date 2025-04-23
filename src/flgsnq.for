      MODULE FLGSNQ

!     Sn quadrature flags

      INTEGER , PARAMETER :: GaussLegendre                = 1
      INTEGER , PARAMETER :: DoubleGaussLegendre          = 2
      INTEGER , PARAMETER :: ChebyshevLegendre            = 3
      INTEGER , PARAMETER :: ChebyshevDoubleLegendre      = 4
      INTEGER , PARAMETER :: ChebyshevTimesLegendre       = 5
      INTEGER , PARAMETER :: ChebyshevTimesDoubleLegendre = 6
      INTEGER , PARAMETER :: LevelSymmetric               = 7
      INTEGER , PARAMETER :: BuiltIn                      = 8
      INTEGER , PARAMETER :: MPNRectangularElements       = 9
      INTEGER , PARAMETER :: MPNTriangularElements        = 10

      INTEGER , PARAMETER :: nsnqw=10
      INTEGER , PARAMETER :: nmpnkw=2

      CHARACTER (LEN=32) , PARAMETER , DIMENSION(nsnqw) :: snqw = (/
     &   'GaussLegendre               ','DoubleGaussLegendre         ',
     &   'ChebyshevLegendre           ','ChebyshevDoubleLegendre     ',
     &   'ChebyshevTimesLegendre      ','ChebyshevTimesDoubleLegendre',
     &   'LevelSymmetric              ','BuiltIn                     ',
     &   'MPNRectangularElements      ','MPNTriangularElements       '/)

      CHARACTER (LEN=32) , PARAMETER, DIMENSION(nmpnkw) :: mpnkw = (/
     &   'NumPolar                    ','NumAzimuthal                '/)


      END MODULE
