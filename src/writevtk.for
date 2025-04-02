      SUBROUTINE VOLVTK(nx,ny,nz,dx,dy,dz,ox,oy,oz,delta3,flx,name)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx,ny,nz
      LOGICAL, INTENT(IN) :: dx,dy,dz
      REAL, INTENT(IN) :: flx(nx,ny,nz)
      REAL, INTENT(IN) :: ox,oy,oz
      REAL, INTENT(IN) :: delta3(3)
      CHARACTER(LEN=10), INTENT(IN) :: name

      INTEGER :: i, j, k, a, b,c

! Écriture du fichier VTK
      OPEN(unit=10, file=name, status="replace", action="write")

! En-tête VTK
      WRITE(10, '(A)') "# vtk DataFile Version 3.0"
      WRITE(10, '(A)') "VTK output from Fortran"
      WRITE(10, '(A)') "ASCII"
      WRITE(10, '(A)') "DATASET STRUCTURED_GRID"

      print *,"A"

      a = nx + 1
      b = ny + 1
      c = nz + 1
      IF (.NOT. dx) THEN
        a=1
      END IF
      IF (.NOT. dy) THEN
        b = 1
      END IF
      IF (.NOT. dz) THEN
        c = 1
      END IF


      WRITE(10, '(A,3I5)') "DIMENSIONS", a, b,c
      WRITE(10, '(A,1I5,A)') "POINTS", a*b*c," float"

      DO k = 0,c-1
        DO j = 0,b-1
          DO i = 0,a-1
            WRITE(10, '(3F10.5)') ox + 1.0*i*delta3(1)/nx,
     &                oy + 1.0*j*delta3(2)/ny, oz + 1.0*k*delta3(3)/nz
          END DO
        END DO
      END DO


      WRITE(10, '(A,1I5)') "CELL_DATA", nx*ny*nz
      WRITE(10, '(A)') "SCALARS CellValues float 1"
      WRITE(10, '(A)') "LOOKUP_TABLE default"

      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            WRITE(10, '(F10.5)', advance='no') flx(i,j,k)
        END DO
        END DO
        WRITE(10, *)  ! Nouvelle ligne
      END DO

      CLOSE(10)
      print *, "Fichier output.vtk écrit avec succès !"

      END SUBROUTINE VOLVTK
