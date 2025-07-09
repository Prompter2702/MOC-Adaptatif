      MODULE OCTREE
      IMPLICIT NONE

      TYPE :: TabPtrChild
        TYPE(CellOctree), POINTER :: ptrChild => NULL()
      END TYPE TabPtrChild
      
      
      TYPE :: CellOctree
        TYPE(TabPtrChild) :: parent
        TYPE(TabPtrChild) :: children(8)

      END TYPE CellOctree


      !-----------------------------------------------------------------

      CONTAINS

      FUNCTION NEW_ROOT_OCTREE() result(root)
        IMPLICIT NONE
        TYPE(CellOctree) :: root
        INTEGER :: i
        root%parent%ptrChild => NULL()
        ! Initialize children pointers to NULL
        DO i = 1, 8
          root%children(i)%ptrChild => NULL()
        END DO
      END FUNCTION NEW_ROOT_OCTREE



      SUBROUTINE NEW_CELL_OCTREE(parent_cell, index, cell)

        IMPLICIT NONE
        TYPE(CellOctree), TARGET, INTENT(INOUT) :: parent_cell
        INTEGER, INTENT(IN) :: index
        TYPE(CellOctree), POINTER :: cell
        INTEGER :: i

        ! Alloue la cellule ici si ce n'est pas déjà fait
        IF (.NOT. ASSOCIATED(cell)) THEN
          ALLOCATE(cell)
        END IF

        cell%parent%ptrChild => parent_cell

        parent_cell%children(index)%ptrChild => cell

        DO i = 1, 8
          cell%children(i)%ptrChild => NULL()
        END DO
      END SUBROUTINE NEW_CELL_OCTREE

      !-----------------------------------------------------------------
      
      SUBROUTINE ADD_CHILDREN(cell)
          IMPLICIT NONE
          TYPE(CellOctree), INTENT(INOUT) :: cell
          TYPE(CellOctree), POINTER :: new_child
          INTEGER :: i
      
          DO i = 1, 8
              IF (.NOT. ASSOCIATED(cell%children(i)%ptrChild)) THEN
                  ALLOCATE(new_child)
                  CALL NEW_CELL_OCTREE(cell, i, new_child)
              END IF
          END DO
      END SUBROUTINE ADD_CHILDREN


      !-----------------------------------------------------------------

      RECURSIVE SUBROUTINE DELETE_CELL(cell)
        IMPLICIT NONE
        TYPE(CellOctree), POINTER :: cell
        INTEGER :: i
        
        ! Recursively delete children first
        DO i = 1, 8
          IF (ASSOCIATED(cell%children(i)%ptrChild)) THEN
            CALL DELETE_CELL(cell%children(i)%ptrChild)
          END IF
        END DO
        
        ! Now delete the current cell
        IF (ASSOCIATED(cell%parent%ptrChild)) THEN
          cell%parent%ptrChild => NULL()  ! Remove pointer to this cell from parent
        END IF
        
        ! Deallocate the cell itself
        DEALLOCATE(cell, STAT=i)
        IF (i /= 0) PRINT *, "Error deallocating cell"

      END SUBROUTINE 

      SUBROUTINE IS_LEAF(cell,isleaf)

      IMPLICIT NONE

      TYPE(CellOctree) :: cell
      LOGICAL, INTENT(OUT) :: isleaf

      isleaf = .FALSE.

      IF (.NOT. ASSOCIATED(cell%children(1)%ptrChild)) THEN
        isleaf = .TRUE.
      ENDIF

      END SUBROUTINE


      !-----------------------------------------------------------------



      END MODULE OCTREE


