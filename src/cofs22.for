      SUBROUTINE COFS22(sx,sy,nc,nbd,sc,si,se,st)


!     Computes signs of collision "sc", incoming "si", escape "se" and
!     transmission "st" coefficients for full linear schemes
!     in one octant defined by its direction cosines
!     "sx" and "sy".


      INTEGER :: sx,sy,nc,nbd
      INTEGER :: sc(nc,nc),si(nc,nbd),se(nbd,nc),st(nbd,nbd)

      INTEGER , PARAMETER , DIMENSION(2) :: vx = (/ 2, 4/)
      INTEGER , PARAMETER , DIMENSION(2) :: vy = (/ 3, 4/)

      INTEGER :: i,j


      CALL SETINT( nc*nc ,1,sc,1)
      CALL SETINT( nc*nbd,1,si,1)
      CALL SETINT( nc*nbd,1,se,1)
      CALL SETINT(nbd*nbd,1,st,1)


      DO i=1,4
      DO j=1,2
         sc(i,vx(j))=sx*sc(i,vx(j))
         sc(i,vy(j))=sy*sc(i,vy(j))
         sc(vy(j),i)=sy*sc(vy(j),i)
         sc(vx(j),i)=sx*sc(vx(j),i)
      ENDDO
      ENDDO

      DO i=1,4
      DO j=1,2
         si(vx(j),i)=sx*si(vx(j),i)
         si(vy(j),i)=sy*si(vy(j),i)
      ENDDO
      ENDDO
      DO i=1,4
         si(i,4)=sx*si(i,4)
         si(i,2)=sy*si(i,2)
      ENDDO

      DO i=1,4
      DO j=1,4
         se(j,i)=si(i,j)
      ENDDO
      ENDDO

      DO i=1,4
         st(i,4)=sx*st(i,4)
         st(4,i)=sx*st(4,i)
         st(i,2)=sy*st(i,2)
         st(2,i)=sy*st(2,i)
      ENDDO

      END
