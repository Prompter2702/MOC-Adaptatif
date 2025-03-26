      SUBROUTINE COFS32(sx,sy,sz,nc,nbd,sc,si,se,st)


!     Computes signs of collision "sc", incoming "si", escape "se" and
!     transmission "st" coefficients for full linear schemes
!     in one octant defined by its direction cosines
!     "sx", "sy" and "sz".


      INTEGER :: sx,sy,sz,nc,nbd
      INTEGER :: sc(nc,nc),si(nc,nbd),se(nbd,nc),st(nbd,nbd)

      INTEGER , PARAMETER , DIMENSION(4) :: vx = (/ 2, 4, 6, 8/)
      INTEGER , PARAMETER , DIMENSION(4) :: vy = (/ 3, 4, 7, 8/)
      INTEGER , PARAMETER , DIMENSION(4) :: vz = (/ 5, 6, 7, 8/)
      INTEGER , PARAMETER , DIMENSION(4) :: bx = (/ 6, 8,10,12/)
      INTEGER , PARAMETER , DIMENSION(4) :: by = (/ 2, 4,11,12/)
      INTEGER , PARAMETER , DIMENSION(4) :: bz = (/ 3, 4, 7, 8/)

      INTEGER :: i,j


      CALL SETINT( nc*nc ,1,sc,1)
      CALL SETINT( nc*nbd,1,si,1)
      CALL SETINT( nc*nbd,1,se,1)
      CALL SETINT(nbd*nbd,1,st,1)


      DO i=1,8
      DO j=1,4
         sc(i,vx(j))=sx*sc(i,vx(j))
         sc(i,vy(j))=sy*sc(i,vy(j))
         sc(i,vz(j))=sz*sc(i,vz(j))
         sc(vy(j),i)=sy*sc(vy(j),i)
         sc(vx(j),i)=sx*sc(vx(j),i)
         sc(vz(j),i)=sz*sc(vz(j),i)
      ENDDO
      ENDDO

      DO i=1,12
      DO j=1,4
         si(vx(j),i)=sx*si(vx(j),i)
         si(vy(j),i)=sy*si(vy(j),i)
         si(vz(j),i)=sz*si(vz(j),i)
      ENDDO
      ENDDO
      DO j=1,4
      DO i=1,8
         si(i,bx(j))=sx*si(i,bx(j))
         si(i,by(j))=sy*si(i,by(j))
         si(i,bz(j))=sz*si(i,bz(j))
      ENDDO
      ENDDO

      DO i=1,8
      DO j=1,12
         se(j,i)=si(i,j)
      ENDDO
      ENDDO

      DO i=1,12
      DO j=1,4
         st(i,bx(j))=sx*st(i,bx(j))
         st(i,by(j))=sy*st(i,by(j))
         st(i,bz(j))=sz*st(i,bz(j))
         st(bx(j),i)=sx*st(bx(j),i)
         st(by(j),i)=sy*st(by(j),i)
         st(bz(j),i)=sz*st(bz(j),i)
      ENDDO
      ENDDO

      END
