      SUBROUTINE COFS21(sx,sy,nc,nbd,sc,si,se,st)


!     Computes signs of collision "sc", incoming "si", escape "se" and
!     transmission "st" coefficients for simplified linear schemes
!     in one octant defined by its direction cosines
!     "sx" and "sy".


      INTEGER :: sx,sy,nc,nbd
      INTEGER :: sc(nc,nc),si(nc,nbd),se(nbd,nc),st(nbd,nbd)

      INTEGER :: i,j


      CALL SETINT( nc*nc ,1,sc,1)
      CALL SETINT( nc*nbd,1,si,1)
      CALL SETINT( nc*nbd,1,se,1)
      CALL SETINT(nbd*nbd,1,st,1)

      DO i=1,3
         sc(i,2)=sx*sc(i,2)
         sc(2,i)=sx*sc(2,i)
         sc(i,3)=sy*sc(i,3)
         sc(3,i)=sy*sc(3,i)
      ENDDO

      DO i=1,4
         si(2,i)=sx*si(2,i)
         si(3,i)=sy*si(3,i)
      ENDDO
      DO i=1,3
         si(i,4)=sx*si(i,4)
         si(i,2)=sy*si(i,2)
      ENDDO

      DO i=1,3
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
