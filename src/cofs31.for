      SUBROUTINE COFS31(sx,sy,sz,nc,nbd,sc,si,se,st)


!     Computes signs of collision "sc", incoming "si", escape "se" and
!     transmission "st" coefficients for simplified linear schemes
!     in one octant defined by its direction cosines
!     "sx", "sy" and "sz".


      INTEGER :: sx,sy,sz,nc,nbd
      INTEGER :: sc(nc,nc),si(nc,nbd),se(nbd,nc),st(nbd,nbd)

      INTEGER :: i,j

      CALL SETINT(nc*nc,1,sc,1)
      CALL SETINT(nc*nbd,1,si,1)
      CALL SETINT(nbd*nc,1,se,1)
      CALL SETINT(nbd*nbd,1,st,1)

      DO i=1,4
         sc(i,2)=sx*sc(i,2)
         sc(2,i)=sx*sc(2,i)
         sc(i,3)=sy*sc(i,3)
         sc(3,i)=sy*sc(3,i)
         sc(i,4)=sz*sc(i,4)
         sc(4,i)=sz*sc(4,i)
      ENDDO

      DO i=1,9
         si(2,i)=sx*si(2,i)
         si(3,i)=sy*si(3,i)
         si(4,i)=sz*si(4,i)
      ENDDO
      DO i=1,4
         si(i,5)=sx*si(i,5)
         si(i,8)=sx*si(i,8)
         si(i,2)=sy*si(i,2)
         si(i,9)=sy*si(i,9)
         si(i,3)=sz*si(i,3)
         si(i,6)=sz*si(i,6)
      ENDDO

      DO i=1,4
      DO j=1,9
         se(j,i)=si(i,j)
      ENDDO
      ENDDO

      DO i=1,9
         st(i,5)=sx*st(i,5)
         st(5,i)=sx*st(5,i)
         st(i,8)=sx*st(i,8)
         st(8,i)=sx*st(8,i)
         st(i,2)=sy*st(i,2)
         st(2,i)=sy*st(2,i)
         st(i,9)=sy*st(i,9)
         st(9,i)=sy*st(9,i)
         st(i,3)=sz*st(i,3)
         st(3,i)=sz*st(3,i)
         st(i,6)=sz*st(i,6)
         st(6,i)=sz*st(6,i)
      ENDDO

      END
