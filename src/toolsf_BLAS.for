

      SUBROUTINE SSCAL(n,a,x,ix)

      DIMENSION  x(*)

      nn=ix*(n-1)+1
      DO i=1,nn,ix
         x(i)=a*x(i)
      ENDDO

      END
      
      
      
      
      SUBROUTINE SAXPY(n,a,x,ix,y,iy)

      DIMENSION  x(*),y(*)

      nn=iy*(n-1)+1
      i=1
      DO j=1,nn,iy
         y(j)=y(j)+a*x(i)
         i=i+ix
      ENDDO

      END



      SUBROUTINE SCOPY(n,x,ix,y,iy)

      DIMENSION  x(*),y(*)

      nn=iy*(n-1)+1
      i=1
      DO j=1,nn,iy
         y(j)=x(i)
         i=i+ix
      ENDDO

      END


      SUBROUTINE DCOPY(n,x,ix,y,iy)

      REAL(KIND=8) :: x(*),y(*)

      nn=iy*(n-1)+1
      i=1
      DO j=1,nn,iy
         y(j)=x(i)
         i=i+ix
      ENDDO

      END



      FUNCTION SDOT(n,x,ix,y,iy)

      DIMENSION         x(*),y(*)
      DOUBLE PRECISION  aux

      aux=0.
      IF(ix.EQ.1)THEN
         IF(iy.EQ.1)THEN
            DO i=1,n
               aux=aux+x(i)*y(i)
            ENDDO
         ELSE
            j=1
            DO i=1,n
               aux=aux+x(i)*y(j)
               j=j+iy
            ENDDO
         ENDIF
      ELSE
         nn=iy*(n-1)+1
         i=1
         DO j=1,nn,iy
            aux=aux+x(i)*y(j)
            i=i+ix
         ENDDO
      ENDIF
      SDOT=aux

      END



      FUNCTION ISAMAX(n,x,inc)

      REAL  x(*)

      nn=inc*(n-1)+1
      ISAMAX=1
      a=ABS(x(1))
      DO i=1,nn,inc
         aa=ABS(x(i))
         IF(a.LT.aa)THEN
            ISAMAX=i
            a=aa
         ENDIF
      ENDDO

      END



      FUNCTION IDAMAX(n,x,inc)

      DOUBLE PRECISION a,x(*)

      nn=inc*(n-1)+1
      IDAMAX=1
      a=ABS(x(1))
      DO i=1,nn,inc
         aa=ABS(x(i))
         IF(a.LT.aa)THEN
            IDAMAX=i
            a=aa
         ENDIF
      ENDDO

      END




      SUBROUTINE DSCAL(n,a,x,ix)

      DOUBLE PRECISION  x(*),a

      nn=ix*(n-1)+1
      DO i=1,nn,ix
         x(i)=a*x(i)
      ENDDO

      END


