      SUBROUTINE SETREA(n,a,x,inc)

      REAL       a,x(*)

      nn=inc*(n-1)+1
      DO i=1,nn,inc
         x(i)=a
      ENDDO

      END
      
      
      
      SUBROUTINE SETDBL(n,a,x,inc)

      REAL(KIND=8)       a,x(*)

      nn=inc*(n-1)+1
      DO i=1,nn,inc
         x(i)=a
      ENDDO

      
      END
      
      
      SUBROUTINE SETINT(n,a,x,inc)

      INTEGER    a,x(*)

      nn=inc*(n-1)+1
      DO i=1,nn,inc
         x(i)=a
      ENDDO

      END



      FUNCTION ISAMIN(n,x,inc)

      REAL  x(*)

      nn=inc*(n-1)+1
      ISAMIN=1
      a=ABS(x(1))
      DO i=1,nn,inc
         aa=ABS(x(i))
         IF(a.GT.aa)THEN
            ISAMIN=i
            a=aa
         ENDIF
      ENDDO

      END



      FUNCTION ICHIEQ(n,x,ix,a)

      INTEGER  x(*),a

      IF(ix.GT.0)THEN
         ii=1
      ELSE
        ii=-ix*n
      ENDIF
      DO ICHIEQ=1,n
         IF(x(ii).EQ.a) RETURN
         ii=ii+ix
      ENDDO
      ICHIEQ=n+1

      END



      FUNCTION ICHINE(n,x,ix,a)

      INTEGER  x(*),a

      IF(ix.GT.0)THEN
         ii=1
      ELSE
        ii=-ix*n
      ENDIF
      DO ICHINE=1,n
         IF(x(ii).NE.a) RETURN
         ii=ii+ix
      ENDDO
      ICHINE=n+1

      END



      FUNCTION ICHDEQ(n,x,ix,a)

      REAL(KIND=8) :: x(*),a

      IF(ix.GT.0)THEN
         ii=1
      ELSE
        ii=-ix*n
      ENDIF
      DO ICHDEQ=1,n
         IF(x(ii).EQ.a) RETURN
         ii=ii+ix
      ENDDO
      ICHDEQ=n+1

      END



      FUNCTION ICHDGE(n,x,ix,a)

         REAL(KIND=8) :: x(*),a
   
         IF(ix.GT.0)THEN
            ii=1
         ELSE
           ii=-ix*n
         ENDIF
         DO ICHDGE=1,n
            IF(x(ii).GE.a) RETURN
            ii=ii+ix
         ENDDO
         ICHDGE=n+1
   
      END



      FUNCTION ICHRGT(n,x,ix,a)

      REAL     x(*),a

      IF(ix.GT.0)THEN
         ii=1
      ELSE
        ii=-ix*n
      ENDIF
      DO ICHRGT=1,n
         IF(x(ii).GT.a) RETURN
         ii=ii+ix
      ENDDO
      ICHRGT=n+1

      END



      FUNCTION ICHCEQ(n,x,ix,a,l1,l2)

      CHARACTER*(*) x(*),a

      IF(ix.GT.0)THEN
         ii=1
      ELSE
         ii=-ix*n
      ENDIF
      DO ICHCEQ=1,n
         IF(x(ii)(l1:l2).EQ.a(l1:l2)) RETURN
         ii=ii+ix
      ENDDO
      ICHCEQ=n+1

      END



      FUNCTION SSUM(n,x,ix)

      DIMENSION        x(*)

      DOUBLE PRECISION aux
      nn=ix*(n-1)+1
      aux=0.
      DO i=1,nn,ix
         aux=aux+x(i)
      ENDDO
      SSUM=aux

      END



      DOUBLE PRECISION FUNCTION DSUM(n,x,ix)

      DOUBLE PRECISION x(*),aux
      nn=ix*(n-1)+1
      aux=0.
      DO i=1,nn,ix
         aux=aux+x(i)
      ENDDO
      DSUM=aux

      END



      SUBROUTINE SVCAL(n,a,x,ix,y,iy)

      DIMENSION  x(*),y(*)

      nn=iy*(n-1)+1
      i=1
      DO j=1,nn,iy
         y(j)=a*x(i)
         i=i+ix
      ENDDO

      END



      SUBROUTINE SPXMY(n,x,y)

      DIMENSION  x(*),y(*)

      DO i=1,n
         y(i)=x(i)-y(i)
      ENDDO

      END



      SUBROUTINE DPXMY(n,x,y)

      DOUBLE PRECISION  x(*),y(*)

      DO i=1,n
         y(i)=x(i)-y(i)
      ENDDO

      END



      SUBROUTINE SPYMX(n,x,y)

      DIMENSION  x(*),y(*)

      DO i=1,n
         y(i)=y(i)-x(i)
      ENDDO

      END
      
      

      SUBROUTINE SPXPY(n,x,y)

      DIMENSION  x(*),y(*)
      
      y(:n)=y(:n)+x(:n)
      
      END



      SUBROUTINE GATHRR(n,x,ind,y)

      REAL       x(*),y(n)
      INTEGER    ind(n)

      DO i=1,n
         y(i)=x(ind(i))
      ENDDO

      END



      SUBROUTINE GATHRC(n,a,x,ind,y)

      REAL       a,x(*),y(n)
      INTEGER    ind(n)

      DO i=1,n
         y(i)=a*x(ind(i))
      ENDDO

      END



      SUBROUTINE SVECT2(n,a1,x1,a2,x2)

      DIMENSION  x1(*),x2(*)

      DO i=1,n
         x1(i)=a1*x1(i)+a2*x2(i)
      ENDDO

      END



      SUBROUTINE SVECT3(n,a1,x1,a2,x2,a3,x3)

      DIMENSION  x1(*),x2(*),x3(*)

      DO i=1,n
         x1(i)=a1*x1(i)+a2*x2(i)+a3*x3(i)
      ENDDO

      END



      SUBROUTINE SVIVR3(n,a,x,y)

      INTEGER    a(*)
      REAL       x(*),y(*)

      DO i=1,n
         y(i)=a(i)*x(i)
      ENDDO

      END



      SUBROUTINE SVIVR4(n,x,y)

      REAL       x(n),y(n)

      y(:n) = x(:n)*y(:n) 

      END

      

      SUBROUTINE SVSVVV(n,x,y,z)

      REAL       x(*),y(*),z(*)

      DO i=1,n
         z(i)=z(i)+x(i)*y(i)
      ENDDO

      END



      SUBROUTINE SVSVVV2(n,x,y,z)

         REAL       x(*)
         REAL(KIND=8) y(*),z(*)
   
         DO i=1,n
            z(i)=z(i)+x(i)*y(i)
         ENDDO
   
      END



      SUBROUTINE SVCAL2(n,a,x)

      REAL       x(*),a

      x(:n)=a*x(:n)

      END



      SUBROUTINE SVCAL3(n,a,x,y)

      REAL       x(*),y(*),a

      y(:n)=y(:n)+a*x(:n)

      END



      SUBROUTINE SVCAL4(n,a,x,y)

      DIMENSION  x(*),y(*)

      y(:n)=a*x(:n)

      END



      SUBROUTINE SPRAD0(n,x,ix,y,iy,z,iz)

      REAL       x(*),y(*),z(*)

      nn=iz*(n-1)+1
      i=1
      j=1
      DO k=1,nn,iz
         z(k)=x(i)*y(j)
         i=i+ix
         j=j+iy
      ENDDO

      END



      SUBROUTINE SPRAD1(n,x,ix,y,iy,z,iz)

      REAL       x(*),y(*),z(*)

      nn=iz*(n-1)+1
      i=1
      j=1
      DO k=1,nn,iz
         z(k)=z(k)+x(i)*y(j)
         i=i+ix
         j=j+iy
      ENDDO

      END
      
      
      
      SUBROUTINE AMUXVE(m,n,a,x,y)
      
      REAL       a(m,*),x(*),y(*)
      
      DO i=1,n
         y(:m) = y(:m) + x(i)*a(:m,i)
      ENDDO 

      END

      
      
      SUBROUTINE AMUXVX(m,n,x,a)
      
      REAL       a(m,*),x(*)
      
      DO i=1,n
         a(:m,i) = x(i)*a(:m,i)
      ENDDO 

      END

      
      
      SUBROUTINE AMUXVZ(m,n,x,a,b)
      
      REAL       b(m,*),a(m,*),x(*)
      
      DO i=1,n
         b(:m,i) = b(:m,i) + x(i)*a(:m,i)
      ENDDO 

      END


      
      DOUBLE PRECISION FUNCTION SDOTDBL(n,x,ix,y,iy)

         REAL(KIND=8) ::  x(*), y(*)
         REAL(KIND=8) ::  aux
         INTEGER :: n,ix,iy

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
         SDOTDBL=aux

      END
