      SUBROUTINE COFSGN(sgnc,sgni,sgne,sgnt,nc,nbd,ndim,mord)

!     Fills the arrays "sgnc", "sgni", "sgne", "sgnt" containing the
!     signs per octant of respectively collision, incoming, escape
!     and transmission coefficients for higher order differencing
!     schemes. "nc" and "nb" are input as number of flux space moments,
!     respectively, in mesh interior and on mesh boundaries.
!     "mord" is input as method order indicator:
!     for 1D, it is equal to the order of polynomial expansion,
!     for 2D & 3D, 0 = diamond, 1 = constant,  2 = pure linear,
!     3 = bilinear,
    
      INTEGER :: sgnc(nc,nc,*),sgni(nc,nbd,*)
      INTEGER :: sgne(nbd,nc,*),sgnt(nbd,nbd,*)
      INTEGER :: nc,nbd,ndim,mord
    
!     Signs of direction cosines per octant.
    
      INTEGER , PARAMETER,DIMENSION(8) :: sx=(/ 1,-1,-1, 1, 1,-1,-1, 1/)
      INTEGER , PARAMETER,DIMENSION(8) :: sy=(/ 1, 1,-1,-1, 1, 1,-1,-1/)
      INTEGER , PARAMETER,DIMENSION(8) :: sz=(/ 1, 1, 1, 1,-1,-1,-1,-1/)
    
      INTEGER :: noct,oct,s,sj,i,j
    
    
!     No sign change for diamond and constant schemes.
    
      IF(nc.EQ.1)THEN  
         noct=2**ndim
         CALL SETINT(noct,1,sgnc,1)
         CALL SETINT(nbd*noct,1,sgni,1)
         CALL SETINT(nbd*noct,1,sgne,1)
         CALL SETINT(nbd*nbd*noct,1,sgnt,1)
         RETURN  
      ENDIF
    
!     Other cases.
    
      SELECT CASE (ndim)
    
      CASE (1) ! 1D
    
      s=1
        DO i=1,nc
            sgni(i,1,1)=1
            sgni(i,1,2)=s
            sgne(1,i,1)=1
            sgne(1,i,2)=s
            sj=1
            DO j=1,nc
               sgnc(i,j,1)=1
               sgnc(i,j,2)=s*sj
               sj=-sj
            ENDDO
                s=-s
            ENDDO
        sgnt(:,:,1:2)=1
    
      CASE (2) ! 2D
    
             SELECT CASE (mord)
             CASE (2) ! Linear
                DO oct=1,4
                   CALL COFS21(sx(oct),sy(oct),nc,nbd,
     &                     sgnc(1,1,oct),sgni(1,1,oct),
     &                     sgne(1,1,oct),sgnt(1,1,oct))
                ENDDO
             CASE (3) ! Bilinear
                DO oct=1,4
                   CALL COFS22(sx(oct),sy(oct),nc,nbd,
     &                     sgnc(1,1,oct),sgni(1,1,oct),
     &                     sgne(1,1,oct),sgnt(1,1,oct))
                ENDDO
             END SELECT
    
          CASE (3) ! 3D
    
             SELECT CASE (mord)
             CASE (2) ! Linear
                DO oct=1,8
                   CALL COFS31(sx(oct),sy(oct),sz(oct),nc,nbd,
     &                     sgnc(1,1,oct),sgni(1,1,oct),
     &                     sgne(1,1,oct),sgnt(1,1,oct))
                ENDDO
             CASE (3) ! Bilinear
                DO oct=1,8
                   CALL COFS32(sx(oct),sy(oct),sz(oct),nc,nbd,
     &                     sgnc(1,1,oct),sgni(1,1,oct),
     &                     sgne(1,1,oct),sgnt(1,1,oct))
                ENDDO
             END SELECT
    
          END SELECT
    
          END
    