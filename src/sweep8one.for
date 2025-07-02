      MODULE  SWEEP8ONE

      IMPLICIT NONE

      CONTAINS
      
      SUBROUTINE ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                      delt,stot,
     &                      sgnc,sgni,sgne,sgnt,
     &                      ccof,icof,ecof,tcof)


!     Computes collision "ccof", incoming "icof", escape "ecof"
!     and transmission "tcof" coefficients for all discrete directions
!     in the first octant of angular space, and for all 
!     energy groups, in 3D Cartesian geometry.

!     Input values are:

!     mu      direction cosines,
!     eta     direction cosines,
!     ksi     direction cosines,
!     ndir    number of angular directions per octant,
!     delt    size of the pixel
!     stot    total cross sections
!     ng      number of groups
!     nn      ng * ndir 

      IMPLICIT NONE 

      INTEGER, PARAMETER :: nc=4, nbd=9

      INTEGER       :: nn, ndir, ng
      REAL (KIND=8) :: mu(ndir),eta(ndir),ksi(ndir)
      REAL          :: delt(3),stot(ng)
      INTEGER       :: sgnc(nc,nc),sgni(nc,nbd),
     &                 sgne(nc,nbd),sgnt(nbd,nbd)
      REAL          :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL          :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)

      INTEGER       :: d,g,cnt
      REAL          :: ox,oy,oz
      INTEGER :: ncol= nc*nc, nesc= nc*nbd, ntrn= nbd*nbd
      REAL , DIMENSION(nc,nc) :: ccof_aux , ccof_bux
      REAL , DIMENSION(nc,nbd) :: icof_aux , icof_bux
      REAL , DIMENSION(nbd,nc) :: ecof_aux , ecof_bux
      REAL , DIMENSION(nbd,nbd) :: tcof_aux , tcof_bux

!   For each direction in the octant and for each energy group 


      cnt = 0 
      DO d=1,ndir
        ox=2*mu(d)/delt(1)
        oy=2*eta(d)/delt(2)
        oz=2*ksi(d)/delt(3)
        DO g=1,ng

!           Coefficient for LL Method of Characteristics (simplified).
            CALL COF3C1(ox,oy,oz,stot(g),
     &                  ccof_aux(1,1),icof_aux(1,1),
     &                  ecof_aux(1,1),tcof_aux(1,1))
            
!      multiply by the sign-matrices to specilize the coefficient
!      to the octant 
            CALL SVIVR3(ncol,sgnc,ccof_aux(1,1),ccof_bux(1,1))
            CALL SVIVR3(nesc,sgni,icof_aux(1,1),icof_bux(1,1))
            CALL SVIVR3(nesc,sgne,ecof_aux(1,1),ecof_bux(1,1))
            CALL SVIVR3(ntrn,sgnt,tcof_aux(1,1),tcof_bux(1,1))
             
            cnt = cnt + 1 
            ccof(cnt,:,:) = ccof_bux(:,:)
            icof(cnt,:,:) = icof_bux(:,:)
            ecof(cnt,:,:) = ecof_bux(:,:)
            tcof(cnt,:,:) = tcof_bux(:,:)
         ENDDO
      ENDDO 

      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE  EIGHT_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                      delt,stot,
     &                      sgnc,sgni,sgne,sgnt,
     &                      ccof,icof,ecof,tcof)

      IMPLICIT NONE 

        
      INTEGER, PARAMETER :: nr = 8, nc=4, nbd=9
      INTEGER       :: nn,ndir, ng
      REAL (KIND=8) :: mu(ndir),eta(ndir),ksi(ndir)
      REAL          :: delt(3),stot(ng,nr)
      INTEGER       :: sgnc(*),sgni(*),sgne(*),sgnt(*)
      REAL          :: ccof(nn,nc,nc,nr),icof(nn,nc,nbd,nr)
      REAL          :: ecof(nn,nbd,nc,nr),tcof(nn,nbd,nbd,nr)

      INTEGER :: r

      DO r=1,nr
        CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                      delt,stot(:,r),
     &                      sgnc,sgni,sgne,sgnt,
     &                      ccof(:,:,:,r),icof(:,:,:,r),
     &                      ecof(:,:,:,r),tcof(:,:,:,r))
      ENDDO

      END SUBROUTINE




      SUBROUTINE ASSIGN_EIGHT_COEF3D(nn,nx,ny,
     &                    order_cell, nmat,zreg,
     &                    sgnc,sgni,sgne,sgnt,
     &                    ccofglb,icofglb,ecofglb,tcofglb,
     &                    ccof,icof,ecof,tcof)  
      
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: nr = 8, nc=4, nbd=9
      INTEGER       :: nn, nmat,nx,ny
      INTEGER       :: sgnc(*),sgni(*),sgne(*),sgnt(*)
      INTEGER       :: order_cell(6,nr)
      REAL          :: ccof(nn,nc,nc,nr),icof(nn,nc,nbd,nr)
      REAL          :: ecof(nn,nbd,nc,nr),tcof(nn,nbd,nbd,nr)
      REAL, INTENT(IN) :: ccofglb(nn ,nc ,nc ,nmat),
     &                    icofglb (nn,nc ,nbd,nmat),
     &                    ecofglb(nn ,nbd,nc ,nmat),
     &                    tcofglb(nn ,nbd,nbd,nmat)

      REAL , DIMENSION(nc,nc)   :: ccof_aux
      REAL , DIMENSION(nc,nbd)  :: icof_aux
      REAL , DIMENSION(nbd,nc)  :: ecof_aux
      REAL , DIMENSION(nbd,nbd) :: tcof_aux

      INTEGER, INTENT(IN) :: zreg(*)
      INTEGER :: ncol= nc*nc, nesc= nc*nbd, ntrn= nbd*nbd
      INTEGER :: i,j,k,mat,n,cnt

      DO cnt=1,nr
        i = order_cell(1,cnt)
        j = order_cell(3,cnt)
        k = order_cell(5,cnt)
        mat = zreg(i + nx*(j-1) + nx*ny*(k-1))
        
        ! Copy the coefficients from the global arrays            
        DO n=1,nn
        !      multiply by the sign-matrices to specilize the coefficient
        !      to the octant 
          CALL SVIVR3(ncol,sgnc,ccofglb(n,:,:,mat),ccof_aux)
          CALL SVIVR3(nesc,sgni,icofglb(n,:,:,mat),icof_aux)
          CALL SVIVR3(nesc,sgne,ecofglb(n,:,:,mat),ecof_aux)
          CALL SVIVR3(ntrn,sgnt,tcofglb(n,:,:,mat),tcof_aux)

          ccof(n,:,:,cnt) = ccof_aux
          icof(n,:,:,cnt) = icof_aux
          ecof(n,:,:,cnt) = ecof_aux
          tcof(n,:,:,cnt) = tcof_aux

        END DO

      END DO

      END SUBROUTINE

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      SUBROUTINE SWEEP_ONEREGION(n,mord,asrc,finc,aflx,fout,
     &                           sgnc,sgni,sgne,sgnt, 
     &                           ccof,icof,ecof,tcof)
!     Spatial sweep of 1 region 
!     "mord"  is method order flag.
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices.
!     "ccof", "icof",
!     "ecof", "tcof" arrays contain respectively collision,
!             surface-to-volume, escape and transmission coefficients.
!     "finc", "fout" incoming/outgoing flux auxiliary scratch arrays.
!     "aflx", "asrc" angular flux and angular source 

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9
      INTEGER      :: n,mord
      REAL         :: aflx(n,nc),asrc(n,nc)
      REAL         :: finc(n,nbd),fout(n,nbd)
      REAL,INTENT(IN) :: ccof(n,nc,nc),icof(n,nc,nbd)
      REAL,INTENT(IN) :: ecof(n,nbd,nc),tcof(n,nbd,nbd)
      REAL         :: ccofaux(n,nc,nc),icofaux(n,nc,nbd)
      REAL         :: ecofaux(n,nbd,nc),tcofaux(n,nbd,nbd)
      INTEGER      :: sgnc(nc,nc),sgni(nc,nbd),
     &                sgne(nc,nbd),sgnt(nbd,nbd)

      INTEGER      :: c,b,g
      INTEGER      :: ncol= nc*nc, nesc= nc*nbd, ntrn= nbd*nbd


      DO g=1,n
        CALL SVIVR3(ncol,sgnc,ccof(g,:,:),ccofaux(g,:,:))
        CALL SVIVR3(nesc,sgni,icof(g,:,:),icofaux(g,:,:))
        CALL SVIVR3(nesc,sgne,ecof(g,:,:),ecofaux(g,:,:))
        CALL SVIVR3(ntrn,sgnt,tcof(g,:,:),tcofaux(g,:,:))
      END DO

      SELECT CASE (mord)

      CASE (1) ! Constant flux scheme.
!              Sweep along x.
       aflx(:,1) =ccofaux(:,1,1)*asrc(:,1)
     &           +icofaux(:,1,1)*finc(:,1)
     &           +icofaux(:,1,2)*finc(:,2)
     &           +icofaux(:,1,3)*finc(:,3)
        fout(:,1)=ecofaux(:,1,1)*asrc(:,1)
     &           +tcofaux(:,1,1)*finc(:,1)
     &           +tcofaux(:,1,2)*finc(:,2)
     &           +tcofaux(:,1,3)*finc(:,3)
        fout(:,2)=ecofaux(:,2,1)*asrc(:,1)
     &           +tcofaux(:,2,1)*finc(:,1)
     &           +tcofaux(:,2,2)*finc(:,2)
     &           +tcofaux(:,2,3)*finc(:,3)
        fout(:,3)=ecofaux(:,3,1)*asrc(:,1)
     &           +tcofaux(:,3,1)*finc(:,1)
     &           +tcofaux(:,3,2)*finc(:,2)
     &           +tcofaux(:,3,3)*finc(:,3)
!    
      CASE (2) ! Linear flux scheme.

      DO c=1,4
        aflx(:,c)=ccofaux(:,c,1)*asrc(:,1)
     &             +ccofaux(:,c,2)*asrc(:,2)
     &             +ccofaux(:,c,3)*asrc(:,3)
     &             +ccofaux(:,c,4)*asrc(:,4)
     &             +icofaux(:,c,1)*finc(:,1)
     &             +icofaux(:,c,2)*finc(:,2)
     &             +icofaux(:,c,3)*finc(:,3)
     &             +icofaux(:,c,4)*finc(:,4)
     &             +icofaux(:,c,5)*finc(:,5)
     &             +icofaux(:,c,6)*finc(:,6)
     &             +icofaux(:,c,7)*finc(:,7)
     &             +icofaux(:,c,8)*finc(:,8)
     &             +icofaux(:,c,9)*finc(:,9)
      ENDDO
      DO b=1,9
        fout(:,b)=ecofaux(:,b,1)*asrc(:,1)
     &           +ecofaux(:,b,2)*asrc(:,2)
     &           +ecofaux(:,b,3)*asrc(:,3)
     &           +ecofaux(:,b,4)*asrc(:,4)
     &           +tcofaux(:,b,1)*finc(:,1)
     &           +tcofaux(:,b,2)*finc(:,2)
     &           +tcofaux(:,b,3)*finc(:,3)
     &           +tcofaux(:,b,4)*finc(:,4)
     &           +tcofaux(:,b,5)*finc(:,5)
     &           +tcofaux(:,b,6)*finc(:,6)
     &           +tcofaux(:,b,7)*finc(:,7)
     &           +tcofaux(:,b,8)*finc(:,8)
     &           +tcofaux(:,b,9)*finc(:,9)
      ENDDO
   

      END SELECT

      END
      
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      
      SUBROUTINE SWEEP_8REGIONS(n,mord,asrc,finc,aflx,fout,
     &                          ccof,icof,ecof,tcof,
     &                          xinc, yinc, zinc)

!     Spatial sweep of 1 region 
!
!     "mord"  is method order flag.
!     "nb"    is nubmer of spatial moments for mesh boundary flux.
!     "nc"    is number of spatial moments for mesh interior flux.
!     "nbd"   is a dimension of the coefficient matrices.
!     "ccof", "icof",
!     "ecof", "tcof" arrays contain respectively collision,
!             surface-to-volume, escape and transmission coefficients.
!     "finc", "fout" incoming/outgoing flux auxiliary scratch arrays.
!     "aflx", "asrc" angular flux and angular source 
         
      IMPLICIT NONE

      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9,
     &                      nr = 8, ns = 4, n2 = 2
      INTEGER,PARAMETER :: xx=1, yy=2, zz=3
      INTEGER,PARAMETER :: nx=2, ny=2
      INTEGER      :: n,mord
      REAL         :: aflx(n,nc,nr),asrc(n,nc,nr)
      REAL         :: finc(n,nb,zz,ns),fout(n,nb,zz,ns),faux(n,nbd)
      REAL         :: ccof(n,nc,nc,nr),icof(n,nc,nbd,nr)
      REAL         :: ecof(n,nbd,nc,nr),tcof(n,nbd,nbd,nr)

      REAL         :: xaux(n,nb),yaux(n,n2,nb),zaux(n,n2,n2,nb)

      INTEGER      :: b,c,r,s,x,y,z
      INTEGER, INTENT(IN) :: xinc, yinc, zinc

      INTEGER :: xout, yout, zout

      xout = 3-xinc
      yout = 3-yinc
      zout = 3-zinc

      SELECT CASE (mord)
!    
      CASE (2) ! Linear flux scheme.
!        Prepare sweep in first z-plane, incoming flux.
        ! fout = finc

         s=0
         DO y=1,2
         DO x=1,2
            s=s+1
            zaux(:,x,y,1)=finc(:,1,zz,s)
            zaux(:,x,y,2)=finc(:,2,zz,s)
            zaux(:,x,y,3)=finc(:,3,zz,s)
         ENDDO
         ENDDO
!        Sweep over z-planes.

        DO z=zinc,zout,zout-zinc
!   Prepare sweep over (x,y) in one z-plane, incomig flux.
            DO x=1,2
               s = (z-1)*nx + x
               yaux(:,x,1)=finc(:,1,yy,s)
               yaux(:,x,2)=finc(:,2,yy,s)
               yaux(:,x,3)=finc(:,3,yy,s)
            ENDDO

!           Sweep over (x,y) in one z-plane.
            DO y=yinc,yout,yout-yinc
!              Prepare sweep along x, incoming flux.
               s=(z-1)*ny+y
               xaux(:,1)=finc(:,1,xx,s)
               xaux(:,2)=finc(:,2,xx,s)
               xaux(:,3)=finc(:,3,xx,s)

!              Sweep along x.
               DO x=xinc,xout,xout-xinc
                  r=nx*(ny*(z-1)+y-1)+x
                  DO c=1,4
                     aflx(:,c,r)=ccof(:,c,1,r)*asrc(:,1,r)
     &                          +ccof(:,c,2,r)*asrc(:,2,r)
     &                          +ccof(:,c,3,r)*asrc(:,3,r)
     &                          +ccof(:,c,4,r)*asrc(:,4,r)
     &                          +icof(:,c,1,r)*xaux(:,1)
     &                          +icof(:,c,2,r)*xaux(:,2)
     &                          +icof(:,c,3,r)*xaux(:,3)
     &                          +icof(:,c,4,r)*yaux(:,x,1)
     &                          +icof(:,c,5,r)*yaux(:,x,2)
     &                          +icof(:,c,6,r)*yaux(:,x,3)
     &                          +icof(:,c,7,r)*zaux(:,x,y,1)
     &                          +icof(:,c,8,r)*zaux(:,x,y,2)
     &                          +icof(:,c,9,r)*zaux(:,x,y,3)
                  ENDDO
                  DO b=1,9
                    faux(:,b)= ecof(:,b,1,r)*asrc(:,1,r)
     &                        +ecof(:,b,2,r)*asrc(:,2,r)
     &                        +ecof(:,b,3,r)*asrc(:,3,r)
     &                        +ecof(:,b,4,r)*asrc(:,4,r)
     &                        +tcof(:,b,1,r)*xaux(:,1)
     &                        +tcof(:,b,2,r)*xaux(:,2)
     &                        +tcof(:,b,3,r)*xaux(:,3)
     &                        +tcof(:,b,4,r)*yaux(:,x,1)
     &                        +tcof(:,b,5,r)*yaux(:,x,2)
     &                        +tcof(:,b,6,r)*yaux(:,x,3)
     &                        +tcof(:,b,7,r)*zaux(:,x,y,1)
     &                        +tcof(:,b,8,r)*zaux(:,x,y,2)
     &                        +tcof(:,b,9,r)*zaux(:,x,y,3)
                  ENDDO

!                 Incomig angular flux for neighbouring cells.
                  xaux(:,1)    =faux(:,1)
                  xaux(:,2)    =faux(:,2)
                  xaux(:,3)    =faux(:,3)
                  yaux(:,x,1)  =faux(:,4)
                  yaux(:,x,2)  =faux(:,5)
                  yaux(:,x,3)  =faux(:,6)
                  zaux(:,x,y,1)=faux(:,7)
                  zaux(:,x,y,2)=faux(:,8)
                  zaux(:,x,y,3)=faux(:,9)

                ENDDO
!              End of sweep along x, outgoing flux.
               fout(:,1,xx,s)=xaux(:,1)
               fout(:,2,xx,s)=xaux(:,2)
               fout(:,3,xx,s)=xaux(:,3)
            ENDDO
!           End of sweep over (x,y) in one z-plane, outgoing flux.
            DO x=1,nx
               s=(z-1)*nx+x
               fout(:,1,yy,s)=yaux(:,x,1)
               fout(:,2,yy,s)=yaux(:,x,2)
               fout(:,3,yy,s)=yaux(:,x,3)
            ENDDO

         ENDDO
!        End of sweep over z, outgoing flux.

         DO y=1,ny
         DO x=1,nx
            s=(y-1)*nx+x
            fout(:,1,zz,s)=zaux(:,x,y,1)
            fout(:,2,zz,s)=zaux(:,x,y,2)
            fout(:,3,zz,s)=zaux(:,x,y,3)
         ENDDO
         ENDDO
   

      END SELECT

      END

      END MODULE SWEEP8ONE



      SUBROUTINE COMPUTE_COEFF_GLB(ng,ndir,nmat,delt3,nb_niv,
     &                      mu,eta,ksi,sigt,
     &                      ccofglb,icofglb,ecofglb,tcofglb)

      IMPLICIT NONE

      INTEGER, PARAMETER :: nc=4, nbd=9
      
      INTEGER, INTENT(IN):: ng,ndir,nmat,nb_niv

      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir)
      REAL, INTENT(IN) :: delt3(3)
      REAL, INTENT(IN) :: sigt(ng,nmat)

      REAL, INTENT(INOUT) :: ccofglb(ng*ndir,nc,nc  ,nmat,nb_niv),
     &                       icofglb(ng*ndir,nc,nbd ,nmat,nb_niv),
     &                       ecofglb(ng*ndir,nbd,nc ,nmat,nb_niv),
     &                       tcofglb(ng*ndir,nbd,nbd,nmat,nb_niv)

      REAL , DIMENSION(nc,nc) :: ccof_aux
      REAL , DIMENSION(nc,nbd) :: icof_aux
      REAL , DIMENSION(nbd,nc) :: ecof_aux
      REAL , DIMENSION(nbd,nbd) :: tcof_aux
      
      INTEGER :: mat,niv,d,g,cnt,c1,c2


      REAL :: ox,oy,oz, deltniv(3)

      
      DO niv=1,nb_niv
        deltniv = delt3/(2**(niv-1))
        DO mat=1,nmat 
            
          cnt = 1
          DO d=1,ndir
            ox=2*mu(d)/deltniv(1)
            oy=2*eta(d)/deltniv(2)
            oz=2*ksi(d)/deltniv(3)

            DO g=1,ng

!           Coefficient for LL Method of Characteristics (simplified).
              CALL COF3C1(ox,oy,oz,sigt(g,mat),
     &                ccofglb(cnt,:,:,mat,niv),icofglb(cnt,:,:,mat,niv),
     &                ecofglb(cnt,:,:,mat,niv),tcofglb(cnt,:,:,mat,niv))

              cnt = cnt+1
        
           ENDDO
         ENDDO

        ENDDO ! mat

      ENDDO ! NIV


      END SUBROUTINE COMPUTE_COEFF_GLB