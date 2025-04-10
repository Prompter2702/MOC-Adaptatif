       SUBROUTINE SWEE3D_ADAPTIVE(
       ! inputs: dimensions 
     &                          nn,ng,nr,nh,nc,nmat,
     &                          nb,nbd,nbfx,nbfy,nbfz,
     &                          nx,ny,nz,
     &                          nani,nhrm,nd,ndir,
       ! inputs: cross section  
     &                          sigt,sigs,
       ! inputs: source and flux moments at previous iteration 
     &                          flxp,srcm,
       ! input/output: boundary flux 
     &                          bflx,bfly,bflz,
       ! input: angular quadrature & spherical harmonics
     &                          mu,eta,ksi,w,pisn,sphr,
       ! input: sing matrices to adapt coefficients to octants 
     &                          sgnc,sgni,sgne,sgnt,
       ! input: auxiliar memory for boundary angular fluxes, angular 
       ! source moments and self-scattering xs 
     &                          finc,fout,tmom,sigg,
       ! input: auxiliar memory angular flux & source 
     &                          aflx,asrc,dsrc,
       ! input: auxiliar memory for the interface angular flux 
     &                          xaux,yaux,zaux,
       ! input: auxiliar memory to drive angular mirror-reflection 
       ! or rotation/translation b.c.
     &                          rdir,
       ! input: pixel-to-medium array 
     &                          zreg, 
       ! input: pixel-to-medium array 
     &                          dira,dirf,
       ! input: logical to add angular source in case of
       ! time-dependent calculation 
     &                          lgki,
       ! output: flux moments 
     &                          flxm,
       ! Delta on each direction of the region
     &                          delt3,
      ! max inner iterations and tolerance for inner iterations
     &                          maxinner, tolinner)

      USE FLGOCT
      USE FLGBCD
      USE FLGINC
      USE SWEEP8ONE
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: noct = 8, n8=8, ns=4
      
      INTEGER :: nn,ng,nr,nh,nc,nx,ny,nz,nmat
      INTEGER :: nb,nbd,nbfx,nbfy,nbfz
      INTEGER :: nani,nhrm,nd,ndir
      REAL    :: flxm(ng,nr,nh,nc)
      REAL    :: bflx(nn,nb,nbfx,2,noct),bfly(nn,nb,nbfy,2,noct),
     &           bflz(nn,nb,nbfz,2,noct)
      REAL    :: sigt(ng,nmat),sigs(ng,0:nani,nmat)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      REAL    :: sphr(nhrm,nd)
      INTEGER :: sgnc(nc,nc,noct),sgni(nc,nbd,noct)
      INTEGER :: sgne(nbd,nc,noct),sgnt(nbd,nbd,noct)
      INTEGER :: zreg(nr)
      REAL    :: aflx(nn,nr,nc,noct),asrc(nn,nr,nc,noct)
      REAL    :: dsrc(nn,nr,nc,*)
      REAL    :: finc(nn,nbd),fout(nn,nbd)
      REAL    :: flxp(ng,nr,nh,nc),srcm(ng,nr,nh,nc),tmom(ng,nr,nh,nc)
      REAL    :: sigg(ng,nr)
      REAL    :: xaux(nn,nb),yaux(nn,nx,nb),zaux(nn,nx,ny,nb)
      INTEGER :: rdir(nd,*),dira(*),dirf(*)
      LOGICAL :: lgki
      REAL, INTENT(IN)     :: delt3(3)

      INTEGER :: xout,yout,zout,fst
      INTEGER :: oc,oct,d,dd,da,c,h,i

      REAL  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL  :: xshom0(ng), xshom1(ng,n8)
      REAL  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL  :: ok 
      INTEGER  :: ii,jj,kk,rr 
      INTEGER , INTENT(IN) :: maxinner
      REAL, INTENT(IN) :: tolinner
      REAL :: errinner
      INTEGER :: cnt,g,r


      errinner = tolinner + 1.0
      flxp = 0.0

      DO WHILE (errinner>tolinner .AND. cnt<maxinner)
    
      cnt = cnt + 1

      flxm = 0.0
      CALL GMOM3D(ng,nani,nh,nc,nr,sigs,flxp,srcm,zreg,sigg,tmom)

!     Loop over octants of angular space in order defined by
!     octant list "olst3D".


      DO oc=1,8
         oct=olst3D(oc)
!        Index of outgoing side.
         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
!        Sweep for all directions in octant "oct".
!        initial direction of the octant 

         da = (oct-1)*ndir+1

!        Directional source.

         CALL GIRSRC(nhrm,ng,ndir,nr,nh,nc,sphr(1,da),
     &               tmom(1,1,1,1),asrc(1,1,1,oct))
         IF(lgki)CALL SPXPY(nr*nc,dsrc(1,1,1,oct),asrc(1,1,1,oct))

!        Sweep.
        ! The boundary entries are computed on the exiting part
        ! in order to work only on the exiting part of the bdfl
        
        bflx(:,:,:,xout,oct) = bflx(:,:,:,xinc(oct),oct)
        bfly(:,:,:,yout,oct) = bfly(:,:,:,yinc(oct),oct)
        bflz(:,:,:,zout,oct) = bflz(:,:,:,zinc(oct),oct)


!     Niveau 0 doit être déjà fait :
        CALL XSSRCHOMO0(nn,ng,nr,nx,ny, ndir,
     &                 1,nx,1,ny,1,nz,
     &                 sigt, zreg, aflx(:,:,:,oct), w,
     &                 asrc(:,:,:,oct), xshom0,
     &                 asrcm0)


        CALL ONE_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                delt3,xshom0,
     &                sgnc(:,:,oct),sgni(:,:,oct),
     &                sgne(:,:,oct),sgnt(:,:,oct),
     &                ccof,icof,ecof,tcof)
      

        CALL MERGEBOUND0(nn,nb,
     &                   1,nx,1,ny,1,nz,
     &                   nx,ny,nz,
     &                   bflx(:,:,:,xinc,oct),
     &                   bfly(:,:,:,yinc,oct),
     &                   bflz(:,:,:,zinc,oct),
     &                   finc0)
        
        CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)
       

      ! RECURIVE SWEEP

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc(oct), yinc(oct), zinc(oct),oct,
     &                             1,nx,1,ny,1,nz,
     &                             0,sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc(:,:,oct),sgni(:,:,oct),
     &                             sgne(:,:,oct),sgnt(:,:,oct),
     &                             zreg, asrc(:,:,:,oct),
     &                             aflx(:,:,:,oct),
     &                             bflx(:,:,:,xout,oct),
     &                             bfly(:,:,:,yout,oct),
     &                             bflz(:,:,:,zout,oct),
     &                             aflx0, aflx1,xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,tolinner,
     &                             errcor,errmul,ok,delt3,ii,jj,kk,rr)  

! Compute the flux moments
         DO c=1,nc
         DO h=1,nh
         fst = 0
         DO d = 1,ndir
            flxm(:ng,:nr,h,c) = flxm(:ng,:nr,h,c) + 
     &        (sphr(h,da)*w(d))*aflx(fst+1:fst+ng,:nr,c,oct)
            fst = fst + ng
         ENDDO
         ENDDO
         ENDDO
         
      ENDDO
!     Boundary conditions: specular reflection, albedo
    
      DO oc=1,8
         oct=olst3D(oc)
!        Index of outgoing side.
         xout=3-xinc(oct)
         yout=3-yinc(oct)
         zout=3-zinc(oct)
         
         CALL BCD2R3(oct,xout,yout,zout,
     &               bflx,bfly,bflz,mu,eta,
     &               ksi,w,pisn,rdir)
     
      ENDDO 

! Compute the inner error
      errinner = 0.0
      DO c=1,nc
      DO h=1,nh
      fst = 0
      DO d = 1,ndir
      DO r = 1,nr
      DO g = 1,ng
        errinner = errinner + (flxm(g,r,h,c)-flxp(g,r,h,c))**2
        flxp(g,r,h,c) = flxm(g,r,h,c)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      print *, errinner
      errinner = sqrt(errinner/(nn*nr*nh*nc))
      
      print *,"Error inner", errinner
      print *, "FIN ITER"

! End inner iterations
      ENDDO

      print *,"Number iteration", cnt
      print *,"Error inner", errinner
      
      END

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      
      RECURSIVE SUBROUTINE RECADA_ONE_OCTANT(nn,ng,ndir,nr,nh,
     &                             nx,ny,nz,
     &                             xinc,yinc,zinc, oct,
     &                             imin,imax,jmin,jmax,kmin,kmax,niv,
     &                             sigt, 
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflx0, aflx1, xshom0, xshom1,
     &                             asrcm0, asrcm1,finc0,finc1,
     &                             fout0,fout1,ccof,icof,ecof,tcof,
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)
     
      USE SWEEP8ONE
      USE SRCCORR 

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: nc = 4, nb = 3, nbd = 9, n8 = 8,
     &                    ns = 4

      INTEGER, INTENT(IN) :: nn,ng,ndir,nr,nh,
     &           imin,imax,jmin,jmax,kmin,kmax,
     &           nx,ny,nz,niv,oct
      INTEGER, INTENT(IN) :: xinc,yinc,zinc
    !   REAL, INTENT(IN)    :: flxm(ng,nr,nh,nc)
      REAL, INTENT(IN)    :: sigt(ng,*)
      REAL(KIND=8), INTENT(IN) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir)
      INTEGER, INTENT(IN) :: sgnc(nc,nc),sgni(nc,nbd)
      INTEGER, INTENT(IN) :: sgne(nbd,nc),sgnt(nbd,nbd)
      INTEGER, INTENT(IN) :: zreg(nr)
      REAL, INTENT(IN) :: tol
      REAL, INTENT(INOUT) :: asrc(nn,nr,nc)
      REAL, INTENT(INOUT) :: aflx(nn,nr,nc)
      REAL, INTENT(INOUT) :: bflx(nn,nb,ny*nz),bfly(nn,nb,nx*nz),
     &           bflz(nn,nb,nx*ny)
     
      REAL, INTENT(INOUT)  :: aflx0(nn,nc), aflx1(nn,nc,n8)
      REAL, INTENT(INOUT)  :: xshom0(ng), xshom1(ng,n8)
      REAL, INTENT(INOUT)  :: asrcm0(ng,ndir,nc), asrcm1(ng,ndir,nc,n8)
      REAL, INTENT(INOUT)  :: finc1(nn,nb,3,ns), fout1(nn,nb,3,ns)
      REAL, INTENT(INOUT)  :: finc0(nn,nb,3), fout0(nn,nb,3)

      REAL, INTENT(INOUT)  :: ccof(nn,nc,nc),icof(nn,nc,nbd)
      REAL, INTENT(INOUT)  :: ecof(nn,nbd,nc),tcof(nn,nbd,nbd)
      REAL(KIND=8),INTENT(INOUT)  :: errcor(ng, ndir),errmul(ng, ndir)
      REAL, INTENT(IN)     :: ccof8(nn,nc,nc,n8),icof8(nn,nc,nbd,n8)
      REAL, INTENT(IN)     :: ecof8(nn,nbd,nc,n8),tcof8(nn,nbd,nbd,n8)
      LOGICAL, INTENT(INOUT)  :: ok 
      REAL, INTENT(IN)        :: delt3(3)
      INTEGER, INTENT(INOUT)  :: i,j,k,r 

      REAL :: aflxtps(nn,nc,n8),asrcmtps(ng,ndir,nc,n8),xshomtps(ng,n8)
      REAL :: fouttps(nn,nb,3,n8)
      REAL :: errtot
      REAL :: ccoftps(nn,nc,nc,n8),icoftps(nn,nc,nbd,n8)
      REAL :: ecoftps(nn,nbd,nc,n8),tcoftps(nn,nbd,nbd,n8)
      
      fout1 = 0.0

! Error check by source-correction estimation
      CALL SRCCOR(ng,ndir,delt3/(2**niv),mu,eta,ksi,asrcm0,
     &            xshom0,errcor)
      errtot = SUM(ABS(errcor))

      IF(imax-imin>0 .AND. jmax-jmin>0 
     & .AND. kmax-kmin>0 .AND.  errtot >tol) THEN
        ok = .FALSE.
      ELSE
        ok = .TRUE.
      ENDIF

      IF (.NOT. ok) THEN 
! If the criterion is not satisfied, compute on lvl 1
        CALL XSSRCHOMO1(nn,ng,nr,nx,ny, ndir,
     &                 imin, imax,jmin,jmax,kmin,kmax,
     &                 sigt, zreg, aflx, w, asrc, xshom1,
     &                 asrcm1)

        CALL MERGEBOUND1(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, finc1)

        CALL EIGHT_COEF3D(nn,ndir,ng,mu,eta,ksi,
     &                    delt3/(2**(niv+1)),xshom1,
     &                    sgnc,sgni,sgne,sgnt,
     &                    ccof8,icof8,ecof8,tcof8)

        CALL SWEEP_8REGIONS(nn,2,asrcm1, finc1, aflx1, fout1,
     &                     ccof8,icof8,ecof8,tcof8, xinc, yinc, zinc)

      ENDIF
      
      IF (ok) THEN
! If the criterion is satisfied, compute on lvl 0 and 
        IF (niv > 0) THEN
            CALL MERGEBOUND0(nn,nb,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx,
     &                   bfly,
     &                   bflz,
     &                   finc0)

            CALL SWEEP_ONEREGION(nn,2,asrcm0, finc0, aflx0, fout0,
     &                          ccof,icof,ecof,tcof)
        END IF 
        
        ! - write the angular flux of node-0 on the fine-mesh flux 
        CALL SPLITAFLX0(nn,nr,nc,nx,ny,
     &                   imin,imax,jmin,jmax,kmin,kmax,
     &                   aflx,aflx0)
!      - write the outgoing angular flux of node-0 
!            on the fine-mesh boundary flux

         CALL SPLITBOUND0(nn, ng,nb,
     &                   imin, imax,jmin,jmax,kmin,kmax,
     &                   nx,ny,nz,
     &                   bflx, bfly, bflz, fout0)

          RETURN

      ELSE IF((imax-imin)>1) THEN
!         15/ si ko : 
!        - call recursively the subroutine for the 8 nodes of level-1

        asrcmtps(:,:,:,1) = asrcm1(:,:,:, xinc+2*(yinc-1) + 4*(zinc-1))  
        aflxtps(:,:,1) =   aflx1(:,:, xinc+2*(yinc-1) + 4*(zinc-1))
        xshomtps(:,1) = xshom1(:,       xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps(:,:,:,2) = asrcm1(:,:,:,3-xinc+2*(yinc-1) + 4*(zinc-1))  
        aflxtps(:,:,2) = aflx1(:,:,  3-xinc+2*(yinc-1) + 4*(zinc-1))
        xshomtps(:,2) = xshom1(:,      3-xinc+2*(yinc-1) + 4*(zinc-1))

        asrcmtps(:,:,:,3)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(zinc-1))  
        aflxtps(:,:,3)= aflx1(:,:,  xinc+2*(2-yinc) + 4*(zinc-1))
        xshomtps(:,3)= xshom1(:,      xinc+2*(2-yinc) + 4*(zinc-1))

        asrcmtps(:,:,:,4)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))  
        aflxtps(:,:,4)= aflx1(:,:,  3-xinc + 2*(2-yinc) + 4*(zinc-1))
        xshomtps(:,4)= xshom1(:,      3-xinc + 2*(2-yinc) + 4*(zinc-1))

        asrcmtps(:,:,:,5)= asrcm1(:,:,:,xinc + 2*(yinc-1) + 4*(2-zinc))  
        aflxtps(:,:,5)= aflx1(:,:,  xinc + 2*(yinc-1) + 4*(2-zinc))
        xshomtps(:,5)= xshom1(:,      xinc + 2*(yinc-1) + 4*(2-zinc))

        asrcmtps(:,:,:,6)=asrcm1(:,:,:,3-xinc + 2*(yinc-1) + 4*(2-zinc))  
        aflxtps(:,:,6)= aflx1(:,:,  3-xinc + 2*(yinc-1) + 4*(2-zinc))
        xshomtps(:,6) = xshom1(:,     3-xinc + 2*(yinc-1) + 4*(2-zinc))  

        asrcmtps(:,:,:,7)= asrcm1(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))  
        aflxtps(:,:,7)= aflx1(:,:,  xinc+2*(2-yinc) + 4*(2-zinc))
        xshomtps(:,7)= xshom1(:,      xinc+2*(2-yinc) + 4*(2-zinc))

        asrcmtps(:,:,:,8)=asrcm1(:,:,:,3-xinc + 2*(2-yinc) + 4*(2-zinc))  
        aflxtps(:,:,8) = aflx1(:,:, 3-xinc + 2*(2-yinc) + 4*(2-zinc))
        xshomtps(:,8)= xshom1(:,      3-xinc + 2*(2-yinc) + 4*(2-zinc))

        
      ccoftps(:,:,:,1) = ccof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      ecoftps(:,:,:,1) = ecof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      tcoftps(:,:,:,1) = tcof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))
      icoftps(:,:,:,1) = icof8(:,:,:, xinc + 2*(yinc-1) + 4*(zinc-1))

      ccoftps(:,:,:,2) = ccof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      ecoftps(:,:,:,2) = ecof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      tcoftps(:,:,:,2) = tcof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))
      icoftps(:,:,:,2) = icof8(:,:,:,3-xinc+ 2*(yinc-1) + 4*(zinc-1))

      ccoftps(:,:,:,3) = ccof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      ecoftps(:,:,:,3) = ecof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      tcoftps(:,:,:,3) = tcof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))
      icoftps(:,:,:,3) = icof8(:,:,:,xinc + 2*(2-yinc) + 4*(zinc-1))

      ccoftps(:,:,:,4) = ccof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      ecoftps(:,:,:,4) = ecof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      tcoftps(:,:,:,4) = tcof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))
      icoftps(:,:,:,4) = icof8(:,:,:,3-xinc + 2*(2-yinc) + 4*(zinc-1))

      ccoftps(:,:,:,5) = ccof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      ecoftps(:,:,:,5) = ecof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      tcoftps(:,:,:,5) = tcof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))
      icoftps(:,:,:,5) = icof8(:,:,:,xinc+2*(yinc-1) + 4*(2-zinc))

      ccoftps(:,:,:,6) = ccof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      ecoftps(:,:,:,6) = ecof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      tcoftps(:,:,:,6) = tcof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))
      icoftps(:,:,:,6) = icof8(:,:,:,3-xinc+2*(yinc-1) + 4*(2-zinc))

      ccoftps(:,:,:,7) = ccof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      ecoftps(:,:,:,7) = ecof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      tcoftps(:,:,:,7) = tcof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))
      icoftps(:,:,:,7) = icof8(:,:,:,xinc+2*(2-yinc) + 4*(2-zinc))

      ccoftps(:,:,:,8) = ccof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      ecoftps(:,:,:,8) = ecof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      tcoftps(:,:,:,8) = tcof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))
      icoftps(:,:,:,8) = icof8(:,:,:,3-xinc+2*(2-yinc) + 4*(2-zinc))


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,1),aflx1,
     &                             xshomtps(:,1),xshom1,
     &                             asrcmtps(:,:,:,1),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,1),icoftps(:,:,:,1),
     &                             ecoftps(:,:,:,1),tcoftps(:,:,:,1),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)
      

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,2),aflx1,
     &                             xshomtps(:,2),xshom1,
     &                             asrcmtps(:,:,:,2),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,2),icoftps(:,:,:,2),
     &                             ecoftps(:,:,:,2),tcoftps(:,:,:,2),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)


        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,3),aflx1,
     &                             xshomtps(:,3),xshom1,
     &                             asrcmtps(:,:,:,3),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,3),icoftps(:,:,:,3),
     &                             ecoftps(:,:,:,3),tcoftps(:,:,:,3),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (zinc-1)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (zinc-1)*(jmax-jmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,4),aflx1,
     &                             xshomtps(:,4),xshom1,
     &                             asrcmtps(:,:,:,4),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,4),icoftps(:,:,:,4),
     &                             ecoftps(:,:,:,4),tcoftps(:,:,:,4),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)


!-----------------------------------------------------------------------

         CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,5),aflx1,
     &                             xshomtps(:,5),xshom1,
     &                             asrcmtps(:,:,:,5),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,5),icoftps(:,:,:,5),
     &                             ecoftps(:,:,:,5),tcoftps(:,:,:,5),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (yinc-1)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (yinc-1)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,6),aflx1,
     &                             xshomtps(:,6),xshom1,
     &                             asrcmtps(:,:,:,6),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,6),icoftps(:,:,:,6),
     &                             ecoftps(:,:,:,6),tcoftps(:,:,:,6),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)  

        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (xinc-1)*(imax-imin+2))/2 ,
     &   (imin + imax + (xinc-1)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,7),aflx1,
     &                             xshomtps(:,7),xshom1,
     &                             asrcmtps(:,:,:,7),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,7),icoftps(:,:,:,7),
     &                             ecoftps(:,:,:,7),tcoftps(:,:,:,7),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)       
        
        CALL RECADA_ONE_OCTANT(nn,ng, ndir,nr,nh,
     &                             nx,ny, nz,
     &                             xinc,yinc,zinc, oct,
     &   (2*imin + (2-xinc)*(imax-imin+2))/2,
     &   (imin + imax + (2-xinc)*(imax-imin))/2,
     &   (2*jmin + (2-yinc)*(jmax-jmin+2))/2,
     &   (jmin + jmax + (2-yinc)*(jmax-jmin))/2,
     &   (2*kmin + (2-zinc)*(kmax-kmin+2))/2,
     &   (kmin + kmax + (2-zinc)*(kmax-kmin))/2,
     &                             niv+1,
     &                             sigt,
     &                             mu,eta,ksi,w,
     &                             sgnc,sgni,sgne,sgnt,
     &                             zreg,asrc,aflx,
     &                             bflx, bfly, bflz,
     &                             aflxtps(:,:,8),aflx1,
     &                             xshomtps(:,8),xshom1,
     &                             asrcmtps(:,:,:,8),asrcm1,finc0,finc1,
     &                             fout0,fout1,
     &                             ccoftps(:,:,:,8),icoftps(:,:,:,8),
     &                             ecoftps(:,:,:,8),tcoftps(:,:,:,8),
     &                             ccof8,icof8,ecof8,tcof8,tol,
     &                             errcor,errmul,ok,delt3,i,j,k,r)       

    !- write the outgoing angular flux of node-0 on the fine-mesh
    ! boundary flux
         RETURN
      
        ELSE 
!   Case where the error is not guaranted but the pixel-lvl is reached
            CALL SPLITAFLX1(nn,nr,nc,nx,ny,nz,
     &                      imin,imax,jmin,jmax,kmin,kmax,
     &                      aflx,aflx1)


            CALL SPLITBOUND1(nn,ng,nb,
     &                      imin, imax,jmin,jmax,kmin,kmax,
     &                      nx,ny,nz,
     &                      bflx, bfly, bflz, fout1)

      ENDIF

      END SUBROUTINE RECADA_ONE_OCTANT