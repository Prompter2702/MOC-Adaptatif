      SUBROUTINE BCD2R2(oct,xout,yout,bflx,bfly,mu,eta,w,pisn,rdir)

!     Boundary conditions on exiting boundary of XY rectangular mesh.
!     The routines called here set the incoming angular flux on the
!     boundary reached at the end of Sn sweep.

      USE FLXCOM
      USE FLGOCT
      USE FLGINC
      USE FLGBCD

      IMPLICIT NONE

      INTEGER  :: oct,xout,yout
      REAL     :: bflx(nb,nbfx,2,*),bfly(nb,nbfy,2,*)
      REAL(KIND=8) :: mu(ndir),eta(ndir),w(ndir),pisn
      INTEGER      :: rdir(nd,*)

!     x-sides.

      SELECT CASE (typc(xout,1))

      CASE (SpecularReflection,AxialSymmetry)
         CALL SPECRF(bflx,nb,nbf(1),ndir,xout,xoct(oct),rdir(1,1))

      CASE (Translation)
         CALL TRANSL(bflx,nb,nbf(1),ndir,3-xout,oct)

      CASE (HalfPiRotation,TranslationOfHalfPi,AxialHalfPiRotation)
         CALL TRANP4(bfly,bflx,nb,nbf(1),ndir,xout,xoc4(oct),rdir(1,4))

      CASE (HalfPiRotationN,TranslationOfHalfPiN)
    !      CALL TRANP4_TST(bfly,bflx,nb,nbf(2),nbf(1),ndir,
    !  &                   xout,xoc4(oct),rdir(1,4),
    !  &                   AuxBc%addx(1),
    !  &                   AuxBc%surx(1),
    !  &                   AuxBc%matx(1,1))

      CASE (IsotropicReflection)
         CALL ALBDRF(1,xout,bflx,mu ,w,pisn,nb,nbf(1),nd,ndir,noct)

      END SELECT

!     y-sides.

      SELECT CASE (typc(yout,2))

      CASE (SpecularReflection,AxialSymmetry)
         CALL SPECRF(bfly,nb,nbf(2),ndir,yout,yoct(oct),rdir(1,2))

      CASE (Translation)
         CALL TRANSL(bfly,nb,nbf(2),ndir,3-yout,oct)

      CASE (HalfPiRotation,TranslationOfHalfPi,AxialHalfPiRotation)
         CALL TRANP4(bflx,bfly,nb,nbf(2),ndir,yout,yoc4(oct),rdir(1,3))

      CASE (HalfPiRotationN,TranslationOfHalfPiN)
    !      CALL TRANP4_TST(bflx,bfly,nb,nbf(1),nbf(2),ndir,
    !  &                   yout,yoc4(oct),rdir(1,3),
    !  &                   AuxBc%addx(2+nbf(2)),
    !  &                   AuxBc%surx(1+AuxBc%ni),
    !  &                   AuxBc%matx(1,1+AuxBc%ni))

      CASE (IsotropicReflection)
         CALL ALBDRF(2,yout,bfly,eta,w,pisn,nb,nbf(2),nd,ndir,noct)

      END SELECT

      END SUBROUTINE
      
!-----------------------------------------------------------------------------------      
!-----------------------------------------------------------------------------------      

      SUBROUTINE BCD2R3(oct,xout,yout,zout,
     &                  bflx,bfly,bflz,mu,eta,ksi,w,pisn,rdir)

!     Boundary conditions on exiting boundary of XYZ rectangular mesh.
!     The routines called here set the incoming angular flux on the
!     boundary reached at the end of Sn sweep.

      USE FLXCOM
      USE FLGOCT
      USE FLGINC
      USE FLGBCD

      IMPLICIT NONE

      INTEGER  :: oct,xout,yout,zout
      REAL     :: bflx(nb,nbfx,2,*),bfly(nb,nbfy,2,*),bflz(nb,nbfz,2,*)
      REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
      INTEGER      :: rdir(nd,*)

      ! x sides 
      SELECT CASE (typc(xout,1))
      
      CASE (SpecularReflection,AxialSymmetry)
         CALL SPECRF(bflx,nb,nbf(1),ndir,xout,xoct(oct),rdir(1,1))
         
      CASE (IsotropicReflection)
         CALL ALBDRF(1,xout,bflx,mu ,w,pisn,nb,nbf(1),nd,ndir,noct)

      CASE (Translation)
         CALL TRANSL(bflx,nb,nbf(1),ndir,3-xout,oct)

      CASE (HalfPiRotation,TranslationOfHalfPi,AxialHalfPiRotation)
         CALL TRANP4(bfly,bflx,nb,nbf(1),ndir,xout,xoc4(oct),rdir(1,5))
         
      CASE (HalfPiRotationN,TranslationOfHalfPiN)
    !      CALL TRANP4_TST(bfly,bflx,nb,nbf(2),nbf(1),ndir,
    !  &                   xout,xoc4(oct),rdir(1,5),
    !  &                   AuxBc%addx(1),
    !  &                   AuxBc%surx(1),
    !  &                   AuxBc%matx(1,1))
      END SELECT

      ! y sides 
      SELECT CASE (typc(yout,2))
      CASE (SpecularReflection,AxialSymmetry)
         CALL SPECRF(bfly,nb,nbf(2),ndir,yout,yoct(oct),rdir(1,2))

      CASE (Translation)
         CALL TRANSL(bfly,nb,nbf(2),ndir,3-yout,oct)

      CASE (HalfPiRotation,TranslationOfHalfPi,AxialHalfPiRotation)
         CALL TRANP4(bflx,bfly,nb,nbf(2),ndir,yout,yoc4(oct),rdir(1,4))
         
      CASE (HalfPiRotationN,TranslationOfHalfPiN)
    !      CALL TRANP4_TST(bflx,bfly,nb,nbf(1),nbf(2),ndir,
    !  &                   yout,yoc4(oct),rdir(1,4),
    !  &                   AuxBc%addx(2+nbf(2)),
    !  &                   AuxBc%surx(1+AuxBc%ni),
    !  &                   AuxBc%matx(1,1+AuxBc%ni))

      CASE (IsotropicReflection)
         CALL ALBDRF(2,yout,bfly,eta,w,pisn,nb,nbf(2),nd,ndir,noct)
      END SELECT

      ! z sides 
      SELECT CASE (typc(zout,3))
      CASE (SpecularReflection,AxialSymmetry)
         CALL SPECRF(bflz,nb,nbf(3),ndir,zout,zoct(oct),rdir(1,3))
      CASE (IsotropicReflection)
         CALL ALBDRF(3,zout,bflz,ksi,w,pisn,nb,nbf(3),nd,ndir,noct)
      END SELECT

      END SUBROUTINE


      SUBROUTINE BCD2R3MPN(oct,xout,yout,zout,bflx,bfly,bflz,
     &                     mu,eta,ksi,w,pisn,rdir)
          
         !     Boundary conditions on exiting boundary of XYZ rectangular mesh.
         !     The routines called here set the incoming angular flux on the
         !     boundary reached at the end of Sn sweep.
         
         USE FLXCOM
         USE FLGOCT
         USE FLGINC
         USE FLGBCD                   
         
         IMPLICIT NONE
         
         INTEGER  :: oct,xout,yout,zout
         REAL     :: bflx(nb_nhc,nbfx,2,*)
         REAL     :: bfly(nb_nhc,nbfy,2,*)
         REAL     :: bflz(nb_nhc,nbfz,2,*)
         REAL(KIND=8) :: mu(ndir),eta(ndir),ksi(ndir),w(ndir),pisn
         INTEGER      :: rdir(nd,*)
         
         !     x-sides.
   
         SELECT CASE (typc(xout,1))
         
         CASE (SpecularReflection)
          CALL SPECRF(bflx,nb_nhc,nbf(1),ndir,xout,xoct(oct),rdir(1,1))
         END SELECT
         
         !     y-sides.
         
         SELECT CASE (typc(yout,2))
   
         CASE (SpecularReflection)
          CALL SPECRF(bfly,nb_nhc,nbf(2),ndir,yout,yoct(oct),rdir(1,2))
         END SELECT
         
         !     z-sides.
         
         SELECT CASE (typc(zout,3))
         
         CASE (SpecularReflection)
          CALL SPECRF(bflz,nb_nhc,nbf(3),ndir,zout,zoct(oct),rdir(1,3))
         END SELECT
          
      END SUBROUTINE
