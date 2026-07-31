#include "cppdefs.h"
      MODULE sediment_biology_mod

#ifdef SEDBIO
      implicit none
      PRIVATE
      PUBLIC :: sediment_bio

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment_bio(ng, Istr, Iend, LBi, UBi,                &
     &                      dtdays,                                   &
#   ifdef MASKING
     &                      rmask,                                     &
#   endif
     &                      temp, DOX, dNO3, dNH4, dPO4, dSiO4,        &
     &                      FPOC, FPON, FPOP, FPSi,                    &
     &                      bPOC, bPON, bPOP, bPSi,                    &
     &                      bNO3, bNH4, bPO4, bSi,                     &
     &                      bUC, bUN, bUP, bUS,                        &
     &                      JNO3, JNH4, JPO4, JSi, SOD,                &
#ifdef CARBON
     &                      dTIC, bTIC, JTIC,                          &
#ifdef TALK_NONCONSERV
     &                      dAlk, bAlk, JAlk                           &
#endif
#   ifdef CACO3
     &                      ,FPCa, bPCa, bUCa, omega_arag              &
#   endif
#   ifdef DIAGNOSTICS_BIO
     &                      ,Jden                                      &
#   endif
#endif
     &                      )
!***********************************************************************
! ORDER OF ARGUMENTS:   temp, oxygen, dNO3, dNH4, dPO4, dSiO4, FPOC, FPON,
!                       FPOP, FPSi, bPOC, bPON, bPOP, bPSi, bNO3, bNH4,
!                       bPO4, bSi, bUC, bUN, bUP, bUS, JNO3, JNH4, JPO4
!                       JSi, SOD
! ADDED ARGUMENTS:      dTIC, bTIC, JTIC, dAlk, bAlk, JAlk

! CALLED IN BIO UMAINE15 HEADER AS FOLLOWS ------------------------------
!# ifdef SEDBIO
!
!       call sediment_bio(ng, Istr, Iend, LBi, UBi,                    &
!     &                      dtdays,                                    &
!#   ifdef MASKING
!     &                      rmask(LBi:UBi,j),                          &
!#   endif
!     &                      Bio(LBi:,1,itemp),Bio(LBi:,1,iOxyg),       &
!     &                      Bio(LBi:,1,iNO3_),Bio(LBi:,1,iNH4_),       &
!     &                      Bio(LBi:,1,iPO4_),Bio(LBi:,1,iSiOH),       &
!     &                      FPOC, FPON, FPOP, FPSi,                    &
!#   ifdef CACO3
!     &                      FPCa,                                      &
!#   endif
!     &                      bPOC(LBi:,j,:), bPON(LBi:,j,:),            &
!     &                      bPOP(LBi:,j,:), bPSi(LBi:,j,:),            &
!     &                      bNO3(LBi:,j), bNH4(LBi:,j),                &
!     &                      bPO4(LBi:,j), bSi(LBi:,j),                 &
!     &                      bUC(LBi:,j,:), bUN(LBi:,j,:),              &
!     &                      bUP(LBi:,j,:), bUS(LBi:,j,:),              &
!     &                      JNO3(LBi:,j), JNH4(LBi:,j),                &
!     &                      JPO4(LBi:,j), JSi(LBi:,j), SOD(LBi:,j)     &
!#  ifdef CARBON
!     &                      Bio(LBi:,1,iTIC_), bTIC(LBi:,j),           &
!     &                      JTIC(LBi:,j),                              &
!#   ifdef TALK_NONCONSERV
!     &                      Bio(LBi:,1,iTAlk), bAlk(LBi:,j),           &
!     &                      JAlk(LBi:,j)                               &
!#   endif
!#    ifdef CACO3
!     &                      ,FPCa, bPCa(LBi:,j,:), bUCa, omega_arag    & 
!#   ifdef DIAGNOSTICS_BIO
!     &                      ,Jden                                      &
!#   endif       
!#  endif
!     &                      )
! -----------------------------------------------------------------------
! ORDER OF ARGUMENTS:   temp, oxygen, dNO3, dNH4, dPO4, dSiO4, FPOC, FPON,
!                       FPOP, FPSi, bPOC, bPON, bPOP, bPSi, bNO3, bNH4,
!                       bPO4, bSi, bUC, bUN, bUP, bUS, JNO3, JNH4, JPO4
!                       JSi, SOD
! ADDED ARGUMENTS:      dTIC, dAlk, bTIC, bAlk, JTIC, JAlk
! CHANGING THE ORDER OF ARGUMEMNTS ON 3/2/23
!                       dTIC, bTIC, JTIC
!                       dAlk, bAlk, JAlk

      USE mod_kinds
# ifdef CACO3
      USE mod_biology, ONLY : bUmax,bUmaxSi,bdep,balpha,bw,btheta_diag,   &
     &                       btheta_Si,bnit,btheta_nit,bdo_nit,bdenit,    &
     &                       btheta_denit, bpsi_n, bdo_c, bao2, bpi,      &
     &                       bpsi_p, bfc, bfn, bfp, bfs, p2n, o2nh,       &
     &                       bUmaxCa, bfca, omega_thresh
# else      
      USE mod_biology, ONLY : bUmax,bUmaxSi,bdep,balpha,bw,btheta_diag,   &
     &                       btheta_Si,bnit,btheta_nit,bdo_nit,bdenit,    &
     &                       btheta_denit, bpsi_n, bdo_c, bao2, bpi,      &
     &                       bpsi_p, bfc, bfn, bfp, bfs, p2n, o2nh
# endif

      USE mod_sedbio, ONLY : nspc

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Istr, Iend
      integer, intent(in) :: LBi, UBi
      real(r8), intent(in) :: dtdays
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:)
#  endif
      real(r8), intent(in) :: temp(LBi:)
      real(r8), intent(in) :: DOX(LBi:)
      real(r8), intent(in) :: dNO3(LBi:)
      real(r8), intent(in) :: dNH4(LBi:)
      real(r8), intent(in) :: dPO4(LBi:)
      real(r8), intent(in) :: dSiO4(LBi:)
      real(r8), intent(in) :: FPOC(LBi:)
      real(r8), intent(in) :: FPON(LBi:)
      real(r8), intent(in) :: FPOP(LBi:)
      real(r8), intent(in) :: FPSi(LBi:)
#  ifdef CACO3
      real(r8), intent(in) :: FPCa(LBi:)
      real(r8), intent(inout) :: bUCa(LBi:,:)
      real(r8), intent(inout) :: bPCa(LBi:,:)
      real(r8), intent(in) :: omega_arag(LBi:)
#  endif
      real(r8), intent(inout) :: bPOC(LBi:,:)
      real(r8), intent(inout) :: bPON(LBi:,:)
      real(r8), intent(inout) :: bPOP(LBi:,:)
      real(r8), intent(inout) :: bPSi(LBi:,:)
      real(r8), intent(inout) :: bNO3(LBi:)
      real(r8), intent(inout) :: bNH4(LBi:)
      real(r8), intent(inout) :: bPO4(LBi:)
      real(r8), intent(inout) :: bSi(LBi:)
      real(r8), intent(inout) :: bUC(LBi:,:)
      real(r8), intent(inout) :: bUN(LBi:,:)
      real(r8), intent(inout) :: bUP(LBi:,:)
      real(r8), intent(inout) :: bUS(LBi:,:)
      real(r8), intent(inout) :: JNO3(LBi:)
      real(r8), intent(inout) :: JNH4(LBi:)
      real(r8), intent(inout) :: JPO4(LBi:)
      real(r8), intent(inout) :: JSi(LBi:)
      real(r8), intent(inout) :: SOD(LBi:)
#  ifdef DIAGNOSTICS_BIO
	real(r8), intent(inout) :: Jden(LBi:)
#  endif
#ifdef CARBON
      real(r8), intent(in) :: dTIC(LBi:)
      real(r8), intent(inout) :: bTIC(LBi:)
      real(r8), intent(inout) :: JTIC(LBi:)
#ifdef TALK_NONCONSERV
      real(r8), intent(in) :: dAlk(LBi:)
      real(r8), intent(inout) :: bAlk(LBi:)
      real(r8), intent(inout) :: JAlk(LBi:)
#endif
#endif

# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi)
#  endif
      real(r8), intent(in) :: temp(LBi:UBi)
      real(r8), intent(in) :: DOX(LBi:UBi)
      real(r8), intent(in) :: dNO3(LBi:UBi)
      real(r8), intent(in) :: dNH4(LBi:UBi)
      real(r8), intent(in) :: dPO4(LBi:UBi)
      real(r8), intent(in) :: dSiO4(LBi:UBi)
      real(r8), intent(in) :: FPOC(LBi:UBi)
      real(r8), intent(in) :: FPON(LBi:UBi)
      real(r8), intent(in) :: FPOP(LBi:UBi)
      real(r8), intent(in) :: FPSi(LBi:UBi)
#  ifdef CACO3
      real(r8), intent(in) :: FPCa(LBi:UBi)
      real(r8), intent(inout) :: bUCa(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPCa(LBi:UBi,nspc)
      real(r8), intent(in) :: omega_arag(LBi:UBi)
#  endif
      real(r8), intent(inout) :: bPOC(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPON(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPOP(LBi:UBi,nspc)
      real(r8), intent(inout) :: bPSi(LBi:UBi,nspc)
      real(r8), intent(inout) :: bNO3(LBi:UBi)
      real(r8), intent(inout) :: bNH4(LBi:UBi)
      real(r8), intent(inout) :: bPO4(LBi:UBi)
      real(r8), intent(inout) :: bSi(LBi:UBi)
      real(r8), intent(inout) :: bUC(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUN(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUP(LBi:UBi,nspc)
      real(r8), intent(inout) :: bUS(LBi:UBi,nspc)
      real(r8), intent(inout) :: JNO3(LBi:UBi)
      real(r8), intent(inout) :: JNH4(LBi:UBi)
      real(r8), intent(inout) :: JPO4(LBi:UBi)
      real(r8), intent(inout) :: JSi(LBi:UBi)
      real(r8), intent(inout) :: SOD(LBi:UBi)
#  ifdef DIAGNOSTICS_BIO
	real(r8), intent(inout) :: Jden(LBi:UBi)
#  endif
#ifdef CARBON
      real(r8), intent(in) :: dTIC(LBi:UBi)
      real(r8), intent(inout) :: bTIC(LBi:UBi)
      real(r8), intent(inout) :: JTIC(LBi:UBi)
#ifdef TALK_NONCONSERV
      real(r8), intent(in) :: dAlk(LBi:UBi)
      real(r8), intent(inout) :: bAlk(LBi:UBi)
      real(r8), intent(inout) :: JAlk(LBi:UBi)
#endif
#endif
# endif

! local variables:

      integer :: i, k
      real(r8) :: Qt, Qtsi
      real(r8), dimension(nspc) :: JPOC
      real(r8), dimension(nspc) :: JPON
      real(r8), dimension(nspc) :: JPOP
      real(r8), dimension(nspc) :: JPSi
# ifdef CACO3
      real(r8), dimension(nspc) :: JPCa
!      real(r8), parameter :: omega_thresh=1.0_r8
# endif
      real(r8) :: stc, CSOD, NSOD, Jnit, Jdenit, JCO2, O2, DOxyg
      real(r8) :: JAlk_sed, JPOC_sum, Cdenit

      real(r8) :: w2, cff, cff1, cff2, cff3, D, B
      real(r8), parameter :: TR = 20.0_r8
      real(r8), parameter :: eps = 1.0e-10_r8
      real(r8), parameter :: o2conv = 32.0e-03_r8
      real(r8), parameter :: ta2po4 = 5.8_r8

      w2=bw(ng)*2.74e-5_r8   ! conversion from cm/yr to m/day

        cff=dtdays/bdep(ng)
        DO i=Istr,Iend
# ifdef DIAGNOSTICS_BIO
          Jden(i)=0.0_r8
# endif
# ifdef MASKING
          IF( rmask(i).gt.0.0) THEN
# endif
! if near-bottom DOX is less than zero, make it zero:
            DOxyg=max(DOX(i),0.0_r8)
!            O2=DOX(i)*o2conv !mmol.m-3 to mg/L
            O2=DOxyg*o2conv !mmol.m-3 to mg/L
            Qt=btheta_diag(ng)**(temp(i)-TR)
            Qtsi=btheta_Si(ng)**(temp(i)-TR)
            DO k=1,nspc
!        update sediment POM concentrations in response to deposition (source term)
                bPOC(i,k)=bPOC(i,k)+FPOC(i)*bfc(ng,k)*cff
                bPON(i,k)=bPON(i,k)+FPON(i)*bfn(ng,k)*cff
                bPOP(i,k)=bPOP(i,k)+FPOP(i)*bfp(ng,k)*cff
                bPSi(i,k)=bPSi(i,k)+FPSi(i)*bfs(ng,k)*cff
# ifdef CACO3
                bPCa(i,k)=bPCa(i,k)+FPCa(i)*bfca(ng,k)*cff
# endif
!        update sediment POM concentrations in response to diagenesis and burial
!        get diagenesis fluxes
                cff1=bUC(i,k)*Qt*dtdays
                cff2=w2*cff
                bPOC(i,k)=bPOC(i,k)/(1.0_r8+cff1+cff2)
                ! diagenesis flux (mmol/m3):
                JPOC(k)=cff1*bPOC(i,k)
                cff1=bUN(i,k)*Qt*dtdays
                bPON(i,k)=bPON(i,k)/(1.0_r8+cff1+cff2)
                ! diagenesis flux:
                JPON(k)=cff1*bPON(i,k)
                cff1=bUP(i,k)*Qt*dtdays
                bPOP(i,k)=bPOP(i,k)/(1.0_r8+cff1+cff2)
                ! diagenesis flux:
                JPOP(k)=cff1*bPOP(i,k)
                cff1=bUS(i,k)*Qtsi*dtdays
                bPSi(i,k)=bPSi(i,k)/(1.0_r8+cff1+cff2)
                ! diagenesis flux:
                JPSi(k)=cff1*bPSi(i,k) 
# ifdef CACO3
                if( omega_arag(i).le.omega_thresh(ng) ) then
                ! dissolution of aragonite
                  cff1=bUCa(i,k)*dtdays
                  bPCa(i,k)=bPCa(i,k)/(1.0_r8+cff1+cff2)
                  ! diagenesis flux:
                  JPCa(k)=cff1*bPCa(i,k)
		    else
		    ! no dissolution
		      JPCa(k)=0.0_r8
		    end if
# endif
            END DO
! update porewater concentrations due to diagenesis fluxes (note that TIC and Alk handled below)
            cff1=0.0_r8
            cff2=0.0_r8
            cff3=0.0_r8
            DO k=1,nspc
                cff1=cff1+JPON(k)
                cff2=cff2+JPOP(k)
                cff3=cff3+JPSi(k)
            END DO
            bNH4(i)=bNH4(i)+cff1
            bPO4(i)=bPO4(i)+cff2
            bSi(i)=bSi(i)+cff3
            
 !         get sum of POC remineralization
            JPOC_sum=0.0_r8
            DO k=1,nspc
               JPOC_sum=JPOC_sum+JPOC(k)
		END DO
!          CSOD needs to be computed prior to denitrification
!            CSOD=0.0_r8
!            DO k=1,nspc
!                CSOD=CSOD+JPOC(k)
!            END DO
            ! CSOD converted to mmol/m^2
!            CSOD=CSOD*(O2/(O2+bdo_c(ng)))*bdep(ng)
		CSOD=JPOC_sum*(O2/(O2+bdo_c(ng)))*bdep(ng)
            ! JCO2: flux of carbon from POC remineralization into sediment pore water
            JCO2 = CSOD
            ! difference between JPOC_sum and CSOD is POC oxidized by oxidants other than O2 (assumed here to be NO3...denitrification)
            Cdenit=JPOC_sum*bdep(ng)-CSOD
            
!        Nitrification (source of NO3, sink of NH4)
            cff1=bnit(ng)*dtdays*(O2/(O2+bdo_nit(ng)))*btheta_nit(ng)**(temp(i)-TR)
            bNH4(i)=bNH4(i)/(1.0_r8+cff1)
!        Nitrification flux (units mmol N/m^3)
            Jnit=cff1*bNH4(i)
            bNO3(i)=bNO3(i)+Jnit

!        Denitrification (sink of NO3, goes to N2 gas, so simply a loss term)
!            cff1=dtdays*bdenit(ng)*(btheta_denit(ng)**(temp(i)-TR))*   &
!     &           exp(-O2/bpsi_n(ng))
! removing O2 dependence:
            cff1=dtdays*bdenit(ng)*(btheta_denit(ng)**(temp(i)-TR))
!            IF( cff1*bNO3(i).lt.0.8_r8*CSOD/bdep(ng) ) THEN
            IF( cff1*bNO3(i).lt.0.8_r8*Cdenit/bdep(ng) ) THEN
             ! treat implicitly
	          bNO3(i)=bNO3(i)/(1.0_r8+cff1)
	          Jdenit=cff1*bNO3(i)
		ELSE
		 ! treat explicitly since in this case Jdenit doesn't depend on bNO3 
!		    Jdenit=0.8_r8*CSOD/bdep(ng)
		    Jdenit=0.8_r8*Cdenit/bdep(ng)
		    bNO3(i)=bNO3(i)-Jdenit
		END IF
#ifdef DIAGNOSTICS_BIO
! convert denitrification to units mmol N/m^2/day
            Jden(i)=Jdenit*bdep(ng)/dtdays
#endif

!           NSOD converted to units of mmol/m^2
            NSOD=bao2(ng)*Jnit*bdep(ng)
!            if( CSOD.lt.0.0_r8.or.NSOD.lt.0.0_r8) then
!                write(*,*)'i=',i,'CSOD=',CSOD,'NSOD=',NSOD
!                write(*,*)'DOX(i)=',DOX(i),'JPOC(1:3)=',JPOC(1:3)
!            end if
            SOD(i)=CSOD+NSOD

            ! Alkalinity flux into pore water
            ! OM degradation to NH4: dTA = + n - p
            D = o2nh(ng)
            B = p2n(ng)
!            JAlk_sed = CSOD * (-(B+1)/D) - NSOD + Jdenit*bdep(ng)
            JAlk_sed = CSOD * ((1 - B * ta2po4)/D) - NSOD + Jdenit*bdep(ng)

!         Update porewater TIC and Alk
#ifdef CARBON
            bTIC(i)=bTIC(i) + JCO2/bdep(ng) + 1.2_r8*Jdenit
# ifdef CACO3
!         Update porewater TIC with contribution from CaCO3 dissolution
            DO k=1,nspc
              bTIC(i)=bTIC(i)+JPCa(k)
		ENDDO
# endif
#ifdef TALK_NONCONSERV
            bAlk(i) = bAlk(i) + JAlk_sed/bdep(ng)
# ifdef CACO3
!         Update porewater Alk with contribution from CaCO3 dissolution
            DO k=1,nspc
              bAlk(i)=bAlk(i)+2.0_r8*JPCa(k)
		ENDDO
# endif
#endif
#endif

!        convert SOD to units of mmol/m2/day
            SOD(i)=SOD(i)/dtdays
!        Sediment-water Nutrient fluxes
!         exchange coefficient (m/day):
            stc=SOD(i)/max(DOxyg,1.0e-6_r8)
            cff1=stc*dtdays/bdep(ng)
!         changes in porewater concentration due to fluxes:
            bNO3(i)=bNO3(i)/(1.0_r8+cff1)+cff1*dNO3(i)/(1.0_r8+cff1)
            bNH4(i)=bNH4(i)/(1.0_r8+cff1)+cff1*dNH4(i)/(1.0_r8+cff1)
            bSi(i)=bSi(i)/(1.0_r8+cff1)+cff1*dSiO4(i)/(1.0_r8+cff1)
            cff=bpi(ng)*exp(-O2/bpsi_p(ng))
            bPO4(i)=bPO4(i)/(1.0_r8+cff*cff1)+cff1*dPO4(i)/(1.0_r8+cff*cff1)
#ifdef CARBON
            bTIC(i)=bTIC(i)/(1.0_r8+cff1)+cff1*dTIC(i)/(1.0_r8+cff1)
# ifdef TALK_NONCONSERV
            bAlk(i)=bAlk(i)/(1.0_r8+cff1)+cff1*dAlk(i)/(1.0_r8+cff1)
# endif
#endif
!          fluxes (multiply by layer thickness (bdep) and divide by dtdays to get mmol/m2/day)
            JNO3(i)=cff1*bdep(ng)*(bNO3(i)-dNO3(i))/dtdays 
            JNH4(i)=cff1*bdep(ng)*(bNH4(i)-dNH4(i))/dtdays 
            JSi(i)=cff1*bdep(ng)*(bSi(i)-dSiO4(i))/dtdays           
            JPO4(i)=cff1*bdep(ng)*(bPO4(i)*cff-dPO4(i))/dtdays
#ifdef CARBON
            JTIC(i) = cff1*bdep(ng)*(bTIC(i) - dTIC(i))/dtdays
# ifdef TALK_NONCONSERV
            JAlk(i) = cff1*bdep(ng)*(bAlk(i) - dAlk(i))/dtdays
# endif
#endif
!         ensure that concentrations are >= 0
            bNO3(i)=max(bNO3(i),0.0_r8)
            bNH4(i)=max(bNH4(i),0.0_r8)
            bPO4(i)=max(bPO4(i),0.0_r8)
            bSi(i)=max(bSi(i),0.0_r8)
#ifdef CARBON
            bTIC(i) = max(bTIC(i),0.0_r8)
#ifdef TALK_NONCONSERV
            bAlk(i) = max(bAlk(i),0.0_r8)
#endif
#endif
            DO k=1,nspc
                bPOC(i,k)=max(bPOC(i,k),0.0_r8)
                bPON(i,k)=max(bPON(i,k),0.0_r8)
                bPOP(i,k)=max(bPOP(i,k),0.0_r8)
                bPSi(i,k)=max(bPSi(i,k),0.0_r8)
            END DO
!         compute decay rates
            DO k=1,nspc
                bUC(i,k)=bUC(i,k)+balpha(ng)*bUmax(ng,k)*bUC(i,k)*dtdays &
     &              -FPOC(i)*bfc(ng,k)*dtdays*bUC(i,k)/bPOC(i,k)/bdep(ng)
                bUN(i,k)=bUN(i,k)+balpha(ng)*bUmax(ng,k)*bUN(i,k)*dtdays &
     &              -FPON(i)*bfn(ng,k)*dtdays*bUN(i,k)/bPON(i,k)/bdep(ng)
                bUP(i,k)=bUP(i,k)+balpha(ng)*bUmax(ng,k)*bUP(i,k)*dtdays &
     &              -FPOP(i)*bfp(ng,k)*dtdays*bUP(i,k)/bPOP(i,k)/bdep(ng)
!         for silicate, decay rates are not bacterially-mediated, so set to max values
                bUS(i,k)=bUmaxSi(ng,k)
                
#ifdef CACO3
!         for CACO3, decay rates are not bacterially-mediated, so set to max values
                bUCa(i,k)=bUmaxCa(ng,k)
#endif
!         apply limits to decay rates and ensure positivity
                bUC(i,k)=max( min(bUC(i,k),bUmax(ng,k)),0.0_r8 )
                bUN(i,k)=max( min(bUN(i,k),bUmax(ng,k)),0.0_r8 )
                bUP(i,k)=max( min(bUP(i,k),bUmax(ng,k)),0.0_r8 )
!         if particulate concentration < eps, then set decay rate to max value
                if( bPOC(i,k).lt.eps) bUC(i,k)=bUmax(ng,k)
                if( bPON(i,k).lt.eps) bUN(i,k)=bUmax(ng,k)
                if( bPOP(i,k).lt.eps) bUP(i,k)=bUmax(ng,k)
            END DO

# ifdef MASKING
          END IF
# endif
        END DO


      RETURN
      END SUBROUTINE sediment_bio
#endif
      END MODULE sediment_biology_mod
