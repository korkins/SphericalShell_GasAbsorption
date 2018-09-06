!******************************************************************************* 
!TRYME: Use <U0> = (U0(IB) + U0(IB+1))/2 in SS
!*******************************************************************************
SUBROUTINE SGLSCAT_TOA_PS(SZG, NSZ, VZG, NVZ, AZA, NAZ, HKM, TAU, NB, NEWLR,   &
                          NL, A1, NLR, I1PS)
!===============================================================================
! PURPOSE:
!   To compute the single scattering radiation on TOA in the pseudospherical
!   approximation.
!
! INPUT:
!   SZG     D(NSZ)        SZA at the ground WRT Z-axis (local normal):
!                               SZG = [0, 90) deg.
!   NSZ     I(1)          Number of SZAs
!   VZG     D(NVZ)        VZA at the ground WRT Z-axis (local normal):
!                               VZG = [0, 90) deg.
!   NVZ     I(1)          Number of VZAs
!   AZA     D(NAZ)        Azimuth in degrees
!   NAZ     I(1)          
!   HKM     D(NB)         Heights of boundaries:
!                               TOA = HKM(1) > HKM(i) > HKM(NB) = BOA 
!   TAU     D(NB)         Profile of OT, TOA to BOA:
!                               0 = TOA = TAU(1) < TAU(i) < TAU(NB) = Tau_total
!   NB      I(1)          Number of boundaries
!   NEWLR   I(NL)         Profile of microlayers in each optical layer
!   NL      I(1)          Number of microlayers, NL=NB-1
!   A1      D(..., NLR)   SSA*Phase_function
!   NLR     I(1)          Number of optical layers
!
! OUTPUT:
!   I1PS      D(NSG, NVG, NAZ)   Primary scattered radiation
!
! TREE:
!   -
!
! COMMENTS:
!   SZA - Solar Zenith Angle
!   VZA - View Zenith Angle
!   OT  - Optical thickness
!   WRT - With Respect To
!   TOA - Top Of Atmosphere
!   BOA - Bottom Of Atmosphere
!   ChF - The Chapman Function
!
!   This subroutine is intended to work for ascending light only.
!
!   See 085.doc for theoretica lderivations.
!
!   *** IT IS UNCLEAR HOW AZIMUTH VARIES WITH HEIGHT ***
!
! REFERENCES:
!   1. Dahlback A & Stamnes K, 1991: Planet.Space.Sci., 39(5), p671, Appendix B
!===============================================================================
!
    IMPLICIT NONE
!
! PARAMETERS
    REAL*8, PARAMETER :: PI =  3.1415926535897932D0 ! The number PI
    REAL*8, PARAMETER :: D2R = PI/180.0D0           ! To convert to radians
    REAL*8, PARAMETER :: R2D = 180.0D0/PI           ! To convert to degrees
    !REAL*8, PARAMETER :: RE = 6356.8D0              ! Earth radius, km
    REAL*8, PARAMETER :: RE = 6371.0D0              ! Earth radius, km: Adams & Kattawar
!
! INPUT VARIABLES
    INTEGER, INTENT(IN) :: NSZ, NVZ, NAZ, NB, NL, NLR
!
! INPUT ARRAYS
    INTEGER, DIMENSION(NL), INTENT(IN) :: NEWLR
    REAL*8, DIMENSION(NSZ), INTENT(IN) :: SZG
    REAL*8, DIMENSION(NVZ), INTENT(IN) :: VZG
    REAL*8, DIMENSION(NAZ), INTENT(IN) :: AZA
    REAL*8, DIMENSION(NB), INTENT(IN) :: HKM, TAU
    REAL*8, DIMENSION(NAZ, NSZ, NVZ, NLR), INTENT(IN) :: A1
!
! OUTPUT ARRAYS
    REAL*8, DIMENSION(NAZ, NSZ, NVZ), INTENT(OUT) :: I1PS
!
! LOCAL VARIABLES
    INTEGER &
        IB,  & ! Loop index over boundaries, IB = 1(TOA):NB(BOA)
        IL,  & ! Loop index over layers, IL = 1:NL
        ILR, & ! Index for optical layer 
        IVZ, & ! Loop index over VZA
        ISZ, & ! Loop index over SZA
        IAZ, & ! Loop index over AZA
        JB     ! Loop index over boundaries JB = 2:IB   
    REAL*8 &
        DVZ, & ! VZR - VZP
        RP,  & ! RE + HKM(IB) for a given IB
        CNU, & ! cos(VZR - VZP), VZR .GE. VZP
        SNU, & ! sin(VZR - VZP), VZR .GE. VZP
        PS2, & ! [RP*sin(SZA)]**2
        PV2, & ! Same as PS2 except for VZA
        U      ! -cos(VZP) < 0
!
! LOCAL ARRAYS
    REAL*8, DIMENSION(NL) :: &
        DTAU, & ! TAU(I+1) - TAU(I)
        EXT,  & ! Extinction in layers, dTau/dH
        R2      ! RP*RP
    REAL*8, DIMENSION(NVZ) :: &
        VZR, & ! VZG in radians
        SVR    ! sin(VZR)
    REAL*8, DIMENSION(NSZ) :: &
        SZR, & ! SZG in radians
        MU0, & ! cos(SZR)
        SMU0   ! sin(SZR)
    REAL*8, DIMENSION(NAZ) :: &
        CAZ, & ! cos(AZA), note AZA must be converted to rad.
        U0,  & ! cos(SZP) > 0
        UU0    ! U0/(U0 - U) for the SS
    REAL*8, DIMENSION(NVZ, NL) :: &
        CHV, & ! ChF for view direction
        VZP    ! VZA=f(HKM, VZG), except for HKM=BOA, in radians
    REAL*8, DIMENSION(NAZ, NSZ, NVZ, NL) :: &
        CHS, & ! ChF for Solar direction
        SZP    ! VZA=f(HKM, geometry), except for HKM=BOA, in radians
!===============================================================================
!
    VZR = VZG*D2R
    SZR = SZG*D2R
!
    SVR = SIN( VZR )
    CAZ = COS( AZA*D2R )
    MU0 = COS( SZR )
    SMU0 = SIN( SZR )
!
    DO IL = 1, NL
        DTAU(IL) = TAU(IL+1) - TAU(IL)
        EXT(IL) =  DTAU(IL)/(HKM(IL) - HKM(IL+1))
        RP = RE + HKM(IL)
        R2(IL) = RP*RP
        DO IVZ = 1, NVZ
            VZP(IVZ, IL) = ASIN( SVR(IVZ)*RE/RP )  ! in rad.
            DVZ = VZR(IVZ) - VZP(IVZ, IL)          ! in rad.
            CNU = COS( DVZ )
            SNU = SIN( DVZ )
            DO ISZ = 1, NSZ
                SZP(:, ISZ, IVZ, IL) = &           ! in rad.
                    ACOS( MU0(ISZ)*CNU - SNU*SMU0(ISZ)*CAZ )
            END DO ! DO ISZ = 1, NSZ
        END DO ! IVZ = 1, NVZ
    END DO ! IL = 1, NL
!
!   ChF for the View direction, CHF. CHF(1) = CHF(TOA) = 0.0
!
!   Loop over boundaries, excluding TOA & BOA
    DO IB = 2, NL
        DO IVZ = 1, NVZ
            PV2 = R2(IB)*SIN( VZP(IVZ, IB) )**2
!           Accumulate ChF: initialize with 0.0
            CHV(IVZ, IB) = 0.0D0
            DO JB = 2, IB
                CHV(IVZ, IB) = CHV(IVZ, IB) + &
                    EXT(JB-1)*(SQRT( R2(JB-1) - PV2 ) - SQRT( R2(JB) - PV2 ))
            END DO ! JB = 2, IB
        END DO ! IVZ = 1, NVZ
    END DO ! IB = 2, NL
!
!   ChF for the View direction, CHF. CHF(1) = CHF(TOA) = 0.0
!
!!!!   Loop over boundaries, excluding TOA & BOA
!!!    DO IB = 2, NL
!!!        DO IVZ = 1, NVZ
!!!            DO ISZ = 1, NSZ
!!!                DO IAZ = 1, NAZ
!!!                    PS2 = R2(IB)*SIN( SZP(IAZ, ISZ, IVZ, IB) )**2
!!!!                   Accumulate ChF-s: initialize with 0.0
!!!                    CHS(IAZ, ISZ, IVZ, IB) = 0.0D0    
!!!                    DO JB = 2, IB
!!!                        CHS(IAZ, ISZ, IVZ, IB) = CHS(IAZ, ISZ, IVZ, IB) + &
!!!                            EXT(JB-1)*(SQRT( R2(JB-1) - PS2 ) - SQRT( R2(JB) - PS2 ))
!!!                    END DO ! JB = 2, IB
!!!                END DO ! IAZ = 1, NAZ
!!!            END DO ! ISZ = 1, NSZ
!!!        END DO ! IVZ = 1, NVZ
!!!    END DO ! IB = 2, NL
!
!   TOA, IL = IB = 1: I1 is initilaized here for subsequent accumulation
    DO IVZ = 1, NVZ
        U = -COS( VZP(IVZ, 1) ) ! IL = 1, MUP = cos(pi-VZP) = -cos(VZP) => MUP < 0
        DO ISZ = 1, NSZ
            U0 = COS( SZP(:, ISZ, IVZ, 1) ) ! IL = 1
            UU0 = U0/( U0 - U )
            I1PS(:, ISZ, IVZ) = A1(:, ISZ, IVZ, 1)*UU0* & ! ILR = 1
                                    (1.0D0 - EXP(DTAU(1)/U/UU0)) ! IL = 1
        END DO ! ISZ = 1, NSZ
    END DO ! IVZ = 1, NVZ
!
!
    DO IL = 2, NL
        ILR = NEWLR(IL) ! Optical layer
        DO IVZ = 1, NVZ
            U = -COS( VZP(IVZ, IL) ) ! IL = 1, MUP = cos(pi-VZP) = -cos(VZP) => MUP < 0
            DO ISZ = 1, NSZ
                U0 = COS( SZP(:, ISZ, IVZ, IL) ) ! IL = 1
                UU0 = U0/(U0 - U)
                I1PS(:, ISZ, IVZ) = I1PS(:, ISZ, IVZ) + &
                                    A1(:, ISZ, IVZ, ILR)*UU0* &
                                    (1.0D0 - EXP( DTAU(IL)/U/UU0 ))* &
                                    !EXP(-CHS(:, ISZ, IVZ, IL)-CHV(IVZ, IL))
                                    EXP( -SUM(DTAU(1:IL-1)) )*EXP(-CHV(IVZ, IL))
            END DO ! ISZ = 1, NSZ
        END DO ! IVZ = 1, NVZ
    END DO ! IL = 2, NL
!
!!!!        DO I = 1, NVZ
!!!!        VZR = VZA(I)*D2R
!!!!        MU = -COS(VZR) ! VZA' = PI - VZA
!!!!        DO J = 1, NSZ
!!!!            MUMU0 = MU0(J)/(MU0(J) - MU) ! > 0
!!!!            TAUE = 0.0D0
!!!!            I1PP(:, J, I) =  MUMU0*A1(:, J, I, 1)*(1.0D0 - EXP(DTAU(1)/MU/MUMU0))
!!!!            DO L = 2, NL
!!!!                TAUE = TAUE + DTAU(L-1) ! Accumulate
!!!!                I1PP(:, J, I) = I1PP(:, J, I) + &
!!!!                                         MUMU0*A1(:, J, I, L)* &
!!!!                                            EXP(TAUE/MU/MUMU0)* &
!!!!                                                (1.0D0 - EXP(DTAU(L)/MU/MUMU0))
!!!!            END DO
!!!!        END DO
!!!!    END DO
!
END SUBROUTINE SGLSCAT_TOA_PS
!===============================================================================
! 12May17 - 
!===============================================================================