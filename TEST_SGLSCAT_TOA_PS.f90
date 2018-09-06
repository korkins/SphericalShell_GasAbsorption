PROGRAM TEST_SGLSCAT_TOA_PS
!
    IMPLICIT NONE
!
    INTEGER, PARAMETER :: &
        NSZ = 10, & !10, & ! Solar
        NVZ = 18, & ! View
        NAZ = 5, & !5,  & ! Azimuth
        NKM = 31, & ! Heights
        NL = NKM-1, & ! Number of layers
        NLR = NL ! Number of optical layers
    REAL*8, PARAMETER :: &
        DF = 0.03D0, & ! Depolarization factor
        PI =  3.1415926535897932D0, & ! The number PI
        D2R = PI/180.0D0, & ! To convert to radians
        R2D = 180.0D0/PI    ! To convert to degrees
!
    INTEGER IX, I, J, K, L ! Loop indices
!
    INTEGER, DIMENSION(NLR) :: NEWLR
    REAL*8, DIMENSION(NSZ) :: MU0, SMU0, SZA, SZR
    REAL*8, DIMENSION(NVZ) :: VZA
    REAL*8, DIMENSION(NAZ) :: AZA, CAZ
    REAL*8, DIMENSION(NKM) :: ZKM, TAUR, TAUG, TAU
    REAL*8, DIMENSION(NL) :: ssa, DTAU
    REAL*8, DIMENSION(NAZ, NSZ, NVZ) :: F11, I1PS, csca, I1PP
    REAL*8, DIMENSION(NAZ, NSZ, NVZ, NL) :: A1
    REAL*8, DIMENSION(NAZ*NSZ*NVZ, 9) :: IPRT
    REAL*8, DIMENSION(NAZ*NSZ*NVZ, 6) :: RERR
!
    REAL*8 D, D1, MU, SMU, MUMU0, SMUMU0, R11, X, VZR, TAUE
!
    SZA(1:NSZ-1) = (/(10.0D0*(I-1), I = 1, NSZ-1)/); SZA(NSZ) = 85.0D0; !SZA(NSZ) = 0.0D0
    VZA(1:NVZ) = (/( 5.0D0*(J-1), J = 1, NVZ)/); !VZA = 80.0D0
    AZA(1:NAZ) = (/(45.0D0*(K-1), K = 1, NAZ)/)
    ZKM(1:NKM) = (/(31.0D0 - L,   L = 1, NKM)/)
!
    TAUR = (/0.0000000000D+00, 1.6635450000D-03, 3.6055130000D-03, &
             5.8734900000D-03, 8.5243005000D-03, 1.1623726500D-02, &
             1.5251039000D-02, 1.9500706500D-02, 2.4482504000D-02, &
             3.0327700500D-02, 3.7191454000D-02, 4.5235166000D-02, &
             5.4644082000D-02, 6.5650877000D-02, 7.8525107000D-02, &
             9.3585277000D-02, 1.1120688200D-01, 1.3182576700D-01, &
             1.5595153700D-01, 1.8416692700D-01, 2.1660994700D-01, &
             2.5332049700D-01, 2.9471597200D-01, 3.4123207200D-01, &
             3.9334779200D-01, 4.5156866200D-01, 5.1642363200D-01, &
             5.8847942200D-01, 6.6833690700D-01, 7.5663307200D-01, &
             8.5404043200D-01/)
!
    TAUG = (/0.000000000000D+00, 4.007070600000D-03, 8.455106700000D-03, &
             1.338439780000D-02, 1.885947130000D-02, 2.492995700000D-02, &
             3.149607080000D-02, 3.841972620000D-02, 4.561358200000D-02, &
             5.280496170000D-02, 5.988451750000D-02, 6.668893110000D-02, &
             7.293175470000D-02, 7.852502340000D-02, 8.337044140000D-02, &
             8.755435270000D-02, 9.126157210000D-02, 9.459432530000D-02, &
             9.764462500000D-02, 1.003355757000D-01, 1.023912443000D-01, &
             1.039089135000D-01, 1.050845035000D-01, 1.060582170700D-01, &
             1.069769744700D-01, 1.078768796700D-01, 1.087999266000D-01, &
             1.097754882700D-01, 1.108479148700D-01, 1.119803700700D-01, &
             1.131312126700D-01/)
!
    TAU = TAUR + TAUG
!
!   Depolarization constants
    D1 = 1.0D0 + 0.5D0*DF
    D = (1.0D0 - DF)/D1
!****************************************************
!   Adams & Kattawar:
!    D = 1.0     ! no depolarization
!    TAUG = 0.0  ! no absorption
!    TAU(1:NKM) = (/(0.25D0*(L-1)/(NKM-1),   L = 1, NKM)/)
!****************************************************
    CAZ = COS(AZA*D2R)
    SZR = SZA*D2R
    MU0 = COS(SZA*D2R)
    SMU0 = SIN(SZA*D2R)
    DO I = 1, NVZ
        VZR = VZA(I)*D2R
        MU = -COS(VZR) ! VZA' = PI - VZA
        SMU = SIN(VZR)
        DO J = 1, NSZ
            MUMU0  = MU*MU0(J)
            SMUMU0 = SMU*SMU0(J) 
            DO K = 1, NAZ
                X = MUMU0 + SMUMU0*CAZ(K)
                R11 = 0.75D0*(1.0D0 + X*X)
                F11(K, J, I) = D*(R11 - 1.0D0) + 1.0D0
                csca(k, j, i) = X
            END DO
        END DO
    END DO
!
    DO L = 1, NL
        NEWLR(L) = L
        DTAU(L) = TAU(L+1) - TAU(L)
        ssa(L) = (TAUR(L+1) - TAUR(L))/DTAU(L)
        !ssa(L) = 1.0D0 ! Kattawar & Adams
        A1(:, :, :, L) = ssa(L)*F11
    END DO
!
    CALL SGLSCAT_TOA_PS(SZA, NSZ, VZA, NVZ, AZA, NAZ, ZKM, TAU, NKM, NEWLR,   &
                        NL, A1, NLR, I1PS) ! Normalization: F_TOA = 4pi
!
!   Unit flux on TOA
    I1PS = 0.25D0*I1PS/PI !
! -----------------------------------------------------------------------------    
! Compute PP & test VS IPRT B2 - SS only: err  0.00% for all geometries -> ok!
    DO I = 1, NVZ
        VZR = VZA(I)*D2R
        MU = -COS(VZR) ! VZA' = PI - VZA
        DO J = 1, NSZ
            MUMU0 = MU0(J)/(MU0(J) - MU) ! > 0
            TAUE = 0.0D0
            I1PP(:, J, I) =  MUMU0*A1(:, J, I, 1)*(1.0D0 - EXP(DTAU(1)/MU/MUMU0))
            DO L = 2, NL
                TAUE = TAUE + DTAU(L-1) ! Accumulate
                I1PP(:, J, I) = I1PP(:, J, I) + &
                                         MUMU0*A1(:, J, I, L)* &
                                            EXP(TAUE/MU/MUMU0)* &
                                                (1.0D0 - EXP(DTAU(L)/MU/MUMU0))
            END DO
        END DO
    END DO
    I1PP =  0.25D0*I1PP/PI !
!------------------------------------------------------------------------------    
!
!   Read benchmark
    CALL READ2D('SS_IQUV.txt', 11, 0, NSZ*NAZ*NVZ, 9, IPRT)
!
!****************************************************
!   Adams & Kattawar:
!    IPRT(:, 6) = 0.0D0
!    I1PS = PI*I1PS
!    I1PP = PI*I1PP
!****************************************************
!
    DO J = 1, NSZ
        DO K = 1, NAZ
            DO I = 1, NVZ
                IX = (J-1)*NAZ*NVZ + (K-1)*NVZ + I
                RERR(IX, 1:3) = IPRT(IX, 2:4) ! Copy SZA, AZA, VZA
                RERR(IX, 4) = 100.0D0*(1.0D0 - I1PS(K, J, I)/IPRT(IX, 6))
                RERR(IX, 5) = I1PS(K, J, I)
                RERR(IX, 6) = I1PP(K, J, I)
            END DO
       END DO
    END DO
!
    CALL PRINT2D(RERR, NAZ*NSZ*NVZ, 6, 'RERR_SS_PS.txt', 14)
!
    WRITE(*, *) 'Done!'
!
END PROGRAM