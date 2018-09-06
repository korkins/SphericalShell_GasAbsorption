SUBROUTINE READ2D(FNAME, L, S, M, N, A)
!===============================================================================
! PURPOSE:
!   To read M-by-N matrix A from a text file FNAME
!
! INPUT:
!   FNAME   A(L)      Input file name
!   L       I(1)      Length of FNAME
!   S       I(1)      Skip S first lines
!   M       I(1)      Number of rows
!   N       I(1)      Number of columns
!
! OUTPUT:
!   A       D(M, N)   Output array
!
! TREE:
!   -
!
! COMMENTS:
!   In each line, TABs and spaces are skipped. Note that A in the file must
!   not contain empty lines (by default they are filled in with zeros).
!
! REFERENCES:
!   -
!===============================================================================
!
! DECLARATION OF VARIABLES:
    IMPLICIT NONE
!
! PARAMETERS:
    INTEGER, PARAMETER :: FNUM = 100 ! Input file number
!
! INPUT SCALAR ARGUMENTS:
    INTEGER, INTENT(IN) :: L, S, M, N
!
! INPUT ARRAY ARGUMENTS:
    CHARACTER(L), INTENT(IN) :: FNAME
!
! OUTPUT ARRAY ARGUMENTS:
    REAL*8, DIMENSION(M, N), INTENT(OUT) :: A
!
! LOCAL SCALARS:
    INTEGER &
        I, &   ! Loop index
        STATUS ! If end-of-file
!===============================================================================
!
    OPEN(FNUM, FILE = FNAME, ACTION = 'READ')
!   Skip S lines
    DO I = 1, S
        READ(FNUM, *)
    END DO
!   Read the data
    DO I = 1, M
        READ(FNUM, *, IOSTAT = STATUS) A(I, :)
        IF (STATUS /= 0) EXIT
    END DO ! I = 1, M
    CLOSE(FNUM)
!
END SUBROUTINE READ2D
!===============================================================================
! 02Oct14 - Cosmetic changes in declaration of variables according to SORD.
!
! 15Mar13 - Created and tested for the following txt file in the project
!           directory:
!
!           A test file. A string after the next one is empty.
!           The last line is  TABed.
!
!           1.0  2.0
!           3.0  4.0
!          -5.0  0.0
!           3.0  1.5
!	            19.0 -6.0
!
!           For the file NOT in the project directory (FNAME = FULL_PATH).
!           Note, the second line with numbers, 3.0 4.0, was deleted to
!           distinguish between the files:
!
!           A test file. A string after the next one is empty.
!           The last line is  TABed.
!
!           1.0  2.0
!          -5.0  0.0
!           3.0  1.5
!	            9.0 -6.0
!===============================================================================