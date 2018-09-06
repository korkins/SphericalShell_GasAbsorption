SUBROUTINE PRINT2D(A, L, M, FNAME, N)
!===============================================================================
! PURPOSE:
!   To print all elements of a 2D array, A, into a text file.
!
! INPUT:
!   A       D(L, M)   Input array
!   L       I(1)      Number of lines in A
!   M       I(1)      Number of columns in A 
!   FNAME   C(N)      Output file name
!   N       I(1)      Length of FNAME
!
! OUTPUT:
!   -
!
! TREE:
!   - 
!
! COMMENTS:
!   Format is fixed by the subroutine. Note 1000 in the FORMAT line.
!
! REFERENCES:
!   -     
!===============================================================================  
! 
! DECLARATION OF VARIABLES
    IMPLICIT NONE
!
! PARAMETERS
    INTEGER, PARAMETER :: FNUM = 100 ! Output file number      
!
! INPUT SCALAR ARGUMENTS
    INTEGER, INTENT(IN) :: L, M, N   
!
! INPUT ARRAY ARGUMENTS
    CHARACTER(N), INTENT(IN) :: FNAME
    REAL*8, DIMENSION(L, M), INTENT(IN) :: A
!
! LOCAL SCALARS
    INTEGER &
        I ! Loop index
!===============================================================================
!
    OPEN(FNUM, FILE = FNAME)
    DO I = 1, L   
        WRITE(FNUM, 10) A(I, :)
    END DO ! I = 1, L
    CLOSE(FNUM)
!
!10  FORMAT(100ES24.12)
!10  FORMAT(4(3X, F8.5), 2(3X, F8.5))
!10  FORMAT(1000E20.8)
!10  FORMAT(6(2X, ES15.8))
10  FORMAT(1X, 1F5.2, 3X, 1F6.2, 3X, 1F5.2, 3X, 1F8.2, 2F12.6)
!	
END SUBROUTINE PRINT2D
!===============================================================================
! 02Oct14 - Cosmetic changes in declaration of variables according to SORD.
!
! 21Dec12 - First created. Tested for 11-by-44 array.
!===============================================================================