MODULE SUBROUTINES

    CONTAINS

    SUBROUTINE BOUNDARY_CONDITIONS(I, P, PN, NI, NJ, ORDER)

        IMPLICIT NONE
        INTEGER :: I, NI, NJ, ORDER
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: P
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: PN

        IF (ORDER.EQ.1) THEN
            PN(I, NJ) = P(I, NJ - 1)

        ELSE IF (ORDER.EQ.2) THEN
            PN(I, NJ) = (4.0 / 3.0) * P(I, NJ - 1) - (1.0 / 3.0) * P(I, NJ - 2)

        END IF

    END SUBROUTINE

    SUBROUTINE CALCULATE_ON_BREAK(I, J, P, PN, NI, NJ, ORDER, DX, DELTA, H2)

        IMPLICIT NONE
        INTEGER :: I, J, NI, NJ, ORDER
        DOUBLE PRECISION :: DX, H2, DELTA
        DOUBLE PRECISION :: RIGHT, LEFT       ! ДОПОЛНИТЕЛЬНЫЕ ПЕРЕМЕННЫЕ
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: PN
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: P

        LEFT = ((H2 + DELTA)**3) / DX
        RIGHT =  H2 / DX

        IF (ORDER == 1) THEN
            PN(I, J) = (RIGHT * P(I + 1, J) + LEFT * P(I - 1, J) + DELTA) / (LEFT + RIGHT)

        ELSE IF (ORDER == 2) THEN
            PN(I, J) = ((RIGHT / 2) * (4.0 * P(I + 1,J) - P(I + 2, J)) - (LEFT / 2) * (-4.0 * P(I - 1, J) + &
            P(I - 2, J)) + DELTA) / ((3.0 / 2) * (LEFT + RIGHT))

        END IF

    END SUBROUTINE

    SUBROUTINE FORCE_CALCULATE(DX, DY, P, F, NI, NJ)

        IMPLICIT NONE
        INTEGER :: I, J, NI, NJ
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: P
        DOUBLE PRECISION :: DX, DY, DS, F, DF

        OPEN (4, FILE='force_residual.txt')

        F = 0
        DS = DX * DY

        DO I = 1, NI - 1
            DO J = 1, NJ - 1

            DF = DS * (P(I, J) + P(I, J + 1) + P(I + 1, J) + P(I + 1, J + 1)) / 4.0
            F = F + DF

            ENDDO
        ENDDO

        WRITE(4,*) F

    END SUBROUTINE FORCE_CALCULATE

    SUBROUTINE OUTPUT_PLT(X, Y, P, H, NI, NJ)

        IMPLICIT NONE
        INTEGER :: I, J, NI, NJ
        DOUBLE PRECISION, DIMENSION(NI) :: X
        DOUBLE PRECISION, DIMENSION(NJ) :: Y
        DOUBLE PRECISION, DIMENSION(NI) :: H
        DOUBLE PRECISION, DIMENSION(NI, NJ) :: P

        OPEN (3, FILE='field.plt')
        WRITE(3,*) 'ZONE I=', NI
        WRITE(3,*) 'J=', NJ

        DO J = 1, NJ
            DO I = 1, NI
                WRITE (3, *) X(I), Y(J), P(I,J), H(I)
            END DO
        END DO

        CLOSE (3)

    END SUBROUTINE OUTPUT_PLT

END MODULE SUBROUTINES