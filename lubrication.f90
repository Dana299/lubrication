PROGRAM LUBRICATION

    USE SUBROUTINES

    IMPLICIT NONE

    INTEGER :: NI, NJ, ORDER        ! ЧИСЛО УЗЛОВ ПО ОСИ X И Y, ПОРЯДОК АППРОКСИМАЦИИ ПРОИЗВОДНЫХ В ГУ
    INTEGER :: S_MAX, S             ! МАКСИМАЛЬНОЕ ЧИСЛО ИТЕРАЦИЙ, ТЕКУЩАЯ ИТЕРАЦИЯ
    INTEGER :: NODE_L1              ! УЗЕЛ, В КОТОРОМ НАХОДИТСЯ СТУПЕНЬКА
    INTEGER :: I, J

    DOUBLE PRECISION :: H1, H2, L1, W, DELTA   ! ГЕОМЕТРИЧЕСКИЕ ХАРАКТЕРИСТИКИ ПОДПЯТНИКА
    DOUBLE PRECISION :: DX, DY, EPS            ! РАЗМЕРЫ ШАГОВ ПО КООРДИНАТЕ, ТОЧНОСТЬ РЕШЕНИЯ
    DOUBLE PRECISION :: D, DP, DP0             ! ПЕРЕМЕННЫЕ ДЛЯ ПРОВЕРКИ КРИТЕРИЯ СХОДИМОСТИ
    DOUBLE PRECISION :: NUMERATOR, DENOMINATOR ! ВСПОМОГ. ПЕРЕМЕННЫЕ ДЛЯ ИТЕР. ПРОЦЕССА
    DOUBLE PRECISION :: F                     ! СИЛА

    DOUBLE PRECISION :: LEFT, RIGHT
    DOUBLE PRECISION :: H_HALF_1, H_HALF_2

    DOUBLE PRECISION, ALLOCATABLE :: P(:,:), PN(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
    DOUBLE PRECISION, ALLOCATABLE :: H(:)


    OPEN (1, FILE='input.txt')

        READ (1,*) H1
        READ (1,*) H2
        READ (1,*) W
        READ (1,*) L1
        READ (1,*) DELTA
        READ (1,*) NI
        READ (1,*) NJ
        READ (1,*) S_MAX
        READ (1,*) EPS
        READ (1,*) ORDER

    CLOSE (1)

    ALLOCATE(X(1:NI))
    ALLOCATE(Y(1:NJ))
    ALLOCATE(H(1:NI))
    ALLOCATE(P(1:NI,1:NJ))
    ALLOCATE(PN(1:NI,1:NJ))

    ! ШАГИ ПО КООРДИНАТАМ
    DX = 1.0 / (NI - 1)
    DY = W / 2 / (NJ - 1)

    NODE_L1 = INT(L1 / DX) + 1  ! УЗЕЛ, В КОТОРОМ НАХОДИТСЯ "СТУПЕНЬКА"

    ! СЕТКА
    DO I = 1, NI
        X(I) = (I - 1) * DX
    ENDDO

    DO J = 1, NJ
        Y(J) = (J - 1) * DY
    ENDDO

    ! ЗАПОЛНЕНИЕ МАССИВА БЕЗРАЗМЕРНЫХ ТОЛЩИН СМАЗОЧНОГО СЛОЯ
    DO I = 1, NI

        IF (X(I).LE.L1) THEN
            H(I) = H1 - (X(I) / L1) * (H1 - H2 - DELTA)

        ELSE IF (X(I).GT.L1) THEN
            H(I) = H2
        END IF

    END DO

    P(:,:) = 0.
    PN(:,:) = 0.

    DP = 1
    DP0 = 1

    ! ИТЕРАЦИОННЫЙ ПРОЦЕСС
    S = 0

    OPEN (2, FILE='residual.plt')

    DO WHILE ((S.LT.S_MAX).AND.((DP / DP0).GT.EPS))

        S = S + 1

        DP = 0.

        DO I = 2, NI - 1
            DO J = 2, NJ - 1

                IF (I.EQ.NODE_L1) THEN         ! РАСЧЕТ НА РАЗРЫВЕ

                    CALL CALCULATE_ON_BREAK(I, J, P, PN, NI, NJ, ORDER, DX, DELTA, H2)

                ELSE
                    IF (X(I).LE.L1) THEN
                        H_HALF_1 = (H_NONDIM(X(I), DELTA, H1, H2, L1) + &
                        H_NONDIM(X(I + 1), DELTA, H1, H2, L1)) / 2

                        H_HALF_2 = (H_NONDIM(X(I), DELTA, H1, H2, L1) + &
                        H_NONDIM(X(I - 1), DELTA, H1, H2, L1)) / 2

                    ELSE IF (X(I).GT.L1) THEN
                        H_HALF_1 = H2
                        H_HALF_2 = (H_NONDIM(X(I), DELTA, H1, H2, L1) + &
                        H_NONDIM(X(I - 1), DELTA, H1, H2, L1)) / 2
                    END IF


                    NUMERATOR = H_HALF_1**3 * P(I + 1, J) / DX**2 + H_HALF_2**3 * P(I - 1, J) / DX**2 + &
                    H(I)**3 * (P(I, J + 1) + P(I, J - 1)) / DY**2 - (H(I + 1) - H(I - 1)) / (2 * DX)

                    DENOMINATOR = (H_HALF_1**3) / DX**2 + (H_HALF_2**3) / DX**2 + &
                    2 * (H(I))**3 / DY**2

                    PN(I, J) = NUMERATOR / DENOMINATOR

                    D = ABS(PN(I, J) - P(I,J))
                    DP = MAX(DP, D)

				END IF

            ENDDO

            CALL BOUNDARY_CONDITIONS(I, P, PN, NI, NJ, ORDER)

            D = ABS(PN(I, J) - P(I,J))
            DP = MAX(DP, D)

        ENDDO

        P(:,:) = PN(:,:)

        CALL FORCE_CALCULATE(DX, DY, P, F, NI, NJ)

        IF (S == 1) THEN
            DP0 = DP
        END IF

        WRITE(2,*) S, DP/DP0, F
        WRITE(*,*) 'Iteration ', S, DP/DP0

    ENDDO

    CLOSE(2)

    CALL OUTPUT_PLT(X, Y, P, NI, NJ, DELTA, H1, H2, L1)

    WRITE(*,*) 'Force = ', F
    WRITE(*,*) 'dx = ', dx
    WRITE(*,*) 'dy = ', dy

END PROGRAM



