! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE VALIDATE_ATLAS_MOD
  USE PARKIND1, ONLY: JPIM, JPRB
  USE CLOUDSC_MPI_MOD

  USE ATLAS_MODULE
  USE ATLAS_FIELDSET_MODULE
  use atlas_functionspace_blockstructuredcolumns_module
  USE, INTRINSIC :: ISO_C_BINDING
  USE EXPAND_MOD, ONLY: LOAD_AND_EXPAND
  USE FILE_IO_MOD, ONLY: INPUT_INITIALIZE, INPUT_FINALIZE

  IMPLICIT NONE

CONTAINS

  SUBROUTINE VALIDATESTATE_ATLAS(FSET, NAME, NLON, NGPTOTG)
    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    CHARACTER(*), INTENT(IN) :: NAME
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG

    CALL VALIDATEVAR_ATLAS(FSET, NAME, NLON, NGPTOTG, "A")
    CALL VALIDATEVAR_ATLAS(FSET, NAME, NLON, NGPTOTG, "Q")
    CALL VALIDATEVAR_ATLAS(FSET, NAME, NLON, NGPTOTG, "T")
    CALL VALIDATEVAR_ATLAS(FSET, NAME, NLON, NGPTOTG, "CLD")
  END SUBROUTINE VALIDATESTATE_ATLAS

  SUBROUTINE VALIDATEVAR_ATLAS(FSET, NAME, NLON, NGPTOTG, STATE_VAR)
    ! Computes and prints errors "in the L1 norm sense"
    TYPE(ATLAS_FIELDSET), INTENT(INOUT) :: FSET
    CHARACTER(*), INTENT(IN) :: NAME
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG
    CHARACTER(*), INTENT(IN), OPTIONAL :: STATE_VAR

    REAL(KIND=JPRB), ALLOCATABLE :: REF_R2(:,:), REF_R3(:,:,:), REF_R4(:,:,:,:)
    REAL(C_DOUBLE), POINTER :: FIELD_R1(:,:), FIELD_R2(:,:,:), FIELD_R3(:,:,:,:)
    TYPE(ATLAS_FUNCTIONSPACE_BLOCKSTRUCTUREDCOLUMNS) :: FSPACE
    TYPE(ATLAS_FIELD) :: FIELD
    INTEGER :: B, BSIZE, JL, JK, JM
    REAL(KIND=JPRB) :: ZMINVAL(1), ZMAX_VAL_ERR(2), ZDIFF, ZSUM_ERR_ABS(2), ZRELERR, ZAVGPGP
    INTEGER :: FRANK, NBLOCKS, NLEV, NGPTOT, NPROMA, VAR_ID, NDIM
    CHARACTER(LEN=256) :: FULLNAME

    IF (PRESENT(STATE_VAR)) THEN
      FULLNAME = NAME//'_'//STATE_VAR
    ELSE
      FULLNAME = NAME
    ENDIF

    FIELD = FSET%FIELD(NAME)
    FRANK = FIELD%RANK()
    FSPACE = FIELD%FUNCTIONSPACE()
    NLEV = FIELD%LEVELS()
    NGPTOT = FSPACE%SIZE()
    NBLOCKS = FSPACE%NBLKS()
    NPROMA = FSPACE%BLOCK_SIZE(1)

    ZMINVAL(1) = +HUGE(ZMINVAL(1))
    ZMAX_VAL_ERR(1) = -HUGE(ZMAX_VAL_ERR(1))
    ZMAX_VAL_ERR(2) = 0.0_JPRB
    ZSUM_ERR_ABS(:) = 0.0_JPRB

    CALL INPUT_INITIALIZE(NAME='reference')
    IF (FRANK == 2) THEN
        CALL LOAD_AND_EXPAND(NAME, REF_R2, NLON, NPROMA, NGPTOT, NBLOCKS, NGPTOTG)
        CALL FIELD%DATA(FIELD_R1)
        !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
        !& REDUCTION(MIN:ZMINVAL, MAX:ZMAX_VAL_ERR, +:ZSUM_ERR_ABS)
        DO B=1, NBLOCKS
          BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
          ZMINVAL(1) = MIN(ZMINVAL(1),MINVAL(FIELD_R1(:,B)))
          ZMAX_VAL_ERR(1) = MAX(ZMAX_VAL_ERR(1),MAXVAL(FIELD_R1(:,B)))
          DO JK=1, bsize
            ! Difference against reference result in one-norm sense
            ZDIFF = ABS(FIELD_R1(JK,B) - REF_R2(JK,B))
            ZMAX_VAL_ERR(2) = MAX(ZMAX_VAL_ERR(2),ZDIFF)
            ZSUM_ERR_ABS(1) = ZSUM_ERR_ABS(1) + ZDIFF
            ZSUM_ERR_ABS(2) = ZSUM_ERR_ABS(2) + ABS(REF_R2(JK,B))
          ENDDO
        END DO
    ELSE IF (FRANK == 3) THEN
        CALL LOAD_AND_EXPAND(NAME, REF_R3, NLON, FIELD%LEVELS(), NPROMA, NGPTOT, NBLOCKS, NGPTOTG)
        CALL FIELD%DATA(FIELD_R2)
        !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
        !& REDUCTION(MIN:ZMINVAL, MAX:ZMAX_VAL_ERR, +:ZSUM_ERR_ABS)
        DO B=1, NBLOCKS
          BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size !! TODO ad other loops
          ZMINVAL(1) = MIN(ZMINVAL(1),MINVAL(FIELD_R2(:,:,B)))
          ZMAX_VAL_ERR(1) = MAX(ZMAX_VAL_ERR(1),MAXVAL(FIELD_R2(:,:,B)))
          DO JL=1, NLEV
            DO JK=1, bsize
              ! Difference against reference result in one-norm sense
              ZDIFF = ABS(FIELD_R2(JK,JL,B) - REF_R3(JK,JL,B))
              ZMAX_VAL_ERR(2) = MAX(ZMAX_VAL_ERR(2),ZDIFF)
              ZSUM_ERR_ABS(1) = ZSUM_ERR_ABS(1) + ZDIFF
              ZSUM_ERR_ABS(2) = ZSUM_ERR_ABS(2) + ABS(REF_R3(JK,JL,B))
            ENDDO
          END DO
        END DO
      ELSE IF (FRANK == 4 .AND. PRESENT(STATE_VAR)) THEN
        CALL FIELD%DATA(FIELD_R3)
        NDIM = FIELD%SHAPE(3) - 3
        IF (STATE_VAR /= 'CLD') THEN
            VAR_ID = 1
            IF (STATE_VAR == 'A') THEN
                VAR_ID = 2
            ENDIF
            IF (STATE_VAR == 'Q') THEN
                VAR_ID = 3
            ENDIF
            CALL LOAD_AND_EXPAND(NAME//'_'//STATE_VAR, REF_R3, NLON, NLEV, NPROMA, NGPTOT, NBLOCKS, NGPTOTG)
            !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
            !& REDUCTION(MIN:ZMINVAL, MAX:ZMAX_VAL_ERR, +:ZSUM_ERR_ABS)
            DO B=1, NBLOCKS
              BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
              ZMINVAL(1) = MIN(ZMINVAL(1),MINVAL(FIELD_R3(:,:,VAR_ID,B)))
              ZMAX_VAL_ERR(1) = MAX(ZMAX_VAL_ERR(1),MAXVAL(FIELD_R3(:,:,VAR_ID,B)))
              DO JL=1, NLEV
                DO JK=1, bsize
                  ! Difference against reference result in one-norm sense
                  ZDIFF = ABS(FIELD_R3(JK,JL,VAR_ID,B) - REF_R3(JK,JL,B))
                  ZMAX_VAL_ERR(2) = MAX(ZMAX_VAL_ERR(2),ZDIFF)
                  ZSUM_ERR_ABS(1) = ZSUM_ERR_ABS(1) + ZDIFF
                  ZSUM_ERR_ABS(2) = ZSUM_ERR_ABS(2) + ABS(REF_R3(JK,JL,B))
                ENDDO
              END DO
            END DO
        ELSE IF (STATE_VAR == 'CLD') THEN
          CALL LOAD_AND_EXPAND(NAME//'_CLD', REF_R4, NLON, NLEV, NDIM, NPROMA, NGPTOT, NBLOCKS, NGPTOTG)
          !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
          !& REDUCTION(MIN:ZMINVAL, MAX:ZMAX_VAL_ERR, +:ZSUM_ERR_ABS)
          DO B=1, NBLOCKS
            BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
            ZMINVAL(1) = MIN(ZMINVAL(1),MINVAL(FIELD_R3(:,:,:,B)))
            ZMAX_VAL_ERR(1) = MAX(ZMAX_VAL_ERR(1),MAXVAL(FIELD_R3(:,:,:,B)))
            DO JM=1, NDIM
              DO JL=1, NLEV
                DO JK=1, BSIZE
                  ! Difference against reference result in one-norm sense
                  ZDIFF = ABS(FIELD_R3(JK,JL,JM,B) - REF_R4(JK,JL,JM,B))
                  ZMAX_VAL_ERR(2) = MAX(ZMAX_VAL_ERR(2),ZDIFF)
                  ZSUM_ERR_ABS(1) = ZSUM_ERR_ABS(1) + ZDIFF
                  ZSUM_ERR_ABS(2) = ZSUM_ERR_ABS(2) + ABS(REF_R4(JK,JL,JM,B))
                ENDDO
              ENDDO
            END DO
          END DO
        ENDIF
      ELSE
        PRINT *, "FIELD RANK NOT SUPPORTED"
        CALL EXIT(1)
    ENDIF
    CALL INPUT_FINALIZE()

    CALL CLOUDSC_MPI_REDUCE_MIN(ZMINVAL, 1, 0)
    CALL CLOUDSC_MPI_REDUCE_MAX(ZMAX_VAL_ERR, 2, 0)
    CALL CLOUDSC_MPI_REDUCE_SUM(ZSUM_ERR_ABS, 2, 0)

    IF (PRESENT(NGPTOTG)) THEN
      ZAVGPGP = ZSUM_ERR_ABS(1) / REAL(NGPTOTG,JPRB)
    ELSE
      ZAVGPGP = ZSUM_ERR_ABS(1) / REAL(NGPTOT,JPRB)
    END IF

    IF (IRANK == 0) THEN
      CALL ERROR_PRINT(FULLNAME, ZMINVAL(1), ZMAX_VAL_ERR(1), ZMAX_VAL_ERR(2), &
        & ZSUM_ERR_ABS(1), ZSUM_ERR_ABS(2), ZAVGPGP, NDIM=FRANK-1)
    END IF
  END SUBROUTINE VALIDATEVAR_ATLAS

  SUBROUTINE VALIDATE(NAME, REF, FIELD, NLON, NLEV, NDIM, NGPTOT, NBLOCKS, NGPTOTG)
    ! Computes and prints errors "in the L1 norm sense"
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(INOUT) :: REF(:,:,:,:), FIELD(:,:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN) :: NLON, NLEV, NDIM, NBLOCKS, NGPTOT
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NGPTOTG

    INTEGER :: B, BSIZE, JL, JK, JM
    REAL(KIND=JPRB) :: ZMINVAL(1), ZMAX_VAL_ERR(2), ZDIFF, ZSUM_ERR_ABS(2), ZRELERR, ZAVGPGP

    ZMINVAL(1) = +HUGE(ZMINVAL(1))
    ZMAX_VAL_ERR(1) = -HUGE(ZMAX_VAL_ERR(1))
    ZMAX_VAL_ERR(2) = 0.0_JPRB
    ZSUM_ERR_ABS(:) = 0.0_JPRB

    !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B, BSIZE) &
    !& REDUCTION(MIN:ZMINVAL, MAX:ZMAX_VAL_ERR, +:ZSUM_ERR_ABS)
    DO B=1, NBLOCKS
      BSIZE = MIN(NLON, NGPTOT - (B-1)*NLON)  ! Field block size
      ZMINVAL(1) = MIN(ZMINVAL(1),MINVAL(FIELD(:,:,:,B)))
      ZMAX_VAL_ERR(1) = MAX(ZMAX_VAL_ERR(1),MAXVAL(FIELD(:,:,:,B)))
      DO JM=1, NDIM
        DO JL=1, NLEV
          DO JK=1, bsize
            ! Difference against reference result in one-norm sense
            ZDIFF = ABS(FIELD(JK,JL,JM,B) - REF(JK,JL,JM,B))
            ZMAX_VAL_ERR(2) = MAX(ZMAX_VAL_ERR(2),ZDIFF)
            ZSUM_ERR_ABS(1) = ZSUM_ERR_ABS(1) + ZDIFF
            ZSUM_ERR_ABS(2) = ZSUM_ERR_ABS(2) + ABS(REF(JK,JL,JM,B))
          END DO
        END DO
      END DO
    END DO

    CALL CLOUDSC_MPI_REDUCE_MIN(ZMINVAL, 1, 0)
    CALL CLOUDSC_MPI_REDUCE_MAX(ZMAX_VAL_ERR, 2, 0)
    CALL CLOUDSC_MPI_REDUCE_SUM(ZSUM_ERR_ABS, 2, 0)

    IF (PRESENT(NGPTOTG)) THEN
      ZAVGPGP = ZSUM_ERR_ABS(1) / REAL(NGPTOTG,JPRB)
    ELSE
      ZAVGPGP = ZSUM_ERR_ABS(1) / REAL(NGPTOT,JPRB)
    END IF

    IF (IRANK == 0) THEN
      CALL ERROR_PRINT(NAME, ZMINVAL(1), ZMAX_VAL_ERR(1), ZMAX_VAL_ERR(2), &
        & ZSUM_ERR_ABS(1), ZSUM_ERR_ABS(2), ZAVGPGP, NDIM=3)
    END IF
  END SUBROUTINE VALIDATE

  SUBROUTINE ERROR_PRINT(NAME, ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP, NDIM)
    ! Print error statistic for a single variable (adapted from diff_mod.F90)
    CHARACTER(*), INTENT(IN) :: NAME
    REAL(KIND=JPRB), INTENT(IN) :: ZMINVAL, ZMAXVAL, ZMAXERR, ZERRSUM, ZSUM, ZAVGPGP
    INTEGER(KIND=JPIM), INTENT(IN) :: NDIM
    REAL(KIND=JPRB) :: zrelerr
    REAL(KIND=JPRB), parameter :: zeps = epsilon(1.0_JPRB)
    INTEGER :: IOPT
    character(len=5) clwarn

    iopt = 0
    if (zerrsum < zeps) then
      zrelerr = 0.0_JPRB
      iopt = 1
    elseif (zsum < zeps) then
      zrelerr = zerrsum/(1.0_JPRB + zsum)
      iopt = 2
    else
      zrelerr = zerrsum/zsum
      iopt = 3
    endif

    !-- If you get 4 exclamation marks next to your error output,
    !   then it is likely that some uninitialized variables exists or
    !   some other screw-up -- watch out this !!!!
    clwarn = ' '
    if (zrelerr > 10.0_JPRB * zeps) clwarn = ' !!!!'
    zrelerr = 100.0_JPRB * zrelerr

    write(*,1000) name,ndim,iopt, &
     & zminval,zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn
1000 format(1X,A20,1X,I1,'D',I1,5(1X,E20.13),A)

  END SUBROUTINE ERROR_PRINT

END MODULE VALIDATE_ATLAS_MOD
