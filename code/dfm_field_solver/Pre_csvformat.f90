! SPDX-FileCopyrightText: 2026 Qiang Sun
! SPDX-License-Identifier: BSD-3-Clause

! Minimal delimited-text writers for release 0.1.0.
! This independent implementation retains only the three entry points used by
! the manuscript drivers.
MODULE Pre_csvformat
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: csv_write_integer,csv_write_dble,csv_write_char

CONTAINS

    SUBROUTINE csv_write_integer(lun,value,advance,separator)
        INTEGER,INTENT(IN) :: lun,value
        LOGICAL,INTENT(IN) :: advance
        CHARACTER(LEN=3),INTENT(IN) :: separator
        CHARACTER(LEN=40) :: buffer

        WRITE(buffer,'(I0)') value
        CALL write_token(lun,TRIM(buffer),advance,separator)
    END SUBROUTINE csv_write_integer

    SUBROUTINE csv_write_dble(lun,value,advance,separator)
        INTEGER,INTENT(IN) :: lun
        DOUBLE PRECISION,INTENT(IN) :: value
        LOGICAL,INTENT(IN) :: advance
        CHARACTER(LEN=3),INTENT(IN) :: separator
        CHARACTER(LEN=40) :: buffer

        WRITE(buffer,'(G20.12)') value
        CALL write_token(lun,TRIM(ADJUSTL(buffer)),advance,separator)
    END SUBROUTINE csv_write_dble

    SUBROUTINE csv_write_char(lun,value,advance,separator)
        INTEGER,INTENT(IN) :: lun
        CHARACTER(LEN=*),INTENT(IN) :: value
        LOGICAL,INTENT(IN) :: advance
        CHARACTER(LEN=3),INTENT(IN) :: separator
        CHARACTER(LEN=:),ALLOCATABLE :: token

        IF (separator == 'cmr') THEN
            token='"'//TRIM(value)//'"'
        ELSE
            token=TRIM(value)
        END IF
        CALL write_token(lun,token,advance,separator)
    END SUBROUTINE csv_write_char

    SUBROUTINE write_token(lun,token,advance,separator)
        INTEGER,INTENT(IN) :: lun
        CHARACTER(LEN=*),INTENT(IN) :: token
        LOGICAL,INTENT(IN) :: advance
        CHARACTER(LEN=3),INTENT(IN) :: separator

        IF (advance) THEN
            WRITE(lun,'(A)') TRIM(token)
        ELSE IF (separator == 'cmr') THEN
            WRITE(lun,'(A)',ADVANCE='NO') TRIM(token)//', '
        ELSE IF (separator == 'spc') THEN
            WRITE(lun,'(A)',ADVANCE='NO') TRIM(token)//' '
        ELSE
            WRITE(*,'(A,A)') 'Unsupported output separator: ',separator
            STOP 82
        END IF
    END SUBROUTINE write_token

END MODULE Pre_csvformat
