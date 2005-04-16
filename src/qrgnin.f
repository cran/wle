      BLOCK DATA added
      COMMON/XYINIT/qinit
      LOGICAL qinit
      DATA qinit/.FALSE./
      END

      LOGICAL FUNCTION qrgnin()
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNIN()
C               Q Random GeNerators INitialized?
C
C     A trivial routine to determine whether or not the random
C     number generator has been initialized.  Returns .TRUE. if
C     it has, else .FALSE.
C
C**********************************************************************
C     .. Local Scalars ..
      LOGICAL qinit
      COMMON/XYINIT/qinit
C     ..
      qrgnin = qinit
      RETURN
      END

      LOGICAL FUNCTION qrgnsn(qvalue)
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNSN( QVALUE )
C               Q Random GeNerators Set whether iNitialized
C
C     Sets state of whether random number generator is initialized
C     to QVALUE.
C     It returns the (meaningless) value .TRUE.
C
C**********************************************************************
      LOGICAL qvalue
      LOGICAL qinit
      COMMON/XYINIT/qinit

      qinit = qvalue
      qrgnsn = .TRUE.
      RETURN

      END
