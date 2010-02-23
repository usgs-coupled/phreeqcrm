!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabase(FILENAME)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: FILENAME
        INTEGER          :: LoadDatabase
#if defined(_WIN32)
        INTERFACE 
          FUNCTION LoadDatabaseF(FILENAME)
            !DEC$ ATTRIBUTES C,REFERENCE::LoadDatabaseF
            !DEC$ ATTRIBUTES ALIAS:'_LoadDatabaseF'::LoadDatabaseF
            CHARACTER(LEN=*) :: FILENAME
            INTEGER(KIND=4) :: LoadDatabaseF
          END FUNCTION LoadDatabaseF
        END INTERFACE
#else
        INTEGER          :: LoadDatabaseF
#endif

        LoadDatabase = LoadDatabaseF(FILENAME)
      END FUNCTION LoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLastError
        IMPLICIT NONE
#if defined(_WIN32)
        INTERFACE 
          SUBROUTINE OutputLastErrorF
            !DEC$ ATTRIBUTES C,REFERENCE::OutputLastErrorF
            !DEC$ ATTRIBUTES ALIAS:'_OutputLastError'::OutputLastErrorF
          END SUBROUTINE OutputLastErrorF
        END INTERFACE
#endif
        CALL OutputLastErrorF
      END SUBROUTINE OutputLastError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AccumulateLine(LINE)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: LINE
        INTEGER          :: AccumulateLine
#if defined(_WIN32)
        INTERFACE 
          FUNCTION AccumulateLineF(LINE)
            !DEC$ ATTRIBUTES C,REFERENCE::AccumulateLineF
            !DEC$ ATTRIBUTES ALIAS:'_AccumulateLineF'::AccumulateLineF
            CHARACTER(LEN=*) :: LINE
            INTEGER(KIND=4)  :: AccumulateLineF
          END FUNCTION AccumulateLineF
        END INTERFACE
#else
        INTEGER          :: AccumulateLineF
#endif
        AccumulateLine = AccumulateLineF(LINE)
      END FUNCTION AccumulateLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION Run(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
        IMPLICIT NONE
        LOGICAL  :: OUTPUT_ON
        LOGICAL  :: ERROR_ON
        LOGICAL  :: LOG_ON
        LOGICAL  :: SELECTED_ON
        INTEGER  :: Run
#if defined(_WIN32)
        INTERFACE 
          FUNCTION RunF(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
            !DEC$ ATTRIBUTES C,REFERENCE::RunF
            !DEC$ ATTRIBUTES ALIAS:'_RunF'::RunF
            LOGICAL(KIND=4) :: OUTPUT_ON
            LOGICAL(KIND=4) :: ERROR_ON
            LOGICAL(KIND=4) :: LOG_ON
            LOGICAL(KIND=4) :: SELECTED_ON
            INTEGER(KIND=4) :: RunF
          END FUNCTION RunF
        END INTERFACE
#else
        INTEGER  :: RunF
#endif
        Run = RunF(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
      END FUNCTION Run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunFile(FILENAME,OUTPUT_ON,ERROR_ON,LOG_ON,SELECTED_ON)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: FILENAME
        LOGICAL          :: OUTPUT_ON
        LOGICAL          :: ERROR_ON
        LOGICAL          :: LOG_ON
        LOGICAL          :: SELECTED_ON
        INTEGER          :: RunFile
#if defined(_WIN32)
        INTERFACE 
          FUNCTION RunFileF(OUT_ON, ERR_ON, LOG_ON, SEL_ON, FILE)
            !DEC$ ATTRIBUTES C,REFERENCE::RunFileF
            !DEC$ ATTRIBUTES ALIAS:'_RunFileF'::RunFileF
            LOGICAL(KIND=4)  :: OUT_ON
            LOGICAL(KIND=4)  :: ERR_ON
            LOGICAL(KIND=4)  :: LOG_ON
            LOGICAL(KIND=4)  :: SEL_ON
            CHARACTER(LEN=*) :: FILE
            INTEGER(KIND=4)  :: RunFileF
          END FUNCTION RunFileF
        END INTERFACE
#else
        INTEGER          :: RunFileF
#endif
        RunFile = RunFileF(OUTPUT_ON, ERROR_ON, LOG_ON,
     &                     SELECTED_ON, FILENAME)
      END FUNCTION RunFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLines
        IMPLICIT NONE
#if defined(_WIN32)
        INTERFACE 
          SUBROUTINE OutputLinesF
            !DEC$ ATTRIBUTES C,REFERENCE::OutputLinesF
            !DEC$ ATTRIBUTES ALIAS:'_OutputLines'::OutputLinesF
          END SUBROUTINE OutputLinesF
        END INTERFACE
#endif
        CALL OutputLinesF
      END SUBROUTINE OutputLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputRowCount()
        IMPLICIT NONE
        INTEGER  :: GetSelectedOutputRowCount
#if defined(_WIN32)
        INTERFACE 
          FUNCTION FRows
            !DEC$ ATTRIBUTES C,REFERENCE::FRows
            !DEC$ ATTRIBUTES ALIAS:'_GetSelectedOutputRowCountF'::FRows
            INTEGER(KIND=4) :: FRows
          END FUNCTION FRows
        END INTERFACE
        GetSelectedOutputRowCount = FRows() - 1
#else
        INTEGER  :: GetSelectedOutputRowCountF
        GetSelectedOutputRowCount = GetSelectedOutputRowCountF() - 1
#endif
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount()
        IMPLICIT NONE
        INTEGER :: GetSelectedOutputColumnCount
#if defined(_WIN32)
        INTERFACE 
          FUNCTION FCols
          !DEC$ ATTRIBUTES C,REFERENCE::FCols
          !DEC$ ATTRIBUTES ALIAS:'_GetSelectedOutputColumnCount'::FCols
            INTEGER(KIND=4) :: FCols
          END FUNCTION FCols
        END INTERFACE
        GetSelectedOutputColumnCount = FCols()
#else
        INTEGER :: GetSelectedOutputColumnCountF
        GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF()
#endif
      END FUNCTION GetSelectedOutputColumnCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputValue(ROW,COL,VTYPE,DVALUE,SVALUE)
        IMPLICIT NONE
        INTEGER          :: ROW
        INTEGER          :: COL
        INTEGER          :: VTYPE
        REAL*8           :: DVALUE
        CHARACTER(LEN=*) :: SVALUE
        INTEGER          :: GetSelectedOutputValue
        INTEGER          :: adjcol
#if defined(_WIN32)
        INTERFACE 
          FUNCTION Get(ROW,COL,VTYPE,DVALUE,SVALUE)
            !DEC$ ATTRIBUTES C,REFERENCE::Get
            !DEC$ ATTRIBUTES ALIAS:'_GetSelectedOutputValueF'::Get
            INTEGER(KIND=4)  :: ROW
            INTEGER(KIND=4)  :: COL
            INTEGER(KIND=4)  :: VTYPE
            REAL(KIND=8)     :: DVALUE
            CHARACTER(LEN=*) :: SVALUE
            INTEGER(KIND=4)  :: Get
          END FUNCTION Get
        END INTERFACE
        adjcol = col - 1
        GetSelectedOutputValue = Get(ROW,adjcol,VTYPE,DVALUE,SVALUE)
#else
        INTEGER          :: GetSelectedOutputValueF
        adjcol = col - 1
        GetSelectedOutputValue = GetSelectedOutputValueF(ROW,
     &                     adjcol,VTYPE,DVALUE,SVALUE)
#endif
      END FUNCTION GetSelectedOutputValue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION System(COMMAND)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: COMMAND
        INTEGER          :: System
        INTEGER          :: SystemF
        System = SystemF(COMMAND)
      END FUNCTION System
