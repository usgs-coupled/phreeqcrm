!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabase(FILENAME)
        CHARACTER(LEN=*) :: FILENAME
        INTEGER(KIND=4)  :: LoadDatabase
        INTERFACE 
          FUNCTION FLoadDatabase(FILENAME)
            !DEC$ ATTRIBUTES C,REFERENCE::FLoadDatabase
            !DEC$ ATTRIBUTES ALIAS:'_LoadDatabaseF'::FLoadDatabase
            CHARACTER(LEN=*) :: FILENAME
            INTEGER(KIND=4) :: FLoadDatabase
          END FUNCTION FLoadDatabase
        END INTERFACE
        LoadDatabase = FLoadDatabase(FILENAME)
      END FUNCTION LoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLastError
        INTERFACE 
          SUBROUTINE FOutputLastError
            !DEC$ ATTRIBUTES C,REFERENCE::FOutputLines
            !DEC$ ATTRIBUTES ALIAS:'_OutputLastError'::FOutputLastError
          END SUBROUTINE FOutputLastError
        END INTERFACE
        CALL FOutputLastError
      END SUBROUTINE OutputLastError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AccumulateLine(LINE)
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: AccumulateLine
        INTERFACE 
          FUNCTION FAccumulate(LINE)
            !DEC$ ATTRIBUTES C,REFERENCE::FAccumulate
            !DEC$ ATTRIBUTES ALIAS:'_AccumulateLineF'::FAccumulate
            CHARACTER(LEN=*) :: LINE
            INTEGER(KIND=4)  :: FAccumulate
          END FUNCTION FAccumulate
        END INTERFACE
        AccumulateLine = FAccumulate(LINE)
      END FUNCTION AccumulateLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION Run(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
        LOGICAL(KIND=4)  :: OUTPUT_ON
        LOGICAL(KIND=4)  :: ERROR_ON
        LOGICAL(KIND=4)  :: LOG_ON
        LOGICAL(KIND=4)  :: SELECTED_ON
        INTEGER(KIND=4)  :: Run
        INTERFACE 
          FUNCTION FRun(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
            !DEC$ ATTRIBUTES C,REFERENCE::FRun
            !DEC$ ATTRIBUTES ALIAS:'_RunF'::FRun
            LOGICAL(KIND=4) :: OUTPUT_ON
            LOGICAL(KIND=4) :: ERROR_ON
            LOGICAL(KIND=4) :: LOG_ON
            LOGICAL(KIND=4) :: SELECTED_ON
            INTEGER(KIND=4) :: FRun
          END FUNCTION FRun
        END INTERFACE
        Run = FRun(OUTPUT_ON, ERROR_ON, LOG_ON, SELECTED_ON)
      END FUNCTION Run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLines
        INTERFACE 
          SUBROUTINE FOutputLines
            !DEC$ ATTRIBUTES C,REFERENCE::FOutputLines
            !DEC$ ATTRIBUTES ALIAS:'_OutputLines'::FOutputLines
          END SUBROUTINE FOutputLines
        END INTERFACE
        CALL FOutputLines
      END SUBROUTINE OutputLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputRowCount
        INTEGER(KIND=4)  :: GetSelectedOutputRowCount
        INTERFACE 
          FUNCTION FRows
            !DEC$ ATTRIBUTES C,REFERENCE::FRows
            !DEC$ ATTRIBUTES ALIAS:'_GetSelectedOutputRowCount'::FRows
            INTEGER(KIND=4) :: FRows
          END FUNCTION FRows
        END INTERFACE
        GetSelectedOutputRowCount = FRows() - 1
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount
        INTEGER(KIND=4)  :: GetSelectedOutputColumnCount
        INTERFACE 
          FUNCTION FCols
          !DEC$ ATTRIBUTES C,REFERENCE::FCols
          !DEC$ ATTRIBUTES ALIAS:'_GetSelectedOutputColumnCount'::FCols
            INTEGER(KIND=4) :: FCols
          END FUNCTION FCols
        END INTERFACE
        GetSelectedOutputColumnCount = FCols()
      END FUNCTION GetSelectedOutputColumnCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputValue(ROW,COL,VTYPE,DVALUE,SVALUE)
        INTEGER(KIND=4)  :: ROW
        INTEGER(KIND=4)  :: COL
        INTEGER(KIND=4)  :: VTYPE
        REAL(KIND=8)     :: DVALUE
        CHARACTER(LEN=*) :: SVALUE
        INTEGER(KIND=4)  :: GetSelectedOutputValue
        INTEGER(KIND=4)  :: adjcol
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
      END FUNCTION GetSelectedOutputValue
