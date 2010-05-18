!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AccumulateLine(ID,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: LINE
        INTEGER(KIND=4)  :: AccumulateLine
        INTEGER(KIND=4)  :: AccumulateLineF
        AccumulateLine = AccumulateLineF(ID,LINE)
      END FUNCTION AccumulateLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ClearAccumulatedLines(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: ClearAccumulatedLines
        INTEGER(KIND=4)  :: ClearAccumulatedLinesF
        ClearAccumulatedLines = ClearAccumulatedLinesF(ID)
      END FUNCTION ClearAccumulatedLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION CreateIPhreeqc()
        IMPLICIT NONE
        INTEGER(KIND=4)  :: CreateIPhreeqc
        INTEGER(KIND=4)  :: CreateIPhreeqcF
        CreateIPhreeqc = CreateIPhreeqcF()
      END FUNCTION CreateIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION DestroyIPhreeqc(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: DestroyIPhreeqc
        INTEGER(KIND=4)  :: DestroyIPhreeqcF
        DestroyIPhreeqc = DestroyIPhreeqcF(ID)
      END FUNCTION DestroyIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponent(ID,N,COMP)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: COMP        
        INTEGER(KIND=4)  :: GetComponent
        INTEGER(KIND=4)  :: GetComponentF
        GetComponent = GetComponentF(ID,N,COMP)
      END FUNCTION GetComponent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponentCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetComponentCount
        INTEGER(KIND=4) :: GetComponentCountF
        GetComponentCount = GetComponentCountF(ID)
      END FUNCTION GetComponentCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetDumpStringLine
        INTEGER(KIND=4)  :: GetDumpStringLineF
        GetDumpStringLine = GetDumpStringLineF(ID,N,LINE)
      END FUNCTION GetDumpStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetDumpStringLineCount
        INTEGER(KIND=4) :: GetDumpStringLineCountF
        GetDumpStringLineCount = GetDumpStringLineCountF(ID)
      END FUNCTION GetDumpStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetDumpFileOn
        INTEGER(KIND=4) :: GetDumpFileOnF
        IF (GetDumpFileOnF(ID).EQ.0) THEN
          GetDumpFileOn = .FALSE.
        ELSE
          GetDumpFileOn = .TRUE.
        ENDIF
      END FUNCTION GetDumpFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpStringOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetDumpStringOn
        INTEGER(KIND=4) :: GetDumpStringOnF
        IF (GetDumpStringOnF(ID).EQ.0) THEN
          GetDumpStringOn = .FALSE.
        ELSE
          GetDumpStringOn = .TRUE.
        ENDIF
      END FUNCTION GetDumpStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetErrorStringLine
        INTEGER(KIND=4)  :: GetErrorStringLineF
        GetErrorStringLine = GetErrorStringLineF(ID,N,LINE)
      END FUNCTION GetErrorStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetErrorStringLineCount
        INTEGER(KIND=4) :: GetErrorStringLineCountF
        GetErrorStringLineCount = GetErrorStringLineCountF(ID)
      END FUNCTION GetErrorStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetErrorFileOn
        INTEGER(KIND=4) :: GetErrorFileOnF
        IF (GetErrorFileOnF(ID).EQ.0) THEN
          GetErrorFileOn = .FALSE.
        ELSE
          GetErrorFileOn = .TRUE.
        ENDIF
      END FUNCTION GetErrorFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetLogFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetLogFileOn
        INTEGER(KIND=4) :: GetLogFileOnF
        IF (GetLogFileOnF(ID).EQ.0) THEN
          GetLogFileOn = .FALSE.
        ELSE
          GetLogFileOn = .TRUE.
        ENDIF
      END FUNCTION GetLogFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetOutputFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetOutputFileOn
        INTEGER(KIND=4) :: GetOutputFileOnF
        IF (GetOutputFileOnF(ID).EQ.0) THEN
          GetOutputFileOn = .FALSE.
        ELSE
          GetOutputFileOn = .TRUE.
        ENDIF
      END FUNCTION GetOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputColumnCount
        INTEGER(KIND=4) :: GetSelectedOutputColumnCountF
        GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF(ID)
      END FUNCTION GetSelectedOutputColumnCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputFileOn(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: GetSelectedOutputFileOn
        INTEGER(KIND=4) :: GetSelectedOutputFileOnF
        IF (GetSelectedOutputFileOnF(ID).EQ.0) THEN
          GetSelectedOutputFileOn = .FALSE.
        ELSE
          GetSelectedOutputFileOn = .TRUE.
        ENDIF
      END FUNCTION GetSelectedOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputRowCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputRowCount
        INTEGER(KIND=4) :: GetSelectedOutputRowCountF
        GetSelectedOutputRowCount = GetSelectedOutputRowCountF(ID)
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputValue(ID,ROW,COL,VTYPE,DVALUE,SVALUE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: ROW
        INTEGER(KIND=4)  :: COL
        INTEGER(KIND=4)  :: VTYPE
        REAL(KIND=8)     :: DVALUE
        CHARACTER(LEN=*) :: SVALUE
        INTEGER(KIND=4)  :: GetSelectedOutputValue
        INTEGER(KIND=4)  :: GetSelectedOutputValueF
        GetSelectedOutputValue = GetSelectedOutputValueF(ID,ROW,
     &                     COL,VTYPE,DVALUE,SVALUE)
      END FUNCTION GetSelectedOutputValue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningStringLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetWarningStringLine
        INTEGER(KIND=4)  :: GetWarningStringLineF
        GetWarningStringLine = GetWarningStringLineF(ID,N,LINE)
      END FUNCTION GetWarningStringLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningStringLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetWarningStringLineCount
        INTEGER(KIND=4) :: GetWarningStringLineCountF
        GetWarningStringLineCount = GetWarningStringLineCountF(ID)
      END FUNCTION GetWarningStringLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabase(ID,FILENAME)
        IMPLICIT NONE
        INTEGER (KIND=4) :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER(KIND=4)  :: LoadDatabase
        INTEGER(KIND=4)  :: LoadDatabaseF
        LoadDatabase = LoadDatabaseF(ID,FILENAME)
      END FUNCTION LoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabaseString(ID,INPUT)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER(KIND=4)  :: LoadDatabaseString
        INTEGER(KIND=4)  :: LoadDatabaseStringF
        LoadDatabaseString = LoadDatabaseStringF(ID,INPUT)
      END FUNCTION LoadDatabaseString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputErrorString(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputErrorStringF(ID)
      END SUBROUTINE OutputErrorString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputAccumulatedLines(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputAccumulatedLinesF(ID)
      END SUBROUTINE OutputAccumulatedLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputWarningString(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputWarningStringF(ID)
      END SUBROUTINE OutputWarningString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetSelectedOutputFileOn(ID,SELECTED_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: SELECTED_ON
        INTEGER(KIND=4) :: SetSelectedOutputFileOn
        INTEGER(KIND=4) :: SetSelOutFileOnF
        SetSelectedOutputFileOn = SetSelOutFileOnF(ID,SELECTED_ON)
      END FUNCTION SetSelectedOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputFileOn(ID,OUTPUT_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: OUTPUT_ON
        INTEGER(KIND=4) :: SetOutputFileOn
        INTEGER(KIND=4) :: SetOutputFileOnF
        SetOutputFileOn = SetOutputFileOnF(ID,OUTPUT_ON)
      END FUNCTION SetOutputFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorFileOn(ID,ERROR_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: ERROR_ON
        INTEGER(KIND=4) :: SetErrorFileOn
        INTEGER(KIND=4) :: SetErrorFileOnF
        SetErrorFileOn = SetErrorFileOnF(ID,ERROR_ON)
      END FUNCTION SetErrorFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogFileOn(ID,LOG_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: LOG_ON
        INTEGER(KIND=4) :: SetLogFileOn
        INTEGER(KIND=4) :: SetLogFileOnF
        SetLogFileOn = SetLogFileOnF(ID,LOG_ON)
      END FUNCTION SetLogFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpFileOn(ID,DUMP_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: DUMP_ON
        INTEGER(KIND=4) :: SetDumpFileOn
        INTEGER(KIND=4) :: SetDumpFileOnF
        SetDumpFileOn = SetDumpFileOnF(ID,DUMP_ON)
      END FUNCTION SetDumpFileOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpStringOn(ID,DUMP_STRING_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: DUMP_STRING_ON
        INTEGER(KIND=4) :: SetDumpStringOn
        INTEGER(KIND=4) :: SetDumpStringOnF
        SetDumpStringOn = SetDumpStringOnF(ID,DUMP_STRING_ON)
      END FUNCTION SetDumpStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunAccumulated(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: RunAccumulated
        INTEGER(KIND=4)  :: RunAccumulatedF
        RunAccumulated = RunAccumulatedF(ID)
      END FUNCTION RunAccumulated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunFile(ID,FILENAME)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER(KIND=4)  :: RunFile
        INTEGER(KIND=4)  :: RunFileF
        RunFile = RunFileF(ID,FILENAME)
      END FUNCTION RunFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunString(ID,INPUT)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER(KIND=4)  :: RunString
        INTEGER(KIND=4)  :: RunStringF
        RunString = RunStringF(ID,INPUT)
      END FUNCTION RunString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION UnLoadDatabase(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: UnLoadDatabase
        INTEGER(KIND=4)  :: UnLoadDatabaseF
        UnLoadDatabase = UnLoadDatabaseF(ID)
      END FUNCTION UnLoadDatabase
