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
      FUNCTION UnLoadDatabase(ID)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: UnLoadDatabase
        INTEGER(KIND=4)  :: UnLoadDatabaseF
        UnLoadDatabase = UnLoadDatabaseF(ID)
      END FUNCTION UnLoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputError(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputErrorF(ID)
      END SUBROUTINE OutputError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputWarning(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputWarningF(ID)
      END SUBROUTINE OutputWarning
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
      FUNCTION SetSelectedOutputOn(ID,SELECTED_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: SELECTED_ON
        INTEGER(KIND=4) :: SetSelectedOutputOn
        INTEGER(KIND=4) :: SetSelectedOutputOnF
        SetSelectedOutputOn = SetSelectedOutputOnF(ID,SELECTED_ON)
      END FUNCTION SetSelectedOutputOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputOn(ID,OUTPUT_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: OUTPUT_ON
        INTEGER(KIND=4) :: SetOutputOn
        INTEGER(KIND=4) :: SetOutputOnF
        SetOutputOn = SetOutputOnF(ID,OUTPUT_ON)
      END FUNCTION SetOutputOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorOn(ID,ERROR_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: ERROR_ON
        INTEGER(KIND=4) :: SetErrorOn
        INTEGER(KIND=4) :: SetErrorOnF
        SetErrorOn = SetErrorOnF(ID,ERROR_ON)
      END FUNCTION SetErrorOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogOn(ID,LOG_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: LOG_ON
        INTEGER(KIND=4) :: SetLogOn
        INTEGER(KIND=4) :: SetLogOnF
        SetLogOn = SetLogOnF(ID,LOG_ON)
      END FUNCTION SetLogOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpOn(ID,DUMP_ON)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
		LOGICAL(KIND=4) :: DUMP_ON
        INTEGER(KIND=4) :: SetDumpOn
        INTEGER(KIND=4) :: SetDumpOnF
        SetDumpOn = SetDumpOnF(ID,DUMP_ON)
      END FUNCTION SetDumpOn
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
      FUNCTION GetDumpLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetDumpLineCount
        INTEGER(KIND=4) :: GetDumpLineCountF
        GetDumpLineCount = GetDumpLineCountF(ID)
      END FUNCTION GetDumpLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetDumpLine
        INTEGER(KIND=4)  :: GetDumpLineF
        GetDumpLine = GetDumpLineF(ID,N,LINE)
      END FUNCTION GetDumpLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetErrorLineCount
        INTEGER(KIND=4) :: GetErrorLineCountF
        GetErrorLineCount = GetErrorLineCountF(ID)
      END FUNCTION GetErrorLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetErrorLine
        INTEGER(KIND=4)  :: GetErrorLineF
        GetErrorLine = GetErrorLineF(ID,N,LINE)
      END FUNCTION GetErrorLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponentCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetComponentCount
        INTEGER(KIND=4) :: GetComponentCountF
        GetComponentCount = GetComponentCountF(ID)
      END FUNCTION GetComponentCount
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
      FUNCTION GetSelectedOutputRowCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputRowCount
        INTEGER(KIND=4) :: GetSelectedOutputRowCountF
        GetSelectedOutputRowCount = GetSelectedOutputRowCountF(ID)
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetSelectedOutputColumnCount
        INTEGER(KIND=4) :: GetSelectedOutputColumnCountF
        GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF(ID)
      END FUNCTION GetSelectedOutputColumnCount
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
      SUBROUTINE OutputLines(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        CALL OutputLinesF(ID)
      END SUBROUTINE OutputLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningLineCount(ID)
        IMPLICIT NONE
        INTEGER(KIND=4) :: ID
        INTEGER(KIND=4) :: GetWarningLineCount
        INTEGER(KIND=4) :: GetWarningLineCountF
        GetWarningLineCount = GetWarningLineCountF(ID)
      END FUNCTION GetWarningLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER(KIND=4)  :: ID
        INTEGER(KIND=4)  :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER(KIND=4)  :: GetWarningLine
        INTEGER(KIND=4)  :: GetWarningLineF
        GetWarningLine = GetWarningLineF(ID,N,LINE)
      END FUNCTION GetWarningLine
