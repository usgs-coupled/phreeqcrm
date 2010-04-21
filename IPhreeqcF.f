!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION CreateIPhreeqc()
        IMPLICIT NONE
        INTEGER          :: CreateIPhreeqc
        INTEGER          :: CreateIPhreeqcF
        CreateIPhreeqc = CreateIPhreeqcF()
      END FUNCTION CreateIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION DestroyIPhreeqc(ID)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: DestroyIPhreeqc
        INTEGER          :: DestroyIPhreeqcF
        DestroyIPhreeqc = DestroyIPhreeqcF(ID)
      END FUNCTION DestroyIPhreeqc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabase(ID,FILENAME)
        IMPLICIT NONE
        INTEGER          :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER          :: LoadDatabase
        INTEGER          :: LoadDatabaseF
        LoadDatabase = LoadDatabaseF(ID,FILENAME)
      END FUNCTION LoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION LoadDatabaseString(ID,INPUT)
        IMPLICIT NONE
        INTEGER          :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER          :: LoadDatabaseString
        INTEGER          :: LoadDatabaseStringF
        LoadDatabaseString = LoadDatabaseStringF(ID,INPUT)
      END FUNCTION LoadDatabaseString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION UnLoadDatabase(ID)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: UnLoadDatabase
        INTEGER          :: UnLoadDatabaseF
        UnLoadDatabase = UnLoadDatabaseF(ID)
      END FUNCTION UnLoadDatabase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLastError(ID)
        IMPLICIT NONE
        INTEGER :: ID
        CALL OutputLastErrorF(ID)
      END SUBROUTINE OutputLastError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLastWarning(ID)
        IMPLICIT NONE
        INTEGER :: ID
        CALL OutputLastWarningF(ID)
      END SUBROUTINE OutputLastWarning
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION AccumulateLine(ID,LINE)
        IMPLICIT NONE
        INTEGER          :: ID
        CHARACTER(LEN=*) :: LINE
        INTEGER          :: AccumulateLine
        INTEGER          :: AccumulateLineF
        AccumulateLine = AccumulateLineF(ID,LINE)
      END FUNCTION AccumulateLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetSelectedOutputOn(ID,SELECTED_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: SELECTED_ON
        INTEGER :: SetSelectedOutputOn
        INTEGER :: SetSelectedOutputOnF
        SetSelectedOutputOn = SetSelectedOutputOnF(ID,SELECTED_ON)
      END FUNCTION SetSelectedOutputOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetOutputOn(ID,OUTPUT_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: OUTPUT_ON
        INTEGER :: SetOutputOn
        INTEGER :: SetOutputOnF
        SetOutputOn = SetOutputOnF(ID,OUTPUT_ON)
      END FUNCTION SetOutputOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetErrorOn(ID,ERROR_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: ERROR_ON
        INTEGER :: SetErrorOn
        INTEGER :: SetErrorOnF
        SetErrorOn = SetErrorOnF(ID,ERROR_ON)
      END FUNCTION SetErrorOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetLogOn(ID,LOG_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: LOG_ON
        INTEGER :: SetLogOn
        INTEGER :: SetLogOnF
        SetLogOn = SetLogOnF(ID,LOG_ON)
      END FUNCTION SetLogOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpOn(ID,DUMP_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: DUMP_ON
        INTEGER :: SetDumpOn
        INTEGER :: SetDumpOnF
        SetDumpOn = SetDumpOnF(ID,DUMP_ON)
      END FUNCTION SetDumpOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SetDumpStringOn(ID,DUMP_STRING_ON)
        IMPLICIT NONE
        INTEGER :: ID
		LOGICAL :: DUMP_STRING_ON
        INTEGER :: SetDumpStringOn
        INTEGER :: SetDumpStringOnF
        SetDumpStringOn = SetDumpStringOnF(ID,DUMP_STRING_ON)
      END FUNCTION SetDumpStringOn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpLineCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetDumpLineCount
        INTEGER :: GetDumpLineCountF
        GetDumpLineCount = GetDumpLineCountF(ID)
      END FUNCTION GetDumpLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetDumpLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER          :: GetDumpLine
        INTEGER          :: GetDumpLineF
        GetDumpLine = GetDumpLineF(ID,N,LINE)
      END FUNCTION GetDumpLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorLineCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetErrorLineCount
        INTEGER :: GetErrorLineCountF
        GetErrorLineCount = GetErrorLineCountF(ID)
      END FUNCTION GetErrorLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetErrorLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER          :: GetErrorLine
        INTEGER          :: GetErrorLineF
        GetErrorLine = GetErrorLineF(ID,N,LINE)
      END FUNCTION GetErrorLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponentCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetComponentCount
        INTEGER :: GetComponentCountF
        GetComponentCount = GetComponentCountF(ID)
      END FUNCTION GetComponentCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetComponent(ID,N,COMP)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: N
        CHARACTER(LEN=*) :: COMP        
        INTEGER          :: GetComponent
        INTEGER          :: GetComponentF
        GetComponent = GetComponentF(ID,N,COMP)
      END FUNCTION GetComponent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunAccumulated(ID)
        IMPLICIT NONE
        INTEGER  :: ID
        INTEGER  :: RunAccumulated
        INTEGER  :: RunAccumulatedF
        RunAccumulated = RunAccumulatedF(ID)
      END FUNCTION RunAccumulated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunFile(ID,FILENAME)
        IMPLICIT NONE
        INTEGER          :: ID
        CHARACTER(LEN=*) :: FILENAME
        INTEGER          :: RunFile
        INTEGER          :: RunFileF
        RunFile = RunFileF(ID,FILENAME)
      END FUNCTION RunFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION RunString(ID,INPUT)
        IMPLICIT NONE
        INTEGER          :: ID
        CHARACTER(LEN=*) :: INPUT
        INTEGER          :: RunString
        INTEGER          :: RunStringF
        RunString = RunStringF(ID,INPUT)
      END FUNCTION RunString
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputRowCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetSelectedOutputRowCount
        INTEGER :: GetSelectedOutputRowCountF
        GetSelectedOutputRowCount = GetSelectedOutputRowCountF(ID)
      END FUNCTION GetSelectedOutputRowCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputColumnCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetSelectedOutputColumnCount
        INTEGER :: GetSelectedOutputColumnCountF
        GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF(ID)
      END FUNCTION GetSelectedOutputColumnCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetSelectedOutputValue(ID,ROW,COL,VTYPE,DVALUE,SVALUE)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: ROW
        INTEGER          :: COL
        INTEGER          :: VTYPE
        REAL*8           :: DVALUE
        CHARACTER(LEN=*) :: SVALUE
        INTEGER          :: GetSelectedOutputValue
        INTEGER          :: GetSelectedOutputValueF
        GetSelectedOutputValue = GetSelectedOutputValueF(ID,ROW,
     &                     COL,VTYPE,DVALUE,SVALUE)
      END FUNCTION GetSelectedOutputValue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE OutputLines(ID)
        IMPLICIT NONE
        INTEGER          :: ID
        CALL OutputLinesF(ID)
      END SUBROUTINE OutputLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningLineCount(ID)
        IMPLICIT NONE
        INTEGER :: ID
        INTEGER :: GetWarningLineCount
        INTEGER :: GetWarningLineCountF
        GetWarningLineCount = GetWarningLineCountF(ID)
      END FUNCTION GetWarningLineCount
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GetWarningLine(ID,N,LINE)
        IMPLICIT NONE
        INTEGER          :: ID
        INTEGER          :: N
        CHARACTER(LEN=*) :: LINE        
        INTEGER          :: GetWarningLine
        INTEGER          :: GetWarningLineF
        GetWarningLine = GetWarningLineF(ID,N,LINE)
      END FUNCTION GetWarningLine
