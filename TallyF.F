!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION build_tally_table()
        IMPLICIT NONE
        INTEGER  :: build_tally_table
#if defined(_WIN32)
        INTERFACE 
          FUNCTION build_tally_tableF()
      !DEC$ ATTRIBUTES C,REFERENCE::build_tally_tableF
      !DEC$ ATTRIBUTES ALIAS:'_build_tally_tableF'::build_tally_tableF
            INTEGER(KIND=4) :: build_tally_tableF
          END FUNCTION build_tally_tableF
        END INTERFACE
#else
        INTEGER          :: build_tally_tableF
#endif

        build_tally_table = build_tally_tableF()
      END FUNCTION build_tally_table


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION get_tally_table_rows_columns(rows, cols)
        IMPLICIT NONE
        INTEGER  :: get_tally_table_rows_columns
        INTEGER  :: rows
        INTEGER  :: cols
#if defined(_WIN32)
        INTERFACE 
          FUNCTION rows_cols(rows,cols)
      !DEC$ ATTRIBUTES C,REFERENCE::rows_cols
      !DEC$ ATTRIBUTES ALIAS:'_get_tally_table_rows_columnsF'::rows_cols
            INTEGER  :: rows
            INTEGER  :: cols
            INTEGER(KIND=4) :: rows_cols
          END FUNCTION rows_cols
        END INTERFACE

        get_tally_table_rows_columns = rows_cols(rows,cols)
#else
        INTEGER get_tally_table_rows_columnsF
        get_tally_table_rows_columns = 
     &        get_tally_table_rows_columnsF(rows,cols)
#endif
      END FUNCTION get_tally_table_rows_columns


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION get_tally_table_column_heading(col,etype,string)
        IMPLICIT NONE
        INTEGER      col
        INTEGER      etype
        CHARACTER(*) string
        INTEGER      adjcol
        INTEGER      get_tally_table_column_heading
#if defined(_WIN32)
        INTERFACE 
          FUNCTION col_head(col,etype,string)
  !DEC$ ATTRIBUTES C,REFERENCE::col_head
  !DEC$ ATTRIBUTES ALIAS:'_get_tally_table_column_headingF'::col_head
            INTEGER      col
            INTEGER      etype
            CHARACTER(*) string
            INTEGER      col_head
          END FUNCTION col_head
        END INTERFACE
        adjcol = col - 1
        get_tally_table_column_heading = col_head(adjcol,etype,string)
#else
        INTEGER get_tally_table_column_headingF
        adjcol = col - 1
        get_tally_table_column_heading =
     &            get_tally_table_column_headingF(adjcol,etype,string)
#endif
      etype = etype + 1
      END FUNCTION get_tally_table_column_heading

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION get_tally_table_row_heading(row,string)
        IMPLICIT NONE
        INTEGER      row
        CHARACTER(*) string
        INTEGER      adjrow
        INTEGER      get_tally_table_row_heading
#if defined(_WIN32)
        INTERFACE 
          FUNCTION row_head(row,string)
  !DEC$ ATTRIBUTES C,REFERENCE::row_head
  !DEC$ ATTRIBUTES ALIAS:'_get_tally_table_row_headingF'::row_head
            INTEGER      row
            CHARACTER(*) string
            INTEGER      row_head
          END FUNCTION row_head
        END INTERFACE
        adjrow = row - 1
        get_tally_table_row_heading = row_head(adjrow,string)
#else
        INTEGER get_tally_table_row_headingF
        adjrow = row - 1
        get_tally_table_row_heading =
     &            get_tally_table_row_headingF(adjrow,string)
#endif
      END FUNCTION get_tally_table_row_heading
