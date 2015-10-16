# files
SET(Python_Files
advect.pqi
advection.py
${CMAKE_SOURCE_DIR}/database/phreeqc.dat
)

# install directory
SET(Python_Dir ${EXAMPLES_DIR}/swig/python)

# install
install(FILES ${Python_Files} DESTINATION ${Python_Dir})
install(DIRECTORY "${CMAKE_BINARY_DIR}/swig/" DESTINATION ${Python_Dir} FILES_MATCHING PATTERN "*.py")
