# files
SET(Java_Files
advect.pqi
advection.java
${CMAKE_SOURCE_DIR}/database/phreeqc.dat
)

# install directory
SET(Java_Dir ${EXAMPLES_DIR}/swig/java)

# install
install(FILES ${Java_Files} DESTINATION ${Java_Dir})
install(DIRECTORY "${CMAKE_BINARY_DIR}/swig/" DESTINATION ${Java_Dir} FILES_MATCHING PATTERN "*.java")

