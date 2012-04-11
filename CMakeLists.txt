cmake_minimum_required(VERSION 2.6)
 
PROJECT(generationDWI)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  MESSAGE(FATAL_ERROR
          "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

SET(src main.cpp)
GENERATECLP(src main.xml)


add_executable(generationDWI  ${src})
target_link_libraries(generationDWI ${ITK_LIBRARIES} vtkCommon ${VTK_LIBRARIES})


