SET(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0")

FIND_PACKAGE(Qt4)
FIND_PACKAGE(CUDA)
FIND_PACKAGE(BOOST REQUIRED COMPONENTS numeric)

MESSAGE(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(QtApp_MOCS )
SET(QtApp_FORMS gui/mainwindow.ui gui/synthesis_dialog.ui)
SET(QtApp_RCCS application.qrc)

set(CMAKE_AUTOMOC ON)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})


SET(COMMON_SOURCES dtts_2dmorph.cpp dtts_ppa.cpp dtts_wavelet.cpp external/noise/noise.cpp  external/maxflow/graph.cpp external/maxflow/maxflow.cpp)
SET(HEADERS dtts_mainwindow.h)

if(CUDA_FOUND)
  SET(SOURCES ${COMMON_SOURCES} dtts_patchsynthesis.cu dtts_image.cu dtts_match.cu dtts_merging.cu external/SpeedIT_Classic/si_classic.cu)
else()
  SET(SOURCES ${COMMON_SOURCES} dtts_patchsynthesis.cpp dtts_image.cpp dtts_match.cpp dtts_merging.cpp external/SpeedIT_Classic/si_classic_seq.cpp)
endif()

if (QT_FOUND)
  SET(BIN_SOURCES dtts_main_gui.cpp dtts_mainwindow.cpp)
else()
  SET(BIN_SOURCES tts_main.cpp)
endif()

if (CUDA_FOUND)
  SET(EXECFILE TerrainSynthesisGPU)
  set(LIBFILE TerrainSynthesisGPU)
else()
  SET(EXECFILE TerrainSynthesisCPU)
  SET(LIBFILE TerrainSynthesisCPU)
endif()

if (QT_FOUND)
  INCLUDE(${QT_USE_FILE})
  ADD_DEFINITIONS(${QT_DEFINITIONS})
  QT4_WRAP_CPP(QtApp_MOC_SRCS ${QtApp_MOCS})
  QT4_WRAP_UI(QtApp_FORMS_HEADERS ${QtApp_FORMS})
  QT4_ADD_RESOURCES(QtApp_RCC_SRCS ${QtApp_RCCS})
  SET(EXECFILE ${EXECFILE}QT)
endif()

ADD_LIBRARY(${LIBFILE} STATIC ${HEADERS} ${SOURCES})

ADD_EXECUTABLE(${EXECFILE} ${HEADERS} ${SOURCES} ${BIN_SOURCES} ${QtApp_FORMS_HEADERS} ${QtApp_MOC_SRCS} ${QtApp_RCC_SRCS})
if (QT_FOUND)
  TARGET_LINK_LIBRARIES(${EXECFILE} ${QT_LIBRARIES} ${QT_QTDECLARATIVE_LIBRARY})
endif()
