project(rootIO)

find_package(ROOT REQUIRED)

option(ENABLE_FASTJET "Build FastJet-related interfaces in rootIO" OFF)

if(ENABLE_FASTJET)
  # FastJet discovery: prefer FASTJET_DIR hint, fallback to common system paths
  find_path(FASTJET_INCLUDE_DIR fastjet/PseudoJet.hh
    HINTS ${FASTJET_DIR}/include /usr/include /usr/local/include
  )
  find_library(FASTJET_LIBRARY fastjet
    HINTS ${FASTJET_DIR}/lib ${FASTJET_DIR}/lib64 /usr/lib64 /usr/lib /usr/local/lib
  )
  if(NOT FASTJET_INCLUDE_DIR OR NOT FASTJET_LIBRARY)
    message(FATAL_ERROR "ENABLE_FASTJET=ON but FastJet not found. Set FASTJET_DIR or install FastJet.")
  endif()
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  $<$<BOOL:${ENABLE_FASTJET}>:${FASTJET_INCLUDE_DIR}>
)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*Interface.h
)

# If FastJet is disabled, drop its interface from sources/headers
if(NOT ENABLE_FASTJET)
  list(FILTER sources EXCLUDE REGEX ".*/fastjetInterface\\.cc$")
  list(FILTER headers EXCLUDE REGEX ".*/fastjetInterface\\.h$")
endif()

include(${ROOT_USE_FILE})
ROOT_GENERATE_DICTIONARY(G__rootIO ${headers} LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)
add_library(rootIO SHARED ${sources} G__rootIO.cxx)
target_include_directories(rootIO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} $<$<BOOL:${ENABLE_FASTJET}>:${FASTJET_INCLUDE_DIR}>)
target_link_libraries(
  rootIO
  $<$<BOOL:${ENABLE_FASTJET}>:${FASTJET_LIBRARY}>
  ${ROOT_LIBRARIES}
)

install(
  TARGETS rootIO
  LIBRARY DESTINATION lib
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/librootIO_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/librootIO.rootmap DESTINATION lib
)
