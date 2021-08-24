
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(fit FitXicZerotoXiPiInvMass.C)
target_link_libraries(fit SignalExtraction ${ROOT_LIBRARIES} EG Minuit)
