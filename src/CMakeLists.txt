set(SOURCES AliHFInvMassFitter.cxx AliHFMassFitter.cxx AliLog.cxx)

string(REPLACE ".cxx" ".h" HEADERS "${SOURCES}")
list(APPEND HEADERS "AliVertexingHFUtils.h")

add_library(SignalExtraction SHARED ${SOURCES} G__SignalExtraction.cxx)
ROOT_GENERATE_DICTIONARY(G__SignalExtraction ${HEADERS} LINKDEF SignalExtractionLinkDef.h)
target_link_libraries(SignalExtraction ${ROOT_LIBRARIES})
