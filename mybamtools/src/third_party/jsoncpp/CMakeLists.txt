# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/third-party/jsoncpp
# ==========================

# add compiler definitions 
add_definitions(-DBAMTOOLS_JSONCPP_LIBRARY) # (for proper exporting of library symbols)

# create jsoncpp library
add_library ( jsoncpp SHARED 
              json_reader.cpp
              json_value.cpp
              json_writer.cpp
            )
add_library ( jsoncppST STATIC 
              json_reader.cpp
              json_value.cpp
              json_writer.cpp
            )

set_target_properties(jsoncppST PROPERTIES OUTPUT_NAME jsoncpp)

# set jsoncpp library properties
set_target_properties( jsoncpp PROPERTIES
                       SOVERSION   1.0.0
                       OUTPUT_NAME jsoncpp
                     )
