 
set(prjname env_test)

file (GLOB headers *.h *.hpp)
source_group("" FILES ${headers})
set(all_headers ${headers} )

file (GLOB sources *.cpp *.c)
source_group("" FILES ${sources})
set(all_sources ${sources} )

add_executable(${prjname} ${all_headers} ${all_sources})

set_target_properties(${prjname} PROPERTIES DEBUG_OUTPUT_NAME ${prjname}${DEBUG_TARGET_SUFFIX})
set_target_properties(${prjname} PROPERTIES RELEASE_OUTPUT_NAME ${prjname})

include_directories(${PROJECTS_SRC_ROOT}/atm)

target_link_libraries(${prjname} debug env${DEBUG_TARGET_SUFFIX} optimized env)



