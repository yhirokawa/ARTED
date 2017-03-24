set(BACKUP_UTIL dump_backup)

set(SOURCES
    dump_backup.f90
    )

add_executable(${BACKUP_UTIL} ${SOURCES})
add_dependencies(${BACKUP_UTIL} ${COMMON_LIB} ${MODULE_LIB})
target_link_libraries(${BACKUP_UTIL} ${COMMON_LIB} ${MODULE_LIB})
