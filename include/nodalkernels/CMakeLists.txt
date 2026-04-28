if(overlap_coupling_FOUND)
    target_sources(${PROJECT_NAME} PRIVATE "CouplingForce.h" "MinimumExample.h")
endif()
target_sources(${PROJECT_NAME} PRIVATE "CylindricalSurfaceNodalContact.h")
