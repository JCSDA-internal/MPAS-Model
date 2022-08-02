
## Source files
# physics/
set(_atm_core_physics_srcs
        mpas_atmphys_camrad_init.F
        mpas_atmphys_constants.F
        mpas_atmphys_control.F
        mpas_atmphys_date_time.F
        mpas_atmphys_driver.F
        mpas_atmphys_driver_cloudiness.F
        mpas_atmphys_driver_convection.F
        mpas_atmphys_driver_gwdo.F
        mpas_atmphys_driver_lsm.F
        mpas_atmphys_driver_microphysics.F
        mpas_atmphys_driver_oml.F
        mpas_atmphys_driver_pbl.F
        mpas_atmphys_driver_radiation_lw.F
        mpas_atmphys_driver_radiation_sw.F
        mpas_atmphys_driver_sfclayer.F
        mpas_atmphys_finalize.F
        mpas_atmphys_functions.F
        mpas_atmphys_init.F
        mpas_atmphys_init_microphysics.F
        mpas_atmphys_interface.F
        mpas_atmphys_landuse.F
        mpas_atmphys_lsm_noahinit.F
        mpas_atmphys_manager.F
        mpas_atmphys_o3climatology.F
        mpas_atmphys_packages.F
        mpas_atmphys_rrtmg_lwinit.F
        mpas_atmphys_rrtmg_swinit.F
        mpas_atmphys_todynamics.F
        mpas_atmphys_update.F
        mpas_atmphys_update_surface.F
        mpas_atmphys_utilities.F
        mpas_atmphys_vars.F)
list(TRANSFORM _atm_core_physics_srcs PREPEND physics/)

# physics/physics_wrf/
set(_atm_core_physics_wrf_srcs
        libmassv.F
        module_bl_gwdo.F
        module_bl_mynn.F
        module_bl_ysu.F
        module_cam_error_function.F
        module_cam_shr_kind_mod.F
        module_cam_support.F
        module_cu_gf.mpas.F
        module_cu_kfeta.F
        module_cu_ntiedtke.F
        module_cu_tiedtke.F
        module_mp_kessler.F
        module_mp_radar.F
        module_mp_thompson.F
        module_mp_thompson_cldfra3.F
        module_mp_wsm6.F
        module_ra_cam.F
        module_ra_cam_support.F
        module_ra_rrtmg_lw.F
        module_ra_rrtmg_sw.F
        module_ra_rrtmg_vinterp.F
        module_sf_bem.F
        module_sf_bep.F
        module_sf_bep_bem.F
        module_sf_mynn.F
        module_sf_noah_seaice.F
        module_sf_noah_seaice_drv.F
        module_sf_noahdrv.F
        module_sf_noahlsm.F
        module_sf_noahlsm_glacial_only.F
        module_sf_oml.F
        module_sf_sfcdiags.F
        module_sf_sfclay.F
        module_sf_urban.F)
list(TRANSFORM _atm_core_physics_wrf_srcs PREPEND physics/physics_wrf/)

# chemistry/
set(_atm_core_chemistry_srcs
        chemistry_opp/o3_gfs.F
        mpas_chemistry_driver.F
        mpas_chemistry_init.F
        mpas_chemistry_interface.F
        mpas_chemistry_manager.F
        mpas_chemistry_opp.F
        mpas_chemistry_packages.F
        mpas_chemistry_todynamics.F)
list(TRANSFORM _atm_core_chemistry_srcs PREPEND chemistry/)


# diagnostics/
set(_atm_core_diagnostics_srcs
        convective_diagnostics.F
        isobaric_diagnostics.F
        mpas_atm_diagnostic_template.F
        mpas_atm_diagnostics_manager.F
        mpas_atm_diagnostics_utils.F
        pv_diagnostics.F
        soundings.F)
list(TRANSFORM _atm_core_diagnostics_srcs PREPEND diagnostics/)

# dynamics/
set(_atm_core_dynamics_srcs
        mpas_atm_boundaries.F
        mpas_atm_iau.F
        mpas_atm_time_integration.F)
list(TRANSFORM _atm_core_dynamics_srcs PREPEND dynamics/)

# utils/
set(_atm_core_utils_srcs
        atmphys_build_tables_thompson.F
        build_tables.F)
list(TRANSFORM _atm_core_utils_srcs PREPEND utils/)

# core_atosphere
set(_atm_core_srcs
        mpas_atm_core.F
        mpas_atm_core_interface.F
        mpas_atm_dimensions.F
        mpas_atm_threading.F)

## Generated includes
set(atm_core_inc
        block_dimension_routines.inc
        core_variables.inc
        define_packages.inc
        domain_variables.inc
        namelist_call.inc
        namelist_defines.inc
        setup_immutable_streams.inc
        structs_and_variables.inc)


add_library(core_atmosphere ${_atm_core_srcs}
                            ${_atm_core_physics_srcs}
                            ${_atm_core_physics_wrf_srcs}
                            ${_atm_core_chemistry_srcs}
                            ${_atm_core_diagnostics_srcs}
                            ${_atm_core_dynamics_srcs})

mpas_core_target(CORE atmosphere TARGET core_atmosphere INCLUDES ${atm_core_inc})

target_compile_definitions(core_atmosphere PRIVATE DO_PHYSICS=1)
target_compile_definitions(core_atmosphere PRIVATE DO_CHEMISTRY=1)
target_compile_definitions(core_atmosphere PRIVATE mpas=1)

#Get physics_wrf tables from MPAS-Data
include(FetchContent)
if(${PROJECT_VERSION} VERSION_GREATER_EQUAL 7.0)
   set(MPAS_DATA_GIT_TAG v${PROJECT_VERSION})
else()
   set(MPAS_DATA_GIT_TAG master)
endif()

FetchContent_Declare(mpas_data
                     GIT_REPOSITORY https://github.com/MPAS-Dev/MPAS-Data.git
                     GIT_TAG ${MPAS_DATA_GIT_TAG}
                     GIT_PROGRESS True
                     GIT_SHALLOW True)
FetchContent_Populate(mpas_data)
set(_physics_wrf_datadir ${mpas_data_SOURCE_DIR}/atmosphere/physics_wrf/files)
file(GLOB _physics_wrf_data RELATIVE ${_physics_wrf_datadir} "${_physics_wrf_datadir}/*")
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/core_atmosphere)
foreach(_datafile IN LISTS _physics_wrf_data)
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${_physics_wrf_datadir}/${_datafile}
                                                                ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/core_atmosphere/${_datafile})
endforeach()
install(DIRECTORY ${_physics_wrf_datadir}/ DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/core_atmosphere)

add_executable(mpas_atmosphere_build_tables ${_atm_core_utils_srcs})
target_link_libraries(mpas_atmosphere_build_tables PUBLIC core_atmosphere)
mpas_fortran_target(mpas_atmosphere_build_tables)
install(TARGETS mpas_atmosphere_build_tables EXPORT ${PROJECT_NAME}ExportsCore
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
