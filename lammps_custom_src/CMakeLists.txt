 add_library(chromatin SHARED
    bond_harmonic_DNA.cpp
    hremd_steve.cpp
    pair_ljlambda.cpp
    pair_aniso.cpp
    chromatinplugin.cpp
    temper_2.cpp
    temper_2_mod.cpp
    RIGID/fix_rigid_nh_small.cpp
    RIGID/fix_rigid_nve_small.cpp
    RIGID/fix_rigid_small.cpp
)

install(TARGETS chromatin EXPORT chromatin_targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
