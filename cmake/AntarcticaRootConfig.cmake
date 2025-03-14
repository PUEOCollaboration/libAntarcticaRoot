include("${CMAKE_CURRENT_LIST_DIR}/AntarcticaRootTargets.cmake")

# CRUCIAL!! This must be here for downstream projects to find RootFftwWrapper transitively through
# AntarcticaRoot. See
#   1) https://stackoverflow.com/questions/59649679/cmake-propagate-dependencies-using-find-package
#   2) https://stackoverflow.com/questions/60856992/if-i-find-package-in-cmakelists-txt-must-i-find-dependency-in-my-installed-conf

include(CMakeFindDependencyMacro)
find_dependency(RootFftwWrapper REQUIRED)

# I guess it doesn't hurt to add checks here for stuff that doesn't propagate
find_dependency(ZLIB REQUIRED) # TODO: figure out if zlib is required
find_dependency(ROOT CONFIG REQUIRED COMPONENTS HistPainter Physics Ged)

