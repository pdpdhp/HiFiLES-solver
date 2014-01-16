TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    src/solver.cpp \
    src/output.cpp \
    src/mpi_inters.cpp \
    src/int_inters.cpp \
    src/inters.cpp \
    src/input.cpp \
    src/HiFiLES.cpp \
    src/global.cpp \
    src/geometry.cpp \
    src/funcs.cpp \
    src/flux.cpp \
    src/eles_tris.cpp \
    src/eles_tets.cpp \
    src/eles_quads.cpp \
    src/eles_pris.cpp \
    src/eles_hexas.cpp \
    src/eles.cpp \
    src/cuda_kernels.cu \
    src/cubature_tri.cpp \
    src/cubature_tet.cpp \
    src/cubature_quad.cpp \
    src/cubature_hexa.cpp \
    src/cubature_1d.cpp \
    src/bdy_inters.cpp \
    src/mesh_deform.cpp \
    src/matrix_structure.cpp \
    src/linear_solvers_structure.cpp \
    src/mesh.cpp \
    src/vector_structure.cpp \
    include/vector_structure.inl \
    include/matrix_structure.inl \
    include/linear_solvers_structure.inl \
    src/matrix_structure_2.cpp

HEADERS += \
    include/util.h \
    include/solver.h \
    include/solution.h \
    include/rename.h \
    include/parmetisbin.h \
    include/output.h \
    include/mpi_inters.h \
    include/macros.h \
    include/int_inters.h \
    include/inters.h \
    include/input.h \
    include/global.h \
    include/geometry.h \
    include/funcs.h \
    include/flux.h \
    include/error.h \
    include/eles_tris.h \
    include/eles_tets.h \
    include/eles_quads.h \
    include/eles_pris.h \
    include/eles_hexas.h \
    include/eles.h \
    include/cuda_kernels.h \
    include/cubature_tri.h \
    include/cubature_tet.h \
    include/cubature_quad.h \
    include/cubature_hexa.h \
    include/cubature_1d.h \
    include/bdy_inters.h \
    include/array.h \
    include/matrix_structure.hpp \
    include/linear_solvers_structure.hpp \
    include/vector_structure.hpp \
    include/mesh.h

OTHER_FILES += \
    data/loc_tri_inter_pts.dat \
    data/loc_tri_alpha_pts.dat \
    data/loc_tet_inter_pts_old.dat \
    data/loc_tet_inter_pts.dat \
    data/loc_tet_alpha_pts.dat \
    data/loc_1d_gauss_pts.dat \
    data/loc_1d_gauss_lobatto_pts.dat \
    data/cubature_tri.dat \
    data/cubature_tet.dat \
    data/cubature_quad.dat \
    data/cubature_hexa.dat \
    data/cubature_1d.dat \
    makefile.in \
    GMSH_Element_Node_Ordering.txt \
    TODO.txt
