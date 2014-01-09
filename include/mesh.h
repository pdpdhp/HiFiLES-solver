/*!
 * \file mesh.h
 * \brief  - Class to control mesh-related activities (motion, adaptation, etc.)
 * \author - Current development: Aerospace Computing Laboratory (ACL) directed
 *                                by Prof. Jameson. (Aero/Astro Dept. Stanford University).
 * \version 1.0.0
 *
 * HiFiLES (High Fidelity Large Eddy Simulation).
 * Copyright (C) 2013 Aerospace Computing Laboratory.
 */

#pragma once

#include <iostream>
#include <sstream>
#include <cmath>

#include "array.h"
#include "eles.h"
#include "global.h"
#include "input.h"
#include "geometry.h"
#include "funcs.h"
#include "error.h"
#include "solution.h"

/*! necessary structure files from SU2 */
#include "matrix_structure.hpp"
#include "linear_solvers_structure.hpp"

#ifdef _TECIO
#include "TECIO.h"
#endif

#ifdef _MPI
#include "mpi.h"
#include "metis.h"
#include "parmetis.h"
#endif

#ifdef _GPU
#include "../include/util.h"
#endif

class mesh
{
public:
	// #### constructors ####
	
	/** default constructor */
	mesh();
	
	/** default destructor */
	~mesh();

	// #### methods ####

	/** peform prescibed mesh motion using linear elasticity method*/
	void deform(struct solution* FlowSol);

    /** create individual-element stiffness matrix - triangles */
    // will I actually need the FlowSol variable for setting up the Stiffnexx Matrix?
    bool set_2D_StiffMat_ele_tri(array<double> &stiffMat_ele,int ele_id, solution *FlowSol);

    /** create individual-element stiffness matrix - triangles */
    bool set_2D_StiffMat_ele_quad(array<double> &stiffMat_ele,int ele_id, solution *FlowSol);

    /** create individual-element stiffness matrix - triangles */
    bool set_2D_StiffMat_ele_tet(array<double> &stiffMat_ele,int ele_id, solution *FlowSol);

    /** create individual-element stiffness matrix - triangles */
    bool set_2D_StiffMat_ele_hex(array<double> &stiffMat_ele,int ele_id, solution *FlowSol);

	// #### members ####
	
    int n_eles, v_verts;

	/** arrays which define the basic mesh geometry */
	array<double> xv;
	array<int> c2v,c2n_v,ctype,bctype_c,ic2icg,iv2ivg;
	array<int> f2c,f2loc_f,c2f,c2e,f2v,f2nv;
	/** global stiffness matrix for linear-elasticity mesh motion */
    array<array<double> > stiff_mat;

    /** enumeration for cell type */
    enum CTYPE {
        TRI = 0,
        QUAD = 1,
        TET = 2,
        WEDGE = 3,
        BRICK = 4
    };

	
private:

	/** Global stiffness matrix for linear elasticity solution */
	CSysMatrix StiffnessMatrix;

};
