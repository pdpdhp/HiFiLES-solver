/*!
 * \file solver.h
 * \brief _____________________________
 * \author - Original code: SD++ developed by Patrice Castonguay, Antony Jameson,
 *                          Peter Vincent, David Williams (alphabetical by surname).
 *         - Current development: Aerospace Computing Laboratory (ACL) directed
 *                                by Prof. Jameson. (Aero/Astro Dept. Stanford University).
 * \version 1.0.0
 *
 * HiFiLES (High Fidelity Large Eddy Simulation).
 * Copyright (C) 2013 Aerospace Computing Laboratory.
 */

#pragma once

#include "array.h"
#include <string>
#include "input.h"
#include "eles.h"
#include "eles_tris.h"
#include "eles_quads.h"
#include "eles_hexas.h"
#include "eles_tets.h"
#include "eles_pris.h"
#include "int_inters.h"
#include "bdy_inters.h"
#include "solution.h"

#ifdef _MPI
#include "mpi.h"
#include "mpi_inters.h"
#endif

/*!
 * \brief Calculate the residual.
 * \param[in] FlowSol - Structure with the entire solution and mesh information.
 */
void CalcResidual(struct solution* FlowSol);

void set_rank_nproc(int in_rank, int in_nproc, struct solution* FlowSol);

/*! get pointer to transformed discontinuous solution at a flux point */
double* get_disu_fpts_ptr(int in_ele_type, int in_ele, int in_field, int n_local_inter, int in_fpt, struct solution* FlowSol);

/*! get pointer to normal continuous transformed inviscid flux at a flux point */
double* get_norm_tconf_fpts_ptr(int in_ele_type, int in_ele, int in_field, int in_local_inter, int in_fpt, struct solution* FlowSol);

/*! get pointer to determinant of jacobian at a flux point */
double* get_detjac_fpts_ptr(int in_ele_type, int in_ele, int in_ele_local_inter, int in_inter_local_fpt, struct solution* FlowSol);

/*! get pointer to magntiude of normal dot inverse of (determinant of jacobian multiplied by jacobian) at a solution point */
double* get_mag_tnorm_dot_inv_detjac_mul_jac_fpts_ptr(int in_ele_type, int in_ele, int in_ele_local_inter, int in_inter_local_fpt, struct solution* FlowSol);

/*! get pointer to normal at a flux point */
double* get_norm_fpts_ptr(int in_ele_type, int in_ele, int in_local_inter, int in_fpt, int in_dim, struct solution* FlowSol);

/*! get pointer to coordinates at a flux point */
double* get_loc_fpts_ptr(int in_ele_type, int in_ele, int in_local_inter, int in_fpt, int in_dim, struct solution* FlowSol);

/*! get pointer to delta of the transformed discontinuous solution at a flux point */
double* get_delta_disu_fpts_ptr(int in_ele_type, int in_ele, int in_field, int n_local_inter, int in_fpt, struct solution* FlowSol);

/*! get pointer to gradient of the discontinuous solution at a flux point */
double* get_grad_disu_fpts_ptr(int in_ele_type, int in_ele, int in_local_inter, int in_field, int in_dim, int in_fpt, struct solution* FlowSol);

/*! get pointer to grid velocity at a flux point */
double* get_vel_fpts_ptr(int in_ele_type, int in_ele, int in_local_inter, int in_fpt, int in_dim, struct solution* FlowSol);

// Initialize the solution in the mesh
void InitSolution(struct solution* FlowSol);

/*! reading a restart file */
void read_restart(int in_file_num, int in_n_files, struct solution* FlowSol);






