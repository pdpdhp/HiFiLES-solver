/*!
 * \file bdy_inters.h
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

#include "inters.h"
#include "bdy_inters.h"
#include "array.h"
#include "solution.h"

struct solution; // forwards declaration

class bdy_inters: public inters
{
public:
	
	// #### constructors ####
	
	// default constructor
	
	bdy_inters();

    // default destructor
 
    ~bdy_inters();

	// #### methods ####
	
	/*! setup inters */
	void setup(int in_n_inters, int in_inter_type, int in_run_type);

    /*! setup array that contains boundary parameters */
    void set_bdy_params();

    /*! Set bdy interface */
    void set_boundary(int in_inter, int bdy_type, int in_ele_type_l, int in_ele_l, int in_local_inter_l, int in_run_type, struct solution* FlowSol);

    /*! Compute right hand side state at boundaries */
    void set_inv_boundary_conditions(int bdy_type, double* u_l, double* u_r, double *gv, double *norm, double *loc, double *bdy_params, int n_dims, int n_fields, double gamma, double R_ref, double time_bound, int equation);
  
    /*! Compute right hand side gradient at boundaries */
    void set_vis_boundary_conditions(int bdy_type, double* u_l, double* u_r, double* grad_u, double *norm, double *loc, double *bdy_params, int n_dims, int n_fields, double gamma, double R_ref, double time_bound, int equation);

	/*! move all from cpu to gpu */
	void mv_all_cpu_gpu(void);

	/*! calculate normal transformed continuous inviscid flux at the flux points on boundaries*/
	void calc_norm_tconinvf_fpts_boundary(double time_bound);

    /*! calculate delta in transformed discontinuous solution at flux points */
    void calc_delta_disu_fpts_boundary(void);
	
    /*! calculate normal transformed continuous viscous flux at the flux points on boundaries*/
    void calc_norm_tconvisf_fpts_boundary(double time_bound);
	
protected:

	// #### members ####

  int max_bdy_params;

  array<int> boundary_type;	
  array<double> bdy_params;

};
