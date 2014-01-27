/*!
 * \file mesh_deform.cpp
 * \brief  - Perform mesh deformation using linear elasticity method
 * \author - Current development: Aerospace Computing Laboratory (ACL) directed
 *                                by Prof. Jameson. (Aero/Astro Dept. Stanford University).
 * \version 1.0.0
 *
 * HiFiLES (High Fidelity Large Eddy Simulation).
 * Copyright (C) 2013 Aerospace Computing Laboratory.
 */

#include "../include/mesh.h"

using namespace std;

mesh::mesh(void)
{
    n_eles = 0;
    n_verts = 0;
    n_dims = 2;
    n_verts_global = 0;
    n_cells_global = 0;
    n_bnds = 0;
    n_faces = 0;
    LinSolIters = 0;
    failedIts = 0;
    min_vol = DBL_MAX;
}

void mesh::deform(struct solution* FlowSol) {
    array<double> stiff_mat_ele;
    int failedIts = 0;
    /*
    if (FlowSol->n_dims == 2) {
        stiff_mat_ele.setup(6,6);
    }else{
        FatalError("3D Mesh motion not implemented yet!");
    }
    */

    int pt_0,pt_1,pt_2,pt_3;
    bool check;

    // Setup stiffness matrices for each individual element,
    // combine all element-level matrices into global matrix
    //stiff_mat.setup(n_eles);
    LinSysSol.Initialize(n_verts,n_dims,0.0); /// should it be n_verts or n_verts_global?
    LinSysRes.Initialize(n_verts,n_dims,0.0);
    StiffnessMatrix.Initialize(*this);

    /*--- Loop over the total number of grid deformation iterations. The surface
    deformation can be divided into increments to help with stability. In
    particular, the linear elasticity equations hold only for small deformations. ---*/
    for (iGridDef_Iter = 0; iGridDef_Iter < run_input.n_deform_iters; iGridDef_Iter++) {

        /*--- Initialize vector and sparse matrix ---*/

        LinSysSol.SetValZero();
        LinSysRes.SetValZero();
        StiffnessMatrix.SetValZero();

        /*--- Compute the stiffness matrix entries for all nodes/elements in the
        mesh. FEA uses a finite element method discretization of the linear
        elasticity equations (transfers element stiffnesses to point-to-point). ---*/
        for (int ic=0; ic<n_eles; ic++) {
            switch(ctype(ic))
            {
            case TRI:
                pt_0 = iv2ivg(c2v(ic,0));
                pt_1 = iv2ivg(c2v(ic,1));
                pt_2 = iv2ivg(c2v(ic,2));
                check = set_2D_StiffMat_ele_tri(stiff_mat_ele,ic);
                add_StiffMat_EleTri(stiff_mat_ele,pt_0,pt_1,pt_2);
                break;
            case QUAD:
                pt_0 = iv2ivg(c2v(ic,0));
                pt_1 = iv2ivg(c2v(ic,1));
                pt_2 = iv2ivg(c2v(ic,2));
                pt_3 = iv2ivg(c2v(ic,3));
                set_2D_StiffMat_ele_quad(stiff_mat_ele,ic);
                add_StiffMat_EleQuad(stiff_mat_ele,pt_0,pt_1,pt_2,pt_3);
                break;
            default:
                FatalError("Element type not yet supported for mesh motion - supported types are tris and quads");
                break;
            }
            if (!check) {
                failedIts++;
                if (failedIts > 5) FatalError("ERROR: negative volumes encountered during mesh motion.");
            }else{
                failedIts=0;
            }
        }

        /*--- Compute the tolerance of the linear solver using MinLength ---*/
        NumError = MinLength * 1E-2;

        /*--- Set the boundary displacements (as prescribed by the design variable
        perturbations controlling the surface shape) as a Dirichlet BC. ---*/
        set_boundary_displacements(FlowSol);

        /*--- Fix the location of any points in the domain, if requested. ---*/
        /*
        if (config->GetHold_GridFixed())
            SetDomainDisplacements(FlowSol);
        */

        /*--- Communicate any prescribed boundary displacements via MPI,
        so that all nodes have the same solution and r.h.s. entries
        across all paritions. ---*/
        /// HELP!!! Need Tom/Francisco to decipher what's being sent & how it's used
        //StiffMatrix.SendReceive_Solution(LinSysSol, FlowSol);
        //StiffMatrix.SendReceive_Solution(LinSysRes, FlowSol);

        /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/
        CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(StiffnessMatrix, FlowSol);
        CPreconditioner* precond      = new CLU_SGSPreconditioner(StiffnessMatrix, FlowSol);
        CSysSolve *system             = new CSysSolve();

        /*--- Solve the linear system ---*/
        LinSolIters = system->FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, run_input.solver_tolerance, 100, false);

        /*--- Deallocate memory needed by the Krylov linear solver ---*/
        delete system;
        delete mat_vec;
        delete precond;

        /*--- Update the grid coordinates and cell volumes using the solution
        of the linear system (usol contains the x, y, z displacements). ---*/
        update_grid_coords(FlowSol);

        /*--- Check for failed deformation (negative volumes). ---*/
        min_vol = check_grid(FlowSol);

        /*
        if (FlowSol->rank == 0) {
            cout << "Non-linear iter.: " << iGridDef_Iter << "/" << config->GetGridDef_Iter()
                 << ". Linear iter.: " << LinSolIters << ". Min vol.: " << MinVol
                 << ". Error: " << run_input.solver_tolerance << "." <<endl;
        }
        */
    }

    set_grid_velocity(FlowSol,run_input.dt);
    /*--- Now that deformation is complete & velocity is set, update the
      'official' vertex coordinates ---*/
    xv = xv_new;

    /*--- Deallocate vectors for the linear system. ---*/
    LinSysSol.~CSysVector();
    LinSysRes.~CSysVector();
    StiffnessMatrix.~CSysMatrix();
}

void mesh::set_grid_velocity(solution* FlowSol, double dt)
{
    // calculate velocity using simple backward-Euler
    for (int i=0; i<n_verts; i++) {
        for (int j=0; j<n_dims; j++) {
            vel_new(i,j) = (xv_new(i,j) - xv(i,j))/dt;
        }
    }

    // Apply velocity to the eles classes at the shape points
    int local_ic;
    array<double> vel(n_dims);

    for (int ic=0; ic<n_eles; ic++) {
        for (int j=0; j<c2n_v(i); j++) {
            for (int dim=0; dim<n_dims; dim++) {
                vel[dim] = vel_new(c2v(ic),dim);
            }
            local_ic = ic2loc_c(ic);
            FlowSol->mesh_eles(ctype(ic))->set_grid_vel_spt(local_ic,j,vel);
        }
    }

    // Interpolate grid vel @ spts to fpts
    /// CREATE MATRIX INTERPOLATION OPERATION
    /*for (int i=0; i<FlowSol->n_ele_types; i++) {
        FlowSol->mesh_eles(i)->set_grid_vel_fpts();
    }*/
}

/*! set individual-element stiffness matrix for a triangle */
bool mesh::set_2D_StiffMat_ele_tri(array<double> &stiffMat_ele, int ele_id)
{
    int iPoint;

    int n_spts = c2n_v(ele_id);

    array<double> pos_spts;
    pos_spts.setup(n_spts,n_dims);

    for (int i=0; i<n_spts; i++) {
        iPoint = c2v(ele_id);
        for (int j=0; j<n_dims; j++) {
            pos_spts(i,j) = xv(iPoint,j);
        }
    }

    // ----------- Create single-element stiffness matrix ---------------
    // Copied from SU2
    unsigned short iDim, iVar, jVar, kVar;
    double B_Matrix[6][12], BT_Matrix[12][6], D_Matrix[6][6], Aux_Matrix[12][6];
    double a[3], b[3], c[3], Area, E, Mu, Lambda;

    /*--- Initialize the element stuffness matrix to zero ---*/
    for (iVar = 0; iVar < 6; iVar++)
        for (jVar = 0; jVar < 6; jVar++)
            stiffMat_ele(iVar,jVar) = 0.0;

    for (iDim = 0; iDim < n_dims; iDim++) {
        a[iDim] = pos_spts(0,iDim)-pos_spts(2,iDim);
        b[iDim] = pos_spts(1,iDim)-pos_spts(2,iDim);
    }

    Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

    if (Area < 0.0) {

        /*--- The initial grid has degenerate elements ---*/
        for (iVar = 0; iVar < 6; iVar++) {
            for (jVar = 0; jVar < 6; jVar++) {
                stiffMat_ele(iVar,jVar) = 0.0;
            }
        }
        return false;
    }else{

        /*--- Each element uses their own stiffness which is inversely
        proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
        is a modification to help allow rigid rotation of elements (see
        "Robust Mesh Deformation using the Linear Elasticity Equations" by
        R. P. Dwight. ---*/

        E = 1.0 / Area * fabs(scale); // implement scale = min_vol() function
        Mu = E;
        Lambda = -E;

        a[0] = 0.5 * (pos_spts(1,0)*pos_spts(2,1)-pos_spts(2,0)*pos_spts(1,1) / Area;
        a[1] = 0.5 * (pos_spts(2,0)*pos_spts(0,1)-pos_spts(0,0)*pos_spts(2,1) / Area;
        a[2] = 0.5 * (pos_spts(0,0)*pos_spts(1,1)-pos_spts(1,0)*pos_spts(0,1) / Area;

        b[0] = 0.5 * (pos_spts(1,1)-pos_spts(2,1) / Area;
        b[1] = 0.5 * (pos_spts(2,1)-pos_spts(0,1) / Area;
        b[2] = 0.5 * (pos_spts(0,1)-pos_spts(1,1) / Area;

        c[0] = 0.5 * (pos_spts(2,0)-pos_spts(1,0) / Area;
        c[1] = 0.5 * (pos_spts(0,0)-pos_spts(2,0) / Area;
        c[2] = 0.5 * (pos_spts(1,0)-pos_spts(0,0) / Area;

        /*--- Compute the B Matrix ---*/
        B_Matrix[0][0] = b[0];  B_Matrix[0][1] = 0.0;   B_Matrix[0][2] = b[1];  B_Matrix[0][3] = 0.0;   B_Matrix[0][4] = b[2];  B_Matrix[0][5] = 0.0;
        B_Matrix[1][0] = 0.0;   B_Matrix[1][1] = c[0];  B_Matrix[1][2] = 0.0;   B_Matrix[1][3] = c[1];  B_Matrix[1][4] = 0.0;   B_Matrix[1][5] = c[2];
        B_Matrix[2][0] = c[0];  B_Matrix[2][1] = b[0];  B_Matrix[2][2] = c[1];  B_Matrix[2][3] = b[1];  B_Matrix[2][4] = c[2];  B_Matrix[2][5] = b[2];

        for (iVar = 0; iVar < 3; iVar++)
            for (jVar = 0; jVar < 6; jVar++)
                BT_Matrix[jVar][iVar] = B_Matrix[iVar][jVar];

        /*--- Compute the D Matrix (for plane strain and 3-D)---*/
        D_Matrix[0][0] = Lambda + 2.0*Mu;		D_Matrix[0][1] = Lambda;            D_Matrix[0][2] = 0.0;
        D_Matrix[1][0] = Lambda;            D_Matrix[1][1] = Lambda + 2.0*Mu;   D_Matrix[1][2] = 0.0;
        D_Matrix[2][0] = 0.0;               D_Matrix[2][1] = 0.0;               D_Matrix[2][2] = Mu;

        /*--- Compute the BT.D Matrix ---*/
        for (iVar = 0; iVar < 6; iVar++) {
            for (jVar = 0; jVar < 3; jVar++) {
                Aux_Matrix[iVar][jVar] = 0.0;
                for (kVar = 0; kVar < 3; kVar++)
                    Aux_Matrix[iVar][jVar] += BT_Matrix[iVar][kVar]*D_Matrix[kVar][jVar];
            }
        }

        /*--- Compute the BT.D.B Matrix (stiffness matrix) ---*/
        for (iVar = 0; iVar < 6; iVar++) {
            for (jVar = 0; jVar < 6; jVar++) {
                stiffMat_ele(iVar,jVar) = 0.0;
                for (kVar = 0; kVar < 3; kVar++)
                    stiffMat_ele(iVar,jVar) += Area * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar];
            }
        }

        return true;
    }
}

/*! set individual-element stiffness matrix for a quadrilateral */
bool mesh::set_2D_StiffMat_ele_quad(array<double> &stiffMat_ele,int ele_id, solution *FlowSol) {
    FatalError("ERROR: Sorry, mesh motion on quads not yet implemented.  :( ");
}

/*!
 * Transform element-defined stiffness matrix into node-base stiffness matrix for inclusion
 * into global stiffness matrix 'StiffMatrix'
 */
void mesh::Add_EleTri_StiffMat(array<double> StiffMatrix_Elem, int id_pt_0,
                                               int id_pt_1, int id_pt_2) {
    unsigned short nVar = n_dims;

    // Transform local -> global node ID
    id_pt_0 = iv2ivg(id_pt_0);
    id_pt_1 = iv2ivg(id_pt_1);
    id_pt_2 = iv2ivg(id_pt_1);

    array<double> StiffMatrix_Node;
    StiffMatrix_Node.setup(nVar,nVar);
    StiffMatrix_Node.initialize_to_zero();

    /*--- Transform the stiffness matrix for the triangular element into the
   contributions for the individual nodes relative to each other. ---*/
    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,1);
    StiffMatrix.AddBlock(id_pt_0, id_pt_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,3);
    StiffMatrix.AddBlock(id_pt_0, id_pt_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,5);
    StiffMatrix.AddBlock(id_pt_0, id_pt_2, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,1);
    StiffMatrix.AddBlock(id_pt_1, id_pt_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,3);
    StiffMatrix.AddBlock(id_pt_1, id_pt_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,5);
    StiffMatrix.AddBlock(id_pt_1, id_pt_2, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,1);
    StiffMatrix.AddBlock(id_pt_2, id_pt_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,3);
    StiffMatrix.AddBlock(id_pt_2, id_pt_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,5);
    StiffnessMatrix.AddBlock(id_pt_2, id_pt_2, StiffMatrix_Node);
}

void mesh::update(solution* FlowSol)
{
    // Update element shape points
    if (FlowSol->rank==0) cout << "Deform: updating element shape points" << endl;

    int ele_type, local_id;
    array<double> pos(FlowSol->n_dims);

    for (int ic=0; ic<FlowSol->num_eles; ic++) {
        ele_type = FlowSol->c2ctype_c(ic);
        local_id = ic2loc_c(ic);
        for (int iv=0; iv<c2n_v(ic); iv++) {
            for (int k=0; k<FlowSol->n_dims; k++) {
                pos(k) = xv_new(c2v(ic,iv),k);
            }
            FlowSol->mesh_eles(ele_type)->set_shape_node(iv,local_id,pos);
        }
    }

    if (FlowSol->rank==0) cout << "Deform: done updating elements' shape" << endl;

    // Update element transforms
    if (FlowSol->rank==0) cout << "Deform: updating element transforms ... " << endl;
    for(int i=0;i<FlowSol->n_ele_types;i++) {
        if (FlowSol->mesh_eles(i)->get_n_eles()!=0) {
            FlowSol->mesh_eles(i)->set_transforms(in_run_type);
        }
    }

    // Set metrics at interface cubpts
    if (FlowSol->rank==0) cout << "Deform: setting element transforms at interface cubature points ... " << endl;
    for(int i=0;i<FlowSol->n_ele_types;i++) {
        if (FlowSol->mesh_eles(i)->get_n_eles()!=0) {
            FlowSol->mesh_eles(i)->set_transforms_inters_cubpts();
        }
    }

    // Set metrics at volume cubpts
    if (FlowSol->rank==0) cout << "Deform: setting element transforms at volume cubature points ... " << endl;
    for(int i=0;i<FlowSol->n_ele_types;i++) {
        if (FlowSol->mesh_eles(i)->get_n_eles()!=0) {
            FlowSol->mesh_eles(i)->set_transforms_vol_cubpts();
        }
    }
}

void mesh::update_grid_coords(void)
{
    unsigned short iDim;
    unsigned long iPoint, total_index;
    double new_coord;

    /*--- Update the grid coordinates using the solution of the linear system
   after grid deformation (LinSysSol contains the x, y, z displacements). ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++)
        for (iDim = 0; iDim < n_dims; iDim++) {
            total_index = iPoint*n_dims + iDim;
            new_coord = xv(iPoint,iDim) + LinSysSol[total_index];
            if (fabs(new_coord) < eps*eps) new_coord = 0.0;
            xv_new(iPoint,iDim) = new_coord;
        }
}

double mesh::check_grid(solution* FlowSol) {
    unsigned short iDim;
    unsigned long iElem, ElemCounter = 0;
    double Area, Volume, MinArea = 1E22, MinVolume = 1E22;
    //double MaxArea = -1E22, MaxVolume = -1E22  // never used
    bool NegVol;

    /*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/

    for (iElem = 0; iElem < n_eles; iElem++) {

        /*--- Triangles ---*/
        if (nDim == 2) {
            double a[3], b[3];

            for (iDim = 0; iDim < n_dims; iDim++) {
                a[iDim] = xv(c2v(iElem,0),iDim)-xv(c2v(iElem,1),iDim);
                b[iDim] = xv(c2v(iElem,1),iDim)-xv(c2v(iElem,2),iDim);
            }

            Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

            //MaxArea = max(MaxArea, Area);
            MinArea = min(MinArea, Area);

            NegVol = (MinArea < 0);
        }

        /*--- Tetrahedra ---*/

        if (nDim == 3) {
            double r1[3], r2[3], r3[3], CrossProduct[3];

            for (iDim = 0; iDim < n_dims; iDim++) {
                r1[iDim] = xv(c2v(iElem,1),iDim) - xv(c2v(iElem,0),iDim);
                r2[iDim] = xv(c2v(iElem,2),iDim) - xv(c2v(iElem,0),iDim);
                r3[iDim] = xv(c2v(iElem,3),iDim) - xv(c2v(iElem,0),iDim);
            }

            CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
            CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
            CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];

            Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;

            //MaxVolume = max(MaxVolume, Volume);
            MinVolume = min(MinVolume, Volume);

            NegVol = (MinVolume < 0);
        }

        if (NegVol) ElemCounter++;
    }

#ifdef MPI
    unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
    double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
    double MinVolume_Local = MinVolume; MinVolume = 0.0;

    MPI::COMM_WORLD.Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
    //MPI::COMM_WORLD.Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI::DOUBLE, MPI::MAX);
    MPI::COMM_WORLD.Allreduce(&MinVolume_Local, &MinVolume, 1, MPI::DOUBLE, MPI::MIN);
#endif
    /*
    if ((ElemCounter != 0) && (FlowSol->rank == MASTER_NODE))
        cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;
    */
    if (n_dims == 2) return MinArea;
    else return MinVolume;
}

/*
void mesh::setup_boundaries()
{
    // have vertex -> bcflag
    // setup (bcflags,vertices)
}
*/

void mesh::set_boundary_displacements(solution *FlowSol)
{
    unsigned short iDim, nDim = FlowSol->n_dims, iBound, axis = 0;
    unsigned long iPoint, total_index, iVertex;
    double *VarCoord, MeanCoord[3], VarIncrement = 1.0;

    /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation equations iteratively with
    successive small deformations. ---*/

    VarIncrement = 1.0/((double)run_input.n_deform_iters);

    /*--- As initialization, set to zero displacements of all the surfaces except the symmetry
     plane and the receive boundaries. ---*/

    for (iBound = 0; iBound < n_bnds; iBound++) {
//        my version: if ((bound_flag(ibound) != SYMMETRY_PLANE) && bound_flag(iBound) != MPI_BOUND)) {
            for (iVertex = 0; iVertex < nBndPts(iBound); iVertex++) {
                /// is ivg needed for this?
                iPoint = iv2ivg(boundPts(iBound)(iVertex));
                for (iDim = 0; iDim < n_dims; iDim++) {
                    total_index = iPoint*n_dims + iDim;
                    LinSysRes[total_index] = 0.0;
                    LinSysSol[total_index] = 0.0;
                    StiffMatrix.DeleteValsRowi(total_index);
                }
            }
//        }
    }

    /*--- Set to zero displacements of the normal component for the symmetry plane condition ---*/
    /*for (iBound = 0; iBound < config->GetnMarker_All(); iBound++) {
        if ((config->GetMarker_All_Boundary(iBound) == SYMMETRY_PLANE) && (nDim == 3)) {

            for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;
            for (iVertex = 0; iVertex < geometry->nVertex[iBound]; iVertex++) {
                iPoint = geometry->vertex[iBound][iVertex]->GetNode();
                VarCoord = geometry->node[iPoint]->GetCoord();
                for (iDim = 0; iDim < nDim; iDim++)
                    MeanCoord[iDim] += VarCoord[iDim]*VarCoord[iDim];
            }
            for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);

            if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
            if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
            if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;

            for (iVertex = 0; iVertex < geometry->nVertex[iBound]; iVertex++) {
                iPoint = geometry->vertex[iBound][iVertex]->GetNode();
                total_index = iPoint*nDim + axis;
                LinSysRes[total_index] = 0.0;
                LinSysSol[total_index] = 0.0;
                StiffMatrix.DeleteValsRowi(total_index);
            }
        }
    }*/

    array<double> VarCoord(n_dims);
    VarCoord(0) = 0.5;
    VarCoord(1) = 0.5;

    /*--- Set the known displacements, note that some points of the moving surfaces
    could be on on the symmetry plane, we should specify DeleteValsRowi again (just in case) ---*/
    for (iBound = 0; iBound < n_bnds; iBound++) {
        /*if (((config->GetMarker_All_Moving(iBound) == YES) && (Kind_SU2 == SU2_CFD)) ||
                ((config->GetMarker_All_DV(iBound) == YES) && (Kind_SU2 == SU2_MDC))) */
        if (bound_flags(iBound) == MOTION_ENABLED) {
            for (iVertex = 0; iVertex < nBndPts(iBound); iVertex++) {
                iPoint = boundPts(iBound)(iVertex)->GetNode();
                // get amount which each point is supposed to move at this time step
                // for now, set to a constant (setup data structure(s) later)
                //VarCoord = geometry->vertex[iBound][iVertex]->GetVarCoord();
                for (iDim = 0; iDim < nDim; iDim++) {
                    total_index = iPoint*nDim + iDim;
                    LinSysRes[total_index] = VarCoord(iDim) * VarIncrement;
                    LinSysSol[total_index] = VarCoord(iDim) * VarIncrement;
                    StiffMatrix.DeleteValsRowi(total_index);
                }
            }
        }
    }
}
