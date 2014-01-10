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

void mesh::deform(struct solution* FlowSol) {
    array<double> stiff_mat_ele;

    if (FlowSol->n_dims == 2) {
        stiff_mat_ele.setup(6,6);
    }else{
        FatalError("3D Mesh motion not implemented yet!");
    }

    int block_i, block_j;
    bool check;

    stiff_mat.setup(n_eles);

    for (int i=0; i<n_eles; i++) {
        switch(ctype(i))
        {
        case TRI:
            set_2D_StiffMat_ele_tri(stiff_mat(i),i,FlowSol);
            break;
        case QUAD:
            set_2D_StiffMat_ele_quad(stiff_mat(i),i,FlowSol);
            break;
        default:
            FatalError("Element type not yet supported for mesh motion - supported types are tris and quads");
            break;
        }
    }

    assemble_stiffness_matrix(); // will use stiff_mat to create StiffnessMatrix
}

bool mesh::set_2D_StiffMat_ele_tri(array<double> &stiffMat_ele,int ele_id, solution *FlowSol) {
    int index;

    int n_spts = c2n_v(ele_id);

    array<double> pos_spts;
    pos_spts.setup(n_spts,n_dims);

    for (int i=0; i<n_spts; i++) {
        index = c2v(ele_id);
        for (int j=0; j<n_dims; j++) {
            pos_spts(i,j) = xv(index,j);
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

        E = 1.0 / Area * fabs(scale);
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

void CVolumetricMovement::AddFEA_StiffMatrix2D(CGeometry *geometry, array<double> StiffMatrix_Elem, unsigned long val_Point_0,
                                               unsigned long val_Point_1, unsigned long val_Point_2) {
    unsigned short iVar, jVar;
    unsigned short nVar = geometry->GetnDim();

    array<double> StiffMatrix_Node;
    StiffMatrix_Node.setup(nVar,nVar);
    StiffMatrix_Node.initialize_to_zero();

    for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
            StiffMatrix_Node[iVar][jVar] = 0.0;


    /*--- Transform the stiffness matrix for the triangular element into the
   contributions for the individual nodes relative to each other. ---*/
    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,1);
    StiffMatrix.AddBlock(val_Point_0, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,3);
    StiffMatrix.AddBlock(val_Point_0, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(0,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(0,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(1,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(1,5);
    StiffMatrix.AddBlock(val_Point_0, val_Point_2, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,1);
    StiffMatrix.AddBlock(val_Point_1, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,3);
    StiffMatrix.AddBlock(val_Point_1, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(2,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(2,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(3,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(3,5);
    StiffMatrix.AddBlock(val_Point_1, val_Point_2, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,0);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,1);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,0);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,1);
    StiffMatrix.AddBlock(val_Point_2, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,2);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,3);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,2);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,3);
    StiffMatrix.AddBlock(val_Point_2, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node(0,0) = StiffMatrix_Elem(4,4);	StiffMatrix_Node(0,1) = StiffMatrix_Elem(4,5);
    StiffMatrix_Node(1,0) = StiffMatrix_Elem(5,4);	StiffMatrix_Node(1,1) = StiffMatrix_Elem(5,5);
    StiffMatrix.AddBlock(val_Point_2, val_Point_2, StiffMatrix_Node);


    /*--- Deallocate memory and exit ---*/
    for (iVar = 0; iVar < nVar; iVar++)
        delete StiffMatrix_Node[iVar];
    delete [] StiffMatrix_Node;

}

/// original try at the top-level structure, just hangin' around for reference
void mesh_deform(struct solution* FlowSol) {
    //for(j=0; j<FlowSol.n_ele_types; j++) FlowSol.mesh_eles(j)->advance_rk11();

    // 1) Build Stiffness Matrix for each element
    //    Will need to limit to triangles for my first implementation
    // 1a) Setup variables

    array<double> stiff_mat_ele;
    if (FlowSol->n_dims == 2) {
        stiff_mat_ele.setup(6,6);
    }else{
        FatalError("3D Mesh motion not implemented yet!");
    }

    /*for(int i=0; i<FlowSol->num_eles; i++) {
        set_2D_StiffMat_ele(stiff_mat_ele,i,FlowSol);
    }*/
	unsigned long i,j;
    bool check;
    int ctype = 0; // just for triangles at the moment
    for (int ele_id=0; ele_id<FlowSol->mesh_eles(ctype)->get_n_eles; ele_id++) {
        check = FlowSol->mesh_eles(ctype)->set_2D_StiffMat_ele(stiff_mat_ele,ele_id);
        StiffnessMatrix.AddBlock(i,j,stiff_mat_ele);
    }

    /*for(int i=0; i<FlowSol->n_ele_types; i++)
        for (int j=0; j<FlowSol->mesh_eles(i)->n_eles; j++)
            set_2D_StiffMat_ele(FlowSol->mesh_eles(i)->setSti*/

    // Build global stiffness matrix from individual elements

    // Iteratively solve linear system

    // Update node positions

    // Re-Set mesh transforms as needed

}

//void set_2D_StiffMat_ele(double **stiffMat_ele, array<double> xv, [individual element / data pertaining to specific element]) {
bool set_2D_StiffMat_ele(array<double> stiffMat_ele,int ele_id, solution *FlowSol) {
    /* need: c2v (cell to index of vertex)
             c2n_v (# of vetices in each cell)
             (coordindates of *physical* vertex) */

    /*
    mesh_file >> id;
    index = index_locate_int(id-1,in_iv2ivg.get_ptr_cpu(),in_n_verts);
        if (index!=-1) { // Vertex belongs to this processor
            for (int m=0;m<FlowSol->n_dims;m++) {
            mesh_file >> pos;
            out_xv(index,m) = pos;
        }
    }
    */

    //n_dims = FlowSol->n_dims;
    //n_spts = FlowSol->ele2n_vert(ele_id);
    //get_ele_local(ele_id,local_id,local_ctype);  // create this function


    // Probably move this up one function (to wrapper function "deform_mesh")
    int global_id, n_spts, n_dims;
    global_id = FlowSol->mesh_eles(i)->get_global_ele(ele_id);
    n_spts = FlowSol->ele2n_vert;

    /*----------------------------------------------------------------------------
    CHANGE OF PLAN -- loop through each ele class (easier to go from local->global
    than from global->local information on each ele)
    Change again - really think this would best be done in each eles_*.cpp class file
    (and length be damned)
    */
    double **pos_spts;
    pos_spts = new double *[n_spts];
    for (int i=0; i<n_spts; i++)
        pos_spts[i] = new double[n_dims];

    for (int i=0; i<FlowSol->n_ele_types; i++)
        for (int j=0; j<FlowSol->mesh_eles(i)->get_n_eles; j++)
            for (int i=0; i<FlowSol->mesh_eles(i)->n_spts_per_ele(j); i++)
                pos_spts[i] = FlowSol->mesh_eles(i)->get_pos_spt(local_id,i);
    //--------------------------------------------------------------------------

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
        a[iDim] = pos_spts[0][iDim]-pos_spts[2][iDim];
        b[iDim] = pos_spts[1][iDim]-pos_spts[2][iDim];
    }

    Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

    if (Area < 0.0) {

        /*--- The initial grid has degenerated elements ---*/
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

        E = 1.0 / Area * fabs(scale);
        Mu = E;
        Lambda = -E;

        a[0] = 0.5 * (pos_spts[1][0]*pos_spts[2][1]-pos_spts[2][0]*pos_spts[1][1]) / Area;
        a[1] = 0.5 * (pos_spts[2][0]*pos_spts[0][1]-pos_spts[0][0]*pos_spts[2][1]) / Area;
        a[2] = 0.5 * (pos_spts[0][0]*pos_spts[1][1]-pos_spts[1][0]*pos_spts[0][1]) / Area;

        b[0] = 0.5 * (pos_spts[1][1]-pos_spts[2][1]) / Area;
        b[1] = 0.5 * (pos_spts[2][1]-pos_spts[0][1]) / Area;
        b[2] = 0.5 * (pos_spts[0][1]-pos_spts[1][1]) / Area;

        c[0] = 0.5 * (pos_spts[2][0]-pos_spts[1][0]) / Area;
        c[1] = 0.5 * (pos_spts[0][0]-pos_spts[2][0]) / Area;
        c[2] = 0.5 * (pos_spts[1][0]-pos_spts[0][0]) / Area;

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

        for (int i=0; i<n_spts; i++)
            delete pos_spts[i];
        delete [] pos_spts;

        return true;
    }
}
