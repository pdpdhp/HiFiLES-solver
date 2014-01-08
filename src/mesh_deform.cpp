// ----------------------------------------------------------------------------
// Helpful SU2 Mesh Deformation Routines
// ----------------------------------------------------------------------------

double CVolumetricMovement::SetFEAMethodContributions_Elem(CGeometry *geometry) {

    unsigned short iVar, iDim;
    unsigned long Point_0, Point_1, Point_2, Point_3, iElem, iEdge, ElemCounter = 0;
    double *Coord_0, *Coord_1, Length, MinLength = 1E10, **StiffMatrix_Elem, Scale;
    double *Edge_Vector = new double [nDim];
    bool RightVol;

    int rank = MASTER_NODE;
#ifndef NO_MPI
    rank = MPI::COMM_WORLD.Get_rank();
#endif

    if (nDim == 2) {
        StiffMatrix_Elem = new double* [6];
        for (iVar = 0; iVar < 6; iVar++)
            StiffMatrix_Elem[iVar] = new double [6];
    }
    if (nDim == 3) {
        StiffMatrix_Elem = new double* [12];
        for (iVar = 0; iVar < 12; iVar++)
            StiffMatrix_Elem[iVar] = new double [12];
    }

    /*--- First, check the minimum edge length in the entire mesh. ---*/

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

        /*--- Points in edge and coordinates ---*/
        Point_0 = geometry->edge[iEdge]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
        Point_1 = geometry->edge[iEdge]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();

        /*--- Compute Edge_Vector ---*/
        Length = 0;
        for (iDim = 0; iDim < nDim; iDim++) {
            Edge_Vector[iDim] = Coord_1[iDim] - Coord_0[iDim];
            Length += Edge_Vector[iDim]*Edge_Vector[iDim];
        }
        Length = sqrt(Length);
        MinLength = min(Length, MinLength);

    }

#ifndef NO_MPI
    double MinLength_Local = MinLength; MinLength = 0.0;
    MPI::COMM_WORLD.Allreduce(&MinLength_Local, &MinLength, 1, MPI::DOUBLE, MPI::MIN);
#endif

    /*--- Second, compute min volume in the entire mesh. ---*/
    Scale = Check_Grid(geometry);

    /*--- Compute contributions from each element by forming the stiffness matrix (FEA) ---*/
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

        if (nDim == 2) {

            /*--- Triangles are loaded directly ---*/
            Point_0 = geometry->elem[iElem]->GetNode(0);
            Point_1 = geometry->elem[iElem]->GetNode(1);
            Point_2 = geometry->elem[iElem]->GetNode(2);
            RightVol = SetFEA_StiffMatrix2D(geometry, StiffMatrix_Elem, Point_0, Point_1, Point_2, Scale);
            AddFEA_StiffMatrix2D(geometry, StiffMatrix_Elem, Point_0, Point_1, Point_2);

        }

        if (nDim == 3) {

            /*--- Tetrahedra are loaded directly ---*/
            Point_0 = geometry->elem[iElem]->GetNode(0);
            Point_1 = geometry->elem[iElem]->GetNode(1);
            Point_2 = geometry->elem[iElem]->GetNode(2);
            Point_3 = geometry->elem[iElem]->GetNode(3);
            RightVol = SetFEA_StiffMatrix3D(geometry, StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3, Scale);
            AddFEA_StiffMatrix3D(geometry, StiffMatrix_Elem, Point_0, Point_1, Point_2, Point_3);

        }

        /*--- Create a list with the degenerated elements ---*/

        if (!RightVol) ElemCounter++;

    }

#ifndef NO_MPI
    unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
    MPI::COMM_WORLD.Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif

    if ((ElemCounter != 0) && (rank == MASTER_NODE))
        cout <<"There are " << ElemCounter << " degenerated elements in the original grid." << endl;

    /*--- Deallocate memory and exit ---*/
    if (nDim == 2) {
        for (iVar = 0; iVar < 6; iVar++)
            delete StiffMatrix_Elem[iVar];
        delete [] StiffMatrix_Elem;
    }
    if (nDim == 3) {
        for (iVar = 0; iVar < 12; iVar++)
            delete StiffMatrix_Elem[iVar];
        delete [] StiffMatrix_Elem;
    }

    delete [] Edge_Vector;
    
    return MinLength;
}



bool CVolumetricMovement::SetFEA_StiffMatrix2D(CGeometry *geometry, double **StiffMatrix_Elem,
                                               unsigned long val_Point_0, unsigned long val_Point_1,
                                               unsigned long val_Point_2, double scale) {
    unsigned short iDim, iVar, jVar, kVar;
    double B_Matrix[6][12], BT_Matrix[12][6], D_Matrix[6][6], Aux_Matrix[12][6];
    double a[3], b[3], c[3], Area, E, Mu, Lambda;

    double *Coord_0 = geometry->node[val_Point_0]->GetCoord();
    double *Coord_1 = geometry->node[val_Point_1]->GetCoord();
    double *Coord_2 = geometry->node[val_Point_2]->GetCoord();

    /*--- Initialize the element stuffness matrix to zero ---*/
    for (iVar = 0; iVar < 6; iVar++)
        for (jVar = 0; jVar < 6; jVar++)
            StiffMatrix_Elem[iVar][jVar] = 0.0;

    for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = Coord_0[iDim]-Coord_2[iDim];
        b[iDim] = Coord_1[iDim]-Coord_2[iDim];
    }

    Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

    if (Area < 0.0) {

        /*--- The initial grid has degenerated elements ---*/

        for (iVar = 0; iVar < 6; iVar++) {
            for (jVar = 0; jVar < 6; jVar++) {
                StiffMatrix_Elem[iVar][jVar] = 0.0;
            }
        }

        return false;

    }
    else {

        /*--- Each element uses their own stiffness which is inversely
     proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
     is a modification to help allow rigid rotation of elements (see
     "Robust Mesh Deformation using the Linear Elasticity Equations" by
     R. P. Dwight. ---*/

        E = 1.0 / Area * fabs(scale);
        Mu = E;
        Lambda = -E;

        a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
        a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
        a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;

        b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
        b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
        b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;

        c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
        c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
        c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;

        /*--- Compute the B Matrix ---*/
        B_Matrix[0][0] = b[0];	B_Matrix[0][1] = 0.0;		B_Matrix[0][2] = b[1];	B_Matrix[0][3] = 0.0;		B_Matrix[0][4] = b[2];	B_Matrix[0][5] = 0.0;
        B_Matrix[1][0] = 0.0;		B_Matrix[1][1] = c[0];	B_Matrix[1][2] = 0.0;		B_Matrix[1][3] = c[1];	B_Matrix[1][4] = 0.0;		B_Matrix[1][5] = c[2];
        B_Matrix[2][0] = c[0];	B_Matrix[2][1] = b[0];	B_Matrix[2][2] = c[1];	B_Matrix[2][3] = b[1];	B_Matrix[2][4] = c[2];	B_Matrix[2][5] = b[2];

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
                StiffMatrix_Elem[iVar][jVar] = 0.0;
                for (kVar = 0; kVar < 3; kVar++)
                    StiffMatrix_Elem[iVar][jVar] += Area * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar];
            }
        }

        return true;

    }

}


void CVolumetricMovement::SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
    unsigned long IterLinSol, iGridDef_Iter;
    double MinLength, NumError, MinVol;

    int rank = MASTER_NODE;
#ifndef NO_MPI
    rank = MPI::COMM_WORLD.Get_rank();
#endif

    /*--- Initialize the number of spatial dimensions, length of the state
   vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

    nDim   = geometry->GetnDim();
    nVar   = geometry->GetnDim();
    nPoint = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();

    /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);

    /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. In
   particular, the linear elasticity equations hold only for small deformations. ---*/

    for (iGridDef_Iter = 0; iGridDef_Iter < config->GetGridDef_Iter(); iGridDef_Iter++) {

        /*--- Initialize vector and sparse matrix ---*/

        LinSysSol.SetValZero();
        LinSysRes.SetValZero();
        StiffMatrix.SetValZero();

        /*--- Compute the stiffness matrix entries for all nodes/elements in the
     mesh. FEA uses a finite element method discretization of the linear
     elasticity equations (transfers element stiffnesses to point-to-point). ---*/

        if (config->GetKind_GridDef_Method() == SPRING) MinLength = SetSpringMethodContributions_Edges(geometry);
        if (config->GetKind_GridDef_Method() == FEA)    MinLength = SetFEAMethodContributions_Elem(geometry);

        /*--- Compute the tolerance of the linear solver using MinLength ---*/

        NumError = MinLength * 1E-2;

        /*--- Set the boundary displacements (as prescribed by the design variable
     perturbations controlling the surface shape) as a Dirichlet BC. ---*/

        SetBoundaryDisplacements(geometry, config);

        /*--- Fix the location of any points in the domain, if requested. ---*/

        if (config->GetHold_GridFixed())
            SetDomainDisplacements(geometry, config);

        /*--- Communicate any prescribed boundary displacements via MPI,
     so that all nodes have the same solution and r.h.s. entries
     across all paritions. ---*/

        StiffMatrix.SendReceive_Solution(LinSysSol, geometry, config);
        StiffMatrix.SendReceive_Solution(LinSysRes, geometry, config);

        /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/

        CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(StiffMatrix, geometry, config);
        CPreconditioner* precond      = new CLU_SGSPreconditioner(StiffMatrix, geometry, config);
        CSysSolve *system             = new CSysSolve();

        /*--- Solve the linear system ---*/

        if (config->GetKind_GridDef_Method() == FEA) IterLinSol = system->FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, 100, false);
        if (config->GetKind_GridDef_Method() == SPRING) IterLinSol = system->ConjugateGradient(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, 100, false);

        /*--- Deallocate memory needed by the Krylov linear solver ---*/

        delete system;
        delete mat_vec;
        delete precond;

        /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/

        UpdateGridCoord(geometry, config);
        if (UpdateGeo)
            UpdateDualGrid(geometry, config);

        /*--- Check for failed deformation (negative volumes). ---*/

        MinVol = Check_Grid(geometry);

        if (rank == MASTER_NODE) {
            cout << "Non-linear iter.: " << iGridDef_Iter << "/" << config->GetGridDef_Iter()
                 << ". Linear iter.: " << IterLinSol << ". Min vol.: " << MinVol
                 << ". Error: " << NumError << "." <<endl;
        }

    }

    /*--- Deallocate vectors for the linear system. ---*/

    LinSysSol.~CSysVector();
    LinSysRes.~CSysVector();
    StiffMatrix.~CSysMatrix();

}

void CVolumetricMovement::AddFEA_StiffMatrix2D(CGeometry *geometry, double **StiffMatrix_Elem, unsigned long val_Point_0,
                                               unsigned long val_Point_1, unsigned long val_Point_2) {
    unsigned short iVar, jVar;
    unsigned short nVar = geometry->GetnDim();

    double **StiffMatrix_Node;
    StiffMatrix_Node = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
        StiffMatrix_Node[iVar] = new double [nVar];

    for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
            StiffMatrix_Node[iVar][jVar] = 0.0;


    /*--- Transform the stiffness matrix for the triangular element into the
   contributions for the individual nodes relative to each other. ---*/
    StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][1];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][1];
    StiffMatrix.AddBlock(val_Point_0, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][3];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][3];
    StiffMatrix.AddBlock(val_Point_0, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[0][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[0][5];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[1][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[1][5];
    StiffMatrix.AddBlock(val_Point_0, val_Point_2, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][1];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][1];
    StiffMatrix.AddBlock(val_Point_1, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][3];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][3];
    StiffMatrix.AddBlock(val_Point_1, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[2][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[2][5];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[3][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[3][5];
    StiffMatrix.AddBlock(val_Point_1, val_Point_2, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][0];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][1];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][0];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][1];
    StiffMatrix.AddBlock(val_Point_2, val_Point_0, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][2];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][3];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][2];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][3];
    StiffMatrix.AddBlock(val_Point_2, val_Point_1, StiffMatrix_Node);

    StiffMatrix_Node[0][0] = StiffMatrix_Elem[4][4];	StiffMatrix_Node[0][1] = StiffMatrix_Elem[4][5];
    StiffMatrix_Node[1][0] = StiffMatrix_Elem[5][4];	StiffMatrix_Node[1][1] = StiffMatrix_Elem[5][5];
    StiffMatrix.AddBlock(val_Point_2, val_Point_2, StiffMatrix_Node);


    /*--- Deallocate memory and exit ---*/
    for (iVar = 0; iVar < nVar; iVar++)
        delete StiffMatrix_Node[iVar];
    delete [] StiffMatrix_Node;

}
