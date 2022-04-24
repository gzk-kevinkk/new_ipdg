/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: ipdg.c,v 1.162 2022/01/27 05:26:24 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

#ifdef PHG_TO_P4EST /* Note: "make USER_CFLAGS=-DPHG_TO_P4EST ..." */
#include "phg-to-p4est.h"
#endif /* defined PHG_TO_P4EST */

#ifndef INTERFACE_TEST /* test adaptive refinement on an interface */
#define INTERFACE_TEST 0
#endif

#if INTERFACE_TEST /* test adaptive refinement on an interface */
#include "phg/xfem.h"

static void ls_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	x -= .5;
	y -= .5;
	z -= .5;
	*value = x * x + y * y + z * z - .25 * .25;
}

static void ls_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad)
{
	x -= .5;
	y -= .5;
	z -= .5;
	grad[0] = x + x;
	grad[1] = y + y;
	grad[2] = z + z;
}
#endif /* TEST_INTERFACE */

/*---------------------------------------------------------------------------
 * Sample IPDG (interior penalty discontinuous Galerkin) code
 * for solving the following Poisson equation:
 *	-\Delta u = f,		in \Omega
 *		u = g_D,	on \Gamma_D
 *	      u_n = g_N\cdot n,	on \Gamma_N
 * Weak form:
 * 		a(u,v) = f(v), \forall v
 * where:
 * 	a(u,v)	= \int_\Omega \grad u . \grad v
 * 		  +
 *		  \sum_{F\in\Gamma_0\cup\Gamma_D} \int_F (
 *			- {\grad u}.n [v]
 *			- \beta [u]{\grad v}.n
 *			+ \gamma_0 p^2 / h_F [u][v]
 *		  )
 *		  +
 *		  \sum_{F\in\Gamma_0\cup\Gamma_N} \int_F (
 *			\gamma_1 h_F / p^2 [\grad u].n [\grad v].n 
 *		  )
 *	f(v)	= \int_\Omega f v
 *		  +
 *		  \sum_{F\in\Gamma_D} \int_F (
 *			g_D (-\beta(\grad v).n + \gamma_0 p^2 / h_F v)
 *		  )
 *		  +
 *		  \sum_{F\in\Gamma_N} \int_F (
 *			g_N.n (v + \gamma_1 h_F / p^2 (\grad v).n)
 *		  )
 * where F: face, h_F: face diameter,
 * 	 \Gamma_0 is the set of interior faces,
 *	 \Gamma = \Gamma_N \cup \Gamma_D = \partial\Omega,
 *
 * Note:
 * 1. "../test/cube_Neumann_on_top.dat" is a test mesh with Neumann boundary.
 * 2. see "../resistivity/resistivity.c" for an example of Robin BC.
 *---------------------------------------------------------------------------*/

/* The IP parameters */
static FLOAT beta = -1.0, gamma0 = 1.0, gamma1 = 0.1;

static BOOLEAN dump_mat = FALSE, dump_rhs = FALSE;
static INT test_case = 0;

#define PI _F(3.14159265358979323846264338327950288)
#define INTERFACE BDRY_USER5
#define PRINT_BD_INFO 0

static FLOAT k_coef = 1.0;

static void materia_coefficient(int reginoid)
{
	switch (reginoid)
	{
	case 0:
		k_coef = 0.5;
		break;
	case 1:
		k_coef = 1.0;
		break;
	}
}

static int
bc_map(int bctype)
{
	// return bctype;
	switch (bctype)
	{
	case 1:
		return DIRICHLET;
	case 2:
		return NEUMANN;
	// case 4:
	// 	return INTERIOR;
	case 5:
		return INTERFACE;
	default:
		return UNDEFINED;
	}
	return DIRICHLET; /* set Dirichlet BC on all boundaries */
}

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	// *value = -3 * (x + y + z);
	*value = Cos(2. * PI * x) * Cos(2. * PI * y) * Cos(2. * PI * z);
}

static void
func_interface_D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	// *value = x + y + z;
	*value = 3 * Cos(2. * PI * x) * Cos(2. * PI * y) * Cos(2. * PI * z);
}

static void
func_grad_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
	values[0] = -3.0;
	values[1] = -3.0;
	values[2] = -3.0;
	// FLOAT Cx = Cos(2. * PI * x), Cy = Cos(2. * PI * y), Cz = Cos(2. * PI * z);
	// values[0] = -2. * PI * Sin(2. * PI * x) * Cy * Cz;
	// values[1] = -2. * PI * Cx * Sin(2. * PI * y) * Cz;
	// values[2] = -2. * PI * Cx * Cy * Sin(2. * PI * z);
}

static void
func_interface_N(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
	// values[0] = 0.0;
	// values[1] = 0.0;
	// values[2] = 0.0;
	FLOAT Cx = Cos(2. * PI * x), Cy = Cos(2. * PI * y), Cz = Cos(2. * PI * z);
	values[0] = -2. * PI * Sin(2. * PI * x) * Cy * Cz;
	values[1] = -2. * PI * Cx * Sin(2. * PI * y) * Cz;
	values[2] = -2. * PI * Cx * Cy * Sin(2. * PI * z);
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	// *value = 0.0;
	func_u(x, y, z, value);
	*value *= 12. * PI * PI;
}

static void
func_f_lambda(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
	if (e->region_mark == 0)
	{
		values[0] *= 2.0;
	}
	else
	{
		values[0] *= 1.0;
	}
	return;
}

#if 1
/* FIXME: linking against the static library libphg.a, or enabling the line
 * below, makes the code noticebly faster (>2x) */
#include "../src/quad-cache.c"
#endif

static void
do_face(SOLVER *solver, DOF *u_h, QCACHE *qc, int Q_gD, int Q_gN, int Q_ifD, int Q_ifN,
		ELEMENT *e, int face)
{
	GRID *g = u_h->g;
	ELEMENT *e1;
	BYTE face1;
	int n, n1, p; /* p = polynomial order */
	INT I, J;
	int i, j;
	FLOAT val, nv[Dim], *rule;
	FLOAT G0, G1, h, a; /* G0 := gamma0*p^2/h, G1 := gamma1*h/p^2 */
	FLOAT k_i, k_j;
	BTYPE bdry_flag;
	FLOAT dir_coef = 1.0;
	if (e->region_mark == 1)
		dir_coef = -1.0;

	n = DofNBas(u_h, e);
	p = DofTypeOrder(u_h, e);
	h = phgGeomGetFaceDiameter(g, e, face);

	rule = phgQuadGetRule2D(g, e, face, 2 * p);
	phgQCSetRule(qc, rule, -1.);
	phgGeomGetFaceOutNormal(g, e, face, nv);
	phgQCSetConstantNormal(qc, nv);

	if ((e1 = phgGetNeighbour(g, e, face)) == NULL)
	{
		/* boundary face */
		bdry_flag = GetFaceBTYPE(g, e, face);
		// phgPrintf("bdry_flag = %d\n", bdry_flag);
		if (bdry_flag != NEUMANN)
			bdry_flag = DIRICHLET;
		n1 = 0;
		face1 = -1;
		G0 = gamma0 * p * p / h;
		G1 = gamma1 * h / (p * (FLOAT)p);
		materia_coefficient(e->region_mark);
		// phgPrintf("k_coef = %lf\n", k_coef);
		// if (e->region_mark == 0)
		// 	phgError(1, "index = %d", e->index);
		/* RHS */
		for (i = 0; i < n; i++)
		{
			I = phgSolverMapE2G(solver, 0, e, i);
			if (bdry_flag == DIRICHLET)
			{
				/* -\beta\int_\Gamma_D g_D (k \grad v).n */
				a = k_coef * phgQCIntegrateFace(
								 qc, e->index, face, Q_gD, PROJ_NONE, 0,
								 qc, e->index, face, Q_GRAD, PROJ_DOT, i);
				val = -beta * a;
				/* G0 \int_\Gamma_D g_D v */
				a = phgQCIntegrateFace(
					qc, e->index, face, Q_gD, PROJ_NONE, 0,
					qc, e->index, face, Q_BAS, PROJ_NONE, i);
				val += G0 * a;
			}
			else
			{
				/* \int_\Gamma_N g_N v */
				val = phgQCIntegrateFace(
					qc, e->index, face, Q_gN, PROJ_DOT, 0,
					qc, e->index, face, Q_BAS, PROJ_NONE, i);
				/* \int_\Gamma_N G1 g_N (k \grad v).n */
				a = k_coef * phgQCIntegrateFace(
								 qc, e->index, face, Q_gN, PROJ_DOT, 0,
								 qc, e->index, face, Q_GRAD, PROJ_DOT, i);
				val += G1 * a;
			}
			phgSolverAddGlobalRHSEntry(solver, I, val);
		}
	}
	else
	{
		/* interior face shared by "e" and "e1" */
		bdry_flag = GetFaceBTYPE(g, e, face);
		// phgPrintf("bdry_flag = %d\n", bdry_flag);
		face1 = phgOppositeFace(g, e, face, e1);
		n1 = DofNBas(u_h, e1);
		i = DofTypeOrder(u_h, e1);
		if (p < i)
			p = i;
		G0 = gamma0 * p * (FLOAT)p / h;
		G1 = gamma1 * h / (p * (FLOAT)p);
	}

	/* loop on {basis funcs in e} \cup {basis funcs in e1} */
	for (i = 0; i < n + n1; i++)
	{
		SHORT face_i, bas_i;
		ELEMENT *e_i;

		if (i < n)
		{
			e_i = e;
			face_i = face;
			bas_i = i;
		}
		else
		{
			e_i = e1;
			face_i = face1;
			bas_i = i - n;
		}
		I = phgSolverMapE2G(solver, 0, e_i, bas_i);
		materia_coefficient(e_i->region_mark);
		k_i = k_coef;
		/* loop on {basis funcs in e} \cup {basis funcs in e1} */
		for (j = 0; j < n + n1; j++)
		{
			SHORT face_j, bas_j;
			ELEMENT *e_j;

			if (j < n)
			{
				e_j = e;
				face_j = face;
				bas_j = j;
			}
			else
			{
				e_j = e1;
				face_j = face1;
				bas_j = j - n;
			}
			J = phgSolverMapE2G(solver, 0, e_j, bas_j);

			/* skip jumps for interior face and continuous element */
			if (DofFESpace(u_h) == FE_H1 && bdry_flag == INTERIOR)
				continue;

			materia_coefficient(e_j->region_mark);
			k_j = k_coef;

			/*-----------------------------------------------------------------
	     * Note: If the normal vector is reversed, the sign of [.] changes,
	     *	      while the sign of {.} is unaffected.
	     *----------------------------------------------------------------*/

			val = 0.0;

			if (bdry_flag != NEUMANN)
			{
				/* -\int {k_j \grad u}.n [v] (i<=>v, j<=>u, n<=>e_j) */
				a = k_j * phgQCIntegrateFace(
							  qc, e_j->index, face_j, Q_GRAD, PROJ_DOT, bas_j,
							  qc, e_i->index, face_i, Q_BAS, PROJ_NONE, bas_i);
				if (bdry_flag & INTERIOR)
					a *= (e_i == e ? 0.5 : -0.5);
				val = -a;

				/* -\int \beta [u]{k_i \grad v}.n (i<=>v, j<=>u, n<=>e_i) */
				a = k_i * phgQCIntegrateFace(
							  qc, e_j->index, face_j, Q_BAS, PROJ_NONE, bas_j,
							  qc, e_i->index, face_i, Q_GRAD, PROJ_DOT, bas_i);
				if (bdry_flag & INTERIOR)
					a *= (e_j == e ? 0.5 : -0.5);
				val += -beta * a;

				/* \int G0 [u][v] (i<=>v, j<=>u, n<=>e) */
				a = phgQCIntegrateFace(
					qc, e_j->index, face_j, Q_BAS, PROJ_NONE, bas_j,
					qc, e_i->index, face_i, Q_BAS, PROJ_NONE, bas_i);
				if ((bdry_flag & INTERIOR) && e_i != e_j)
					a = -a;
				val += G0 * a;
			}

			if (bdry_flag != DIRICHLET)
			{
				/* \int G1 [k_j \grad u].n [k_i \grad v].n (i<=>v, j<=>u, n<=>e_i) */
				a = k_j * k_i * phgQCIntegrateFace(qc, e_j->index, face_j, Q_GRAD, PROJ_DOT, bas_j, qc, e_i->index, face_i, Q_GRAD, PROJ_DOT, bas_i);
				if ((bdry_flag & INTERIOR) && e_i != e_j)
					a = -a;
				val += G1 * a;
			}

			phgSolverAddGlobalMatrixEntry(solver, I, J, val);
		}

		if (bdry_flag & INTERFACE)
		{
			val = 0.0;

			/* -\beta\int_\Gamma if_D {(k \grad v).n} */
			{
				a = k_i * phgQCIntegrateFace(
								 qc, e_i->index, face_i, Q_ifD, PROJ_NONE, 0,
								 qc, e_i->index, face_i, Q_GRAD, PROJ_DOT, bas_i);
				a *= 0.5;
				val = -beta * a;
			}

			/* G0 \int_\Gamma if_D [v] */
			{
				a = phgQCIntegrateFace(
					qc, e_i->index, face_i, Q_ifD, PROJ_NONE, 0,
					qc, e_i->index, face_i, Q_BAS, PROJ_NONE, bas_i);
				if (e != e_i)
					a = -a;
				val += G0 * a;
			}

			/* \int_\Gamma if_N {v} */
			{
				// a = phgQCIntegrateFace(
				// 	qc, e_i->index, face_i, Q_ifN, PROJ_DOT, 0,
				// 	qc, e_i->index, face_i, Q_BAS, PROJ_NONE, bas_i);
				// a *= 0.5;
				// val += a;
			}

			/* \int_\Gamma G1 if_N [(k \grad v).n] */
			{
				// a = k_i * phgQCIntegrateFace(
				// 				 qc, e_i->index, face_i, Q_ifN, PROJ_DOT, 0,
				// 				 qc, e_i->index, face_i, Q_GRAD, PROJ_DOT, bas_i);
				// if (e != e_i)
				// 	a = -a;
				// val += G1 * a;
			}
			phgSolverAddGlobalRHSEntry(solver, I, dir_coef * val);
		}
	}

	phgFree(rule);
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
	GRID *g = u_h->g;
	ELEMENT *e;
	QCACHE *qc;
	int Q_f, Q_gD, Q_gN, Q_ifD, Q_ifN;

	assert(u_h->dim == 1);

	qc = phgQCNew(QD_PHG, u_h);
	Q_f = phgQCAddFEFunction(qc, f_h);
	// Q_f = phgQCAddXYZFunction(qc, func_f, 1);
	// Q_f1 = phgQCAddXYZFunction(qc, func_f1, 1);
	// Q_f2 = phgQCAddXYZFunction(qc, func_f2, 1);
	Q_gD = phgQCAddXYZFunction(qc, func_u, 1);
	Q_gN = phgQCAddXYZFunction(qc, func_grad_u, 3);
	Q_ifD = phgQCAddXYZFunction(qc, func_interface_D, 1);
	Q_ifN = phgQCAddXYZFunction(qc, func_interface_N, 3);

#if USE_OMP
#pragma omp parallel for private(e)
#endif /* USE_OMP */
	ForAllElementsBegin(g, e)
	{
		materia_coefficient(e->region_mark);
		int N = DofNBas(u_h, e);
		INT I, J, eno = e->index;
		int i, j;
		FLOAT val, *rule;
		/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
		rule = phgQuadGetRule3D(g, e, 2 * DofTypeOrder(u_h, e) /*- 2*/);
		phgQCSetRule(qc, rule, -1.);
		for (i = 0; i < N; i++)
		{
			I = phgSolverMapE2G(solver, 0, e, i);
			for (j = 0; j <= i; j++)
			{
				J = phgSolverMapE2G(solver, 0, e, j);
				/* \int_T \grad u . \grad v */
				val = k_coef * phgQCIntegrate(qc, eno, Q_GRAD, j, qc, eno, Q_GRAD, i);
				phgSolverAddGlobalMatrixEntry(solver, I, J, val);
				if (i != j)
					phgSolverAddGlobalMatrixEntry(solver, J, I, val);
			}
			/* \int_T fv */
			val = phgQCIntegrate(qc, eno, Q_BAS, i, qc, eno, Q_f, 0);
			phgSolverAddGlobalRHSEntry(solver, I, val);
		}
		phgFree(rule);

		/* loop on the faces of the elements */
		for (i = 0; i < NFace; i++)
		{
			ELEMENT *e1 = phgGetNeighbour(g, e, i);
			if (e1 != NULL)
			{
				/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
				if (e->generation < e1->generation)
					continue;
				if (e->generation == e1->generation &&
					GlobalElement(g, e->index) > GlobalElement(g, e1->index))
					continue; /* process each interior face just once */
			}
			do_face(solver, u_h, qc, Q_gD, Q_gN, Q_ifD, Q_ifN, e, i);
		}
	}
	ForAllElementsEnd

		phgQCFree(&qc);
}

int main(int argc, char *argv[])
{
	// char *fn = "../p4est/cube.mesh"; /* hexahedral mesh */
	char *fn = "cube.medit"; /* tetrahedron mesh */
	INT refine0 = 0, refine = 0, refine_step = 1,
		periodicity = 0 /* {1|2|4} (1=x, 2=y, 4=z) */;
	GRID *g;
	DOF *u_h, *u_old, *f_h, *error, *gerror, *gu_h;
	DOF *u_acc, *u_temp, **u_olds;
	SOLVER *solver;
	size_t mem_peak;
	double t, L2norm, H1norm, L2err, H1err;
	BOOLEAN vtk_flag = FALSE, rel_err = FALSE, show_cond = FALSE;
	int level, total_level;
	int index_it = 0;
	int temp_it;
#ifndef PHG_TO_P4EST
	BOOLEAN hp_test = FALSE;
	phgOptionsRegisterNoArg("-hp_test", "Test h-p DOF", &hp_test);
#endif /* PHG_TO_P4EST */

#if INTERFACE_TEST
	XFEM_INFO *xi = phgXFEMInit(NULL, NULL, 0, NULL, 0);
#endif /* INTERFACE_TEST */

	phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
	phgOptionsRegisterInt("-periodicity", "Set periodicity", &periodicity);
	phgOptionsRegisterInt("-refine0", "Initial refinement levels", &refine0);
	phgOptionsRegisterInt("-refine", "Repeated refinement levels", &refine);
	phgOptionsRegisterInt("-refine_step", "Refinement step", &refine_step);
	phgOptionsRegisterFloat("-beta", "The parameter beta", &beta);
	phgOptionsRegisterFloat("-gamma0", "The parameter gamma0", &gamma0);
	phgOptionsRegisterFloat("-gamma1", "The parameter gamma1", &gamma1);
	phgOptionsRegisterInt("-test_case", "Test case", &test_case);

	phgOptionsRegisterNoArg("-rel", "Print relative errors", &rel_err);
	phgOptionsRegisterNoArg("-dump_matrix", "Write matrix to 'A.m'", &dump_mat);
	phgOptionsRegisterNoArg("-dump_rhs", "Write RHS to 'b.m'", &dump_rhs);
	phgOptionsRegisterNoArg("-vtk", "Output 'ipdg.vtk'", &vtk_flag);
	phgOptionsRegisterNoArg("-show_cond", "Estimate condition number",
							&show_cond);

	phgOptionsPreset("-dof_type DG1 -solver gmres");

#if USE_MUMPS
	phgOptionsPreset("-solver mumps -mumps_symmetry unsym");
#endif
#if USE_HYPRE
	/* BoomerAMG seems to work pretty well */
	phgOptionsPreset("-solver hypre -hypre_solver gmres -hypre_pc boomeramg -solver_maxit 2000 -hypre_amg_relax_type l1-gs -solver_rtol 1e-15 -solver_btol 1e-12");
#endif

	phgInit(&argc, &argv);

	g = phgNewGrid(-1);
	// phgSetPeriodicity(g, periodicity);
	phgImportSetBdryMapFunc(bc_map);

	if (!phgImport(g, fn, TRUE))
		phgError(1, "can't read file \"%s\".\n", fn);

	t = phgGetTime(NULL);
	total_level = 0;
	while (refine0 > 0)
	{
		level = refine0 > refine_step ? refine_step : refine0;
		phgRefineAllElements(g, level);
		phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
		refine0 -= level;
		total_level += level;
	}

	/* The discrete solution */
#ifndef PHG_TO_P4EST
	if (hp_test)
	{
		ELEMENT *e;
		HP_TYPE *hp;
		hp = phgHPNew(g, DOF_DEFAULT->fe_space == FE_H1 ? HP_HB : HP_DG);
		ForAllElements(g, e)
			e->hp_order = DOF_DEFAULT->order + GlobalElement(g, e->index) % 3;
		phgHPSetup(hp, FALSE);
		u_h = phgHPDofNew(g, hp, 1, "u_h", DofInterpolation);
		phgHPFree(&hp);
	}
	else
#endif /* PHG_TO_P4EST */
	{
		u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
	}
	phgDofSetDirichletBoundaryMask(u_h, 0);
	phgDofSetDataByValue(u_h, 0.0);
	u_olds = (DOF **)phgAlloc(sizeof(DOF *) * (refine / refine_step + 1));

	while (TRUE)
	{
#if INTERFACE_TEST
		XFEM_INFO *xi = NULL;
		xi = phgXFEMInit(g, ls_func, 2, ls_grad, 2 * u_h->type->order + 3);
		phgXFEMFree(&xi);
#endif /* INTERFACE_TEST */
		phgSetupHalo(g, HALO_FACE);
		phgPrintf("********** Level %d, %d proc%s, %" dFMT " elem%s, LIF %0.2lf, refine time: %0.4lg\n",
				  total_level, g->nprocs, g->nprocs > 1 ? "s" : "",
				  g->nleaf_global, g->nleaf_global > 1 ? "s" : "",
				  (double)g->lif, phgGetTime(NULL) - t);

#if 0
#warning ****************** using 0 as initial solution!
	phgDofSetDataByValue(u_h, 0.0);
#endif

#if 0
	/* Test DOF processing in phgSetupHalo and phgRemoveHalo */
	phgRemoveHalo(g);
	phgSetupHalo(g, HALO_FACE);
#endif

		phgPrintf("Building linear equations: ");
		t = phgGetTime(NULL);
		solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
		/* RHS function */
		f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h", func_f);
		phgDofSetDataByLambdaFunction(f_h, func_f_lambda);
		build_linear_system(solver, u_h, f_h);
		phgDofFree(&f_h);
		phgPrintf("%" dFMT " DOF, build time: %0.4lg\n", DofDataCountGlobal(u_h),
				  phgGetTime(NULL) - t);
		if (show_cond)
			phgPrintf("    Condition number: %0.2le\n",
					  phgMatConditionNumber(solver->mat));
		if (dump_mat)
		{
			phgMatDumpMATLAB(solver->mat, "A", "A.m");
			phgPrintf("    coefficient matrix saved in \"A.m\".\n");
		}
		if (dump_rhs)
		{
			phgVecDumpMATLAB(solver->rhs, "b", "b.m");
			phgPrintf("    RHS saved in \"b.m\".\n");
		}
		phgPrintf("Solving linear equations with %s:\n",
				  solver->oem_solver->name);
		t = phgGetTime(NULL);
		phgSolverSolve(solver, TRUE, u_h, NULL);
		phgMemoryUsage(g, &mem_peak);
		phgPrintf("    nits: %d; residual: %lg; peak mem: %0.2lfMB, "
				  "solve time: %0.4lg\n\n",
				  solver->nits,
				  (double)solver->residual, mem_peak / (1024.0 * 1024.0),
				  phgGetTime(NULL) - t);
		phgSolverDestroy(&solver);
		u_olds[index_it] = phgDofNew(g, u_h->type, 1, "u_old", DofInterpolation);
		phgDofCopy(u_h, &u_olds[index_it], NULL, NULL);
		// phgPrintf("u_h: L2 = %0.5le; u_olds[%d]: L2 = %0.5le;\n\n ",
		// 		  phgDofNormL2(u_h), index_it, phgDofNormL2(u_olds[index_it]));
		index_it++;

#if PRINT_BD_INFO
		ELEMENT *ee;
		int _i;
		int _index = 0;
		BTYPE _bdry_flag;
		ForAllElements(g, ee)
		{
			phgPrintf("index %d:\n", _index++);
			for (_i = 0; _i < NFace; _i++)
			{
				_bdry_flag = GetFaceBTYPE(g, ee, _i);
				phgPrintf("face %d, flag = %d\n", _i, _bdry_flag);
			}
			// phgInfo(-1, "ee->coord = (%lf, %lf, %lf)\n", g->verts[ee->verts[0]][0] , \ 
			// 		g->verts[ee->verts[0]][1] , g->verts[ee->verts[0]][2]);
		}
		return 0;
#endif

		t = phgGetTime(NULL);

		if (refine <= 0)
			break;

		level = refine > refine_step ? refine_step : refine;
		phgRefineAllElements(g, level);
		phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
		refine -= level;
		total_level += level;
	}

	u_acc = phgDofNew(g, DOF_DG3, 1, "u_acc", DofInterpolation);
	phgDofCopy(u_h, &u_acc, NULL, NULL);
	{
		phgPrintf("Building linear equations: ");
		t = phgGetTime(NULL);
		solver = phgSolverCreate(SOLVER_DEFAULT, u_acc, NULL);
		/* RHS function */
		f_h = phgDofNew(g, u_acc->type, 1, "f_h", func_f);
		phgDofSetDataByLambdaFunction(f_h, func_f_lambda);
		build_linear_system(solver, u_acc, f_h);
		phgDofFree(&f_h);
		phgPrintf("%" dFMT " DOF, build time: %0.4lg\n", DofDataCountGlobal(u_acc),
				  phgGetTime(NULL) - t);
		if (show_cond)
			phgPrintf("    Condition number: %0.2le\n",
					  phgMatConditionNumber(solver->mat));
		phgPrintf("Solving linear equations with %s:\n",
				  solver->oem_solver->name);
		t = phgGetTime(NULL);
		phgSolverSolve(solver, TRUE, u_acc, NULL);
		phgMemoryUsage(g, &mem_peak);
		phgPrintf("    nits: %d; residual: %lg; peak mem: %0.2lfMB, "
				  "solve time: %0.4lg\n",
				  solver->nits,
				  (double)solver->residual, mem_peak / (1024.0 * 1024.0),
				  phgGetTime(NULL) - t);
		phgSolverDestroy(&solver);
		t = phgGetTime(NULL);
	}

	temp_it = 0;
	while (temp_it < index_it)
	{
		u_temp = phgDofCopy(u_acc, NULL, NULL, "u_temp");
		// phgPrintf("u_temp: L2 = %0.5le; u_olds[%d]: L2 = %0.5le;\n\n ",
		// 		  phgDofNormL2(u_temp), temp_it, phgDofNormL2(u_olds[temp_it]));
		phgDofAXPY(-1.0, u_olds[temp_it], &u_temp);
		gerror = phgDofGradient(u_temp, NULL, NULL, NULL);
		L2err = phgDofNormL2(u_temp);
		H1err = phgDofNormL2(gerror);
		H1err += L2err;
		phgPrintf("index = %d\n", temp_it);
		phgPrintf("%s errors: L2 = %0.5le; H1 = %0.5le;\n\n ",
				  "Absolute", L2err, H1err);
		phgDofFree(&u_olds[temp_it]);
		phgDofFree(&u_temp);
		phgDofFree(&gerror);
		temp_it++;
	}

	phgDofFree(&u_h);
	phgDofFree(&u_acc);
	phgFreeGrid(&g);

	phgFinalize();

	return 0;
}
