
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <unistd.h>
#include <time.h>

#include <complex>
using std::complex;

#include <list>

#include "wgms3d.h"
#include "mgp.h"
#include "stencil.h"
#include "sparse.h"
#include "pml.h"
#include "diffops.h"
#include "fortran_interface.h"

#define MMIN(X,Y) ((X)<(Y)?(X):(Y))
#define MMAX(X,Y) ((X)>(Y)?(X):(Y))

#define Z0 (4e-7*M_PI*299792458) /* free-space impedance */

#define PML_POWER 2

/* the following is just to store the command-line arguments. */
struct pml_spec {
    char where;
    int numcells;
    double sigmamult;
};

wgms3d::wgms3d (void)
{
    nir = -1; niz = -1; fcsize = -1;
    memset(bconds, 0, sizeof(bconds));
    memset(bcondsym, 0, sizeof(bcondsym));

    complex_calculation = false;
    arpack_evec = NULL;
    arpack_eval = NULL;
    Ht = NULL;
    active_mode = -1;

    mgp = NULL;
    epsis = NULL;
    diffops = NULL;

    derived_field_tmp = NULL;
    matrix_Ht_to_Erho = NULL;
    matrix_Ht_to_Ez = NULL;
    matrix_Ht_to_Ephi = NULL;
    matrix_Ht_to_Hphi = NULL;
    vector_Hz_to_Erho = NULL;
    vector_Hrho_to_Ez = NULL;
}

wgms3d::~wgms3d (void)
{
    if(complex_calculation) {
	delete[] (std::complex<double>*)arpack_evec;
    } else {
	delete[] (double*)arpack_evec;
    }
    delete[] arpack_eval;
    delete[] Ht;

    delete mgp;
    delete[] epsis;
    delete diffops;

    delete[] derived_field_tmp;
    delete matrix_Ht_to_Erho;
    delete matrix_Ht_to_Ez;
    delete matrix_Ht_to_Ephi;
    delete matrix_Ht_to_Hphi;
    delete[] vector_Hz_to_Erho;
    delete[] vector_Hrho_to_Ez;
}

void
wgms3d::set_geometry (const char *mgp_filename)
{
    mgp = new MGP(mgp_filename);
}

/* -- GRID + GHOST-POINT HANDLING ---------------------------------------------------------- */

/* The following function replaces references to points 'from' with
 * references to points 'to'. */

template <class T>
static void
handle_ghost_points_sub (T *m,
			 const int *from,
			 const int *to,
			 int Htan,
			 int Hnorm,
			 int sym,
			 int num)
{
    double htan_symtype, hnorm_symtype;

    if(sym == 0) {
	/* electric wall: Hnorm = 0 ("antisymmetric") */
	htan_symtype = +1.0;
	hnorm_symtype = -1.0;
    } else {
	/* magnetic wall: Htan = 0 */
	htan_symtype = -1.0;
	hnorm_symtype = +1.0;
    }

    for(int k = 0; k < num; k++) {
	int s_from = from[k];
	int s_to   = to[k];
	m[s_to + Htan]    += htan_symtype  * m[s_from + Htan];
	m[s_to + Hnorm]   += hnorm_symtype * m[s_from + Hnorm];
	m[s_from + Htan]   = 0.0;
	m[s_from + Hnorm]  = 0.0;
    }
}

#if NUM_GHOST_POINTS == 1
void
wgms3d::handle_ghost_points (complex<double> *m,
			     int ri,
			     int zi)
{
    static const int Wpoints[3] = { 6, 7, 8 };
    static const int Vpoints[3] = { 5, 0, 1 };
    static const int Epoints[3] = { 4, 3, 2 };
    static const int Npoints[3] = { 8, 1, 2 };
    static const int Hpoints[3] = { 7, 0, 3 };
    static const int Spoints[3] = { 6, 5, 4 };

    int Htan, Hnorm;

    /* Vertical boundaries: tangential H is Hy, normal H is Hx */
    Htan = NSP; Hnorm = 0;
    if(ri == 1) {
	handle_ghost_points_sub(m, Wpoints, bcondsym[0] ? Vpoints : Epoints,
				Htan, Hnorm, bconds[0], 3);
    }
    if(ri == nir) {
	handle_ghost_points_sub(m, Epoints, bcondsym[1] ? Vpoints : Wpoints,
				Htan, Hnorm, bconds[1], 3);
    }

    /* Horizontal boundaries: tangential H is Hx, normal H is Hy */
    Htan = 0; Hnorm = NSP;
    if(zi == 1) {
	handle_ghost_points_sub(m, Spoints, bcondsym[3] ? Hpoints : Npoints,
				Htan, Hnorm, bconds[3], 3);
    }
    if(zi == niz) {
	handle_ghost_points_sub(m, Npoints, bcondsym[2] ? Hpoints : Spoints,
				Htan, Hnorm, bconds[2], 3);
    }
}
#endif

void
wgms3d::get_retain_list_sub (int ri,
			     int rinc,
			     int zi,
			     int zinc,
			     int n,
			     int x_or_y)
{
    int offset = ri + zi*nir;
    if(x_or_y == 1) {
	offset += fcsize;
    }
    while(n--) {
	retain_list[offset] = -1;
	offset += rinc + zinc*nir;
    }
}

void
wgms3d::prepare_retain_list (void)
{
    int i;

    assert(fcsize > 0);
    retain_list.resize(2*fcsize, 0);

    /* Eliminate boundary points with Dirichlet BCs. */
    if(bcondsym[0] == 0) {
	/* W wall: Hnorm = Hx, Htan = Hy. */
	get_retain_list_sub(0, 0, 0, 1, niz, bconds[0] == 1);
    }
    if(bcondsym[1] == 0) {
	/* E wall: Hnorm = Hx, Htan = Hy. */
	get_retain_list_sub(nir-1, 0, 0, 1, niz, bconds[1] == 1);
    }

    if(bcondsym[2] == 0) {
	/* N wall: Hnorm = Hy, Htan = Hx. */
	get_retain_list_sub(0, 1, niz-1, 0, nir, bconds[2] == 0);
    }
    if(bcondsym[3] == 0) {
	/* S wall: Hnorm = Hy, Htan = Hx. */
	get_retain_list_sub(0, 1, 0, 0, nir, bconds[3] == 0);
    }

    /* Eliminate points for semi-vectorial calculation. */
    if(sp.fd_mode == FD_MODE_SEMI_VECTORIAL_HZ || sp.fd_mode == FD_MODE_SCALAR) {
	/* only Hz field retained. */
	for(i = 0; i < fcsize; i++) {
	    retain_list[i] = -1;
	}
    } else if(sp.fd_mode == FD_MODE_SEMI_VECTORIAL_HR) {
	/* only Hrho field retained. */
	for(i = 0; i < fcsize; i++) {
	    retain_list[fcsize + i] = -1;
	}
    }

    /* Now give the retained points new numbers. */
    for(i = 0, number_of_unknowns = 0; i < 2*fcsize; i++) {
	if(retain_list[i] == 0) {
	    retain_list[i] = number_of_unknowns++;
	}
    }

    std::cout << "Eliminated "
	      << 2*fcsize - number_of_unknowns
	      << " unknowns with Dirichlet BCs." << std::endl;
}

sparse_matrix<complex<double> > *
wgms3d::matrix_remove_zero_points (sparse_matrix<complex<double> > *in)
{
    sparse_matrix<complex<double> > *out
	= new sparse_matrix<complex<double> >(number_of_unknowns);
    unsigned int inrow, j;

    in->order(1);
    for(inrow = 0; inrow < in->n; inrow++) {
	int outrow = retain_list[inrow];
	if(outrow >= 0) {
	    for(j = in->indextable[inrow]; j < in->indextable[inrow+1]; j++) {
		int incol = in->entries[j].j;
		int outcol = retain_list[incol];
		if(outcol >= 0) {
		    out->add_entry(outrow, outcol, in->entries[j].v);
		}
	    }
	}
    }

    return out;
}

static void
grid_add_ghost_points (std::vector<double> &grid)
{
    double delta;

    assert(grid.size() >= 2);

    delta = grid[1] - grid.front();
    grid.insert(grid.begin(), grid.front() - delta);

    delta = grid.back() - grid[grid.size()-2];
    grid.push_back(grid.back() + delta);
}

void
wgms3d::set_grid (std::vector<double> *rho_grid,
		  std::vector<double> *z_grid)
{
    _rhos = *rho_grid;
    _zs = *z_grid;

    nir = _rhos.size();
    niz = _zs.size();
    fcsize = nir * niz; /* number of values stored for a single field component */

    grid_add_ghost_points(_rhos);
    grid_add_ghost_points(_zs);
}

/* Copy complex-conjugate eigenvector returned by DNEUPD to a complex
 * array, optionally conjugating the vector in the process. */
static void
cpfield (complex<double> *to,
	 double *dneupd_vector,
	 int realpartoffset,
	 int imagpartoffset,
	 int imagpartmult,
         int fcsize)
{
    int i;

    if(!dneupd_vector)
	return;

    for(i = 0; i < fcsize; i++) {
	to[i] = complex<double>(dneupd_vector[i+realpartoffset],
				imagpartmult*dneupd_vector[i+imagpartoffset]);
    }
}

static void
get_complex_ht_sub (int modetype,
                    complex<double> *Ht,
		    double *evec,
		    int evecsize,
		    int n)
{
    switch(modetype) {
    case 0:
	/* first of complex-conjugate pair:
	 * Hxy0[0:fcsize-1] is real part,
	 * Hxy0[evecsize:evecsize+fcsize-1] is imaginary
	 * part */
	cpfield(Ht, evec, 0, evecsize, 1, n);
	break;
    case 1:
	/* second of pair:
	 * Hxy0[-evecsize:-evecsize+fcsize-1] is real
	 * part, Hxy0[0:fcsize-1] is imaginary part (to be
	 * conjugated!) */
	cpfield(Ht, evec, -evecsize, 0, -1, n);
	break;
    case -1:
	/* Real eigenvalue. Just copy over the field. There's
	 * a real part only. */
	int j;
	for(j = 0; j < n; j++) {
	    Ht[j] = evec[j];
	}
	break;
    default:
	assert(1 == 0);
	break;
    }
}

/* Import field of the specified mode into 'Ht', taking care of the
 * ARPACK storage details as well as the unknowns which were
 * previously eliminated from the system since they were known to be
 * zero in advance. */
void
wgms3d::activate_mode (wgms3d_mode *which)
{
    if(which->number == active_mode) {
	return;
    }

    active_mode = which->number;

    if(Ht == NULL) {
	Ht = new std::complex<double>[2*fcsize];
	memset(Ht, 0, 2*fcsize*sizeof(std::complex<double>));
    }

    /* How the eigenvectors are stored by ARPACK depends on whether
     * the system matrix was purely real or a general complex
     * matrix. */

    void *src;
    if(complex_calculation) {
	/* complex calculation */
	src = (std::complex<double>*)arpack_evec + (active_mode * number_of_unknowns);
    } else {
	/* real calculation */
	src = (double*)arpack_evec + (active_mode * number_of_unknowns);
    }

    /* From the eigenvector of unknowns, build the transverse H field
     * taking into account the unknowns that have been deleted from the
     * original user-specified grid. */

    int start, end = -1;
    while(end < 2*fcsize - 1) {
	/* search start of next block to copy */
	for(start = end + 1; start < 2*fcsize; start++) {
	    if(retain_list[start] >= 0) {
		break;
	    }
	}
	if(start == 2*fcsize) {
	    /* there's no next block. */
	    break;
	}

	/* search end of the block */
	for(end = start + 1; end < 2*fcsize; end++) {
	    if(retain_list[end] < 0) {
		break;
	    }
	}

	/* copy data from start to end-1 */
	if(complex_calculation) {
	    /* Complex calculation. Copy over eigenvector (transverse H field)
	     * to Hx/Hy arrays. */
	    int inc = 1;
	    COPY(end-start, (complex<double>*)src + retain_list[start], inc,
		            Ht + start, inc);
	} else {
	    /* Real calculation. Take care of eigenvector storage format. */
	    get_complex_ht_sub(which->modetype,
			       Ht + start, (double*)src + retain_list[start],
			       number_of_unknowns, end - start);
	}
    }
}

const std::complex<double> *
wgms3d::get_stretched_rhos (void)
{
    return &stretched_rhos[1];
}

const std::complex<double> *
wgms3d::get_stretched_zs (void)
{
    return &stretched_zs[1];
}

/* -- CODE FOR SETTING UP THE FINITE-DIFFERENCE SYSTEM MATRIX ---------------------------------- */

void
wgms3d::get_epsis (const std::complex<double> **epsis,
		   int *ldepsis)
{
    *epsis = this->epsis;
    *ldepsis = _rhos.size();
}

void
wgms3d::add_matrix_entries (sparse_matrix<complex<double> > *A,
			    int to,      /* # of equation where we have to add the entries */
			    int Poffset, /* index of central stencil point */
			    const complex<double> *m)
{
    A->add_entry(to, Poffset,       *m++); /* P  */
    A->add_entry(to, Poffset+nir,   *m++); /* N  */
    A->add_entry(to, Poffset+nir+1, *m++); /* NE */
    A->add_entry(to, Poffset    +1, *m++); /* E  */
    A->add_entry(to, Poffset-nir+1, *m++); /* SE */
    A->add_entry(to, Poffset-nir,   *m++); /* S  */
    A->add_entry(to, Poffset-nir-1, *m++); /* SW */
    A->add_entry(to, Poffset    -1, *m++); /* W  */
    A->add_entry(to, Poffset+nir-1, *m); /* NW */
}

void
wgms3d::add_pml (char where,
		 int numcells,
		 double sigmamult)
{
    pml_spec newpml;
    newpml.where = where;
    newpml.numcells = numcells;
    newpml.sigmamult = sigmamult;
    pml_specs.push_back(newpml);
}

static void
init_pml_arrays_sub (const PML *pml,
		     const std::vector<double> &x,
		     int offset,
		     int nx,
		     /* store (1/s), (s'/s^3) and (\tilde\rho) here: */
		     std::vector<complex<double>> &s1,
		     std::vector<complex<double>> &s2,
		     std::vector<complex<double>> &sx)
{
    int i;
    complex<double> s, sp;

    for(i = offset; i <= offset + nx; i++) {
	s = pml->get_s(x[i]);
	sp = pml->get_sprime(x[i]);
	sx[i] = pml->get_stretched_x(x[i]);
	s1[i] = 1.0 / s;
	s2[i] = sp / (s * s * s);
    }
}

sparse_matrix<complex<double> > *
wgms3d::initmatrix (void)
{
    unsigned int i, j, k;

    /* ------------------- Setup PML configuration ----------------- */

    std::vector<complex<double>> s1(_rhos.size(), 1.0);
    std::vector<complex<double>> s2(_rhos.size(), 0.0);
    std::vector<complex<double>> t1(_zs.size(), 1.0);
    std::vector<complex<double>> t2(_zs.size(), 0.0);

    stretched_rhos.resize(_rhos.size());
    for(i = 0; i < _rhos.size(); i++) {
	stretched_rhos[i] = _rhos[i];
    }
    stretched_zs.resize(_zs.size());
    for(i = 0; i < _zs.size(); i++) {
	stretched_zs[i] = _zs[i];
    }

    PML *pml;
    for(pml_spec &p : pml_specs) {
	switch(p.where) {
	case 'n':
	    pml = &sp.pml[0];
	    pml->init(_zs[_zs.size()-1-p.numcells], 1, PML_POWER);
	    pml->set_optimal_strength(_zs[_zs.size()-1-p.numcells+1] - _zs[_zs.size()-1-p.numcells],
				      _zs[_zs.size()-1] - _zs[_zs.size()-1-p.numcells],
				      sp.k0, p.sigmamult);
	    init_pml_arrays_sub(pml, _zs, _zs.size()-1-p.numcells, p.numcells, t1, t2, stretched_zs);
	    break;
	case 'e':
	    pml = &sp.pml[1];
	    pml->init(_rhos[_rhos.size()-1-p.numcells], 1, PML_POWER);
	    pml->set_optimal_strength(_rhos[_rhos.size()-1-p.numcells+1] - _rhos[_rhos.size()-1-p.numcells],
				      _rhos[_rhos.size()-1] - _rhos[_rhos.size()-1-p.numcells],
				      sp.k0, p.sigmamult);
	    init_pml_arrays_sub(pml, _rhos, _rhos.size()-1-p.numcells, p.numcells, s1, s2, stretched_rhos);
	    break;
	case 's':
	    pml = &sp.pml[2];
	    pml->init(_zs[p.numcells], -1, PML_POWER);
	    pml->set_optimal_strength(_zs[p.numcells] - _zs[p.numcells-1],
				      _zs[p.numcells] - _zs[0],
				      sp.k0, p.sigmamult);
	    init_pml_arrays_sub(pml, _zs, 0, p.numcells, t1, t2, stretched_zs);
	    break;
	case 'w':
	    pml = &sp.pml[3];
	    pml->init(_rhos[p.numcells], -1, PML_POWER);
	    pml->set_optimal_strength(_rhos[p.numcells] - _rhos[p.numcells-1],
				      _rhos[p.numcells] - _rhos[0],
				      sp.k0, p.sigmamult);
	    init_pml_arrays_sub(pml, _rhos, 0, p.numcells, s1, s2, stretched_rhos);
	    break;
	}
    }

    /* -------------------------------------------------------------------------- */

    epsis = mgp->get_epsis_on_grid(_rhos, _zs);
    ldepsis = _rhos.size();

    nmax = 1.0;

    std::cout << "Setting up FD system matrix (initial dimension = "
	      << 2*fcsize << ")... " << std::endl;

    sparse_matrix<complex<double> > *A0 = new sparse_matrix<complex<double> >(2*fcsize);

    diffops = new Diffops(mgp, &sp);

    for(j = NUM_GHOST_POINTS; j < _zs.size() - NUM_GHOST_POINTS; j++) {
	const double zp = _zs[j];
	const double n = _zs[j+1]-zp;
	const double s = zp-_zs[j-1];

	for(i = NUM_GHOST_POINTS; i < _rhos.size() - NUM_GHOST_POINTS; i++) {
	    const double rp = _rhos[i];
	    const double e = _rhos[i+1]-rp;
	    const double w = rp-_rhos[i-1];

	    /* number of current grid point */
	    /* The current grid point has the number gpn = '(j-1)*nir
	     * + (i-1)'. The FD equation for the Hrho field component is
	     * stored in matrix row gpn, while the equation for the Hz
	     * field is stored in matrix row gpn+fcsize. */
	    const int gpn = (j-NUM_GHOST_POINTS)*nir + (i-NUM_GHOST_POINTS);

	    complex<double> *M0;

	    complex<double> epsp = epsis[j*ldepsis + i];
	    complex<double> complex_n = sqrt(epsp);
	    if(re(complex_n) > nmax) {
		nmax = re(complex_n);
	    }

	    if(debugmgp) {
		std::cout << "* (" << i << "," << j << "): n = "
			  << complex_n << std::endl;
	    }

	    bool standard;

	    const direction dirs[] = DIRS;
	    M0 = diffops->calculate_diffop(rp, zp, epsp, dirs);
	    if(!M0) {
		/* No interfaces found in this stencil. Use explicit
		 * expressions for FD weights => less "numerical
		 * noise", less non-zero matrix entries. */
		M0 = diffops->get_standard_diffop(n, e, s, w);
		standard = true;
	    } else {
		standard = false;
	    }

	    if(s1[i] != 1.0 || t1[j] != 1.0) {
		/* We're inside a PML. Replace the derivatives by the
		 * "PML-stretched" versions, so that in the remainder
		 * of the program we can simply pretend there is no
		 * PML (using the complex stretched \rho in the
		 * differential equations, though) -- everything else
		 * is absorbed in these matrices. */
		if(standard) {
		    /* the M0 returned by get_standard_diffop() must
		     * not be modified, so we make a backup here. */
		    complex<double> *M0orig = M0;
		    M0 = new complex<double>[(2*NDO) * (2*NSP)];
		    memcpy(M0, M0orig, (2*NDO) * (2*NSP) * sizeof(*M0));
		    standard = false;
		    /* (The memory just allocated for M0 will be freed
		     * in the Diffops destructor.) */
		}
		for(k = 0; k < 2; k++) {
		    /* h^X_\rho\rho  ==>  (1/s^2) h^X_\rho\rho - (s'/s^3) h^X_\rho */
		    SCAL(2*NSP, s1[i]*s1[i], M0 + 2 + k*NDO, 2*NDO);
		    AXPY(2*NSP, -s2[i],      M0 + 0 + k*NDO, 2*NDO,
			                     M0 + 2 + k*NDO, 2*NDO);
		    /* h^X_zz        ==>  (1/t^2) h^X_zz - (t'/t^3) h^X_z */
		    SCAL(2*NSP, t1[j]*t1[j], M0 + 4 + k*NDO, 2*NDO);
		    AXPY(2*NSP, -t2[j],      M0 + 1 + k*NDO, 2*NDO,
			                     M0 + 4 + k*NDO, 2*NDO);
		    /* h^X_\rho      ==>  (1/s) h^X_\rho */
		    SCAL(2*NSP, s1[i],       M0 + 0 + k*NDO, 2*NDO);
		    /* h^X_z         ==>  (1/t) h^X_z */
		    SCAL(2*NSP, t1[j],       M0 + 1 + k*NDO, 2*NDO);
		    /* h^X_\rho z    ==>  (1/(s*t)) h^X_\rho z */
		    SCAL(2*NSP, s1[i]*t1[j], M0 + 3 + k*NDO, 2*NDO);
		}
	    }

	    /* If not standard, save the derivatives for later use in
	     * export_mode(). */
	    if(!standard) {
		diffops->store_diffops(M0, i, j);
	    }

	    /* Now continue to setup system matrix. */
	    const complex<double> onepxc = 1.0 + stretched_rhos[i] * sp.c;
	    const complex<double> onepxc2 = onepxc * onepxc;
	    complex<double> coeffs1[2*NSP];
	    complex<double> coeffs2[2*NSP];
	    memset(coeffs1, 0, sizeof(coeffs1));
	    memset(coeffs2, 0, sizeof(coeffs2));

	    /* Scalar Helmholtz operator: */
	    AXPY(2*NSP, onepxc2, M0 + 2, 2*NDO, coeffs1, 1);
	    AXPY(2*NSP, onepxc2, M0 + 4, 2*NDO, coeffs1, 1);
	    coeffs1[0]   += onepxc2*sp.k0*sp.k0*epsp;

	    AXPY(2*NSP, onepxc2, M0 + NDO+2, 2*NDO, coeffs2, 1);
	    AXPY(2*NSP, onepxc2, M0 + NDO+4, 2*NDO, coeffs2, 1);
	    coeffs2[NSP] += onepxc2*sp.k0*sp.k0*epsp;

	    if(sp.c != 0.0) {
		/* --- Additional curvature terms --- */
		AXPY(2*NSP, 3*sp.c*onepxc, M0 + 0,     2*NDO, coeffs1, 1);
		AXPY(2*NSP, 2*sp.c*onepxc, M0 + NDO+1, 2*NDO, coeffs1, 1);
		coeffs1[0] += sp.c*sp.c;
		AXPY(2*NSP, sp.c*onepxc,   M0 + NDO+0, 2*NDO, coeffs2, 1);
	    }

	    /* Make sure that references to ghost points are replaced
	     * by references to real grid points (either inner grid
	     * points or boundary points). */
	    handle_ghost_points(coeffs1, i, j);
	    handle_ghost_points(coeffs2, i, j);

	    add_matrix_entries(A0, gpn, gpn, coeffs1);
	    add_matrix_entries(A0, gpn, gpn+fcsize, coeffs1 + NSP);
	    add_matrix_entries(A0, gpn+fcsize, gpn, coeffs2);
	    add_matrix_entries(A0, gpn+fcsize, gpn+fcsize, coeffs2 + NSP);
	}
    }

    std::cout << "Stored "
	      << diffops->get_num_stored_diffops() << "/" << fcsize
	      << " non-standard diffops (~"
	      << (diffops->get_num_stored_diffops()*(2*NSP)*(2*NDO)*sizeof(complex<double>))/(1<<20)
	      << "MB)." << std::endl;

    return A0;
}

static inline complex<double>
reduce_complexdouble (const complex<double> &x,
		      const complex<double> *dummy)
{
    return x;
}

static inline double
reduce_complexdouble (const complex<double> &x,
		      const double *dummy)
{
    return x.real();
}

static void
GSTRS (trans_t a, SuperMatrix *b, SuperMatrix *c, int *d, int *e,
       SuperMatrix *B, SuperLUStat_t *f, int *g)
{
    if(B->Dtype == SLU_D) {
	dgstrs(a,b,c,d,e,B,f,g);
    } else if(B->Dtype == SLU_Z) {
	zgstrs(a,b,c,d,e,B,f,g);
    }
}

static void
GSTRF (superlu_options_t *a, SuperMatrix *AC,
       int d, int e, int *f, void *g, int h, int *i, int *j, 
       SuperMatrix *k, SuperMatrix *l, SuperLUStat_t *m, int *n)
{
    if(AC->Dtype == SLU_D) {
	dgstrf(a,AC,d,e,f,g,h,i,j,k,l,m,n);
    } else if(AC->Dtype == SLU_Z) {
	zgstrf(a,AC,d,e,f,g,h,i,j,k,l,m,n);
    }
}

static int     
QuerySpace (SuperMatrix *L, SuperMatrix *U, mem_usage_t *m)
{
    if(L->Dtype == SLU_D) {
	dQuerySpace(L, U, m);
    } else if(L->Dtype == SLU_Z) {
	zQuerySpace(L, U, m);
    }
    return 0;
}

/* superlu(): convert the sparse FD system matrix (always in complex
 * double format) A to SuperLU format (in the format specified by T). */

template <class T>
static bool
superlu (sparse_matrix<complex<double> > *A,
	 SuperMatrix *L,
	 SuperMatrix *U,
	 int **perm_c,
	 int **perm_r,
	 const T *calctype)
{
    SuperLUStat_t stat;
    superlu_options_t options;
    SuperMatrix A_slu;
    int permc_spec;
    int nnz = A->length;
    T *a = NULL;
    int *asub = NULL, *xa = NULL;
    int ap = 0;
    int *etree = NULL;
    SuperMatrix AC;
    int panel_size;
    int relax;
    int lwork = 0;
    int info;
    unsigned int i, j;
    mem_usage_t mem_usage;

    a = new T[nnz];
    asub = new int[nnz];
    xa = new int[A->n+1];
    etree = new int[A->n];

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    /* Convert system matrix to SuperLU format */
    A->order(2);
    ap = 0;
    for(j = 0; j < A->n; j++) {
	xa[j] = ap;
	for(i = A->indextable[j]; i < A->indextable[j+1]; i++) {
	    a[ap] = reduce_complexdouble(A->entries[i].v, calctype);
	    asub[ap] = A->entries[i].i;
	    ap++;
	}
    }
    xa[j] = ap;

    set_default_options(&options);
    options.ColPerm = COLAMD;

    Create_CompCol_Matrix(&A_slu, A->n, A->n, ap, a, asub, xa, SLU_NC, SLU_GE);

    *perm_c = new int[A->n];
    permc_spec = options.ColPerm;
    if (permc_spec != MY_PERMC && options.Fact == DOFACT)
	get_perm_c(permc_spec, &A_slu, *perm_c);

    StatInit(&stat);

    sp_preorder(&options, &A_slu, *perm_c, etree, &AC);

    *perm_r = new int[A->n];

    GSTRF(&options, &AC, relax, panel_size,
	  etree, NULL, lwork, *perm_c, *perm_r, L, U, &stat, &info);
    if(info) {
	std::cerr << "SuperLU Xgstrf failed with INFO = " << info << std::endl;
	exit(1);
    }

    Destroy_SuperMatrix_Store(&A_slu);
    Destroy_CompCol_Permuted(&AC);

    QuerySpace(L, U, &mem_usage);
    std::cout << " (~" << (int)(mem_usage.total_needed/(1<<20)) << "MB)" << std::endl;

    delete[] etree;
    delete[] xa;
    delete[] asub;
    delete[] a;

    return true;
}

template <class T>
static bool
eigensolve (SuperMatrix *L,
	    SuperMatrix *U,
	    double sigma,
	    int *perm_c,
	    int *perm_r,
	    T **evec, /* Store pointer to eigenvector array here (on success) */
	    complex<double> **eval, /* Store pointer to eigenvalue  array here (on success) */
	    int nev) // in eigs.m: 'k'
{
    int ido = 0;
    char bmat = 'I';
    int n = L->nrow; // dimension of system matrix
    char which[] = "LM";
    double tol = 0.0;
    T *resid = new T[2*n]; // n
    int ncv = MMIN(MMAX(2*nev+1,20),n); // from Matlab's eigs.m, there: 'p'
    int ldv = n;
    int iparam[11];
    int ipntr[15];
    T *workd = NULL; 
    int lworkl = 3*ncv*ncv + 6*ncv;
    T *workl = new T[lworkl];
    int info = 0;
    bool rc = false;
    T sigma2 = sigma;

    logical rvec = 1;
    char howmny = 'A';
    logical *select = new logical[ncv];
    T *workev = new T[3*ncv];
    double *rwork = new double[ncv];

    /* SuperLU variables */
    SuperLUStat_t stat;
    int slu_info;

    if(ncv < nev + 2 || ncv > n) {
	std::cerr << "Too many eigenvalues requested." << std::endl;
	exit(1);
    }

    memset(iparam, 0, sizeof(iparam));
    iparam[0] = 1;
    iparam[2] = 10000;
    iparam[6] = 3; // shift-invert mode, M = id

    std::cout << "Eigensolving using ARPACK (nev=" << nev << ", ncv=" << ncv << ")..." << std::endl;

    workd = new T[3*n];
    *evec = new T[ldv*ncv];
    *eval = new complex<double>[nev+1];

    StatInit(&stat);

    while(1) {
	NAUPD(&ido,
	      &bmat, &n, which, &nev, &tol, resid, &ncv, *evec, &ldv,
	      iparam, ipntr, workd, workl, &lworkl, rwork, &info);

	if(ido == -1 || ido == 1) {
	    /* Solve linear system with (A - \sigma I) */
	    T *rhs = workd + ipntr[0] - 1;
	    T *dst = workd + ipntr[1] - 1;

	    trans_t trans = NOTRANS;
	    NCformat Bstore;
	    SuperMatrix B;

	    Bstore.nnz = n;
	    Bstore.nzval = dst;
	    Bstore.rowind = NULL;
	    Bstore.colptr = NULL;
	    B.Stype = SLU_DN;
	    B.Dtype = get_slu_type(rhs);
	    B.Mtype = SLU_GE;
	    B.nrow = n;
	    B.ncol = 1;
	    B.Store = &Bstore;

	    memcpy(dst, rhs, n * sizeof(T));
	    GSTRS(trans, L, U, perm_c, perm_r, &B, &stat, &slu_info);
	    if(slu_info) {
		std::cerr << "SuperLU Xgstrs failed with INFO = "
			  << slu_info << std::endl;
		exit(1);
	    }
	} else {
	    break;
	}
    }

    if(ido != 99 || info != 0) {
	std::cerr << "ARPACK Xnaupd finished with IDO = "
		  << ido << ", INFO = " << info << std::endl;
	exit(1);
    }

    std::cout << "Eigencalculation finished successfully (niter="
	      << iparam[2] << ", nconv=" << iparam[4] << std::endl;

    NEUPD(&rvec, &howmny, select, *eval, *evec, &ldv,
	  sigma2, workev,
	  &bmat, &n, which, &nev, &tol, resid, &ncv, *evec, &ldv,
	  iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    if(info != 0) {
	std::cout << "Error with dneupd, INFO = " << info << std::endl;
	exit(1);
    }

    rc = true;

    delete[] resid;
    delete[] workl;
    delete[] select;
    delete[] workev;
    delete[] rwork;
    delete[] workd;

    return rc;
}

char
wgms3d_mode::estimate_polarization (void)
{
    int i;
    double sum_Hrho = 0.0, sum_Hz = 0.0;

    wg->activate_mode(this);
    const std::complex<double> *H = wg->Ht;

    for(i = 0; i < wg->fcsize; i++) {
	sum_Hrho += abs(*H++);
    }
    for(i = 0; i < wg->fcsize; i++) {
	sum_Hz += abs(*H++);
    }

    if(sum_Hrho > 2*sum_Hz) {
	return 'V';
    } else if(sum_Hz > 2*sum_Hrho) {
	return 'H';
    } else {
	return '?';
    }
}

void
wgms3d_mode::adjust_phase (void)
{
    int i, j;
    double avgreal = 0.0;
    double avgimag = 0.0;
    int n = 0;

    wg->activate_mode(this);
    std::complex<double> *Hrho = wg->Ht;
    std::complex<double> *Hz = Hrho + wg->fcsize;
    int ldHt = wg->nir;

    for(j = 0; j < wg->niz; j++) {
	double zz = wg->_zs[j+NUM_GHOST_POINTS];
	for(i = 0; i < wg->nir; i++) {
	    double rr = wg->_rhos[i+NUM_GHOST_POINTS];
	    if(wg->mgp->is_in_core(rr,zz)) {
		avgreal += re(Hrho[j*ldHt + i]) + re(Hz[j*ldHt + i]);
		avgimag += im(Hrho[j*ldHt + i]) + im(Hz[j*ldHt + i]);
		n++;
	    }
	}
    }

    if(n == 0) {
	return;
    }

    double phase = atan2(avgimag / n, avgreal / n);
    const std::complex<double> factor = exp(std::complex<double>(0.0, -phase));
    SCAL(wg->fcsize, factor, Hrho, 1);
    SCAL(wg->fcsize, factor, Hz, 1);
}

const std::complex<double> *
wgms3d_mode::get_hr (void)
{
    wg->activate_mode(this);
    return wg->Ht;
}

const std::complex<double> *
wgms3d_mode::get_hz (void)
{
    wg->activate_mode(this);
    return wg->Ht + wg->fcsize;
}

void
wgms3d::init_derivation (int which)
{
    static const complex<double> jay(0,1);

    if(which & 1) {
	if(matrix_Ht_to_Erho != NULL && vector_Hz_to_Erho != NULL) {
	    which &= ~1;
	} else {
	    matrix_Ht_to_Erho = new sparse_matrix<complex<double> >(fcsize,2*fcsize);
	    vector_Hz_to_Erho = new complex<double>[fcsize];
	}
    }

    if(which & 2) {
	if(matrix_Ht_to_Ez != NULL && vector_Hrho_to_Ez != NULL) {
	    which &= ~2;
	} else {
	    matrix_Ht_to_Ez = new sparse_matrix<complex<double> >(fcsize,2*fcsize);
	    vector_Hrho_to_Ez = new complex<double>[fcsize];
	}
    }

    if(which & 4) {
	if(matrix_Ht_to_Ephi != NULL) {
	    which &= ~4;
	} else {
	    matrix_Ht_to_Ephi = new sparse_matrix<complex<double> >(fcsize,2*fcsize);
	}
    }

    if(which & 8) {
	if(matrix_Ht_to_Hphi != NULL) {
	    which &= ~8;
	} else {
	    matrix_Ht_to_Hphi = new sparse_matrix<complex<double> >(fcsize,2*fcsize);
	}
    }

    if(!which) {
	return;
    }

    for(unsigned j = NUM_GHOST_POINTS; j < _zs.size() - NUM_GHOST_POINTS; j++) {
	for(unsigned i = NUM_GHOST_POINTS; i < _rhos.size() - NUM_GHOST_POINTS; i++) {
	    const int gpn = (j-NUM_GHOST_POINTS)*nir + (i-NUM_GHOST_POINTS);
	    complex<double> eps = epsis[j*ldepsis + i];
	    complex<double> *M0 = diffops->get_diffops(_rhos, _zs, i, j);
	    const complex<double> onepxc = 1.0 + stretched_rhos[i]*sp.c;
	    const complex<double> Zkn2 = Z0 / (sp.k0 * eps);
	    complex<double> scale;
	    complex<double> coeffs[2*NSP];

	    if(which & 1) {
		vector_Hz_to_Erho[gpn] = -Zkn2 / onepxc; /* h^z */
		memset(coeffs, 0, sizeof(coeffs));
		scale = +Zkn2 * onepxc;
		AXPY(2*NSP, scale, M0 + 3,       2*NDO, coeffs, 1); /* h^\rho_{\rho z} */
		AXPY(2*NSP, scale, M0 + NDO + 4, 2*NDO, coeffs, 1); /* h^z_{z z} */
		scale = +Zkn2 * sp.c;
		AXPY(2*NSP, scale, M0 + 1,       2*NDO, coeffs, 1); /* h^\rho_z */
		handle_ghost_points(coeffs, i, j);
		add_matrix_entries(matrix_Ht_to_Erho, gpn, gpn, coeffs + 0);
		add_matrix_entries(matrix_Ht_to_Erho, gpn, gpn+fcsize, coeffs + NSP);
	    }

	    if(which & 2) {
		vector_Hrho_to_Ez[gpn] = Zkn2 / onepxc; /* h^\rho */
		memset(coeffs, 0, sizeof(coeffs));
		coeffs[0] = -Zkn2*sp.c*sp.c / onepxc; /* h^\rho */
		scale = -Zkn2 * onepxc;
		AXPY(2*NSP, scale, M0 + 2,       2*NDO, coeffs, 1); /* h^\rho_{\rho\rho} */
		AXPY(2*NSP, scale, M0 + NDO + 3, 2*NDO, coeffs, 1); /* h^z_{\rho z} */
		scale = -3.0 * sp.c * Zkn2;
		AXPY(2*NSP, scale, M0 + 0,       2*NDO, coeffs, 1); /* h^\rho_\rho */
		scale = -2.0 * sp.c * Zkn2;
		AXPY(2*NSP, scale, M0 + NDO + 1, 2*NDO, coeffs, 1); /* h^z_z */
		handle_ghost_points(coeffs, i, j);
		add_matrix_entries(matrix_Ht_to_Ez, gpn, gpn, coeffs + 0);
		add_matrix_entries(matrix_Ht_to_Ez, gpn, gpn+fcsize, coeffs + NSP);
	    }

	    if(which & 4) {
		memset(coeffs, 0, sizeof(coeffs));
		scale = jay * Zkn2;
		AXPY(2*NSP, scale, M0 + 1,       2*NDO, coeffs, 1); /* h^\rho_z */
		scale *= -1.0;
		AXPY(2*NSP, scale, M0 + NDO + 0, 2*NDO, coeffs, 1); /* h^z_\rho */
		handle_ghost_points(coeffs, i, j);
		add_matrix_entries(matrix_Ht_to_Ephi, gpn, gpn, coeffs + 0);
		add_matrix_entries(matrix_Ht_to_Ephi, gpn, gpn+fcsize, coeffs + NSP);
	    }

	    if(which & 8) {
		memset(coeffs, 0, sizeof(coeffs));
		coeffs[0] = jay * sp.c; /* h^\rho */
		scale = jay * onepxc;
		AXPY(2*NSP, scale, M0 + 0,       2*NDO, coeffs, 1); /* h^\rho_\rho */
		AXPY(2*NSP, scale, M0 + NDO + 1, 2*NDO, coeffs, 1); /* h^z_z */
		handle_ghost_points(coeffs, i, j);
		add_matrix_entries(matrix_Ht_to_Hphi, gpn, gpn, coeffs + 0);
		add_matrix_entries(matrix_Ht_to_Hphi, gpn, gpn+fcsize, coeffs + NSP);
	    }
	}
    }

    if(derived_field_tmp == NULL) {
	derived_field_tmp = new std::complex<double>[fcsize];
    }
}

static void
AXPY2 (int n,
       complex<double> alpha,
       complex<double> *x1,
       complex<double> *x2,
       complex<double> *y)
{
    while(n--) {
	*y++ += alpha * (*x1++) * (*x2++);
    }
}

const std::complex<double> *
wgms3d::get_er (wgms3d_mode *which)
{
    init_derivation(1);
    matrix_Ht_to_Erho->vecmult(derived_field_tmp, Ht);
    SCAL(fcsize, 1.0/which->beta, derived_field_tmp, 1);
    AXPY2(fcsize, which->beta, vector_Hz_to_Erho, Ht + fcsize, derived_field_tmp);
    return derived_field_tmp;
}

const std::complex<double> *
wgms3d::get_ez (wgms3d_mode *which)
{
    init_derivation(2);
    matrix_Ht_to_Ez->vecmult(derived_field_tmp, Ht);
    SCAL(fcsize, 1.0/which->beta, derived_field_tmp, 1);
    AXPY2(fcsize, which->beta, vector_Hrho_to_Ez, Ht, derived_field_tmp);
    return derived_field_tmp;
}

const std::complex<double> *
wgms3d::get_ep (wgms3d_mode *which)
{
    init_derivation(4);
    matrix_Ht_to_Ephi->vecmult(derived_field_tmp, Ht);
    return derived_field_tmp;
}

const std::complex<double> *
wgms3d::get_hp (wgms3d_mode *which)
{
    init_derivation(8);
    matrix_Ht_to_Hphi->vecmult(derived_field_tmp, Ht);
    SCAL(fcsize, 1.0/which->beta, derived_field_tmp, 1);
    return derived_field_tmp;
}

const std::complex<double> *
wgms3d_mode::get_er (void)
{
    wg->activate_mode(this);
    return wg->get_er(this);
}

const std::complex<double> *
wgms3d_mode::get_ez (void)
{
    wg->activate_mode(this);
    return wg->get_ez(this);
}

const std::complex<double> *
wgms3d_mode::get_ep (void)
{
    wg->activate_mode(this);
    return wg->get_ep(this);
}

const std::complex<double> *
wgms3d_mode::get_hp (void)
{
    wg->activate_mode(this);
    return wg->get_hp(this);
}

bool
wgms3d::calculate_modes (int number_of_modes,
			 double near_this_effective_index)
{
    sparse_matrix<complex<double> > *A, *A0;
    double sigma;
    bool result;
    int modetype;

    /* Prepare the finite-difference system matrix. */
    A0 = initmatrix();

    /* Prepare a list of those field points which are known to be zero
     * in advance.  */
    prepare_retain_list();

    /* Remove those unknowns from the system. (We know the field must
     * be zero there, so this adds unnecessary information to the
     * system, which can result in spurious eigenvalues.) */
    A = matrix_remove_zero_points(A0);
    delete A0;

    std::cout << "Final matrix dimension is "
	      << A->n << "; "
	      << A->length << " non-zero entries." << std::endl;

    /* Prepare for ARPACK's shift-and-invert mode. */
    if(near_this_effective_index < 0.0) {
	near_this_effective_index = nmax;
    }
    std::cout << "Searching for modes near n_eff = "
	      << near_this_effective_index << "." << std::endl;
    sigma = sq(sp.k0*near_this_effective_index);
    A->add_constdiag(-sigma);

    /* Use SuperLU/ARPACK functions optimized for real matrices if
     * possible (faster, less memory consumption). */
    complex_calculation = !A->is_real();
    SuperMatrix L, U;
    int *perm_c, *perm_r;
    
    if(complex_calculation) {
	const complex<double> dummy = 0.0;
	std::cout << "Factorizing FD matrix using (complex) SuperLU..." << std::flush;
	result = superlu(A, &L, &U, &perm_c, &perm_r, &dummy);
    } else {
	const double dummy = 1.0;
	std::cout << "Factorizing FD matrix using (real) SuperLU..." << std::flush;
	result = superlu(A, &L, &U, &perm_c, &perm_r, &dummy);
    }
    delete A;
    if(!result) {
	std::cerr << "Couldn't LU-factorize system matrix." << std::endl;
	goto ende_nosuperlu;
    }

    if(complex_calculation) {
	complex<double> **evec = (complex<double>**)&arpack_evec;
	result = eigensolve(&L, &U, sigma, perm_c, perm_r,
			    evec, &arpack_eval, number_of_modes);
    } else {
	double **evec = (double**)&arpack_evec;
	result = eigensolve(&L, &U, sigma, perm_c, perm_r,
			    evec, &arpack_eval, number_of_modes);
    }
    if(!result) {
	std::cerr << "Couldn't perform eigensolve." << std::endl;
	goto ende_noarpack;
    }

    modes.resize(number_of_modes);

    modetype = 0;
    for(int i = 0; i < number_of_modes; i++) {
	wgms3d_mode &mode = modes[i];

	complex<double> beta2 = arpack_eval[i]; /* This is beta^2 */
 	complex<double> migamma = sqrt(beta2); /* sqrt(beta^2) = (-i * gamma) */

	/* Make sure we have modes with non-negative real part of the effective index
	 * = "forward-travelling waves". */
	if(re(migamma) < 0) {
	    migamma = -migamma; // simply use the other square root...
	}

	mode.wg = this;
	mode.number = i;
	mode.beta = migamma;

	/* the following is meaningful only for the real calculation
	 * modes. this is to categorize the modes as:
	 *  -1 = real beta^2
	 *   0 = complex beta^2, first of the complex-conjugate pair
	 *  +1 = complex beta^2, second of the complex-conjugate pair
	 *
	 * This info is needed in activate_mode().
	 */
	if(im(beta2) == 0.0) {
	    mode.modetype = -1;
	} else {
	    mode.modetype = modetype++;
	    if(modetype == 2) {
		modetype = 0;
	    }
	}
    }

    result = true;

  ende_noarpack:
    delete[] perm_c;
    delete[] perm_r;
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

  ende_nosuperlu:
    return result;
}
