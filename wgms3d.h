
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

#ifndef _WGMS3D_H
#define _WGMS3D_H

#include <vector>

#include "mgp.h"
#include "pml.h"
#include "sparse.h"
#include "fortran_interface.h"
#include "diffops.h"

enum FD_MODE {
    FD_MODE_FULL_VECTORIAL,
    FD_MODE_SEMI_VECTORIAL_HZ, /* assume vertical (z-directed) magnetic field */
    FD_MODE_SEMI_VECTORIAL_HR, /* assume horizontal (rho-directed) magnetic field */
    FD_MODE_SCALAR
};

extern int debugwgms3d, debugmgp;

class wgms3d_mode;
struct pml_spec;

struct wgms3d_simulation_parameters {
    int fd_mode;
    bool use_five_point_standard;
    double k0; /* free-space wave number */
    double c; /* Waveguide curvature c = 1/R */
    PML pml[4]; /* pml_north, pml_east, pml_south, pml_west */

    wgms3d_simulation_parameters (void) {
	fd_mode = FD_MODE_FULL_VECTORIAL;
	use_five_point_standard = false;
	k0 = 2 * M_PI / 1.55e-6;
	c = 0;
    }
};

class wgms3d {
    friend class wgms3d_mode;

  private:

    /* Basic simulation parameters: */
    wgms3d_simulation_parameters sp;

    /* This is just to temporally store what the user specifies using
     * add_pml() until the system matrix is set up. These data are
     * then converted to the final PML objects in 'sp'. */
    std::vector<pml_spec> pml_specs;

    /* Order: left, right, top, bottom.
       0 = Electric wall
       1 = Magnetic wall. */
    int bconds[4];

    /* 0 = grid point lies right on wall,
       1 = two grid points lie symmetrically around wall
       (for pseudo-2D mode) */
    int bcondsym[4];

    /* nir = 'n'umber of 'i'nner 'r'ho-grid points = all points that are
     * not ghost points. */
    int nir, niz;

    /* number of values stored for a single field component: (=nir*niz) */
    int fcsize;
    
    double nmax; /* maximum refractive index occurring in structure */

    /* grid including ghost points: */
    std::vector<double> _rhos, _zs;
    /* _rhos and _zs with complex stretching: */
    std::vector<std::complex<double>> stretched_rhos, stretched_zs;

    /* waveguide geometry: */
    MGP *mgp;

    /* refractive-index distribution, computed from mgp: */
    std::complex<double> *epsis;
    int ldepsis;

    /* for computing and storing finite-difference approximations for
     * differential operators: */
    class Diffops *diffops;

    /* 'retain_list' is an array of size 2*nir*niz, one entry for each
     * discretization point. If retain_list[i] < 0, then point i is
     * known to be zero in advance (Dirichlet BC, or SV/scalar
     * approximation).  Otherwise, retain_list[i] contains a
     * non-negative integer. The non-negative entries of retain_list[]
     * are consecutive integers starting at zero. 'number_of_unknowns'
     * is the number of non-negative entries in retain_list. */
    std::vector<int> retain_list;
    int number_of_unknowns;

    bool complex_calculation;
    void *arpack_evec;
    std::complex<double> *arpack_eval;

    /* 'Ht' stores the transverse H field of the mode specified by
       'active_mode'.  If active_mode < 0, Ht is invalid. */
    std::complex<double> *Ht;
    int active_mode;
    void activate_mode (wgms3d_mode *which);

    std::complex<double> *derived_field_tmp;
    const std::complex<double> * get_er(wgms3d_mode *which);
    const std::complex<double> * get_ez(wgms3d_mode *which);
    const std::complex<double> * get_ep(wgms3d_mode *which);
    const std::complex<double> * get_hp(wgms3d_mode *which);
    sparse_matrix<std::complex<double> > *matrix_Ht_to_Erho;
    sparse_matrix<std::complex<double> > *matrix_Ht_to_Ez;
    sparse_matrix<std::complex<double> > *matrix_Ht_to_Ephi;
    sparse_matrix<std::complex<double> > *matrix_Ht_to_Hphi;
    std::complex<double> *vector_Hz_to_Erho;
    std::complex<double> *vector_Hrho_to_Ez;

    /* -------- functions ---------- */

    sparse_matrix<std::complex<double> > * initmatrix (void);

    void add_matrix_entries (sparse_matrix<std::complex<double> > *A,
			     int to,
			     int Poffset,
			     const std::complex<double> *m);

    /* In a finite-difference expression, replace references to ghost
     * points with references to inner or boundary grid points, applying
     * the proper symmetry conditions. */
    void handle_ghost_points (std::complex<double> *m,
			      int xi,
			      int yi);

    void get_retain_list_sub (int ri,
			      int rinc,
			      int zi,
			      int zinc,
			      int n,
			      int x_or_y);

    void prepare_retain_list (void);
    
    sparse_matrix<std::complex<double> > *
	matrix_remove_zero_points (sparse_matrix<std::complex<double> > *in);

  public:

    wgms3d (void);

    ~wgms3d (void);

    void set_fd_mode (int mode) {
	sp.fd_mode = mode;
    }

    int get_fd_mode (void) {
	return sp.fd_mode;
    }

    void set_geometry (const char *mgp_filename);

    bool is_core_layer_defined (void) {
	assert(mgp != NULL);
	return mgp->is_core_layer_defined();
    }

    void set_grid (std::vector<double> *rho_grid,
		   std::vector<double> *z_grid);

    void set_wavelength (double lambda) {
	sp.k0 = 2.0 * M_PI / lambda;
    }

    void add_pml (char where,
		  int numcells,
		  double sigmamult);

    void set_curvature (double curvature) {
	sp.c = curvature;
    }

    void set_walltype (int which_wall,
		       int type) {
	assert(which_wall >= 0 && which_wall <= 3);
	assert(type == 0 || type == 1);

	bconds[which_wall] = type;
    }

    void set_five_point_standard (bool to) {
	sp.use_five_point_standard = to;
    }

    /* if specified effective index < 0, then use maximum index
     * occurring in structure. */
    bool calculate_modes (int number_of_modes,
			  double near_this_effective_index);

    std::vector<wgms3d_mode> modes;

    const std::complex<double> * get_stretched_rhos();
    const std::complex<double> * get_stretched_zs();
    void get_epsis(const std::complex<double> **epsis, int *ldepsis);

    /* init_derivation() prepares the conversion matrices required for
     * the calculation of the derived fields (Er, Ez, Ep, Hp) from the
     * transverse H field. The argument 'which' is a bitmask that
     * specifies which of the matrices should be prepared. */
    void init_derivation (int which);
};

class wgms3d_mode {
    friend class wgms3d;

  private:
    wgms3d *wg;
    int number;
    int modetype;

  public:
    std::complex<double> beta;

    std::complex<double> get_neff (void) {
	return beta / wg->sp.k0;
    }

    /* Power loss in dB/unit length */
    double get_alpha_per_uol (void) {
	return im(beta) * 20.0 / M_LN10;
    }

    /* Power loss in dB/90° */
    double get_alpha_per_90deg (void) {
	if(wg->sp.c == 0.0) {
	    return 0.0;
	} else {
	    return get_alpha_per_uol() * M_PI / 2.0 * (1.0 / wg->sp.c);
	}
    }

    char estimate_polarization (void);

    /* Add phase factor to the entire mode field such that its phase
     * averaged (in some sense) over the core region (defined by the 'C'
     * directive in the geometry file) is zero. This does not have any
     * physical meaning and is only for convenience when visualizing the
     * mode fields. */
    void adjust_phase (void);

    /* The following functions return pointers to the requested field
     * data. The data is only valid until the next call of one these
     * functions. If more persistence is needed, a copy should be made
     * by the caller. Do not delete[] or free this result. */

    const std::complex<double> * get_hr(void);
    const std::complex<double> * get_hz(void);
    const std::complex<double> * get_hp(void);
    const std::complex<double> * get_er(void);
    const std::complex<double> * get_ez(void);
    const std::complex<double> * get_ep(void);
};

#endif
