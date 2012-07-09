
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

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

#include "wgms3d.h"
#include "fortran_interface.h"

int debugwgms3d = 0;
int debugmgp = 0;

static void
parse_grid_spec (double *&x,
		 int &n,
		 const char *spec)
{
    double p1, p2;
    int i, ret;
    FILE *f;

    if(3 == sscanf(spec, "%lf:%d:%lf", &p1, &n, &p2)) {
	x = new double[n];
	for(i = 0; i < n; i++) {
	    x[i] = p1 + ((p2 - p1) * i) / (n - 1);
	}
    } else {
	f = fopen(spec, "r");
	if(!f) {
	    cerr << "Can't open " << spec << " for reading." << endl;
	    exit(1);
	}

	ret = fscanf(f, "%lf", &p1);
	if(ret == EOF) {
	    cerr << "Error reading " << spec << endl;
	    exit(1);
	}
	n = int(p1);
	x = new double[n];
	for(i = 0; i < n; i++) {
	    ret = fscanf(f, "%lf", &p1);
	    if(ret == EOF) {
		cerr << "Error reading " << spec << endl;
		exit(1);
	    }
	    x[i] = p1;
	}
	fclose(f);
    }
}

static void
write_field (const std::complex<double> *p,
	     int size,
	     const char *fn,
	     int n,
	     std::complex<double> beta)
{
    char buf[256];
    FILE *f;

    sprintf(buf, "%s-%02d.bin", fn, n);
    f = fopen(buf, "w");
    if(!f) {
	std::cerr << "Can't open " << buf << " for writing." << std::endl;
	exit(1);
    }
    fwrite(&beta, sizeof(beta), 1, f);
    fwrite(p, sizeof(p[0]), size, f);
    fclose(f);
}
int
main (int argc,
      char **argv)
{
    int reqargs = 0;
    int num_modes = 2;
    double neff = -1.0;

    double *r0 = NULL, *z0 = NULL;
    int nr, nz;

    double curvature = 0;
    double lambda = 1.55e-6;
    int numpmls = 0;
    int output = 1;
    int write_epsis = 0;
    int write_erho = 0;
    int write_ez = 0;
    int write_ephi = 0;
    int write_hphi = 0;

    cout << "* wgms3d version " VERSION << " *" << endl;

    if(sizeof(int) != 4) {
	cerr << "sizeof(int) != 4!" << endl;
	exit(1);
    }

    wgms3d wg;

    while (1) {
	int oo;
	oo = getopt(argc, argv, "l:n:g:c:s:R:hM:eEFGH5dU:V:P:uvp");
	if(oo == -1)
	    break;

	switch(oo) {
	case 'h':
	    goto usage;
	    break;

	case 'l':
	    lambda = atof(optarg);
	    break;

	case 'n':
	    num_modes = atoi(optarg);
	    break;

	case 'g':
	    wg.set_geometry(optarg);
	    reqargs |= 1;
	    break;

	case 'c':
	    curvature = atof(optarg);
	    break;

	case 'R':
	    curvature = 1.0 / atof(optarg);
	    break;

	case 's':
	    neff = atof(optarg);
	    break;

	case 'M':
	    switch(optarg[0]) {
	    case 'w': wg.set_walltype(0, 1); break;
	    case 'n': wg.set_walltype(2, 1); break;
	    case 'e': wg.set_walltype(1, 1); break;
	    case 's': wg.set_walltype(3, 1); break;
	    }
	    break;

	case '5':
	    wg.set_five_point_standard(true);
	    break;

	case 'd':
	    output = 0;
	    break;

	case 'U':
	    parse_grid_spec(r0, nr, optarg);
	    reqargs |= 2;
	    break;

	case 'V':
	    parse_grid_spec(z0, nz, optarg);
	    reqargs |= 4;
	    break;

	case 'P': {
	    char where;
	    int numcells;
	    double sigmamult;

	    if(3 == sscanf(optarg, "%c:%d:%lf", &where, &numcells, &sigmamult)) {
		if(where == 'n' || where == 'e' || where == 's' || where == 'w') {
		    if(numcells >= 2) {
			wg.add_pml(where, numcells, sigmamult);
			numpmls++;
			break;
		    }
		}
	    }

	    std::cerr << "Malformed PML specification '" << optarg << "'" << std::endl;
	    exit(1);
	    break;
	}

	case 'u':
	    wg.set_fd_mode(FD_MODE_SEMI_VECTORIAL_HZ);
	    break;

	case 'v':
	    wg.set_fd_mode(FD_MODE_SEMI_VECTORIAL_HR);
	    break;

	case 'p':
	    wg.set_fd_mode(FD_MODE_SCALAR);
	    break;

	case 'e':
	    write_epsis = true;
	    break;

	case 'E':
	    write_erho = true;
	    break;

	case 'F':
	    write_ez = true;
	    break;

	case 'G':
	    write_ephi = true;
	    break;

	case 'H':
	    write_hphi = true;
	    break;

	default:
	    exit(1);
	    break;
	}

    }

    if(reqargs != 7) {
      usage:
	printf("Usage:\n");
	printf(" -l <lambda>      Wavelength (1.55e-6)\n");
	printf(" -g <filename>    MGP geometry file (required)\n");
	printf(" -U <desc>        rho-axis grid (required, see below)\n");
	printf(" -V <desc>         z-axis  grid (required, see below)\n");
	printf(" -P <pmlspec>     Add a perfectly matched layer\n");
	printf(" -c <curvature>   Curvature (0)\n");
	printf(" -R <radius>      Radius of curvature (infty)\n");
	printf(" -n <n>           Number of modes to calculate (2)\n");
	printf(" -s <neff>        Return modes near this n_eff (n_max)\n");
	printf(" -M <n|e|s|w>     Use `magnetic' instead of `electric' wall\n");
	printf(" -d               Disable generation of output files\n");
	printf(" -u               Calculate semi-vectorial rho-polarized modes\n");
	printf(" -v               Calculate semi-vectorial z-polarized modes\n");
	printf(" -p               Calculate scalar modes\n");
	printf(" -e               Also export epsis.bin\n");
	printf(" -E               Also export E^rho field\n");
	printf(" -F               Also export E^z field\n");
	printf(" -G               Also export E^phi field\n");
	printf(" -H               Also export H^phi field\n");
	printf(" -5               Use 5-point formulas in hom. regions\n");
	printf("\n");
	printf(" Grid specification <desc> is either:\n");
	printf("    p1:n:p2       N grid points between P1 and P2\n");
	printf(" or:\n");
	printf("    <filename>    N + 1 lines, N on first line\n");
	printf("\n");
	printf(" PML specification <pmlspec> is:\n");
	printf("    <n|e|s|w>:num_cells:sigma_mult\n");
	printf("    To start, set num_cells=10, sigma_mult=1.\n");
	exit(1);
    }

    wg.set_grid(new std::vector<double>(r0, r0+nr),
		new std::vector<double>(z0, z0+nz));
    const int fieldsize = nr * nz;

    wg.set_curvature(curvature);
    std::cout << "Curvature = " << curvature << "/UOL (Radius of curvature = "
	      << 1.0 / curvature << "UOL)" << std::endl;

    if(wg.get_fd_mode() == FD_MODE_SCALAR) {
	/* Some sanity checks for the scalar mode. */
	if(write_erho || write_ez || write_ephi || write_hphi) {
	    std::cerr << "Exportation of derived fields (Er, Ez, Ep, Hp) in scalar mode not meaningful." << std::endl;
	    exit(1);
	}
	if(curvature != 0.0 || numpmls != 0) {
	    std::cerr << std::endl;
	    std::cerr << "WARNING: Scalar mode in conjunction with PMLs or non-zero" << std::endl;
	    std::cerr << "         waveguide curvature not yet verified. Use at your own risk." << std::endl;
	    std::cerr << std::endl;
	}
    }

    wg.set_wavelength(lambda);
    std::cout << "Wavelength = " << lambda << "UOL" << std::endl;
    
    time_t time0 = time(NULL);

    if(!wg.calculate_modes(num_modes, neff)) {
	std::cerr << "Error calculating modes." << std::endl;
	exit(1);
    }

    /* Export results ---------------------------------------------------------- */

    for(int i = 0; i < num_modes; i++) {
	wgms3d_mode &mode = wg.modes[i];
	std::complex<double> neff = mode.get_neff();
	const std::complex<double> *field;

	printf("EV %3d: n_eff = %1.16f + i% .16e\n"
	       "        alpha =% .2edB/UOL [%.2edB/90deg], pol = '%c'\n",
	       i, re(neff), im(neff),
	       mode.get_alpha_per_uol(),
	       mode.get_alpha_per_90deg(),
	       mode.estimate_polarization());

	if(output) {
	    if(wg.is_core_layer_defined()) {
		/* Multiply field with a phase factor such that it is mostly
		 * real in the core region (only meaningful for lossy / leaky
		 * / curved waveguides). */
		mode.adjust_phase();
	    }

	    /* Export transverse H field */
	    field = mode.get_hz();
	    if(wg.get_fd_mode() == FD_MODE_SCALAR) {
		write_field(field, fieldsize, "sc", i, mode.beta);
	    } else {
		write_field(field, fieldsize, "hz", i, mode.beta);
		
		field = mode.get_hr();
		write_field(field, fieldsize, "hr", i, mode.beta);

		/* Compute and export derived fields. */

		/* Speed up computation a little bit by precalculating
		 * the conversion matrices all in one turn: */
		if(i == 0) {
		    wg.init_derivation(   (write_erho ? 1 : 0)
				        | (write_ez   ? 2 : 0)
				        | (write_ephi ? 4 : 0)
				        | (write_hphi ? 8 : 0) );
		}

		if(write_erho) {
		    field = mode.get_er();
		    write_field(field, fieldsize, "er", i, mode.beta);
		}
		if(write_ez) {
		    field = mode.get_ez();
		    write_field(field, fieldsize, "ez", i, mode.beta);
		}
		if(write_ephi) {
		    field = mode.get_ep();
		    write_field(field, fieldsize, "ep", i, mode.beta);
		}
		if(write_hphi) {
		    field = mode.get_hp();
		    write_field(field, fieldsize, "hp", i, mode.beta);
		}
	    }
	}
    }

    if(output) {
	FILE *f;
	const std::complex<double> *data;
	int ldepsis;
	int i, j;

	f = fopen("r.txt", "w");
	if(!f) {
	    std::cerr << "Can't open r.txt for writing." << std::endl;
	    exit(1);
	}
	data = wg.get_stretched_rhos();
	for(i = 0; i < nr; i++) {
	    std::complex<double> point = data[i];
	    fprintf(f, "%.16e %.16e\n", re(point), im(point));
	}
	fclose(f);

	f = fopen("z.txt", "w");
	if(!f) {
	    std::cerr << "Can't open z.txt for writing." << std::endl;
	    exit(1);
	}
	data = wg.get_stretched_zs();
	for(i = 0; i < nz; i++) {
	    std::complex<double> point = data[i];
	    fprintf(f, "%.16e %.16e\n", re(point), im(point));
	}
	fclose(f);

	if(write_epsis) {
	    f = fopen("epsis.bin", "w");
	    if(!f) {
		std::cerr << "Can't open epsis.bin for writing." << std::endl;
		exit(1);
	    }
	    wg.get_epsis(&data, &ldepsis);
	    for(j = 0; j < nz; j++) {
		fwrite(data + j*ldepsis, sizeof(*data), nr, f);
	    }
	    fclose(f);
	}
    }

    /* ------------------------------------------------------------------------- */

    time0 = time(NULL) - time0;
    int min = time0 / 60;
    int sec = time0 - 60*min;
    printf("Total walltime (min:sec) = %d:%02d.\n", min, sec);

    return 0;
}
