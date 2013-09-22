
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2013  Michael Krause <m.krause@tu-harburg.de>

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

#ifndef WGMS3D_COMPLEX_FUNCTIONS_H
#define WGMS3D_COMPLEX_FUNCTIONS_H

/** \file
 *
 * Some helper functions for dealing with complex numbers.
 */

#include <complex>
using std::complex;

template <class T> T sq (T x) {
    return x * x;
}

template <class T> T re (T &x) {
    return x;
}
    
template <class T> T re (complex<T> &x) {
    return x.real();
}
    
template <class T> T im (T &x) {
    return 0.0;
}

template <class T> T im (complex<T> &x) {
    return x.imag();
}
    
#endif // WGMS3D_COMPLEX_FUNCTIONS_H
