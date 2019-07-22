/***************************************************************************
 *   Copyright (C) 2018,2019 by Marek Sawerwain                            *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
 *   Copyright (C) 2005 -- 2016 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Part of the Quantum Computing Simulator:                              *
 *   http://code.google.com/p/qcs                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

%define DOCSTRING_MODULE
"Quantum Computing Simulator (QCS) port of ANSI C/C++/FORTRAN library of quantum
computations routines for Python and other languages supported by SWIG."
%enddef

%module (docstring=DOCSTRING_MODULE) qcs

%{

#include <complex.h>

#include <math.h>
#include <time.h>
#include <stdarg.h>

#ifdef PYTHON_SCRIPT
#include <Python.h>
#endif


/*
#include "qcs_matrix.h"
#include "qcs_complex.h"
#include "qcs_qubit.h"
#include "qcs_quantum_reg.h"
#include "qcs_gates.h"
#include "qcs_qubit_gates.h"
#include "qcs_info.h"
*/


#include "qcs_quantum_register.h"
#include "qcs_info.h"

#define QCS_INTERFACE_VER "0.2.1"
static const char* _QCS_I_CompileSystem=": compilation date (" __DATE__ " "  __TIME__")";

%}

%include "qcs_quantum_register.h"

#ifdef PYTHON_SCRIPT
%init %{
    initialize_qcs_core_library();
%}
#endif

%extend QuantumRegister {

	%feature("autodoc", "QuantumRegister()") QuantumRegister();
	QuantumRegister()
	{
	    QuantumRegister *qr = NULL;

	    qr = qcs_new_quantum_register( 2 );

	    return qr;
	}

	%feature("autodoc", "QuantumRegister(int size)") QuantumRegister(int qr_size);
	QuantumRegister(int qr_size)
	{
	    QuantumRegister *qr = NULL;

	    qr = qcs_new_quantum_register( qr_size );

	    return qr;
	}

	void Reset()
	{
	    qcs_quantum_register_reset( self );
	}

	void X(int t)
	{
	    qcs_quantum_register_pauli_x(self, t);
	}

	void Y(int t)
	{
	}

	void Z(int t)
	{
	}

	void Had(int t)
	{
	}

	void CNot(int c, int t)
	{
	}

	int M(int t)
	{
	}

	%feature("autodoc", "Noop()") Noop();
	void Noop()
	{
	    return;
	}

	%feature("autodoc", "Pr()");
	void Pr()
	{
	    if(self->mode == USE_STATE_VECTOR_QUBIT)
		qcs_quantum_register_print_bin( self );
	}

	%feature("autodoc", "PrFull()");
	void PrFull()
	{
	    if(self->mode == USE_STATE_VECTOR_QUBIT)
		qcs_quantum_register_print_bin_full( self );
	}

	~QuantumRegister()
	{
	}
};


