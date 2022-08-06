/***************************************************************************
 *   Copyright (C) 2018, 2019 by Marek Sawerwain                            *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
 *   Copyright (C) 2005 -- 2016 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *   Part of the Quantum Computing Simulator:                              *
 *   https://github.com/qMSUZ/QCS                                          *
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
#include "Python.h"
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
	PySys_WriteStdout(
	"/---\\ /---\\ /---\\ \n"
	"|   | |     | 		   Quantum\n"
	"|   | |     \\---\\     Computing System\n"
	"|   | |         | 	\n"
	"\\-\\-/ \\---/ \\---/ \n"
	"   \\\n"
	);

	PySys_WriteStdout("%s\n", version());
	PySys_WriteStdout("%s\n", compile_system());
	PySys_WriteStdout("%s\n", compilator_name());
	qcs_core_library_initialization();
	PySys_WriteStdout("+ qcs gate's cache initialised\n");
	PySys_WriteStdout("+ QCS for Python interface version: %s\n", QCS_INTERFACE_VER);
	PySys_WriteStdout("%s\n", _QCS_I_CompileSystem);
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

	%feature("autodoc", "QuantumRegister()") ~QuantumRegister();
	~QuantumRegister()
	{
		qcs_delete_quantum_register( $self );
	}

	%feature("autodoc", "Reset()");
	void Reset()
	{
	    qcs_quantum_register_reset( $self );
	}

	// method for state set

	%feature("autodoc", "SetGHZState()");
	void SetGHZState()
	{
		qcs_quantum_register_set_ghz_state( $self );
	}

	%feature("autodoc", "SetGHZ01State()");
	void SetGHZ01State()
	{
		qcs_quantum_register_set_ghz_state( $self );
	}
 
	// single gate operations

	%feature("autodoc", "X(int t)");
	void X(int t)
	{
		qcs_quantum_register_pauli_x_gate( $self, t );
	}

	%feature("autodoc", "Y(int t)");
	void Y(int t)
	{
		qcs_quantum_register_pauli_y_gate( $self, t );
	}

    %feature("autodoc", "Z(int t)");
	void Z(int t)
	{
		qcs_quantum_register_pauli_z_gate( $self, t );
	}

	%feature("autodoc", "MXRot90N(int t)");
	void MXRot90N(int t)
	{
		qcs_quantum_register_mx_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "MYRot90N(int t)");
	void MYRot90N(int t)
	{
		qcs_quantum_register_my_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "MZRot90N(int t)");
	void MZRot90N(int t)
	{
		qcs_quantum_register_mz_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "XRot90N(int t)");
	void XRot90N(int t)
	{
		qcs_quantum_register_x_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "YRot90N(int t)");
	void YRot90N(int t)
	{
		qcs_quantum_register_y_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "ZRot90N(int t)");
	void ZRot90N(int t)
	{
		qcs_quantum_register_z_rot90_n_gate( $self, t);
	}

	%feature("autodoc", "Had(int t)");
	void Had(int t)
	{
		qcs_quantum_register_had_n_gate( $self, t );
	}

	%feature("autodoc", "HadN(int t)");
	void HadN(int t)
	{
		qcs_quantum_register_had_n_gate( $self, t );
	}

	%feature("autodoc", "HadAll()");
	void HadAll()
	{
		qcs_quantum_register_had_gate_for_whole_register( $self );
	}

	%feature("autodoc", "SquareRootOfNotN(int t)");
    void SquareRootOfNotN(int t)
	{
		qcs_quantum_register_square_root_not_n_gate($self, t);
	}

	%feature("autodoc", "SqrtOfNotN(int t)");
	void SqrtOfNotN(int t)
	{
		qcs_quantum_register_square_root_not_n_gate($self, t);
	}

	// two qubit/qudit gate operations

	%feature("autodoc", "CNot(int c, int t)");
	void CNot(int c, int t)
	{
		qcs_quantum_register_cnot( $self, c, t);
	}

	// measure operation

	%feature("autodoc", "Measure()");
	int Measure()
	{
#ifdef PYTHON_SCRIPT
		PySys_WriteStdout("Function unimplemented, yet!\n");
#endif
		return -1;
	}

	%feature("autodoc", "M(int t)");
	int M(int t)
	{
#ifdef PYTHON_SCRIPT
		PySys_WriteStdout("Function unimplemented, yet!\n");
#endif
	return -1;
	}

	%feature("autodoc", "MeasureN(int _from, int _to)");
	int MeasureN(int _from, int _to)
	{
		return qcs_quantum_register_measure_from_to( $self, _from, _to);
	}

    %feature("autodoc", "MeasureOneQubit(int t) -> int") MeasureOneQubit(int t);
	int MeasureOneQubit(int t)
	{
		return qcs_quantum_register_measure_one_qubit( $self, t);
	}

	%feature("autodoc", "ProbeQubitStdBase(int i) -> [p0,p1]") ProbeQubitStdBase(int i, float *p0, float *p1);
	%apply float *OUTPUT { float *p0, float *p1 };
	void ProbeQubitStdBase(int i, float *p0, float *p1)
	{
		float _tmp_p0, _tmp_p1;
		qcs_quantum_register_probe_one_qubit_in_std_base( $self, i, &_tmp_p0, &_tmp_p1);
		
		*p0 = _tmp_p0;
		*p1 = _tmp_p1;
	}

	%feature("autodoc", "Noop()") Noop();
	void Noop()
	{
	    return;
	}

	%feature("autodoc", "Pr()");
	void Pr()
	{
	    if( $self->mode == USE_STATE_VECTOR_QUBIT )
			qcs_quantum_register_print_bin( $self );
	}

	%feature("autodoc", "PrSqr()");
	void PrSqr()
	{
	    if( $self->mode == USE_STATE_VECTOR_QUBIT )
			qcs_quantum_register_print_bin_sqr( $self );
	}


	%feature("autodoc", "PrFull()");
	void PrFull()
	{
	    if( $self->mode == USE_STATE_VECTOR_QUBIT )
			qcs_quantum_register_print_bin_full( $self );
	}


	%feature("autodoc", "PrFullSqr()");
	void PrFullSqr()
	{
	    if( $self->mode == USE_STATE_VECTOR_QUBIT )
			qcs_quantum_register_print_bin_full_sqr( $self );
	}

	%feature("autodoc", "PrAsMatlab()");
	void PrAsMatlab()
	{
		qcs_quantum_register_print_bin_in_matlab_format( $self );
	}

	%feature("autodoc", "PrAsMathematica()");
	void PrAsMathematica()
	{
#ifdef PYTHON_SCRIPT
		PySys_WriteStdout("Function unimplemented, yet!\n");
#endif
	}
};


