/***************************************************************************
 *   Copyright (C) 2005 -- 2013 by Marek Sawerwain                         *
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


#include "qcs_matrix.h"
#include "qcs_complex.h"
#include "qcs_qubit.h"
#include "qcs_qudit.h"
#include "qcs_quantum_reg.h"
#include "qcs_gates.h"
#include "qcs_qubit_gates.h"
#include "qcs_qudit_gates.h"
#include "qcs_info.h"
#include "qcs_chp.h"
#include "qcs_quarray.h"
#include "qcs_qwalk_1d.h"
#include "qcs_qwalk_2d.h"

#define QCS_INTERFACE_VER "0.2.0"
static const char* _QCS_I_CompileSystem=": compilation date (" __DATE__ " "  __TIME__")";

%}

%rename (Fidelity) qcs_fidelity;
%rename (SquareFidelity) qcs_square_of_fidelity;
%rename (SuperFidelity) qcs_super_fidelity;
%rename (TraceDistance) qcs_trace_distance;
%rename (HSDistance) qcs_hilbert_schmidt_distance;
%rename (BuresMetric) qcs_bures_metric;
%rename (AngleMetric) qcs_angle_metric;
%rename (SineMetric) qcs_sine_metric;

%rename (BasicSum) qcs_basic_sum;
%rename (TensorProduct) qcs_tensor_product;

%rename (QUBIT_SYMBOLIC) USE_SYMBOLIC_STATE_VECTOR_QUBIT;
%rename (QUDIT_SYMBOLIC) USE_SYMBOLIC_STATE_VECTOR_QUDIT;
%rename (DENSITY_MATRIX) USE_DENSITY_MATRIX;

%rename (QW1D_Line_Lattice)      QWalk1D_Line_Lattice;
%rename (QW1D_Segment_Lattice)   QWalk1D_Segment_Lattice;
%rename (QW1D_Cycle_Lattice)     QWalk1D_Cycle_Lattice;
%rename (QW2D_Diagonal_Lattice)  QWalk2D_Diagonal_Lattice;
%rename (QW2D_Natural_Lattice)   QWalk2D_Natural_Lattice;
%rename (QW2D_Hexagonal_Lattice) QWalk2D_Hexagonal_Lattice;

%rename (SchmidtDecompositionsCompare) qcs_compare_schmidt_decompositions;

%rename (MixingTwoDensityMatrix) qcs_mixing_two_density_matrix;

%rename (CreateUnitaryRandomRealMatrix) qcs_create_unitary_random_real_matrix;
%rename (CreateUnitaryRandomMatrix) qcs_create_unitary_random_matrix;

%rename (CreateDensityRandomRealMatrix) qcs_create_density_random_real_matrix;
%rename (CreateDensityRandomMatrix) qcs_create_density_random_matrix;

%rename (CreateHermiteanRandomRealMatrix) qcs_create_hermitian_random_real_matrix;
%rename (CreateHermiteanRandomMatrix) qcs_create_hermitian_random_matrix;

%rename (DestroyOperator) qcs_create_destroy_operator;
%rename (CreateOperator) qcs_create_create_operator;

%rename (ReCreateQuantumRegBasedOnGivenDensityMatrix) qcs_recreate_quantum_reg_based_on_given_density_matrix;

%rename (MatrixPauliX) qcs_create_matrix_for_qubit_gate_x;
%rename (MatrixPauliY) qcs_create_matrix_for_qubit_gate_y;
%rename (MatrixPauliZ) qcs_create_matrix_for_qubit_gate_z;
%rename (MatrixHadamard) qcs_create_matrix_for_qubit_gate_hadamard;
%rename (MatrixCNotQubit) cnot_qubit_syntesis_u_matrix_one_control_one_target;
%rename (MatrixCNotQudit) cnot_qudit_syntesis_one_control_one_target;

%rename (CreateUxGateForQGames) qcs_create_matrix_for_ux_1q_gate_float_arg;
%rename (CreateEGateForQGames) qcs_create_matrix_for_e_2q_gate_float_arg;

%rename (CreateXYPerfectSpinTransferUnitrayOp) qcs_create_matrix_of_unitary_operation_of_xy_spin_perfect_transfer_float_arg;

%rename (DisplayTwoRegs) qcs_display_2_quantum_regs;
%rename (DisplayThreeRegs) qcs_display_3_quantum_regs;
%rename (DisplayFourRegs) qcs_display_4_quantum_regs;

%rename (Bin2Dec) qcs_bin2dec;

%rename (CreateThetaOp) get_general_theta_matrix;
%rename (CreateBetaOp) get_general_beta_matrix;
%rename (CreateEtaOp) get_general_eta_matrix;

%include "complex.i"
%include "typemaps.i"

%include "qcs_matrix.h"
%include "qcs_complex.h"
%include "qcs_qubit.h"
%include "qcs_qubit_gates.h"
%include "qcs_qudit.h"
%include "qcs_qudit_gates.h"
%include "qcs_quantum_reg.h"
%include "qcs_info.h"
%include "qcs_quarray.h"
%include "qcs_qwalk_1d.h"
%include "qcs_qwalk_2d.h"


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
 initialize_qcs_core_library();
 PySys_WriteStdout("+ qcs gate's cache initialised\n");
 PySys_WriteStdout("+ QCS for Python interface version: %s\n", QCS_INTERFACE_VER);
 PySys_WriteStdout("%s\n", _QCS_I_CompileSystem);
%}
#endif

#ifdef PYTHON_SCRIPT
%pythoncode %{

def Arange(_s, _e, _step):
	return qcs_create_matrix_arange_operation_with_float_args(_s, _e, _step)

def Linspace(_s, _e, _n, endpoint = True):
	if endpoint == True:
		return qcs_create_matrix_linspace_operation_with_endpoint_with_float_args(_s, _e, _n)
	else:
		return qcs_create_matrix_linspace_operation_without_endpoint_with_float_args(_s, _e, _n)

%}
#endif

%extend QubitBaseDesc {

	%feature("autodoc", "QubitBaseDesc()") QubitBaseDesc();
	QubitBaseDesc()
	{
	    tf_qcs_qubit_base_desc *base;

	    base=malloc(sizeof(tf_qcs_qubit_base_desc));

		qcs_make_std_base(base);
	    
	    return base;
	}

	%feature("autodoc", "QubitBaseDesc(int base_type)") QubitBaseDesc(int base_type);
	QubitBaseDesc(int base_type)
	{
		tf_qcs_qubit_base_desc *base;

		base = malloc( sizeof(tf_qcs_qubit_base_desc) );
		
		base->gate_s0 = NULL;
		base->gate_s1 = NULL;

		if ( base_type == STD_BASE )
		{
			qcs_make_std_base(base);
		}
		
		if ( base_type == HADAMARD_BASE )
		{
			qcs_make_hadamard_base(base);
		}
    	
		if ( base_type == M_BASE )
    	{
    		qcs_make_m_base(base, 0);
		}
    	
		if ( base_type == PAULI_X )
    	{
			qcs_make_pauli_x_base(base);
		}
    	
		if ( base_type == PAULI_Y )
    	{
			qcs_make_pauli_y_base(base);
		}
    	
		if ( base_type == PAULI_Z )
    	{
			qcs_make_pauli_z_base(base);
		}
	
		return base;
	}

	%feature("autodoc", "~QubitBaseDesc()") ~QubitBaseDesc();
	~QubitBaseDesc()
	{
    	qcs_delete_base(self);
    	free(self);
	}

	%feature("autodoc", "SetMBase(float theta)") SetMBase(float theta);
	void SetMBase(float theta)
	{
		qcs_make_m_base(self, theta);
	}

	%feature("autodoc", "SetBBase(float theta)") SetBBase(float theta);
	void SetBBase(float theta)
	{
		qcs_make_b_base(self, theta);
	}
		
	%feature("autodoc", "Pr()") Pr();
	void Pr()
	{
    	qcs_base_print( self );
	}
}

%extend Complex {
    Complex()
    {
        tf_qcs_complex *n=qcs_create_complex( 0.0, 0.0); 
        
        n->re=0;
        n->im=0;
        
        return n;
    }

    Complex(float re, float im)
    {
        tf_qcs_complex *n=qcs_create_complex_float_arg( re, im); 
        
        n->re=re;
        n->im=im;
        
        return n;
    }

  %feature("autodoc", "operator + for float value") __add__(float a);
  Complex __add__(float a)
  {
    Complex t;

    t.re=self->re;
    t.im=self->im;
    
    t.re += a;
    t.im += a;

    return t;
  }

  %feature("autodoc", "__sub__(float a)") __sub__(float a);
  Complex __sub__(float a)
  {
    Complex t;

    t.re=self->re;
    t.im =self->im;

    t.re -= a;
    t.im -= a;
    
    return t;
  }

  %feature("autodoc", "__mul__(float a)") __mul__(float a);
  Complex __mul__(float a)
  {
    Complex t;

    t.re=self->re;
    t.im=self->im;
    
    t.re *= a;
    t.im *= a;

    return t;
  }

  %feature("autodoc", "__div__(float a)") __div__(float a);
  Complex __div__(float a)
  {
    Complex t;

    t.re=self->re;
    t.im=self->im;
    
    t.re /= a;
    t.im /= a;

    return t;
  }

  %feature("autodoc", "operator + for two complex values") __add__(Complex a);
  Complex __add__(Complex a)
  {
    Complex t;

    t.re=0;
    t.im=0;
    
    qcs_complex_add(self, &a, &t);
    
    return t;
  }

  %feature("autodoc", "__sub__(Complex a)") __sub__(Complex a);
  Complex __sub__(Complex a)
  {
    Complex t;

    t.re=0;
    t.im=0;
    
    qcs_complex_sub(self, &a, &t);
    
    return t;
  }

  %feature("autodoc", "__mul__(Complex a)") __mul__(Complex a);
  Complex __mul__(Complex a)
  {
    Complex t;

    t.re=0;
    t.im=0;
    
    qcs_complex_mul(self, &a, &t);
    
    return t;
  }

  %feature("autodoc", "__div__(Complex a)") __div__(Complex a);
  Complex __div__(Complex a)
  {
    Complex t;

    t.re=0;
    t.im=0;
    
    qcs_complex_div(self, &a, &t);
    
    return t;
  }

  %feature("autodoc", "Inv()") Inv();
  Complex Inv()
  {
    Complex t;
    
    qcs_inv_complex(self, &t);
    
    return t;
  }

  %feature("autodoc", "Conj()") Conj();
  Complex Conj()
  {
    Complex t;
    
    qcs_conj_complex(self, &t);
    
    return t;
  }

  %feature("autodoc", "Exp()") Exp();
  Complex Exp()
  {
    Complex t;
    
    qcs_exp_complex(self, &t);
    
    return t;
  }

  %feature("autodoc", "Mod()") Mod();
  float Mod()
  {
    float t;
    
    qcs_mod_complex(self, &t);
    
    return t;
  }

  %feature("autodoc", "Re()") Re();
  float Re()
  {
	return self->re;
  }
  
  %feature("autodoc", "Im()") Im();
  float Im()
  {
	return self->im;
  }
    
  
  %feature("autodoc", "Pr()") Pr();
  void Pr() 
  {
    qcs_print_complex(self);
  }

  %feature("autodoc", "__str__() ") __str__() ;
  char *__str__() 
  {
    static char temp[512];

    sprintf(temp, "%f+%fi", self->re, self->im);

    return temp; 
  }
  
}; /* end of %extend complex */

%extend SchmidtDecomposition {

	SchmidtDecomposition()
	{
		tf_qcs_schmidt_decomposition *sd;
		
		sd = qcs_create_schmidt_decomposition_empty();
		
		return sd;
	}

	~SchmidtDecomposition()
	{
		qcs_delete_schmidt_decomposition(self);
	}
	
	void PrBase1()
	{
		qcs_print_matrix(self->base1);
	}
	
	void PrBase2()
	{
		qcs_print_matrix(self->base2);
	}
	
	void PrCoeff()
	{
		qcs_print_matrix(self->schmidt_coeff);
	}	
}

%extend SpectralDecomposition {

	SpectralDecomposition()
	{
		tf_qcs_spectral_decomposition *sd;
		
		sd = qcs_create_spectral_decomposition();
		
		return sd;
	}

	~SpectralDecomposition()
	{
		qcs_delete_spectral_decomposition(self);
	}
	
	void PrBase()
	{
		qcs_print_matrix(self->eigenvectors);
	}
	
	void PrEigenvalues()
	{
		qcs_print_matrix(self->eigenvalues);
	}	
}

%extend SVDDecomposition {

	SVDDecomposition()
	{
		tf_qcs_svd_decomposition *sd;
		
		sd = qcs_create_svd_decomposition();
		
		return sd;
	}

	~SVDDecomposition()
	{
		qcs_delete_svd_decomposition(self);
	}
	
	void PrU()
	{
		qcs_print_matrix(self->U);
	}

	void PrS()
	{
		qcs_print_matrix(self->S);
	}

	void PrV()
	{
		qcs_print_matrix(self->V);
	}
	
}

%extend Matrix {
	%feature("autodoc", "Matrix() -- create empty matrix without rows and cols") Matrix();
	Matrix()
	{
		tf_qcs_matrix *tqm;

		tqm=qcs_create_empty_matrix();
	
#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
		PySys_WriteStdout("Matrix()\n");
#endif    
		return tqm;
	}

	%feature("autodoc", "Matrix(int rows, int cols)") Matrix(int _rows, int _cols);
	Matrix(int _rows, int _cols)
	{
		tf_qcs_matrix *tqm;

		tqm=qcs_create_matrix(_rows, _cols);
#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
		PySys_WriteStdout("Matrix(int, int)\n");
#endif    
		return tqm;
	}

	%feature("autodoc", "~Matrix()") ~Matrix();
	~Matrix()
	{
#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
		PySys_WriteStdout("~Matrix()\n");
#endif    
		qcs_delete_matrix(self);
	}

	%feature("autodoc", "__add__(Matrix a)") __add__(Matrix a);
	Matrix __add__(Matrix a)
	{
	  	tf_qcs_matrix *tqm;

		tqm=qcs_create_matrix(a.rows, a.cols);
		
		qcs_add_matrix(&a, self, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__sub__(Matrix a)") __sub__(Matrix a);
	Matrix __sub__(Matrix a)
	{
	  	tf_qcs_matrix *tqm;

		tqm=qcs_create_matrix(a.rows, a.cols);
		
		qcs_sub_matrix(&a, self, tqm);
		
		return *tqm;
	}

	%feature("autodoc", "__mul__(Matrix a)") __mul__(Matrix a);
	Matrix __mul__(Matrix a)
	{
	  	tf_qcs_matrix *tqm;

		tqm=qcs_create_matrix(a.rows, a.cols);
		
		qcs_mul_matrix(&a, self, tqm);
		
		return *tqm;
	}

	/*
		mul complex scalar and matrix
	*/
	
	%feature("autodoc", "__rmul__(complex a)") __rmul__(complex a);
	Matrix __rmul__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__mul__(complex a)") __mul__(complex a);
	Matrix __mul__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	/*
		mul float scalar and matrix
	*/
	
	%feature("autodoc", "__rmul__(float a)") __rmul__(float a);
	Matrix __rmul__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__mul__(float a)") __mul__(float a);
	Matrix __mul__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		mul int scalar and matrix
	*/
	
	
	%feature("autodoc", "__rmul__(int a)") __rmul__(int a);
	Matrix __rmul__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__mul__(int a)") __mul__(int a);
	Matrix __mul__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_mul_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		div complex scalar and matrix
	*/
	
	%feature("autodoc", "__rdiv__(complex a)") __rdiv__(complex a);
	Matrix __rdiv__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__div__(complex a)") __div__(complex a);
	Matrix __div__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	
	/*
		div float scalar and matrix
	*/
	
	%feature("autodoc", "__rdiv__(float a)") __rdiv__(float a);
	Matrix __rdiv__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__div__(float a)") __div__(float a);
	Matrix __div__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
	div int scalar and matrix
	*/
	
	%feature("autodoc", "__rdiv__(int a)") __rdiv__(int a);
	Matrix __rdiv__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__div__(int a)") __div__(int a);
	Matrix __div__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_div_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}	

	/*
		add complex scalar and matrix
	*/
	
	%feature("autodoc", "__add__(complex a)") __add__(complex a);
	Matrix __add__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	%feature("autodoc", "__radd__(complex a)") __radd__(complex a);
	Matrix __radd__(complex  a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re = __real__ a;
		cmplx.im = __imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		add float scalar and matrix
	*/
	
	%feature("autodoc", "__radd__(float a)") __radd__(float a);
	Matrix __radd__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__add__(float a)") __add__(float a);
	Matrix __add__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		add int scalar and matrix
	*/	
	
	%feature("autodoc", "__radd__(int a)") __radd__(int a);
	Matrix __radd__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__add__(int a)") __add__(int a);
	Matrix __add__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_add_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		sub complex scalar and matrix
	*/
	
	%feature("autodoc", "__rsub__(complex a)") __rsub__(complex a);
	Matrix __rsub__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=__real__ a;
		cmplx.im=__imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__sub__(complex a)") __sub__(complex a);
	Matrix __sub__(complex a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=__real__ a;
		cmplx.im=__imag__ a;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	
	/*
		sub float scalar and matrix
	*/
	
	%feature("autodoc", "__rsub__(float a)") __rsub__(float a);
	Matrix __rsub__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__sub__(float a)") __sub__(float a);
	Matrix __sub__(float a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}

	/*
		sub int scalar and matrix
	*/	
	
	%feature("autodoc", "__rsub__(int a)") __rsub__(int a);
	Matrix __rsub__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}
	
	%feature("autodoc", "__sub__(int a)") __sub__(int a);
	Matrix __sub__(int a)
	{
	  	tf_qcs_matrix *tqm;
		tf_qcs_complex cmplx;

		cmplx.re=a;
		cmplx.im=0;
		
		tqm=qcs_create_matrix(self->rows, self->cols);
		
		qcs_sub_scalar_matrix(self, &cmplx, tqm);
		
		return *tqm;
	}	
	
	%feature("autodoc", "Calc_D_dot_DT()") Calc_D_dot_DT();
	Matrix Calc_D_dot_DT()
	{
		tf_qcs_matrix *tqm;
		
		tqm = qcs_calculate_d_dot_dt_matrix(self);
		
		return *tqm;
	}

	%feature("autodoc", "Chop()") Chop();
	void Chop()
	{
		qcs_chop_matrix( self );
	}
	
	void SetQF(int q, int f)
	{
		self->q=q;
		self->freedom_level=f;
	}

	%apply int *OUTPUT { int *q, int *f };
	void GetQF(int *q, int *f)
	{
		*q = self->q;
		*f = self->freedom_level;
	}
	
	%feature("autodoc", "At(int r, int c, Complex cnum)") At(int r, int c, Complex cnum);
	void At(int r, int c, Complex cnum)
	{
		qcs_set_cell_at_matrix_complex(self, r, c, &cnum);
	}

	%feature("autodoc", "AtDirect(int r, int c, float re, float im)") AtDirect(int r, int c, float re, float im);
	void AtDirect(int r, int c, float re, float im)
	{
		qcs_set_cell_at_matrix_direct(self, r, c, re, im);
	}

	%feature("autodoc", "AddNoise(float n)") AddNoise(float n);
	void AddNoise(float n)
	{
		qcs_add_noise_to_matrix(self, n);
	}

	%feature("autodoc", "Pr()") Pr();
	void Pr() 
	{
		qcs_print_matrix(self);
	}

	%feature("autodoc", "PrMatlab()") PrMatlab();
	void PrMatlab()
	{
		qcs_print_matrix_in_matlab_format(self);
	}

	%feature("autodoc", "PrMathematica()") PrMathematica();
	void PrMathematica()
	{	
		qcs_print_matrix_in_mathematica_format(self);
	}

	%feature("autodoc", "PrE()") PrE();
	void PrE()
	{
		qcs_print_eigenvalue_of_matrix(self);
	}

	%feature("autodoc", "SaveToFlatTextFile(char *fname)") SaveToFlatTextFile(char *fname);
	void SaveToFlatTextFile(char *fname)
	{
		qcs_save_matrix_as_flat_text_file_for_gnuplot(self, fname);
	}
	
	%feature("autodoc", "MakeEV()") MakeEV();
	void MakeEV()
	{
		qcs_make_eigenvectors_of_matrix(self);
	}

	%feature("autodoc", "Entropy()") Entropy();
	float Entropy()
	{
		return qcs_entropy_of_matrix(self);
	}

	%feature("autodoc", "LinearEntropy()") LinearEntropy();
	float LinearEntropy()
	{
		return qcs_linear_entropy_of_matrix(self);
	}

	%feature("autodoc", "Negativity()") Negativity();
	float Negativity()
	{
		return qcs_negativity_of_matrix(self);
	}
	
    %feature("autodoc", "GetCol(int c)") GetCol(int c);
    Matrix GetCol(int c)
    {
        return *qcs_get_column_from_matrix( self, c);
    }

    %feature("autodoc", "GetRow(int c)") GetRow(int r);
    Matrix GetRow(int r)
    {
        return *qcs_get_row_from_matrix( self, r);
    }
    
	%feature("autodoc", "PTraceQubit1(int n)") PTraceQubit1(int n);
	Matrix PTraceQubit1(int n)
	{
		return *qcs_create_partial_trace_matrix_1_qubit(self, n);
	}

	%feature("autodoc", "PTraceQudit1(int n)") PTraceQudit1(int n);
	Matrix PTraceQudit1(int n)
	{
		return *qcs_create_partial_trace_matrix_1_qudit(self, n);
	}
	
	%feature("autodoc", "PTraceQubitN(int _from, int _to)") PTraceQubitN(int _from, int _to);
	Matrix PTraceQubitN(int _from, int _to)
	{
		return *qcs_create_partial_trace_matrix_n_qubit(self, _from, _to);
	}

	%feature("autodoc", "Sqrt()") Sqrt();
	void Sqrt()
	{
		qcs_square_root_of_operator_matrix_self( self );
	}	

	%feature("autodoc", "Exp()") Exp();
	void Exp()
	{
		qcs_exp_of_matrix( self );
	}

	%feature("autodoc", "Inv()") Inv();
	void Inv()
	{
		qcs_inv_matrix( self );
	}
	
	%feature("autodoc", "Transpose()") Transpose();
	void Transpose()
	{
		qcs_transpose_matrix(self);
	}

	%feature("autodoc", "PTransposeQubit(int n)") PTransposeQubit(int n);
	void PTransposeQubit(int n)
	{
		qcs_partial_transpose_matrix(self, n);
	}	

	%feature("autodoc", "PTransposeQudit(int n)") PTransposeQudit(int n);
	void PTransposeQudit(int n)
	{
		qcs_partial_transpose_matrix_qudit(self, n);
	}	

	
	%feature("autodoc", "Tr()") Tr();
	Complex Tr()
	{
		return qcs_trace_matrix(self);
	}	  

	%feature("autodoc", "TrSqr()") TrSqr();
	Complex TrSqr()
	{
		return qcs_trace_square_matrix(self);
	}
	
	%feature("autodoc", "Eye()") Eye();
	void Eye()
	{
		qcs_eye_matrix(self);
	}

	%feature("autodoc", "EyeP(float p)") EyeP(float p);
	void EyeP(float p)
	{
		qcs_eye_matrix_with_param(self, p);
	}
	
	%feature("autodoc", "Zero()") Zero();
	void Zero()
	{
		qcs_zero_matrix(self);
	}
	
	%feature("autodoc", "IsPure() -> int") IsPure();
	int IsPure()
	{
		return qcs_matrix_ispure(self);
	}
	
	%feature("autodoc", "SpectralDecomposition() -- create spectral decomposition") SpectralDecomposition();
	SpectralDecomposition SpectralDecomposition()
	{
		SpectralDecomposition *sd;
		tf_qcs_matrix *eigenvectors=NULL, *eigenvalues=NULL;
		
		sd = qcs_create_spectral_decomposition();
		
		eigenvectors = qcs_create_matrix(self->rows, self->cols);
		eigenvalues = qcs_create_matrix(self->rows, 1);

		qcs_spectral_decompose_of_matrix(self, eigenvalues, eigenvectors);
		
		sd->eigenvalues =  eigenvalues;
		sd->eigenvectors = eigenvectors;
		
		return *sd;
	}

	%feature("autodoc", "SVDDecomposition() -- create svd decomposition") SVDDecomposition();
	SVDDecomposition SVDDecomposition()
	{
		SVDDecomposition *sd;
		
		sd = qcs_create_svd_decomposition();
		
        sd->U = qcs_create_empty_matrix();
        sd->S = qcs_create_empty_matrix();
        sd->V = qcs_create_empty_matrix();
		
		qcs_svd_decompose_of_matrix(self, sd->S, sd->U, sd->V);
			
		return *sd;
	}
    
    
	%feature("autodoc", "CreateW0() -- create 9x9 W_0 witnesses") CreateW0();
	void CreateW0()
	{
		qcs_delete_matrix(self);
	
		self = qcs_create_w0_9x9_matrix();
	}

	%feature("autodoc", "CreateMaximallyMixedState(int x) -- ") CreateMaximallyMixedState(int x);  
	void CreateMaximallyMixedState(int x)
	{
		qcs_delete_matrix(self);
	
		self = qcs_create_maximally_mixed_state(x);
	}

	%feature("autodoc", "CreateMaximallyMixedState(int x, int y) -- ") CreateMaximallyMixedState(int x, int y);  
	void CreateMaximallyMixedState(int x, int y)
	{
		qcs_delete_matrix(self);
	
		self = qcs_create_maximally_mixed_state(x);
	}

	%feature("autodoc", "CreateMaximallyMixedStateP() -- ") CreateMaximallyMixedStateP(int x, float g);  
	void CreateMaximallyMixedStateP(int x, float g, int q, int f)
	{
		qcs_delete_matrix(self);
	
		self = qcs_create_maximally_mixed_state_with_param(x, g);
		self->q=q;
		self->freedom_level=f;
	}
	
	%feature("autodoc", "CreateHa9x9State(float g) -- ") CreateHa9x9State(float g);  
	void CreateHa9x9State(float g)
	{
		qcs_delete_matrix(self);
		
		self = qcs_create_ha_9x9_state(g);
	}

	%feature("autodoc", "CreateHa9x9StateP(float g, float p) -- ") CreateHa9x9StateP(float g, float p);  
	void CreateHa9x9StateP(float g, float p)
	{
		qcs_delete_matrix(self);
		
		self = qcs_create_ha_9x9_state_with_param(g, p);
	}
	
	%feature("autodoc", "CreateHorodecky9x9State(float a)") CreateHorodecky9x9State(float a);
	void CreateHorodecky9x9State(float a)
	{
		qcs_delete_matrix(self);
		
		self = qcs_create_horodecky_9x9_state(a);
	}

	%feature("autodoc", "CreateHorodecky4x4State(float p, float a, float b)") CreateHorodecky4x4State(float p, float a, float b);
	void CreateHorodecky4x4State(float p, float a, float b)
	{
		qcs_delete_matrix(self);
		
		self = qcs_create_horodecki_4x4_state(p, a, b);
	}

	%feature("autodoc", "CreateMatrixForState(int q, int f)") CreateMatrixForState(int q, int f);
	void CreateMatrixForState(int q, int f)
	{
		qcs_delete_matrix(self);
		
		self = qcs_create_matrix_for_quantum_state(q, f);
	}
	
	
	%feature("autodoc", "PPTest() -> float") PPTest();
	float PPTest()
	{
		return qcs_ppt_criterion( self );
	}

	%feature("autodoc", "CCNRTest() -> float") CCNRTest();
	float CCNRTest()
	{
		return qcs_ccnr_criterion( self );
	}
	
	int NonZeros()
	{
		return qcs_non_zero_elements_of_matrix( self );
	}
};


%extend Qubit {

	%feature("autodoc", "Qubit()") Qubit();
	Qubit()
	{
      	tf_qcs_qubit *q;
	
  	    q=qcs_new_qubit();

	    return q;
	}

	%feature("autodoc", "~Qubit()") ~Qubit();
	~Qubit()
	{
	    qcs_delete_qubit(self);
	}

	%feature("autodoc", "Id()") Id();
	void Id()
	{
	    tf_qcs_qubit t;

	    qcs_1q_id_gate(self, &t);

	    self->alpha=t.alpha;
	    self->beta=t.beta;
	}
  
	%feature("autodoc", "Not()") Not();
	void Not()
	{
		tf_qcs_qubit t;

		qcs_1q_not_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
	}

   %feature("autodoc", "PauliX()") PauliX();
   void PauliX()
   {
      tf_qcs_qubit t;

      qcs_1q_pauli_x_gate(self, &t);

      self->alpha=t.alpha;
      self->beta=t.beta;
   }

   %feature("autodoc", "PauliY()") PauliY();
   void PauliY() 
   {
      tf_qcs_qubit t;

      qcs_1q_pauli_y_gate(self, &t);

      self->alpha=t.alpha;
      self->beta=t.beta;
   }

   %feature("autodoc", "PauliZ()") PauliZ();
   void PauliZ()
   {
      tf_qcs_qubit t;

      qcs_1q_pauli_z_gate(self, &t);

      self->alpha=t.alpha;
      self->beta=t.beta;
   }

  %feature("autodoc", "Had()") Had();
  void Had()
  {
      tf_qcs_qubit t;

      qcs_1q_hadamard_gate(self, &t);

      self->alpha=t.alpha;
      self->beta=t.beta;
  }

  %feature("autodoc", "SquareRoot()") SquareRoot();
  void SquareRoot()
  {
      tf_qcs_qubit t;

      qcs_1q_square_root_gate(self, &t);

      self->alpha=t.alpha;
      self->beta=t.beta;
  }
	
  %feature("autodoc", "RotAlpha(float k)") RotAlpha(float k);
  void RotAlpha(float k)
  {
		tf_qcs_qubit t;
		
		qcs_1q_rotate_alpha_gate(self, k, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
  }	

  %feature("autodoc", "RotTheta(float theta)") RotTheta(float theta);
  void RotTheta(float theta)
  {
		tf_qcs_qubit t;
		
		qcs_1q_rotate_theta_gate(self, theta, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
  }	
  
  %feature("autodoc", "XRot90() ") XRot90();
  void XRot90() 
  {
		tf_qcs_qubit t;
		
		qcs_1q_x_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
  }
	
  %feature("autodoc", "YRot90()") YRot90();
  void YRot90()
  {
		tf_qcs_qubit t;
		
		qcs_1q_y_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
	}

  %feature("autodoc", "ZRot90()") ZRot90();
  void ZRot90()
  {
		tf_qcs_qubit t;
		
		qcs_1q_z_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
		
  }
	
  %feature("autodoc", "MXRot90()") MXRot90();
  void MXRot90()
  {
		tf_qcs_qubit t;
		
		qcs_1q_minus_x_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
 	}
	
  %feature("autodoc", "MYRot90()") MYRot90();
  void MYRot90()
  {
		tf_qcs_qubit t;
		
		qcs_1q_minus_y_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
	}

  %feature("autodoc", "MZRot90()") MZRot90();
  void MZRot90()
  {
		tf_qcs_qubit t;
		
		qcs_1q_minus_z_rotate90_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
	}
	
  %feature("autodoc", "T()") T();
  void T() 
  {
		tf_qcs_qubit t;
		
		qcs_1q_t_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
  }
	
  %feature("autodoc", "S()") S();
  void S() 
  {
		tf_qcs_qubit t;
		
		qcs_1q_s_gate(self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
  }
	
  %feature("autodoc", "Measure()") Measure();
  int Measure() 
  {
    return qcs_qubit_measure(self);    
  } 
  
  %feature("autodoc", "MiscGate(Matrix m)") MiscGate(Matrix m);
  void MiscGate(Matrix m) 
  {
    tf_qcs_qubit t;
    Matrix tmp_m;
      
    qcs_1q_arbitrary_gate(&tmp_m, self, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;
		
  }

  %feature("autodoc", "PhaseNK(float k)") PhaseNK(float k);
  void PhaseNK(float k) 
  {
		tf_qcs_qubit t;
		
		qcs_1q_phase_gate(self, k, &t);

		self->alpha=t.alpha;
		self->beta=t.beta;		
  }

	void SetRandomState()
	{
		qcs_set_random_real_state_qubit(self);
	}
  
  %feature("autodoc", "SetState(int n)") SetState(int n);
  void SetState(int n)
  {
	if( n == 0 )
		qcs_set_state0_qubit(self);
	else
		qcs_set_state1_qubit(self);
  }

  
  %feature("autodoc", "Set0()") Set0();
  void Set0()
  {
	qcs_set_state0_qubit(self);
  }

  %feature("autodoc", "Set1()") Set1();
  void Set1()
  {
	qcs_set_state1_qubit(self);
  }

  %feature("autodoc", "SetPer(float _per)") SetPer(float _per);
  void SetPer(float _per)
  {
  		qcs_set_superposition_state_per(self, _per);
  }

  %feature("autodoc", "SetThetaPsi(int _theta, int _psi)") SetThetaPsi(int _theta, int _psi);
  void SetThetaPsi(int _theta, int _psi)
  {
  		qcs_set_superposition_state_theta_psi(self, _theta, _psi);
  }

  %feature("autodoc", "Pr()") Pr();
  void Pr()
  {
		qcs_printf_qubit(self);
  }

  %feature("autodoc", "PrSqr()") PrSqr();
  void PrSqr() 
  {
		qcs_printf_qubit_sqr(self);
  }

  %feature("autodoc", "PrVec()") PrVec();
  void PrVec() 
  {	
		qcs_printf_qubit_vec(self);
  }

  %feature("autodoc", "PrVecSqr()") PrVecSqr();
  void PrVecSqr() 
  {
		qcs_printf_qubit_vec_sqr(self);
  }

  %feature("autodoc", "") ;
  char *__str__() 
  {
		static char temp[512];
	  
		sprintf(temp, "(%1.4f+%1.4fi|0> + %1.4f+%1.4fi|1>\n", 
			 self->alpha.re, self->alpha.im,
			 self->beta.re, self->beta.im);
    
		return temp; 
  }
  
}; /* end of %extend qubit */

%extend Qudit {

	%feature("autodoc", "Qudit()") Qudit();
	Qudit()
	{
      	tf_qcs_qudit *q;
	
  	    q=qcs_new_qudit( 2 );

	    return q;
	}

	%feature("autodoc", "Qudit(int level)") Qudit(int level);
	Qudit(int level)
	{
      	tf_qcs_qudit *q;
	
  	    q=qcs_new_qudit( level );

	    return q;
	}	
	
	%feature("autodoc", "Reset()") Reset();
	void Reset()
	{
		qcs_set_state_n_qudit(self, 0);
	}

	%feature("autodoc", "SetState(int s)") SetState(int s);
	void SetState(int s)
	{
		qcs_set_state_n_qudit(self, s);
	}
	
	%feature("autodoc", "SetRandomState()") SetRandomState();
	void SetRandomState()
	{
		qcs_set_random_real_state_qudit( self );
	}

	%feature("autodoc", "Pr()") Pr();
	void Pr()
	{
		qcs_printf_qudit(self);
	}

	%feature("autodoc", "PrSqr()") PrSqr();
	void PrSqr()
	{
		qcs_printf_qudit_sqr(self);
	}

	%feature("autodoc", "~Qudit()") ~Qudit();
	~Qudit()
	{
	    qcs_delete_qudit(self);
	}

}; /* end of %extend qudit */

/* %extend QubitArray {

	%feature("autodoc", "QubitArray()") QubitArray();
	QubitArray()
	{
		tf_qcs_qubit_array *qa;
	   
	   	qa=qcs_new_qubit_array(2);
	   
		return qa;
	}
	
	%feature("autodoc", "QubitArray(int i)") QubitArray(int i);
	QubitArray(int i)
	{
		tf_qcs_qubit_array *qa;
	   
	   	qa=qcs_new_qubit_array(i);
	   
		return qa;
	}
	
	%feature("autodoc", "~QubitArray()") ~QubitArray();
    ~QubitArray()
  	{
  		qcs_delete_qubit_array(self); 
	}
	
	%feature("autodoc", "SetQubitN(int n, Qubit q)") SetQubitN(int n, Qubit q);
	void SetQubitN(int n, Qubit q)
	{
		qcs_set_qubit_array_n(self, n, &q);
	}
	
	%feature("autodoc", "GetQubitN(int n)") GetQubitN(int n);
	Qubit GetQubitN(int n)
	{
		return *qcs_get_qubit_array_n(self, n);
	}
	
};

%extend QuditArray {

	%feature("autodoc", "QuditArray()") QuditArray();
	QuditArray()
	{
		tf_qcs_qudit_array *qa;
	   
	   	qa=qcs_new_qudit_array(2, 2);
	   
		return qa;
	}
	
	%feature("autodoc", "QuditArray(int i, int level)") QuditArray(int i, int level);
	QuditArray(int i, int level)
	{
		tf_qcs_qudit_array *qa;
	   
	   	qa=qcs_new_qudit_array(i, level);
	   
		return qa;
	}

	%feature("autodoc", "~QuditArray()") ~QuditArray();
    ~QuditArray() 
  	{
  		qcs_delete_qudit_array(self); 
	}
	
	%feature("autodoc", "SetQuditN(int n, Qudit q)") SetQuditN(int n, Qudit q);
	void SetQuditN(int n, Qudit q)
	{
		qcs_set_qudit_array_n(self, n, &q);
	}
	
	%feature("autodoc", "GetQuditN(int n)") GetQuditN(int n);
	Qudit GetQuditN(int n)
	{
		return *qcs_get_qudit_array_n(self, n);
	}

}; */

%extend QuArray {

	%feature("autodoc", "QuArray()") QuArray();
	QuArray()
	{
		ts_qcs_quarray *qua;
		
		qua = qcs_quarray_qubit_create( 2 );
		
		return qua;
	}

	%feature("autodoc", "QuArray(int n)") QuArray(int n);
	QuArray(int n)
	{
		ts_qcs_quarray *qua;
		
		qua = qcs_quarray_qubit_create(n);
	}

	%feature("autodoc", "QuArray(int n, int freedom_level)") QuArray(int n, int freedom_level);
	QuArray(int n, int freedom_level)
	{
		ts_qcs_quarray *qua;
		
		qua = qcs_quarray_qudit_create(n, freedom_level);
		
		return qua;
	}
	
/*	QuArray(int n, int type, int freedom_level)
	{
		ts_qcs_quarray *qua;
		
		switch(type)
		{
			case ARRAY_QUBIT:
				qua = qcs_quarray_qubit_create(n);
			break;

			case ARRAY_QUDIT:
				qua = qcs_quarray_qudit_create(n, freedom_level);
			break;

			case ARRAY_QUBIT_SYMBOLIC:
				qua = qcs_quarray_qubit_symbolic_create(n);
			break;

			case ARRAY_QUDIT_SYMBOLIC:
				qua = NULL;
			break;
 		}
		
		return qua;
	}
*/
	
	%feature("autodoc", "SetN(int n, QuUnit q)") SetN(int n, void* q_ptr);
	void SetN(int n, void* q_ptr)
	{
		qcs_set_quarray_n(self, n, q_ptr);
	}

	%feature("autodoc", "GetN(int n) -> void*") GetN(int n);
	void *GetN(int n)
	{
		return qcs_get_quarray_n(self, n);
	}

	%feature("autodoc", "GetQubitN(int n) -> Qubit") GetQubitN(int n);
	Qubit GetQubitN(int n)
	{
		Qubit *a;
		
		a = (Qubit*)qcs_get_quarray_n(self, n);
		
		return *a;
	}

	%feature("autodoc", "GetQuditN(int n)") GetQuditN(int n);
	Qudit GetQuditN(int n)
	{
		Qudit *a;
		
		a = (Qudit*)qcs_get_quarray_n(self, n);
		
		return *a;
	}
	
	%feature("autodoc", "~QuArray()") ~QuArray();
	~QuArray()
	{
		qcs_quarray_destroy(self);
	}	

}; /* end of %extend quarray */

%extend QWalk1D {

	QWalk1D()
	{
		tf_qwalk_1d *qw;
		
		qw = qcs_create_qwalk_1d( get_hadamard_gate(), QWalk1D_Segment_Lattice, 100);
		qw->rbl_prob=0.00f;
		qw->experiments=1;
		
		return qw;
	}
	
	QWalk1D(int type, int steps)
	{
		tf_qwalk_1d *qw;
		
		qw = qcs_create_qwalk_1d( get_hadamard_gate(), type, steps);
		qw->rbl_prob=0.00f;
		qw->experiments=1;
		
		return qw;
	}

	QWalk1D(int _type, int _steps, int _experiments)
	{
		tf_qwalk_1d *qw;
		
		qw = qcs_create_qwalk_1d( get_hadamard_gate(), _type, _steps);
		qw->rbl_prob=0.00f;
		qw->experiments=_experiments;
		
		return qw;
	}

	QWalk1D(int _type, int _steps, int _experiments, float _rbl_prob)
	{
		tf_qwalk_1d *qw;
		
		qw = qcs_create_qwalk_1d( get_hadamard_gate(), _type, _steps);
		qw->rbl_prob=_rbl_prob;
		qw->experiments=_experiments;
		
		return qw;
	}	
	
	void SetCoin(Matrix coin)
	{
		qcs_qwalk_1d_set_coin(self, &coin);
	}
	
	void Reset()
	{
		qcs_qwalk_1d_reset( self );
	}
	
	void SetHadamardState()
	{
		qcs_qwalk_1d_set_hadamard_state( self );
	}
	
	void RandomBrokenLink()
	{
		qcs_qwalk_1d_random_broken_link(self);
	}
	
	void OneIteration()
	{
		qcs_qwalk_1d_initialise_random_broken_links( self );
        qcs_qwalk_1d_random_broken_link( self );
        qcs_qwalk_1d_one_iteration( self );
		
		self->iteration++;
	}
	
	void DoExperiments()
	{
		int e,i;
	    for( e=0; e<self->experiments; e++ )
		{
			qcs_qwalk_1d_set_hadamard_state( self );

			self->iteration=0;
			for( i=0; i<self->steps; i++ )
			{
				qcs_qwalk_1d_initialise_random_broken_links( self );
				qcs_qwalk_1d_random_broken_link( self );

				qcs_qwalk_1d_one_iteration( self );

				self->iteration++;

			}
			qcs_qwalk_1d_update_probability( self );
		}
	}
	
	void SaveAverageProbatility(char *fname)
	{
		qcs_qwalk_1d_save_average_probatility( self, fname);
	}

	void SaveActualProbatility(char *fname)
	{
		qcs_qwalk_1d_save_actual_probatility( self, fname);
	}
	
	float GetActualProbability(int n)
	{
		return qcs_qwalk_1d_get_actual_probability(self, n);
	}
	
	float GetAverageProbability(int n)
	{
		return qcs_qwalk_1d_get_average_probability(self, n);
	}
	
	Complex GetElementFromOrgA(int x, int y)
	{
		return *qcs_qwalk_1d_get_org_a_element(self, x, y);
	}
	
	void SetElementInOrgA(int x, int y, float re, float im)
	{
		qcs_qwalk_1d_set_org_a_element(self, x, y, re, im);
	}

	int GetElementFromRBL(int x, int y)
	{
		return qcs_get_cell_at_int_matrix(self->RBL, x, y);
	}
	
	void SetElementInRBL(int x, int y, int r)
	{
		qcs_set_cell_at_int_matrix(self->RBL, x, y, r);
	}
	
	float GetRBLProg()
	{
		return self->rbl_prob;
	}
	
	void SetRBLProg(float f)
	{
		self->rbl_prob = f;
	}
	
	void PrOrgA()
	{
		qcs_print_matrix( self->org_A );
	}
	
	void PrRBL()
	{
		qcs_print_int_matrix( self->RBL );
	}
	
	~QWalk1D()
	{
		qcs_destroy_qwalk_1d( self );
	}	
	
};

%extend QWalk2D {

	QWalk2D()
	{
	}
	
	QWalk2D(int step, int type)
	{
	}
	
	void Reset()
	{
	}
	
};


%extend QuantumReg {

  %feature("autodoc", "QuantumReg()") QuantumReg();
  QuantumReg()
  {
    QuantumReg* m;
    m=qcs_new_quantum_reg(2);

#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
    PySys_WriteStdout("QubitReg() new 2 qubits\n");
#endif
    return m;
  }

  %feature("autodoc", "QuantumReg(Matrix a)") QuantumReg(Matrix a);
  QuantumReg(Matrix a)
  {
    QuantumReg* m;
	
	m = qcs_new_quantum_reg_qubit_from_density_matrix( &a );
	
	return m;
  }
  
  %feature("autodoc", "QuantumReg(int i)") QuantumReg(int i) ;
  QuantumReg(int i) 
  {
    QuantumReg* m;
    m=qcs_new_quantum_reg(i);

#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
    PySys_WriteStdout("QuantumReg(int i) new %d qubits \n", i);
#endif
    return m;
  }

  %feature("autodoc", "QuantumReg(int i, int freedom_level)") QuantumReg(int i, int freedom_level);
  QuantumReg(int i, int freedom_level)
  {
    QuantumReg *m = NULL;

	switch(freedom_level)
	{
		case USE_SYMBOLIC_STATE_VECTOR_QUBIT:
		{
			m = qcs_new_quantum_reg_symbolic_qubit(i);
#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
				PySys_WriteStdout("QuantumReg(int i) new %d qubits in symbolic mode\n", i);
#endif
		} break;
		case USE_DENSITY_MATRIX:
		{
			m = qcs_new_quantum_reg_density_matrix(i);
#if defined(PYTHON_SCRIPT) && defined( __qcs_core_library_debug_mode__ )
				PySys_WriteStdout("QuantumReg(int i) new %d qubits in density matrix mode\n", i);
#endif
		} break;
		default:
			m=qcs_new_quantum_reg_qudit(i, freedom_level);
	}
	return m;
  }
  
  %feature("autodoc", "QuantumReg(int i, int freedom_level, int k)") QuantumReg(int i, int freedom_level, int k);
  QuantumReg(int i, int freedom_level, int k)
  {
    QuantumReg* m;
    
    if(k == NO_STATE_VECTOR)
    {
      m=qcs_new_quantum_reg_no_state_vec(i);
    }
    if(k == USE_STATE_VECTOR_QUBIT)
    {
      m=qcs_new_quantum_reg(i);
    }

    if(k == USE_DENSITY_MATRIX)
    {
      m=qcs_new_quantum_reg_density_matrix(i);
    }
    
    if(k == USE_GRAPH_STATE_DESC)
    {
      m=qcs_new_quantum_reg_graph_state_desc(i);
    }

    if(k == PQC_MODE)
    {
      m=qcs_new_quantum_reg_pqc_mode(i);
    }

    if(k == CHP_MODE)
    {
		m = qcs_new_quantum_reg_chp_mode(i);
    }

    if(k == ONEWAY_MODEL)
    {
		m = qcs_new_quantum_reg_oneway_model(i, i);
    }
    
	if( k == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
	{
		m = qcs_new_quantum_reg_symbolic_qudit(i, freedom_level);
	}
	
    // printf("QubitReg(int i, int k) new %d %d\n", i, k);

    return m;
  }

  void Substitution(char *name, float value)
  {
	qcs_quantum_reg_substitution(self, name, value);
  }

  void ClearAmplitude( int state)
  {
	qcs_quantum_reg_clear_amplitude(self, state);
  }

#ifdef PYTHON_SCRIPT
%pythoncode %{
	def ClearAmplitudes( self, list ):
		for i in list:
			print( i )
%}
#endif
  
	%feature("autodoc", "GetQF() -> [q,f]") GetQF(int *q, int *f);  
	%apply float *OUTPUT { int *q, int *f };
	void GetQF(int *q, int *f)
	{
		int qq, ff;
		qcs_quantum_reg_get_qf(self, &qq, &ff);
		
		*q=qq;
		*f=ff;
	}
  
  %feature("autodoc", "ResetError()") ResetError();
  void ResetError()
  {
	qcs_quantum_reg_error_reset(self);
  }

  %feature("autodoc", "ErrorLevel()") ErrorLevel();
  int ErrorLevel()
  {
  	return qcs_quantum_reg_error_level(self);
  }

  %feature("autodoc", "~QuantumReg()") ~QuantumReg();
  ~QuantumReg()
  {
    qcs_delete_quantum_reg(self);
  }

  %feature("autodoc", "Info()") Info();
  void Info()
  {
#ifdef PYTHON_SCRIPT
    PySys_WriteStdout("Use %d qubits\nMode:", self->size);
#endif
    switch(self->mode)
    {
          case NO_STATE_VECTOR:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tno state vector\n");
#endif
          } break;
          case USE_STATE_VECTOR_QUBIT:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tstate vector\n");
#endif
          } break;
          case USE_DENSITY_MATRIX:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tdensity matrix\n");
#endif
          } break;
          case USE_GRAPH_STATE_DESC:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tgraph state\n");
#endif
          } break;
          case PQC_MODE:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tpqc mode\n");
#endif
          } break;
          case CHP_MODE:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tchp mode\n");
#endif
          } break;
          case ONEWAY_MODEL:
          {
#ifdef PYTHON_SCRIPT
            PySys_WriteStdout("\tone way quantum computation mode\n");
#endif
          } break;
    };
  }
  
  %feature("autodoc", "GetVecStateN(int i) -> Complex") GetVecStateN(int i);
  Complex GetVecStateN(int i)
  {
    return self->vec_state[i];
  }

  %feature("autodoc", "__getitem__(int i) -> Complex") __getitem__(int i);
  Complex __getitem__(int i)
  {
	return self->vec_state[i];
  }
  
  %feature("autodoc", "SetVecStateN(int i, float re, float im)") SetVecStateN(int i, float re, float im);
  void SetVecStateN(int i, float re, float im)
  {
    (self->vec_state+i)->re=re;
    (self->vec_state+i)->im=im;
  }

  %feature("autodoc", "GetQubitN(int i) -> Qubit") GetQubitN(int i);
  Qubit GetQubitN(int i)
  {
   // return *self->qubits[i];
   return *qcs_quantum_reg_get_qubit_n(self, i);
  }

  %feature("autodoc", "Reset()") Reset();
  void Reset()
  {
    qcs_quantum_reg_reset(self);
  }

  %feature("autodoc", "SetGHZState()") SetGHZState();
  void SetGHZState()
  {
	qcs_quantum_reg_ghzstate(self);
  }

  %feature("autodoc", "SetGHZState()") SetWState();
  void SetWState()
  {
	qcs_quantum_reg_wstate(self);
  }
  
  %feature("autodoc", "ResetSymbolic()") ResetSymbolic();
  void ResetSymbolic()
  {
    qcs_quantum_reg_reset_symbolic(self);
  }
  
  %feature("autodoc", "FillZero()") FillZero();
  void FillZero()
  {
    qcs_quantum_reg_fill_zero(self);
  }
  
  %feature("autodoc", "Cleanup()") Cleanup();
  void Cleanup()
  {
	qcs_quantum_reg_cleanup( self );
  }
  
  %feature("autodoc", "SetN(int i)") SetN(int i);
  void SetN(int i)
  {
    qcs_quantum_reg_set_state_dec(self, i);
  }

  %feature("autodoc", "SetKet(char *_str)") SetKet(char *_str);
  void SetKet(char *_str)
  {
    qcs_quantum_reg_set_state_bin(self, _str);
  }

  void SetInitialStateForShorAlg(int b, int n, int i)
  {
	qcs_quantum_reg_factorization_initial_state(self, b, n, i);
  }
 
#ifdef PYTHON_SCRIPT
%pythoncode %{
def SetMultistate(self, states):
    import math
    self.FillZero()
    ampl=1.0/math.sqrt(len(states))
    for i in states:        
        self.SetVecStateN(i,ampl,0)
%}
#endif

  %feature("autodoc", "SetFromQuArray(QuArray qa)") SetFromQuArray(QuArray qa);
  void SetFromQuArray(QuArray qa)
  {
  	qcs_quantum_reg_set_state_from_quarray(self, &qa);
  }
  
  %feature("autodoc", "Noop()") Noop();
  void Noop()
  {
	return;
  }

	%feature("autodoc", "RandGateRealN(int i)") RandGateRealN(int i);
	void RandGateRealN(int i)
	{
		qcs_quantum_reg_random_one_gate_n(self, i);
	}
  
	%feature("autodoc", "NotN(int i)") NotN(int i);
	void NotN(int i)
	{
		qcs_quantum_reg_not_n(self, i);
	}

  %feature("autodoc", "NotN_Conj(int i)") NotN_Conj(int i);
  void NotN_Conj(int i)
  {
    qcs_quantum_reg_not_n_conj(self, i);
  }

  %feature("autodoc", "PauliX(int i)") PauliX(int i);
  void PauliX(int i)
  {
    qcs_quantum_reg_pauli_x_n(self, i);
  }
  
  %feature("autodoc", "PauliY(int i)") PauliY(int i);
  void PauliY(int i)
  {
    qcs_quantum_reg_pauli_y_n(self, i);
  }

  %feature("autodoc", "PauliZ(int i)") PauliZ(int i);
  void PauliZ(int i)
  {
    qcs_quantum_reg_pauli_z_n(self, i);
  }

  %feature("autodoc", "GateT(int i)") GateT(int i);
  void GateT(int i)
  {
  	qcs_quantum_reg_t_n(self, i);
  }

  %feature("autodoc", "GateS(int i)") GateS(int i);
  void GateS(int i)
  {
  	qcs_quantum_reg_s_n(self, i);
  }	

  %feature("autodoc", "GateSAdj(int i)") GateSAdj(int i);
  void GateSAdj(int i)
  {
  	qcs_quantum_reg_s_adj_n(self, i);
  }	
   
  %feature("autodoc", "GateV(int i)") GateV(int i);
  void GateV(int i)
  {
  	qcs_quantum_reg_v_n(self, i);
  }	

  %feature("autodoc", "HadN(int i)") HadN(int i);
  void HadN(int i)
  {
    qcs_quantum_reg_had_n(self, i);
  }

  %feature("autodoc", "HadN_Conj(int i)") HadN_Conj(int i);
  void HadN_Conj(int i)
  {
    qcs_quantum_reg_had_n_conj(self, i);
  }

  %feature("autodoc", "HadD_N(int i)") HadD_N(int i);
  void HadD_N(int i)
  {
    qcs_quantum_reg_had_d_n(self, i);
  }

  %feature("autodoc", "HadD_N_Conj(int i)") HadD_N_Conj(int i);
  void HadD_N_Conj(int i)
  {
    qcs_quantum_reg_had_d_n(self, i);
  }
  
  %feature("autodoc", "SquareRootN(int i)") SquareRootN(int i);
  void SquareRootOfNotN(int i)
  {
    qcs_quantum_reg_square_root_not_n(self, i);
  }

  %feature("autodoc", "Phase(int i)") Phase(int i);
  void PhaseN(int i) 
  {
    qcs_quantum_reg_phase_n(self, i);
  }

  %feature("autodoc", "PhaseFN(int i)") PhaseFN(int i);
  void PhaseFN(int i)
  {
    qcs_quantum_reg_phase_f_n(self, i);
  }

  %feature("autodoc", "RotAlphaN(int i, float alpha)") RotAlphaN(int i, float alpha);
  void RotAlphaN(int i, float alpha)
  {
  	qcs_quantum_reg_rotate_alpha_n(self, alpha, i);
  }	

  %feature("autodoc", "RotThetaN(int i, float theta)") RotThetaN(int i, float theta);
  void RotThetaN(int i, float theta)
  {
  	qcs_quantum_reg_rotate_theta_n( self, theta, i );
  }	
  
  %feature("autodoc", "XRot90N(int i)") XRot90N(int i);
  void XRot90N(int i)
  {
    qcs_quantum_reg_x_rot90_n(self, i);
  }
  
  %feature("autodoc", "YRot90N(int i)") YRot90N(int i);
  void YRot90N(int i)
  {
    qcs_quantum_reg_y_rot90_n(self, i);
  }
  
  %feature("autodoc", "ZRot90N(int i)") ZRot90N(int i);
  void ZRot90N(int i)
  {
    qcs_quantum_reg_z_rot90_n(self, i);
  }

  %feature("autodoc", "MXRot90N(int i)") MXRot90N(int i);
  void MXRot90N(int i)
  {
    qcs_quantum_reg_mx_rot90_n(self, i);
  }
  
  %feature("autodoc", "MYRot90N(int i)") MYRot90N(int i);
  void MYRot90N(int i)
  {
    qcs_quantum_reg_my_rot90_n(self, i);
  }
  
  %feature("autodoc", "MZRot90N(int i)") MZRot90N(int i);
  void MZRot90N(int i)
  {
    qcs_quantum_reg_mz_rot90_n(self, i);
  }

  %feature("autodoc", "Had()") Had();
  void Had()
  {
    qcs_quantum_reg_had(self);
  }
 
  %feature("autodoc", "ShiftUp()") ShiftUp();
  void ShiftUp()
  {
    qcs_quantum_reg_shift_up(self);
  }

  %feature("autodoc", "ShiftDown()") ShiftDown();
  void ShiftDown()
  {
    qcs_quantum_reg_shift_down(self);
  }

  %feature("autodoc", "ArbitrarySingleGate(int i, Matrix gate)") ArbitrarySingleGate(int i, Matrix gate);
  void ArbitrarySingleGate(int i, Matrix gate)
  {
  	qcs_arbitrary_single_gate(self, &gate, i);
  }

  %feature("autodoc", "ArbitrarySingleOperation(int i, Matrix gate)") ArbitrarySingleOperation(int i, Matrix gate);
  void ArbitrarySingleOperation(int i, Matrix gate)
  {
  	qcs_arbitrary_single_gate(self, &gate, i);
  }

  %feature("autodoc", "AmplitudeDamping(int t, float p)") AmplitudeDamping(int t, float p);
  void AmplitudeDamping(int t, float p)
  {
	qcs_quantum_reg_apply_amplitude_damping(self, t, p);
  }

  %feature("autodoc", "PhaseDamping(int t, float p)") PhaseDamping(int t, float p);
  void PhaseDamping(int t, float p)
  {
	qcs_quantum_reg_apply_phase_damping(self, t, p);
  }
  
  %feature("autodoc", "PhaseDampingV2(int t, float p)") PhaseDampingV2(int t, float p);
  void PhaseDampingV2(int t, float p)
  {
	qcs_quantum_reg_apply_phase_damping_v2(self, t, p);
  }
  
  %feature("autodoc", "FullyCorrelatedPhaseFlip(float p)") FullyCorrelatedPhaseFlip(float p);
  void FullyCorrelatedPhaseFlip(float p)
  {
	qcs_quantum_reg_apply_fully_correlated_phase_flip(self, p);
  }
  
  %feature("autodoc", "ArbitraryOneQubitGate(int i, Matrix gate)") ArbitraryOneQubitGate(int i, Matrix gate);
  void ArbitraryOneQubitGate(int i, Matrix gate)
  {
  	qcs_qubit_arbitrary_one_qubit_gate(self, &gate, i);
  }

  %feature("autodoc", "Arbitrary2QubitGate(int c1, int t, Matrix gate)") Arbitrary2QubitGate(int c1, int t, Matrix gate);
  void Arbitrary2QubitGate(int c1, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_2_qubit_gate_one_control(self, c1, t, &gate);
  }

  %feature("autodoc", "Arbitrary3QubitGate(int c1, int c2, int t, Matrix gate)") Arbitrary3QubitGate(int c1, int c2, int t, Matrix gate);
  void Arbitrary3QubitGate(int c1, int c2, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_3_qubit_gate_one_control(self, c1, c2, t, &gate);
  }

  %feature("autodoc", "Arbitrary4QubitGate(int c1, int c2, int c3, int t, Matrix gate)") Arbitrary4QubitGate(int c1, int c2, int c3, int t, Matrix gate);
  void Arbitrary4QubitGate(int c1, int c2, int c3, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_4_qubit_gate_one_control(self, c1, c2, c3, t, &gate);
  }

  %feature("autodoc", "Arbitrary5QubitGate(int c1, int c2, int c3, int c4, int t, Matrix gate)") Arbitrary5QubitGate(int c1, int c2, int c3, int c4, int t, Matrix gate);
  void Arbitrary5QubitGate(int c1, int c2, int c3, int c4, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_5_qubit_gate_one_control(self, c1, c2, c3, c4, t, &gate);
  }

  %feature("autodoc", "Arbitrary6QubitGate(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate)") Arbitrary6QubitGate(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate);
  void Arbitrary6QubitGate(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_6_qubit_gate_one_control(self, c1, c2, c3, c4, c5, t, &gate);
  }


  %feature("autodoc", "Arbitrary2QubitGateZero(int c1, int t, Matrix gate)") Arbitrary2QubitGateZero(int c1, int t, Matrix gate);
  void Arbitrary2QubitGateZero(int c1, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_2_qubit_gate_zero_control(self, c1, t, &gate);
  }

  %feature("autodoc", "Arbitrary3QubitGateZero(int c1, int c2, int t, Matrix gate)") Arbitrary3QubitGateZero(int c1, int c2, int t, Matrix gate);
  void Arbitrary3QubitGateZero(int c1, int c2, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_3_qubit_gate_zero_control(self, c1, c2, t, &gate);
  }

  %feature("autodoc", "Arbitrary4QubitGateZero(int c1, int c2, int c3, int t, Matrix gate)") Arbitrary4QubitGateZero(int c1, int c2, int c3, int t, Matrix gate);
  void Arbitrary4QubitGateZero(int c1, int c2, int c3, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_4_qubit_gate_zero_control(self, c1, c2, c3, t, &gate);
  }

  %feature("autodoc", "Arbitrary5QubitGateZero(int c1, int c2, int c3, int c4, int t, Matrix gate)") Arbitrary5QubitGateZero(int c1, int c2, int c3, int c4, int t, Matrix gate);
  void Arbitrary5QubitGateZero(int c1, int c2, int c3, int c4, int t, Matrix gate)
  {
  	qcs_qubit_arbitrary_5_qubit_gate_zero_control(self, c1, c2, c3, c4, t, &gate);
  }

  %feature("autodoc", "Arbitrary6QubitGateZero(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate)") Arbitrary6QubitGateZero(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate);
  void Arbitrary6QubitGateZero(int c1, int c2, int c3, int c4, int c5, int t, Matrix gate)
  {
		qcs_qubit_arbitrary_6_qubit_gate_zero_control(self, c1, c2, c3, c4, c5, t, &gate);
  }

  %feature("autodoc", "CNot(int c1, int t)") CNot(int c1, int t);
  void CNot(int c1, int t)
  {
		qcs_quantum_reg_cnot(self, c1, t);
  }

  %feature("autodoc", "CNot_Conj(int c1, int t)") CNot_Conj(int c1, int t);
  void CNot_Conj(int c1, int t)
  {
		qcs_quantum_reg_cnot_conj(self, c1, t);
  }

  %feature("autodoc", "CNot1(int c1, int t)") CNot1(int c1, int t);
  void CNot1(int c1, int t)
  {
		qcs_quantum_reg_cnot1(self, c1, t);
  }

  %feature("autodoc", "CNot2(int c1, int c2, int t)") CNot2(int c1, int c2, int t);
  void CNot2(int c1, int c2, int t)
  {
		qcs_quantum_reg_cnot2(self, c1, c2, t);
  }

  %feature("autodoc", "CNot3(int c1, int c2, int c3, int t)") CNot3(int c1, int c2, int c3, int t);
  void CNot3(int c1, int c2, int c3, int t)
  {
		qcs_quantum_reg_cnot3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CNot4(int c1, int c2, int c3, int c4, int t)") CNot4(int c1, int c2, int c3, int c4, int t);
  void CNot4(int c1, int c2, int c3, int c4, int t)
  {
		qcs_quantum_reg_cnot4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CNot5(int c1, int c2, int c3, int c4, int c5, int t)") CNot5(int c1, int c2, int c3, int c4, int c5, int t);
  void CNot5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cnot5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CSquareRootOfNot(int c1, int t)") CSquareRootOfNot(int c1, int t);
  void CSquareRootOfNot(int c1, int t)
  {
	qcs_quantum_reg_csquare_root_of_not(self, c1, t);
  }
  
  %feature("autodoc", "CPauliX(int c1, int t)") CPauliX(int c1, int t);
  void CPauliX(int c1, int t)
  {
    qcs_quantum_reg_cpauli_x(self, c1, t);
  }

  %feature("autodoc", "CPauliX1(int c1, int t)") CPauliX1(int c1, int t);
  void CPauliX1(int c1, int t)
  {
    qcs_quantum_reg_cpauli_x1(self, c1, t);
  }

  %feature("autodoc", "CPauliX2(int c1, int c2, int t)") CPauliX2(int c1, int c2, int t);
  void CPauliX2(int c1, int c2, int t)
  {
    qcs_quantum_reg_cpauli_x2(self, c1, c2, t);
  }

  %feature("autodoc", "CPauliX3(int c1, int c2, int c3, int t)") CPauliX3(int c1, int c2, int c3, int t);
  void CPauliX3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cpauli_x3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CPauliX4(int c1, int c2, int c3, int c4, int t)") CPauliX4(int c1, int c2, int c3, int c4, int t);
  void CPauliX4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cpauli_x4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CPauliX5(int c1, int c2, int c3, int c4, int c5, int t)") CPauliX5(int c1, int c2, int c3, int c4, int c5, int t);
  void CPauliX5(int c1, int c2, int c3, int c4, int c5, int t)
  { 
    qcs_quantum_reg_cpauli_x5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CPauliY(int c1, int t)") CPauliY(int c1, int t);
  void CPauliY(int c1, int t)
  {
    qcs_quantum_reg_cpauli_y(self, c1, t);
  }

  %feature("autodoc", "CPauliY1(int c1, int t)") CPauliY1(int c1, int t);
  void CPauliY1(int c1, int t)
  {
    qcs_quantum_reg_cpauli_y1(self, c1, t);
  }

  %feature("autodoc", "CPauliY2(int c1, int c2, int t)") CPauliY2(int c1, int c2, int t);
  void CPauliY2(int c1, int c2, int t)
  {
    qcs_quantum_reg_cpauli_y2(self, c1, c2, t);
  }

  %feature("autodoc", "CPauliY3(int c1, int c2, int c3, int t)") CPauliY3(int c1, int c2, int c3, int t);
  void CPauliY3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cpauli_y3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CPauliY4(int c1, int c2, int c3, int c4, int t)") CPauliY4(int c1, int c2, int c3, int c4, int t);
  void CPauliY4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cpauli_y4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CPauliY5(int c1, int c2, int c3, int c4, int c5, int t)") CPauliY5(int c1, int c2, int c3, int c4, int c5, int t);
  void CPauliY5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cpauli_y5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CPauliZ(int c1, int t)") CPauliZ(int c1, int t);
  void CPauliZ(int c1, int t)
  {
    qcs_quantum_reg_cpauli_z(self, c1, t);
  }

  %feature("autodoc", "CPauliZ1(int c1, int t)") CPauliZ1(int c1, int t);
  void CPauliZ1(int c1, int t)
  {
    qcs_quantum_reg_cpauli_z1(self, c1, t);
  }

  %feature("autodoc", "CPauliZ2(int c1, int c2, int t)") CPauliZ2(int c1, int c2, int t);
  void CPauliZ2(int c1, int c2, int t)
  {
    qcs_quantum_reg_cpauli_z2(self, c1, c2, t);
  }

  %feature("autodoc", "CPauliZ3(int c1, int c2, int c3, int t)") CPauliZ3(int c1, int c2, int c3, int t);
  void CPauliZ3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cpauli_z3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CPauliZ4(int c1, int c2, int c3, int c4, int t)") CPauliZ4(int c1, int c2, int c3, int c4, int t);
  void CPauliZ4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cpauli_z4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CPauliZ5(int c1, int c2, int c3, int c4, int c5, int t)") CPauliZ5(int c1, int c2, int c3, int c4, int c5, int t);
  void CPauliZ5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cpauli_z5(self, c1, c2, c3, c4, c5, t);
  }
  
  %feature("autodoc", "CGateV(int c1, int t)") CGateV(int c1, int t);
  void CGateV(int c1, int t)
  {
    qcs_quantum_reg_cgate_v1(self, c1, t);
  }

  %feature("autodoc", "CGateV1(int c1, int t)") CGateV1(int c1, int t);
  void CGateV1(int c1, int t)
  {
    qcs_quantum_reg_cgate_v1(self, c1, t);
  }

  %feature("autodoc", "CGateV2(int c1, int c2, int t)") CGateV2(int c1, int c2, int t);
  void CGateV2(int c1, int c2, int t)
  {
    qcs_quantum_reg_cgate_v2(self, c1, c2, t);
  }

  %feature("autodoc", "CGateV3(int c1, int c2, int c3, int t)") CGateV3(int c1, int c2, int c3, int t);
  void CGateV3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cgate_v3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CGateV4(int c1, int c2, int c3, int c4, int t)") CGateV4(int c1, int c2, int c3, int c4, int t);
  void CGateV4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cgate_v4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CGateV5(int c1, int c2, int c3, int c4, int c5, int t)") CGateV5(int c1, int c2, int c3, int c4, int c5, int t);
  void CGateV5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cgate_v5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CGateS(int c1, int t)") CGateS(int c1, int t);
  void CGateS(int c1, int t)
  {
    qcs_quantum_reg_cgate_s1(self, c1, t);
  }

  %feature("autodoc", "CGateS1(int c1, int t)") CGateS1(int c1, int t);
  void CGateS1(int c1, int t)
  {
    qcs_quantum_reg_cgate_s1(self, c1, t);
  }

  %feature("autodoc", "CGateS2(int c1, int c2, int t)") CGateS2(int c1, int c2, int t);
  void CGateS2(int c1, int c2, int t)
  {
    qcs_quantum_reg_cgate_s2(self, c1, c2, t);
  }

  %feature("autodoc", "CGateS3(int c1, int c2, int c3, int t)") CGateS3(int c1, int c2, int c3, int t);
  void CGateS3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cgate_s3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CGateS4(int c1, int c2, int c3, int c4, int t)") CGateS4(int c1, int c2, int c3, int c4, int t);
  void CGateS4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cgate_s4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CGateS5(int c1, int c2, int c3, int c4, int c5, int t)") CGateS5(int c1, int c2, int c3, int c4, int c5, int t);
  void CGateS5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cgate_s5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CHad(int c1, int t)") CHad(int c1, int t);
  void CHad(int c1, int t)
  {
    qcs_quantum_reg_chad(self, c1, t);
  }

  %feature("autodoc", "CHad1(int c1, int t)") CHad1(int c1, int t);
  void CHad1(int c1, int t)
  {
    qcs_quantum_reg_chad1(self, c1, t);
  }

  %feature("autodoc", "CHad2(int c1, int c2, int t)") CHad2(int c1, int c2, int t);
  void CHad2(int c1, int c2, int t)
  {
    qcs_quantum_reg_chad2(self, c1, c2, t);
  }

  %feature("autodoc", "CHad3(int c1, int c2, int c3, int t)") CHad3(int c1, int c2, int c3, int t);
  void CHad3(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_chad3(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CHad4(int c1, int c2, int c3, int c4, int t)") CHad4(int c1, int c2, int c3, int c4, int t);
  void CHad4(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_chad4(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CHad5(int c1, int c2, int c3, int c4, int c5, int t)") CHad5(int c1, int c2, int c3, int c4, int c5, int t);
  void CHad5(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_chad5(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CNot_zero(int c1, int t)") CNot_zero(int c1, int t);
  void CNot_zero(int c1, int t)
  {
    qcs_quantum_reg_cnot_zero(self, c1, t);
  }

  %feature("autodoc", "CNot1_zero(int c1, int t)") CNot1_zero(int c1, int t);
  void CNot1_zero(int c1, int t)
  {
    qcs_quantum_reg_cnot_zero(self, c1, t);
  }

  %feature("autodoc", "CNot2_zero(int c1, int c2, int t)") CNot2_zero(int c1, int c2, int t);
  void CNot2_zero(int c1, int c2, int t)
  {
    qcs_quantum_reg_cnot2_zero(self, c1, c2, t);
  }

  %feature("autodoc", "CNot3_zero(int c1, int c2, int c3, int t)") CNot3_zero(int c1, int c2, int c3, int t);
  void CNot3_zero(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_cnot3_zero(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CNot4_zero(int c1, int c2, int c3, int c4, int t)") CNot4_zero(int c1, int c2, int c3, int c4, int t);
  void CNot4_zero(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_cnot4_zero(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CNot5_zero(int c1, int c2, int c3, int c4, int c5, int t)") CNot5_zero(int c1, int c2, int c3, int c4, int c5, int t);
  void CNot5_zero(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_cnot5_zero(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CHad_zero(int c1, int t)") CHad_zero(int c1, int t);
  void CHad_zero(int c1, int t)
  {
    qcs_quantum_reg_chad_zero(self, c1, t);
  }

  %feature("autodoc", "CHad1_zero(int c1, int t)") CHad1_zero(int c1, int t);
  void CHad1_zero(int c1, int t)
  {
    qcs_quantum_reg_chad1_zero(self, c1, t);
  }

  %feature("autodoc", "CHad2_zero(int c1, int c2, int t)") CHad2_zero(int c1, int c2, int t);
  void CHad2_zero(int c1, int c2, int t)
  {
    qcs_quantum_reg_chad2_zero(self, c1, c2, t);
  }

  %feature("autodoc", "CHad3_zero(int c1, int c2, int c3, int t)") CHad3_zero(int c1, int c2, int c3, int t);
  void CHad3_zero(int c1, int c2, int c3, int t)
  {
    qcs_quantum_reg_chad3_zero(self, c1, c2, c3, t);
  }

  %feature("autodoc", "CHad4_zero(int c1, int c2, int c3, int c4, int t)") CHad4_zero(int c1, int c2, int c3, int c4, int t);
  void CHad4_zero(int c1, int c2, int c3, int c4, int t)
  {
    qcs_quantum_reg_chad4_zero(self, c1, c2, c3, c4, t);
  }

  %feature("autodoc", "CHad5_zero(int c1, int c2, int c3, int c4, int c5, int t)") CHad5_zero(int c1, int c2, int c3, int c4, int c5, int t);
  void CHad5_zero(int c1, int c2, int c3, int c4, int c5, int t)
  {
    qcs_quantum_reg_chad5_zero(self, c1, c2, c3, c4, c5, t);
  }

  %feature("autodoc", "CNot2_any(int c1, int c2, int t, char b1, char b2)") CNot2_any(int c1, int c2, int t, char b1, char b2);
  void CNot2_any(int c1, int c2, int t, char b1, char b2)
  {
    qcs_quantum_reg_cnot2_any(self, c1, c2, t, b1, b2);
  }
  
  %feature("autodoc", "CRot45(int c, int t)") CRot45(int c, int t);
  void CRot45(int c, int t)
  {
    qcs_quantum_reg_crot45(self, c, t);
  }

  %feature("autodoc", "CRot90(int c, int t)") CRot90(int c, int t);
  void CRot90(int c, int t)
  {
    qcs_quantum_reg_crot90(self, c, t);
  }
  
  %feature("autodoc", "CRotAlpha(int c, int t, float alpha)") CRotAlpha(int c, int t, float alpha);
  void CRotAlpha(int c, int t, float alpha)
  {
    qcs_quantum_reg_crot_alpha(self, c, t, alpha);
  }

  %feature("autodoc", "CRotTheta(int c1, int t, float theta)") CRotTheta(int c1, int t, float theta);
  void CRotTheta(int c1, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta(self, c1, t, theta);
  }
  
 %feature("autodoc", "CRotTheta1(int c1, int t, float theta)") CRotTheta1(int c1, int t, float theta);
  void CRotTheta1(int c1, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta1(self, c1, t, theta);
  }

  %feature("autodoc", "CRotTheta2(int c1, int c2, int t, float theta)") CRotTheta2(int c1, int c2, int t, float theta);
  void CRotTheta2(int c1, int c2, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta2(self, c1, c2, t, theta);
  }

  %feature("autodoc", "CRotTheta3(int c1, int c2, int c3, int t, float theta)") CRotTheta3(int c1, int c2, int c3, int t, float theta);
  void CRotTheta3(int c1, int c2, int c3, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta3(self, c1, c2, c3, t, theta);
  }

  %feature("autodoc", "CRotTheta4(int c1, int c2, int c3, int c4, int t, float theta)") CRotTheta4(int c1, int c2, int c3, int c4, int t, float theta);
  void CRotTheta4(int c1, int c2, int c3, int c4, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta4(self, c1, c2, c3, c4, t, theta);
  }

  %feature("autodoc", "CRotTheta5(int c1, int c2, int c3, int c4, int c5, int t, float theta)") CRotTheta5(int c1, int c2, int c3, int c4, int c5, int t, float theta);
  void CRotTheta5(int c1, int c2, int c3, int c4, int c5, int t, float theta)
  {
  	qcs_quantum_reg_crot_theta5(self, c1, c2, c3, c4, c5, t, theta);
  }
 
  %feature("autodoc", "CPhaseF(int c1, int t)") CPhaseF(int c1, int t);
  void CPhaseF(int c1, int t)
  {
  	qcs_quantum_reg_cphase_f_n(self, c1, t);
  }

  %feature("autodoc", "CPhaseF1(int c1, int t)") CPhaseF1(int c1, int t);
  void CPhaseF1(int c1, int t)
  {
  	qcs_quantum_reg_cphase_f_n1(self, c1, t);
  }

  %feature("autodoc", "CPhaseF2(int c1, int c2, int t)") CPhaseF2(int c1, int c2, int t);
  void CPhaseF2(int c1, int c2, int t)
  {
  	qcs_quantum_reg_cphase_f_n2(self, c1, c2, t);
  }

  %feature("autodoc", "CPhaseF3(int c1, int c2, int c3, int t)") CPhaseF3(int c1, int c2, int c3, int t);
  void CPhaseF3(int c1, int c2, int c3, int t)
  {
  	qcs_quantum_reg_cphase_f_n3( self, c1, c2, c3, t );
  }
  
  %feature("autodoc", "CPhaseF4(int c1, int c2, int c3, int c4, int t)") CPhaseF4(int c1, int c2, int c3, int c4, int t);
  void CPhaseF4(int c1, int c2, int c3, int c4, int t)
  {
  	qcs_quantum_reg_cphase_f_n4( self, c1, c2, c3, c4, t );
  }

  %feature("autodoc", "CPhaseF5(int c1, int c2, int c3, int c4, int c5, int t)") CPhaseF5(int c1, int c2, int c3, int c4, int c5, int t);
  void CPhaseF5(int c1, int c2, int c3, int c4, int c5, int t)
  {
  	qcs_quantum_reg_cphase_f_n5( self, c1, c2, c3, c4, c5, t );
  }
  
  %feature("autodoc", "void Swap(int a, int b)") Swap(int a, int b);
  void Swap(int a, int b)
  {
    qcs_quantum_reg_cnot(self, a, b);
    qcs_quantum_reg_cnot(self, b, a);
    qcs_quantum_reg_cnot(self, a, b);
  }
 
  %feature("autodoc", "Toffoli(int c1, int c2, int t)") Toffoli(int c1, int c2, int t);
  void Toffoli(int c1, int c2, int t)
  {
    qcs_quantum_reg_cnot2(self, c1, c2, t);
  }
  
  %feature("autodoc", "Fredkin(int a, int b, int c)") Fredkin(int a, int b, int c);
  void Fredkin(int a, int b, int c)
  {
  }

  %feature("autodoc", "MatrixApply(Matrix a)") MatrixApply(Matrix a);
  void MatrixApply(Matrix a)
  {
		qcs_quantum_reg_mult_by_matrix(self, &a);
  }	
   
  %feature("autodoc", "Pr()") Pr();
  void Pr()
  {
    if(self->mode == CHP_MODE)
    	qcs_chp_PrintState32(self->chp_state);
    	
	if(self->mode == USE_STATE_VECTOR_QUBIT || self->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
		qcs_quantum_reg_print_bin(self);
		  
	if(self->mode == USE_STATE_VECTOR_QUDIT || self->mode == USE_SYMBOLIC_STATE_VECTOR_QUDIT)
		qcs_quantum_reg_print_d_base(self);
	
	if(self->mode == USE_DENSITY_MATRIX)
		qcs_print_matrix( self->dms->state );	
  }

  %feature("autodoc", "PrHorz()") PrHorz();
  void PrHorz()
  {    	
	if(self->mode == USE_STATE_VECTOR_QUBIT || self->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT)
		qcs_quantum_reg_print_bin_horz(self);		  
  }  
  
  %feature("autodoc", "PrFull()") PrFull();
  void PrFull()
  {
	if(self->mode == USE_STATE_VECTOR_QUBIT || self->mode == USE_SYMBOLIC_STATE_VECTOR_QUBIT )
    	qcs_quantum_reg_print_bin_full(self);  
	if(self->mode == USE_STATE_VECTOR_QUDIT)
		qcs_quantum_reg_print_d_base_full(self);

	if(self->mode == USE_DENSITY_MATRIX)
		qcs_print_matrix( self->dms->state );
  }
  
  %feature("autodoc", "PrDenMat()") PrDenMat();
  void PrDenMat()
  {
  	qcs_quantum_reg_print_den_matrix(self);
  }

  %feature("autodoc", "GenDenMat()") GenDenMat();
  Matrix GenDenMat()
  {
		return *qcs_quantum_reg_generate_density_matrix(self);
  }

  %feature("autodoc", "SchmidtDecomposition( int size_x, int size_y)") SchmidtDecomposition( int size_x, int size_y);
  SchmidtDecomposition SchmidtDecomposition( int size_x, int size_y)
  {
		tf_qcs_matrix *_tmp_base1, *_tmp_base2, *_tmp_schmidt_coeff;
		tf_qcs_complex *v;
		SchmidtDecomposition *sd;
		int cc, i;
		
		sd = qcs_create_schmidt_decomposition_empty();
		
		_tmp_base1 = qcs_create_empty_matrix();
		_tmp_base2 = qcs_create_empty_matrix();
		_tmp_schmidt_coeff = qcs_create_empty_matrix();
  
		qcs_quantum_reg_schmidt_decompositon( self, size_x, size_y, _tmp_base1, _tmp_base2, _tmp_schmidt_coeff);
		
	    sd->base1 = _tmp_base1;
		sd->base2 = _tmp_base2;
		sd->schmidt_coeff = _tmp_schmidt_coeff;
		
		cc=0;
        for(i=0;i<_tmp_schmidt_coeff->cols;i++)
        {
            v=qcs_get_cell_at_matrix_complex( _tmp_schmidt_coeff, 0, i);
			if (v->re!=0.0 || v->im!=0.0) cc++;
        }

		
		if (cc > 1) 
			sd->is_entangled = 1;
		else
			sd->is_entangled = 0;
			
		sd->schmidt_coeff_number = cc;
		
		return *sd;
  }
  
  %feature("autodoc", "GroverChangeSignGen(int eigenstate)") GroverChangeSignGen(int eigenstate);
  Matrix GroverChangeSignGen(int eigenstate)
  {
		return *change_sign_matrix_syntesis(self->size, eigenstate);
  } 

  %feature("autodoc", "GroverChangeSignGen2(int eigenstate1, int eigenstate2)") GroverChangeSignGen2(int eigenstate1, int eigenstate2);
  Matrix GroverChangeSignGen2(int eigenstate1, int eigenstate2)
  {
	return *change_sign_matrix_syntesis2(self->size, eigenstate1, eigenstate2);
  } 
   
  %feature("autodoc", "GroverTurnAroundMean()") GroverTurnAroundMean();
  Matrix GroverTurnAroundMean()
  {
	return *inverse_around_mean_matrix_syntesis(self->size);
  }
  
  %feature("autodoc", "CorrMatForTeleport(int n, int m)") CorrMatForTeleport(int n, int m);  
  Matrix CorrMatForTeleport(int n, int m)
  {
	if(self->mode == USE_STATE_VECTOR_QUBIT)
		return *qudit_d_level_teleportation_correction_matrix( n, m, 2 );

	if(self->mode == USE_STATE_VECTOR_QUDIT)
		return *qudit_d_level_teleportation_correction_matrix( n, m, self->freedom_level );
  }

  %feature("autodoc", "RemoveQubit(int t)") RemoveQubit(int t);
  QuantumReg RemoveQubit(int t)
  {
	return *qcs_quantum_reg_remove_one_qubit(self, t);
  }
  
  %feature("autodoc", "PrTr()") PrTr();  
  void PrTr()
  {
  	qcs_quantum_reg_print_trace(self);
  }

  %feature("autodoc", "PrPurity()") PrPurity();
  void PrPurity()
  {
  	qcs_quantum_reg_print_purity(self);
  }
	
	%feature("autodoc", "PrDetail()") PrDetail();
	void PrDetail()
	{
		qcs_quantum_reg_print_detail(self);
	}
	
  %feature("autodoc", "PrSqr()") PrSqr();
  void PrSqr()
  {
    qcs_quantum_reg_print_bin_sqr(self);  
  }

  %feature("autodoc", "PrBinDec()") PrBinDec();
  void PrBinDec()
  {
    qcs_quantum_reg_print_bin_and_dec(self);
  }

  %feature("autodoc", "PrBinDecSqr()") PrBinDecSqr();
  void PrBinDecSqr()
  {
    qcs_quantum_reg_print_bin_and_dec_sqr(self);
  }

  %feature("autodoc", "SaveToFileFull(char *fname)") SaveToFileFull(char *fname);
  void SaveToFileFull(char *fname)
  {
    qcs_quantum_reg_save_to_file_full(self, fname);
  }

  %feature("autodoc", "SaveToFileWithRemarkFull(char *fname, char *remark_txt)") SaveToFileWithRemarkFull(char *fname, char *remark_txt);
  void SaveToFileWithRemarkFull(char *fname, char *remark_txt)
  {
    qcs_quantum_reg_save_to_file_with_remark_full(self, fname, remark_txt);
  }

  %feature("autodoc", "SaveToFile(char *fname)") SaveToFile(char *fname);
  void SaveToFile(char *fname)
  {
    qcs_quantum_reg_save_to_file(self, fname);
  }

  %feature("autodoc", "SaveToFileWithRemark(char *fname, char *remark_txt)") SaveToFileWithRemark(char *fname, char *remark_txt);
  void SaveToFileWithRemark(char *fname, char *remark_txt)
  {
    qcs_quantum_reg_save_to_file_with_remark(self, fname, remark_txt);
  }

  %feature("autodoc", "PrKet()") PrKet();
  void PrKet()
  {
    qcs_quantum_reg_qubits_ket_print_measure(self);  
  }
  
  %feature("autodoc", "PrBin()") PrBin();
  void PrBin()
  {
     qcs_quantum_reg_print_bin(self);  
  }

  %feature("autodoc", "PrDec()") PrDec();
  void PrDec()
  {
      qcs_quantum_reg_print_dec(self);  
  }

  %feature("autodoc", "ProbeQubitStdBase(int i) -> [p0,p1]") ProbeQubitStdBase(int i, float *p0, float *p1);  
  %apply float *OUTPUT { float *p0, float *p1 };
  void ProbeQubitStdBase(int i, float *p0, float *p1)
  {
		float _tmp_p0, _tmp_p1;
		qcs_quantum_reg_probe_one_qubit_in_std_base( self, i, &_tmp_p0, &_tmp_p1);
		
		*p0 = _tmp_p0;
		*p1 = _tmp_p1;
  }
  
  %feature("autodoc", "Measure() -> int") Measure();    
  int Measure()
  {	
#ifdef PYTHON_SCRIPT
		PySys_WriteStdout("Please, write me!\n");
#endif
	return -1;
  }
  
	%feature("autodoc", "MeasureN(int _from, int _to) -> int") MeasureN(int _from, int _to);
	int MeasureN(int _from, int _to)
	{
		return qcs_quantum_reg_measure_from_to(self, _from, _to);
	}

    %feature("autodoc", "MeasureOneQubit(int k) -> int") MeasureOneQubit(int k);
	int MeasureOneQubit(int k)
	{
		return qcs_quantum_reg_measure_one_qubit(self, k);
	}

    %feature("autodoc", "MeasureOneQubitInStdBase(int k) -> int") MeasureOneQubitInStdBase(int k);
	int MeasureOneQubitInStdBase(int k)
	{
		return qcs_quantum_reg_measure_one_qubit_in_std_base(self, k);
	}

	%feature("autodoc", "MeasureOneQubitInBBase(int k, float a) -> int") MeasureOneQubitInBBase(int k, float a);
	int MeasureOneQubitInBBase(int k, float a)
	{
		int s;
		tf_qcs_qubit_base_desc base;
		
		qcs_zero_base( &base );
		
	    qcs_make_b_base(&base, a);
		s=qcs_quantum_reg_measure_one_qubit_in_base(self, k, &base);
		
		return s;
	}
	
	%feature("autodoc", "MeasureOneQubitForce(int k, int force_value) -> int") MeasureOneQubitForce(int k, int force_value);
	int MeasureOneQubitForce(int k, int force_value)
	{
		return qcs_quantum_reg_measure_one_qubit_force(self, k, force_value);
	}
	
	%feature("autodoc", "MeasureOneQudit(int k) -> int") MeasureOneQudit(int k);
	int MeasureOneQudit(int k)
	{
		return qcs_quantum_reg_measure_one_qudit(self, k);
	}

	%feature("autodoc", "MeasureOneQuditForce(int k, int force_value) -> int") MeasureOneQuditForce(int k, int force_value);
	int MeasureOneQuditForce(int k, int force_value)
	{
		return qcs_quantum_reg_measure_one_qudit_force(self, k, force_value);
	}

	%feature("autodoc", "MeasureOneQubitInBase(int k, QubitBaseDesc base) -> int") MeasureOneQubitInBase(int k, QubitBaseDesc base);
	int MeasureOneQubitInBase(int k, QubitBaseDesc base)
	{
		return qcs_quantum_reg_measure_one_qubit_in_base(self, k, &base);
	}

	%feature("autodoc", "MemAllocateBytes() -> int") MemAllocateBytes();
	unsigned int MemAllocateBytes()
	{
		return qcs_qubit_vec_state_size_in_bytes(self);
	}

}
/* end of %extend QuantumReg */
