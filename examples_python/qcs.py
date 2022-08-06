# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

"""
Quantum Computing Simulator (QCS) port of ANSI C/C++/FORTRAN library of quantum
computations routines for Python and other languages supported by SWIG.
"""

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _qcs
else:
    import _qcs

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


USE_NO_STATE_VECTOR = _qcs.USE_NO_STATE_VECTOR
USE_STATE_VECTOR_QUBIT = _qcs.USE_STATE_VECTOR_QUBIT
USE_STATE_VECTOR_QUDIT = _qcs.USE_STATE_VECTOR_QUDIT
USE_DENSITY_MATRIX = _qcs.USE_DENSITY_MATRIX
USE_GRAPH_STATE_DESC = _qcs.USE_GRAPH_STATE_DESC
USE_CHP_MODE = _qcs.USE_CHP_MODE
USE_ONEWAY_MODEL = _qcs.USE_ONEWAY_MODEL
USE_PQC_MODE = _qcs.USE_PQC_MODE
USE_STATE_VECTOR_MULTI_QUBITSQUDITS = _qcs.USE_STATE_VECTOR_MULTI_QUBITSQUDITS
USE_SYMBOLIC_STATE_VECTOR_QUBIT = _qcs.USE_SYMBOLIC_STATE_VECTOR_QUBIT
USE_SYMBOLIC_STATE_VECTOR_QUDIT = _qcs.USE_SYMBOLIC_STATE_VECTOR_QUDIT
ERROR_BAD_QUBIT_NUMBER = _qcs.ERROR_BAD_QUBIT_NUMBER
class QuantumRegister(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    n = property(_qcs.QuantumRegister_n_get, _qcs.QuantumRegister_n_set)
    vec_state_size = property(_qcs.QuantumRegister_vec_state_size_get, _qcs.QuantumRegister_vec_state_size_set)
    fdl = property(_qcs.QuantumRegister_fdl_get, _qcs.QuantumRegister_fdl_set)
    mode = property(_qcs.QuantumRegister_mode_get, _qcs.QuantumRegister_mode_set)
    el = property(_qcs.QuantumRegister_el_get, _qcs.QuantumRegister_el_set)
    dim_of_qs = property(_qcs.QuantumRegister_dim_of_qs_get, _qcs.QuantumRegister_dim_of_qs_set)
    vs = property(_qcs.QuantumRegister_vs_get, _qcs.QuantumRegister_vs_set)

    def __init__(self, *args):
        r"""QuantumRegister()QuantumRegister(int size)"""
        _qcs.QuantumRegister_swiginit(self, _qcs.new_QuantumRegister(*args))
    __swig_destroy__ = _qcs.delete_QuantumRegister

    def Reset(self):
        r"""Reset()"""
        return _qcs.QuantumRegister_Reset(self)

    def SetGHZState(self):
        r"""SetGHZState()"""
        return _qcs.QuantumRegister_SetGHZState(self)

    def SetGHZ01State(self):
        r"""SetGHZ01State()"""
        return _qcs.QuantumRegister_SetGHZ01State(self)

    def X(self, t):
        r"""X(int t)"""
        return _qcs.QuantumRegister_X(self, t)

    def Y(self, t):
        r"""Y(int t)"""
        return _qcs.QuantumRegister_Y(self, t)

    def Z(self, t):
        r"""Z(int t)"""
        return _qcs.QuantumRegister_Z(self, t)

    def MXRot90N(self, t):
        r"""MXRot90N(int t)"""
        return _qcs.QuantumRegister_MXRot90N(self, t)

    def MYRot90N(self, t):
        r"""MYRot90N(int t)"""
        return _qcs.QuantumRegister_MYRot90N(self, t)

    def MZRot90N(self, t):
        r"""MZRot90N(int t)"""
        return _qcs.QuantumRegister_MZRot90N(self, t)

    def XRot90N(self, t):
        r"""XRot90N(int t)"""
        return _qcs.QuantumRegister_XRot90N(self, t)

    def YRot90N(self, t):
        r"""YRot90N(int t)"""
        return _qcs.QuantumRegister_YRot90N(self, t)

    def ZRot90N(self, t):
        r"""ZRot90N(int t)"""
        return _qcs.QuantumRegister_ZRot90N(self, t)

    def Had(self, t):
        r"""Had(int t)"""
        return _qcs.QuantumRegister_Had(self, t)

    def HadN(self, t):
        r"""HadN(int t)"""
        return _qcs.QuantumRegister_HadN(self, t)

    def HadAll(self):
        r"""HadAll()"""
        return _qcs.QuantumRegister_HadAll(self)

    def SquareRootOfNotN(self, t):
        r"""SquareRootOfNotN(int t)"""
        return _qcs.QuantumRegister_SquareRootOfNotN(self, t)

    def SqrtOfNotN(self, t):
        r"""SqrtOfNotN(int t)"""
        return _qcs.QuantumRegister_SqrtOfNotN(self, t)

    def CNot(self, c, t):
        r"""CNot(int c, int t)"""
        return _qcs.QuantumRegister_CNot(self, c, t)

    def Measure(self):
        r"""Measure()"""
        return _qcs.QuantumRegister_Measure(self)

    def M(self, t):
        r"""M(int t)"""
        return _qcs.QuantumRegister_M(self, t)

    def MeasureN(self, _from, _to):
        r"""MeasureN(int _from, int _to)"""
        return _qcs.QuantumRegister_MeasureN(self, _from, _to)

    def MeasureOneQubit(self, t):
        r"""MeasureOneQubit(int t) -> int"""
        return _qcs.QuantumRegister_MeasureOneQubit(self, t)

    def ProbeQubitStdBase(self, i):
        r"""ProbeQubitStdBase(int i) -> [p0,p1]"""
        return _qcs.QuantumRegister_ProbeQubitStdBase(self, i)

    def Noop(self):
        r"""Noop()"""
        return _qcs.QuantumRegister_Noop(self)

    def Pr(self):
        r"""Pr()"""
        return _qcs.QuantumRegister_Pr(self)

    def PrSqr(self):
        r"""PrSqr()"""
        return _qcs.QuantumRegister_PrSqr(self)

    def PrFull(self):
        r"""PrFull()"""
        return _qcs.QuantumRegister_PrFull(self)

    def PrFullSqr(self):
        r"""PrFullSqr()"""
        return _qcs.QuantumRegister_PrFullSqr(self)

    def PrAsMatlab(self):
        r"""PrAsMatlab()"""
        return _qcs.QuantumRegister_PrAsMatlab(self)

    def PrAsMathematica(self):
        r"""PrAsMathematica()"""
        return _qcs.QuantumRegister_PrAsMathematica(self)

# Register QuantumRegister in _qcs:
_qcs.QuantumRegister_swigregister(QuantumRegister)


def qcs_new_quantum_register(size):
    return _qcs.qcs_new_quantum_register(size)

def qcs_delete_quantum_register(q_reg):
    return _qcs.qcs_delete_quantum_register(q_reg)

def qcs_quantum_register_reset_error_level(q_reg, v):
    return _qcs.qcs_quantum_register_reset_error_level(q_reg, v)

def qcs_quantum_register_set_error_level(q_reg, v):
    return _qcs.qcs_quantum_register_set_error_level(q_reg, v)

def qcs_quantum_register_get_error_level(q_reg):
    return _qcs.qcs_quantum_register_get_error_level(q_reg)

def qcs_quantum_register_reset(q_reg):
    return _qcs.qcs_quantum_register_reset(q_reg)

def qcs_quantum_register_set_state_dec(q_reg, n):
    return _qcs.qcs_quantum_register_set_state_dec(q_reg, n)

def qcs_quantum_register_set_state_bin(q_reg, state_desc):
    return _qcs.qcs_quantum_register_set_state_bin(q_reg, state_desc)

def qcs_quantum_register_print_bin(q_reg):
    return _qcs.qcs_quantum_register_print_bin(q_reg)

def qcs_quantum_register_print_bin_in_matlab_format(q_reg):
    return _qcs.qcs_quantum_register_print_bin_in_matlab_format(q_reg)

def qcs_quantum_register_print_bin_sqr(q_reg):
    return _qcs.qcs_quantum_register_print_bin_sqr(q_reg)

def qcs_quantum_register_print_bin_full(q_reg):
    return _qcs.qcs_quantum_register_print_bin_full(q_reg)

def qcs_quantum_register_print_bin_full_sqr(q_reg):
    return _qcs.qcs_quantum_register_print_bin_full_sqr(q_reg)

def qcs_quantum_register_print_bin_with_prefix(q_reg, prefix):
    return _qcs.qcs_quantum_register_print_bin_with_prefix(q_reg, prefix)

def qcs_quantum_register_print_dec(q_reg):
    return _qcs.qcs_quantum_register_print_dec(q_reg)

def qcs_quantum_register_fill_zero(q_reg):
    return _qcs.qcs_quantum_register_fill_zero(q_reg)

def qcs_quantum_register_set_ghz_state(q_reg):
    return _qcs.qcs_quantum_register_set_ghz_state(q_reg)

def applied_1q_gate_to_quantum_register(q_reg, t, u):
    return _qcs.applied_1q_gate_to_quantum_register(q_reg, t, u)

def applied_2q_gate_to_quantum_register_one_control(q_reg, c1, t, u):
    return _qcs.applied_2q_gate_to_quantum_register_one_control(q_reg, c1, t, u)

def qcs_quantum_register_pauli_x_gate(q_reg, i):
    return _qcs.qcs_quantum_register_pauli_x_gate(q_reg, i)

def qcs_quantum_register_pauli_y_gate(q_reg, i):
    return _qcs.qcs_quantum_register_pauli_y_gate(q_reg, i)

def qcs_quantum_register_pauli_z_gate(q_reg, i):
    return _qcs.qcs_quantum_register_pauli_z_gate(q_reg, i)

def qcs_quantum_register_had_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_had_n_gate(q_reg, i)

def qcs_quantum_register_had_n_conj_gate(q_reg, i):
    return _qcs.qcs_quantum_register_had_n_conj_gate(q_reg, i)

def qcs_quantum_register_had_gate_for_whole_register(q_reg):
    return _qcs.qcs_quantum_register_had_gate_for_whole_register(q_reg)

def qcs_quantum_register_square_root_not_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_square_root_not_n_gate(q_reg, i)

def qcs_quantum_register_x_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_x_rot90_n_gate(q_reg, i)

def qcs_quantum_register_y_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_y_rot90_n_gate(q_reg, i)

def qcs_quantum_register_z_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_z_rot90_n_gate(q_reg, i)

def qcs_quantum_register_mx_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_mx_rot90_n_gate(q_reg, i)

def qcs_quantum_register_my_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_my_rot90_n_gate(q_reg, i)

def qcs_quantum_register_mz_rot90_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_mz_rot90_n_gate(q_reg, i)

def qcs_quantum_register_rotate_alpha_n_gate(q_reg, k, i):
    return _qcs.qcs_quantum_register_rotate_alpha_n_gate(q_reg, k, i)

def qcs_quantum_register_rotate_theta_n_gate(q_reg, theta, i):
    return _qcs.qcs_quantum_register_rotate_theta_n_gate(q_reg, theta, i)

def qcs_quantum_register_t_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_t_n_gate(q_reg, i)

def qcs_quantum_register_v_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_v_n_gate(q_reg, i)

def qcs_quantum_register_s_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_s_n_gate(q_reg, i)

def qcs_quantum_register_s_adj_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_s_adj_n_gate(q_reg, i)

def qcs_quantum_register_phase_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_phase_n_gate(q_reg, i)

def qcs_quantum_register_phase_f_n_gate(q_reg, i):
    return _qcs.qcs_quantum_register_phase_f_n_gate(q_reg, i)

def qcs_qubit_arbitrary_one_qubit_gate(q_reg, gate, i):
    return _qcs.qcs_qubit_arbitrary_one_qubit_gate(q_reg, gate, i)

def qcs_arbitrary_single_gate(q_reg, m, i):
    return _qcs.qcs_arbitrary_single_gate(q_reg, m, i)

def qcs_quantum_register_cnot(q_reg, c1, t):
    return _qcs.qcs_quantum_register_cnot(q_reg, c1, t)

def qcs_quantum_register_cnot_conj(q_reg, c1, t):
    return _qcs.qcs_quantum_register_cnot_conj(q_reg, c1, t)

def qcs_quantum_register_swap_gate(q_reg, a, b):
    return _qcs.qcs_quantum_register_swap_gate(q_reg, a, b)

def qcs_quantum_reg_fredkin_gate(q_reg, a, b, c):
    return _qcs.qcs_quantum_reg_fredkin_gate(q_reg, a, b, c)

def qcs_quantum_register_measure_one_qubit(q_reg, k):
    return _qcs.qcs_quantum_register_measure_one_qubit(q_reg, k)

def qcs_quantum_register_measure_one_qubit_in_std_base(q_reg, k):
    return _qcs.qcs_quantum_register_measure_one_qubit_in_std_base(q_reg, k)

def qcs_quantum_register_measure_one_qubit_in_std_base_force(q_reg, k, force_result):
    return _qcs.qcs_quantum_register_measure_one_qubit_in_std_base_force(q_reg, k, force_result)

def qcs_quantum_register_measure_from_to(q_reg, q_from, q_to):
    return _qcs.qcs_quantum_register_measure_from_to(q_reg, q_from, q_to)

def qcs_quantum_register_probe_one_qubit_in_std_base(q_reg, t, out_value_0, out_value_1):
    return _qcs.qcs_quantum_register_probe_one_qubit_in_std_base(q_reg, t, out_value_0, out_value_1)

def qcs_quantum_register_generate_density_matrix(q_reg):
    return _qcs.qcs_quantum_register_generate_density_matrix(q_reg)


