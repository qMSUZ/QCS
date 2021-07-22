/***************************************************************************
 *   Copyright (C) 2018,2019 by Marek Sawerwain                            *
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

#include "qcs.h"
#include "qcs_info.h"
#include "qcs_qubit_gates.h"

#ifdef __qcs_core_library_debug_mode__
#define MEMWATCH
#define MEMWATCH_STDIO
//#include "memwatch/memwatch.h"
#endif

#ifdef __qcs_core_library_debug_mode__
static char* _Version="Quantum Computing System v"QCS_VERSION" __debug_mode__\nRelease name: "QCS_RELEASE_NAME;
#else
static char* _Version="Quantum Computing System v"QCS_VERSION"\nRelease name: "QCS_RELEASE_NAME;
#endif

int _amplitude_disp_type,
    no_state_vector_max_qubits,
    use_state_vector_max_qubits,
    use_graph_state_desc_max_qubits,
    pqc_mode_max_qubits,
    chp_mode_max_qubits;

static const char* _CompileSystem=": compilation date (" __DATE__ " "  __TIME__")";

DYNAMIC_LIB_DECORATION char* version()
{
	return (char*)_Version;
}

DYNAMIC_LIB_DECORATION char* compile_system()
{
  return (char*)_CompileSystem;
}

DYNAMIC_LIB_DECORATION char* compilator_name()
{
       static char _CompilatorName[512];
       _CompilatorName[0]=0;
       
#if defined(__GNUC__) && defined(__MINGW32_VERSION) && defined(_WIN32)
       sprintf(&_CompilatorName[0], ": target WIN32, compiled by GCC %s from MinGW v%d.%d" ,
       __VERSION__, __MINGW32_MAJOR_VERSION, __MINGW32_MINOR_VERSION);
#endif

#if defined(__GNUC__) && defined(__MINGW64_VERSION_STR) && defined(_WIN32)
       sprintf(&_CompilatorName[0], ": target WIN32, compiled by GCC %s from MinGW64 v%s" ,
       __VERSION__, __MINGW64_VERSION_STR);
#endif

#if defined(__GNUC__) && defined(__linux__)
       sprintf(&_CompilatorName[0], ": target Linux, compiled by GCC %s" , __VERSION__);
#endif

  return _CompilatorName;
}

DYNAMIC_LIB_DECORATION void amplitude_display_format(int _amplitude_disp_type)
{
  printf("amplitude_display_format: please, write me!!!!\n");
}

DYNAMIC_LIB_DECORATION void info()
{
    printf("no state vector max qubits: %d\n",no_state_vector_max_qubits);
    printf("use state vector max qubits: %d\n",use_state_vector_max_qubits);
    printf("use graph state desc max qubits: %d\n",use_graph_state_desc_max_qubits);
    printf("pqc mode max qubits: %d\n",pqc_mode_max_qubits);
    printf("chp mode max qubits: %d\n",chp_mode_max_qubits);
}

DYNAMIC_LIB_DECORATION void qcs_core_library_initialization()
{
    srand(time(NULL));
    rand();
    qcs_mt_init_by_int((unsigned)time(NULL));

    RandomInitialise(1802, 9373);

    qcs_build_qubit_gates_matrix();

    no_state_vector_max_qubits=1000;
    use_state_vector_max_qubits=32;
    use_graph_state_desc_max_qubits=50;
    pqc_mode_max_qubits=10000;
    chp_mode_max_qubits=10000;
}

DYNAMIC_LIB_DECORATION void qcs_core_library_deinitialization()
{
    qcs_destroy_qubit_gates_matrix();

    no_state_vector_max_qubits=0;
    use_state_vector_max_qubits=0;
    use_graph_state_desc_max_qubits=0;
    pqc_mode_max_qubits=0;
    chp_mode_max_qubits=0;
}
