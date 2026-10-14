#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# /***************************************************************************
# *   Copyright (C) 2022 -- 2026 by Marek Sawerwain                         *
# *                                  <M.Sawerwain@gmail.com>                *
# *                                  <M.Sawerwain@issi.uz.zgora.pl>         *
# *                                                                         *
# *                              by Joanna Wiśniewska                       *
# *                                  <Joanna.Wisniewska@wat.edu.pl>         *
# *                                                                         *
# *   Part of the Quantum Distance Classifier:                              *
# *         https://github.com/qMSUZ/QDCLIB                                 *
# *                                                                         *
# *   Licensed under the EUPL-1.2-or-later, see LICENSE file.               *
# *                                                                         *
# *   Licensed under the EUPL, Version 1.2 or - as soon they will be        *
# *   approved by the European Commission - subsequent versions of the      *
# *   EUPL (the "Licence");                                                 *
# *                                                                         *
# *   You may not use this work except in compliance with the Licence.      *
# *   You may obtain a copy of the Licence at:                              *
# *                                                                         *
# *   https://joinup.ec.europa.eu/software/page/eupl                        *
# *                                                                         *
# *   Unless required by applicable law or agreed to in writing,            *
# *   software distributed under the Licence is distributed on an           *
# *   "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,          *
# *   either express or implied. See the Licence for the specific           *
# *   language governing permissions and limitations under the Licence.     *
# *                                                                         *
# ***************************************************************************/

import datetime
import logging

import colorlog
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit import version as qversion
from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_bloch_multivector
from qiskit_aqt_provider import AQTProvider
from qiskit_aqt_provider.primitives import AQTSampler


def log_write_empty_line():
    for handler in logger.handlers:
        handler.stream.write("\n")
        handler.flush()


def log_close(logger_object):
    handlers = logger_object.handlers[:]
    for handler in handlers:
        handler.close()
        logger_object.removeHandler(handler)


def pr_state(state):
    print(f"0 000 {state[0]}")
    print(f"1 001 {state[1]}")
    print(f"2 010 {state[2]}")
    print(f"3 011 {state[3]}")
    print(f"4 100 {state[4]}")
    print(f"5 101 {state[5]}")
    print(f"6 110 {state[6]}")
    print(f"7 111 {state[7]}")


n_shots = 200

#
# logger creation
#

logger = logging.getLogger("Ex2-Switch-Three-Qubits-ColorLogger")
logger.setLevel(logging.DEBUG)

console_output_handler = colorlog.StreamHandler()
console_output_formatter = colorlog.ColoredFormatter(
    "%(log_color)s%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    log_colors={
        "DEBUG": "cyan",
        "INFO": "green",
        "WARNING": "yellow",
        "ERROR": "red",
        "CRITICAL": "red,bg_white",
    },
)
console_output_handler.setFormatter(console_output_formatter)

datetimeval = datetime.datetime.now().strftime("-D%Y-%m-%d-T%H-%M-%S")
filename=f"ex2-switch-three-qubits-{datetimeval}.log"
file_handler = logging.FileHandler( filename, encoding="utf-8")
file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(file_formatter)


logger.addHandler(console_output_handler)
logger.addHandler(file_handler)

#
# data load
#

data = np.load("list-of-qubits.npz")
params = data["params"]
qubits = data["qubits"]

all_probes = params.shape[0]
#all_probes = 100


logger.warning("begin ... ")
logger.info("-" * 60)
logger.info("ex2-switch-on-the-three-qubits example")
   


# log_write_empty_line()
logger.info("-" * 60)
logger.info("executing on AQT")
logger.info(f"qiskit version: {qversion.get_version_info()}")

# simulation
logger.warning("simulation begin ...")
logger.warning("simulation with probability amplitudes")
idx=0
for idx in range(all_probes):
    rqc = QuantumCircuit(3)
    rqc.h(0)
    rqc.rx(params[idx, 0], 0)
    rqc.ry(params[idx, 1], 0)
    rqc.rz(params[idx, 2], 0)
    
    rqc.x(2)
    
    rqc.cx(0, 1)
    rqc.ccx(2, 1, 0)
    rqc.cx(0, 1)
    
    state = Statevector(rqc)
    
    #pr_state(state)
    alpha0=state[4]
    alpha1=state[6]
    error_value = np.abs(alpha0 - qubits[idx][0]) + np.abs(alpha1 - qubits[idx][1])
    tmp_str=f"error on {idx} is {error_value}"
    logger.info(tmp_str)

logger.warning("simulation with probabilities")

idx=0
for idx in range(all_probes):
    rqc = QuantumCircuit(3)
    rqc.h(0)
    rqc.rx(params[idx, 0], 0)
    rqc.ry(params[idx, 1], 0)
    rqc.rz(params[idx, 2], 0)
    
    rqc.x(2)
    
    rqc.cx(0, 1)
    rqc.ccx(2, 1, 0)
    rqc.cx(0, 1)
    
    result = Statevector(rqc)
    quasi_dist = result.probabilities()[4], result.probabilities()[6]
    
    #pr_state(state)
    error_value = (quasi_dist[0] - np.abs(qubits[idx][0])**2 ) + (quasi_dist[1] - np.abs(qubits[idx][1])**2)
    tmp_str = f"error on {idx} is {error_value}"
    logger.info(tmp_str)


logger.warning("... simulation end")

# real hardware

logger.warning("real hardware begin ...")

provider = AQTProvider("ACCESS_TOKEN")
backend = provider.get_backend("offline_simulator_no_noise")

# transpiled_rqc = transpile(rqc, backend)
# print(dict(transpiled_rqc.count_ops()))
# transpiled_rqc.draw()
# transpiled_rqc.draw('latex', filename='./file.eps')

# Instantiate a sampler on the execution backend.
sampler = AQTSampler(backend)

# Optional: set the transpiler's optimization level.
# Optimization level 3 typically provides the best results.
sampler.set_transpile_options(optimization_level=3)

# Sample the circuit on the execution backend.

idx=0
key_errors=0
for idx in range(all_probes):
    rqc = QuantumCircuit(3)
    rqc.h(0)
    rqc.rx(params[idx, 0], 0)
    rqc.ry(params[idx, 1], 0)
    rqc.rz(params[idx, 2], 0)
    
    rqc.x(2)
    
    rqc.cx(0, 1)
    rqc.ccx(2, 1, 0)
    rqc.cx(0, 1)
    
    rqc.measure_all()

    result = sampler.run(rqc, shots=n_shots).result()

    quasi_dist = result.quasi_dists[0]
    
    #pr_state(state)
    if 4 in quasi_dist and 6 in quasi_dist:
        error_value = (quasi_dist[4] - np.abs(qubits[idx][0])**2 ) + (quasi_dist[6] - np.abs(qubits[idx][1])**2)
        tmp_str = f"error on {idx} is {error_value} qubits: {qubits[idx]}"
        logger.info(tmp_str)
    else:
        tmp_str = f"key error on {idx}:  {quasi_dist} {qubits[idx]} "
        logger.info(tmp_str)
        key_errors=key_errors+1
        
        
tmp_str=f"total key errors: {key_errors}"
logger.info(tmp_str)

logger.warning("... hardware end")
log_close(logger)
