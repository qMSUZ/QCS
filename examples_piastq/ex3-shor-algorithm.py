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

import math

import colorlog
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit import version as qversion
from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_bloch_multivector

from qiskit.primitives import StatevectorSampler

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


def gcd( m, n ):
	a=m
	b=n
	while b!=0:
		a,b = b, math.fmod(a,b)
	return a

n_shots = 200

#
# logger creation
#

logger = logging.getLogger("Ex3-Shor-Algorithm-ColorLogger")
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
filename=f"ex3-shor-algorithm-seven-qubits-{datetimeval}.log"
file_handler = logging.FileHandler( filename, encoding="utf-8")
file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(file_formatter)


logger.addHandler(console_output_handler)
logger.addHandler(file_handler)


logger.warning("begin ... ")
logger.info("-" * 60)
logger.info("ex3-shor-algorithm-seven-qubits-example")
       


# log_write_empty_line()
logger.info("-" * 60)
logger.info("executing on AQT")
logger.info(f"qiskit version: {qversion.get_version_info()}")

# simulation
logger.warning("simulation begin ...")

n=15
x=7

r1=3
r2=4

y_classical_reg = ClassicalRegister(3, "y_classical_reg")
qreg7 = QuantumRegister(7)
rqc = QuantumCircuit(qreg7, y_classical_reg)

rqc.h(0)
rqc.h(1)
rqc.h(2)

# Power Mod
rqc.cx(2,4)
rqc.cx(2,5)
rqc.cx(3,5)
rqc.ccx(1,5,3)
rqc.cx(3,5)
rqc.cx(6,4)
rqc.ccx(1,4,6)
rqc.cx(6,4)

# IQFT
rqc.h(0)
rqc.crx(math.pi/2.0, 2, 1)
rqc.h(1)
rqc.crx(math.pi/4.0, 2, 0)
rqc.crx(math.pi/2.0, 1, 0)
rqc.h(2)

# swap
#0 2
rqc.cx(0,2)
rqc.cx(2,0)
rqc.cx(0,2)

rqc.measure([0,1,2], y_classical_reg)

sampler = StatevectorSampler()  
job = sampler.run([rqc], shots=n_shots)
result=job.result()[0]


sim_result =["", "", "", "", "", "", "", ""]
hard_result=["", "", "", "", "", "", "", ""]

for rslt in result.data.y_classical_reg.get_counts():
    decimal = int(rslt, 2)
    prob=result.data.y_classical_reg.get_counts()[rslt] / n_shots
    y = decimal    
    i=gcd(y, (2 ** r1))
    r=(2 ** r1) / i

    p1=gcd((x ** (r/2)) - 1, n)
    p2=gcd((x ** (r/2)) + 1, n)

    if int(p1*p2) == int(n):
        tmp_str=f"{rslt} bin {decimal:>3}(dec) prob={prob} {p1} * {p2}={p1*p2}"
        sim_result[decimal]=prob
        logger.info(tmp_str)
        sim_result[decimal]=tmp_str

    else:
        tmp_str=f"{rslt} bin {decimal:>3}(dec) prob={prob}"
        logger.info(tmp_str)
        sim_result[decimal]=tmp_str


logger.warning("simulation end")

# real hardware

logger.warning("real hardware begin ...")

provider = AQTProvider("ACCESS_TOKEN")
backend = provider.get_backend("offline_simulator_no_noise")

# transpiled_rqc = transpile(rqc, backend)
# transpiled_rqc.draw()

# Instantiate a sampler on the execution backend.
sampler = AQTSampler(backend)

# Optional: set the transpiler's optimization level.
# Optimization level 3 typically provides the best results.
sampler.set_transpile_options(optimization_level=3)

# Sample the circuit on the execution backend.

n=15
x=7

r1=3
r2=4

y_classical_reg = ClassicalRegister(3, "y_classical_reg")
qreg7 = QuantumRegister(7)
rqc = QuantumCircuit(qreg7, y_classical_reg)

rqc.h(0)
rqc.h(1)
rqc.h(2)

# Power Mod
rqc.cx(2,4)
rqc.cx(2,5)
rqc.cx(3,5)
rqc.ccx(1,5,3)
rqc.cx(3,5)
rqc.cx(6,4)
rqc.ccx(1,4,6)
rqc.cx(6,4)

# IQFT
rqc.h(0)
rqc.crx(math.pi/2.0, 2, 1)
rqc.h(1)
rqc.crx(math.pi/4.0, 2, 0)
rqc.crx(math.pi/2.0, 1, 0)
rqc.h(2)

# swap
#0 2
rqc.cx(0,2)
rqc.cx(2,0)
rqc.cx(0,2)

rqc.measure([0,1,2], y_classical_reg)

 
job = sampler.run([rqc], shots=n_shots)
result=job.result()

for rslt in result.quasi_dists[0]:
    decimal = int(rslt)
    binary=bin(decimal)[2:].zfill(3)
    prob=result.quasi_dists[0][rslt]
    y = decimal    
    i=gcd(y, (2 ** r1))
    r=(2 ** r1) / i

    p1=gcd((x ** (r/2)) - 1, n)
    p2=gcd((x ** (r/2)) + 1, n)

    #print ("r=", r)
    #print ("p1=", p1)
    #print ("p2=", p2)
    #print ("p1*p2=", p1*p2)
    if int(p1*p2) == int(n):
        tmp_str=f"{binary:>3} bin {decimal:>3}(dec) prob={prob} {p1} * {p2}={p1*p2}"
        logger.info(tmp_str)
        hard_result[decimal]=tmp_str

    else:
        tmp_str=f"{binary:>3} bin {decimal:>3}(dec) prob={prob}"
        logger.info(tmp_str)
        hard_result[decimal]=tmp_str

logger.warning("... hardware end")


tmp_str=f"                         *** Final results ***"
logger.info(tmp_str)
tmp_str="Simulation" + " " * 35 + "| Hardware Piast-Q"
logger.info(tmp_str)
for i in range(8):
    tmp_str = sim_result[i].ljust(44) + " | " + hard_result[i].ljust(44)
    logger.info(tmp_str)

log_close(logger)
