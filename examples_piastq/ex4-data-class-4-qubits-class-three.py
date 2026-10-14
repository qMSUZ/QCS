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

import time
import datetime
import logging

import math

import colorlog
import matplotlib.pyplot as plt
import numpy as np

from pcss_qapi import AuthorizationService
from pcss_qapi.aqt.provider import PCSS_AQTProvider

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit import version as qversion
from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_bloch_multivector

from qiskit.primitives import StatevectorSampler

from qiskit_algorithms.optimizers import ADAM, COBYLA, SPSA

from qiskit_aqt_provider import AQTProvider
from qiskit_aqt_provider.primitives import AQTSampler

class bcol:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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

def parity(bitstring):
    """Returns 1 if parity of `bitstring` is even, otherwise 0."""
    hamming_weight = sum(int(k) for k in list(bitstring))
    return (hamming_weight+1) % 2


def sim_circuit_execute(probe, params):
    theta = params
    
    qreg4 = QuantumRegister(4)
    y_classical_reg = ClassicalRegister(4, "y_classical_reg")
    rqc = QuantumCircuit(qreg4, y_classical_reg)

    # input state
    
    rqc.ry(probe[0], 0)
    rqc.ry(probe[1], 1)
    rqc.ry(probe[2], 2)
    rqc.ry(probe[3], 3)

    rqc.h(0)
    rqc.h(1)
    rqc.h(2)
    rqc.h(3)

    # parametric circuit

    rqc.ry(theta[0], 0)
    rqc.ry(theta[1], 1)
    rqc.ry(theta[2], 2)
    rqc.ry(theta[3], 3)
    rqc.cx(0, 1)
    rqc.cx(2, 3)
    rqc.ry(theta[4], 0)
    rqc.ry(theta[5], 1)
    rqc.ry(theta[6], 2)
    rqc.ry(theta[7], 3)
    rqc.cx(0, 1)
    rqc.cx(2, 3)
    rqc.ry(theta[8], 1)
    rqc.ry(theta[9], 3)
    rqc.cx(1, 3)
    rqc.ry(theta[10], 2)
    rqc.ry(theta[11], 3)
    rqc.cx(2, 3)

    rqc.barrier()

    rqc.measure([0,1,2,3], y_classical_reg)

    #rqc.draw()

    sampler = StatevectorSampler()  
    job = sampler.run([rqc], shots=n_shots)
    result=job.result()[0]

    output = np.zeros(16)

    for rslt in result.data.y_classical_reg.get_counts():
        decimal = int(rslt, 2)
        prob=result.data.y_classical_reg.get_counts()[rslt] / n_shots
        output[decimal] = prob

    return output

def hardware_circuit_execute(sampler, probe, params):
    theta = params
    
    qreg4 = QuantumRegister(4)
    y_classical_reg = ClassicalRegister(4, "y_classical_reg")
    rqc = QuantumCircuit(qreg4, y_classical_reg)

    # input state
    
    rqc.ry(probe[0], 0)
    rqc.ry(probe[1], 1)
    rqc.ry(probe[2], 2)
    rqc.ry(probe[3], 3)

    rqc.h(0)
    rqc.h(1)
    rqc.h(2)
    rqc.h(3)

    # parametric circuit

    rqc.ry(theta[0], 0)
    rqc.ry(theta[1], 1)
    rqc.ry(theta[2], 2)
    rqc.ry(theta[3], 3)
    rqc.cx(0, 1)
    rqc.cx(2, 3)
    rqc.ry(theta[4], 0)
    rqc.ry(theta[5], 1)
    rqc.ry(theta[6], 2)
    rqc.ry(theta[7], 3)
    rqc.cx(0, 1)
    rqc.cx(2, 3)
    rqc.ry(theta[8], 1)
    rqc.ry(theta[9], 3)
    rqc.cx(1, 3)
    rqc.ry(theta[10], 2)
    rqc.ry(theta[11], 3)
    rqc.cx(2, 3)

    rqc.barrier()

    rqc.measure([0,1,2,3], y_classical_reg)

    #rqc.draw()

    job = sampler.run([rqc], shots=n_shots)
    result=job.result()

    output = np.zeros(16)

    #for rslt in result.data.y_classical_reg.get_counts():
    for rslt in result.quasi_dists[0]:
        decimal = int(rslt)
        #prob=result.data.y_classical_reg.get_counts()[rslt] / n_shots
        prob=result.quasi_dists[0][rslt]
        output[decimal] = prob

    return output
    
def objective_function_class_one(params):
    """Cost function of circuit parameters on training data.
    The optimizer will attempt to minimize this."""
    global n_idx
    
    # for class one
    cost=0
    for i in range(50):
        output = sim_circuit_execute(data[i], params)
        cost = cost + (1.0 - np.dot(one_class, output)) 
    
    #print(f"objective_function[step={n_idx}] cost: {cost}", end='\r' )
    print(f"objective_function[step={n_idx}] cost: {cost}" )
   
    costs.append( cost )
    evaluations.append( n_idx )
    n_idx = n_idx + 1
    
    return cost

def objective_function_class_two(params):
    """Cost function of circuit parameters on training data.
    The optimizer will attempt to minimize this."""
    global n_idx
    
    # for class one
    cost=0
    for i in range(50,100):
        output = sim_circuit_execute(data[i], params)
        cost = cost + (1.0 - np.dot(two_class, output)) 
    
    #print(f"objective_function[step={n_idx}] cost: {cost}", end='\r' )
    print(f"objective_function[step={n_idx}] cost: {cost}" )
   
    costs.append( cost )
    evaluations.append( n_idx )
    n_idx = n_idx + 1
    
    return cost


def objective_function_class_three(params):
    """Cost function of circuit parameters on training data.
    The optimizer will attempt to minimize this."""
    global n_idx
    
    # for class one
    cost=0
    for i in range(100, 150):
        output = sim_circuit_execute(data[i], params)
        cost = cost + (1.0 - np.dot(three_class, output)) 
    
    #print(f"objective_function[step={n_idx}] cost: {cost}", end='\r' )
    print(f"objective_function[step={n_idx}] cost: {cost}" )
   
    costs.append( cost )
    evaluations.append( n_idx )
    n_idx = n_idx + 1
    
    return cost

#
# logger creation
#

logger = logging.getLogger("Ex4-data-class-4-qubits")
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
filename=f"ex4-data-class-4-qubits-{datetimeval}.log"
file_handler = logging.FileHandler( filename, encoding="utf-8")
file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(file_formatter)


logger.addHandler(console_output_handler)
logger.addHandler(file_handler)


logger.info("begin ... ")
logger.info("-" * 60)
logger.info("ex4-data-class-4-qubits-example")
       


# log_write_empty_line()
logger.info("-" * 60)
logger.info("executing on AQT")
logger.info(f"qiskit version: {qversion.get_version_info()}")

# simulation

n_idx = 0
n_shots = 200
num_of_angles = 12


logger.info("begin data loading ...")

data_iris = np.load("iris.npz")

data = data_iris["data"][:,0:4]
labels = data_iris["data"][:,4]

data1 = data[data_iris["data"][:,4]==1]
data2 = data[data_iris["data"][:,4]==2]
data3 = data[data_iris["data"][:,4]==3]

logger.info("... end of data loading")

probe_idx = 0 

costs = [ ]
evaluations = [ ]

one_class   = np.zeros(16)
two_class   = np.zeros(16)
three_class = np.zeros(16)

base_idx=0
for k in range(5):
    one_class[base_idx + 0] = 1
    two_class[base_idx + 1] = 1
    three_class[base_idx + 2] = 1
    base_idx = base_idx + 3
    base_idx = base_idx % 15

one_class=one_class/np.linalg.norm(one_class)
two_class=two_class/np.linalg.norm(two_class)
three_class=three_class/np.linalg.norm(three_class)

#  0 0000   1
#  1 0001   2 
#  2 0010   3
#  3 0011   1
#  4 0100   2 
#  5 0101   3
#  6 0110   1
#  7 0111   2
#  8 1000   3
#  9 1001   1
# 10 1010   2
# 11 1011   3
# 12 1100   1 
# 13 1101   2
# 14 1110   3
# 15 1111   1

#theta = [ 0 ] * num_parameters
theta = np.random.rand( num_of_angles )
initial_point = np.zeros( num_of_angles )

objective_values = []
max_iters = 1024


#output = sim_circuit_execute(data[0], initial_point)
#total_prob=output[:].sum()
#print(output)
#print("total_prob", total_prob)

logger.info("training begin ...")

costs = [ ]
evaluations = [ ]
n_idx = 0
optimizer = COBYLA( maxiter = max_iters )
start = time.time()
result = optimizer.minimize(objective_function_class_three, initial_point)
elapsed = time.time() - start
params_for_class_T = result.x
tmp_str=f"class one training time = {elapsed} sec."
logger.info(tmp_str)

# costs = [ ]
# evaluations = [ ]
# n_idx = 0
# optimizer = COBYLA( maxiter = max_iters )
# start = time.time()
# result = optimizer.minimize(objective_function_class_two, initial_point)
# elapsed = time.time() - start
# params_for_class_two = result.x
# tmp_str=f"class two training time = {elapsed} sec."
# logger.warning(tmp_str)

# costs = [ ]
# evaluations = [ ]
# n_idx = 0
# optimizer = COBYLA( maxiter = max_iters )
# start = time.time()
# result = optimizer.minimize(objective_function_class_three, initial_point)
# elapsed = time.time() - start
# params_for_class_three = result.x
# tmp_str=f"class three training time = {elapsed} sec."
# logger.warning(tmp_str)

logger.info("... training end")

good_labels_one=0
good_labels_two=0
good_labels_three=0

good_labels=0
bad_labels=0

logger.info("simulation probes classification...")
for idx in range(150):
    output = sim_circuit_execute(data[idx], params_for_class_T)
    #output = sim_circuit_execute(data[idx], params_for_class_two)
    # output = sim_circuit_execute(data[idx], params_for_class_three)
    f_c1=np.dot(output, one_class)
    f_c2=np.dot(output, two_class)
    f_c3=np.dot(output, three_class)
    label=labels[idx]
    propose_label=-1

    if f_c1>0.29 and f_c2<0.2 and f_c3<0.2:
        propose_label=1
    if f_c1<0.2 and f_c2>0.29 and f_c3<0.2:
        propose_label=2
    if f_c1<0.2 and f_c2<0.2 and f_c3>0.29:
        propose_label=3

    if int(label) == propose_label and int(label) == 1:
        good_labels_one=good_labels_one+1
    if int(label) == propose_label and int(label) == 2:
        good_labels_two=good_labels_two+1
    if int(label) == propose_label and int(label) == 3:
        good_labels_three=good_labels_three+1

    if int(label) == propose_label:
        good_labels=good_labels+1
        #tmp_str=bcol.OKGREEN  + f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label} OK" + bcol.ENDC
        tmp_str=f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label} OK"
        logger.warning(tmp_str)
    else:
        bad_labels=bad_labels + 1
        #tmp_str=bcol.FAIL  + f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label}" + bcol.ENDC
        tmp_str=f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label}"
        logger.error(tmp_str)


tmp_str = f"good_labels={good_labels:3}, bad_labels={bad_labels:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_one={good_labels_one:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_two={good_labels_two:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_three={good_labels_three:3}"
logger.info(tmp_str)

logger.info("... simulation probes classification end")



logger.info("begin real hardware probes classification...")

provider = AQTProvider("ACCESS_TOKEN")
backend = provider.get_backend("offline_simulator_no_noise")

#AuthorizationService.login("TOKEN")
#provider = PCSS_AQTProvider()
#backend = provider.get_direct_access_backend()


# transpiled_rqc = transpile(rqc, backend)
# transpiled_rqc.draw()

# Instantiate a sampler on the execution backend.
sampler = AQTSampler(backend)

# Optional: set the transpiler's optimization level.
# Optimization level 3 typically provides the best results.
sampler.set_transpile_options(optimization_level=3)

good_labels_one=0
good_labels_two=0
good_labels_three=0

good_labels=0
bad_labels=0


for idx in range(150):
    output = hardware_circuit_execute(sampler, data[idx], params_for_class_T)
    # output = sim_circuit_execute(data[idx], params_for_class_two)
    # output = sim_circuit_execute(data[idx], params_for_class_three)
    f_c1=np.dot(output, one_class)
    f_c2=np.dot(output, two_class)
    f_c3=np.dot(output, three_class)
    label=labels[idx]
    propose_label=-1

    if f_c1>0.29 and f_c2<0.2 and f_c3<0.2:
        propose_label=1
    if f_c1<0.2 and f_c2>0.29 and f_c3<0.2:
        propose_label=2
    if f_c1<0.2 and f_c2<0.2 and f_c3>0.29:
        propose_label=3

    if int(label) == propose_label and int(label) == 1:
        good_labels_one=good_labels_one+1
    if int(label) == propose_label and int(label) == 2:
        good_labels_two=good_labels_two+1
    if int(label) == propose_label and int(label) == 3:
        good_labels_three=good_labels_three+1

    if int(label) == propose_label:
        good_labels=good_labels+1
        #tmp_str=bcol.OKGREEN  + f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label} OK" + bcol.ENDC
        tmp_str=f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label} OK"
        logger.warning(tmp_str)
    else:
        bad_labels=bad_labels + 1
        #tmp_str=bcol.FAIL  + f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label}" + bcol.ENDC
        tmp_str=f"probe label={label}: dots: {f_c1:5.3f} {f_c2:5.3f} {f_c3:5.3f} propose_label={propose_label}"
        logger.error(tmp_str)


tmp_str = f"good_labels={good_labels:3}, bad_labels={bad_labels:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_one={good_labels_one:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_two={good_labels_two:3}"
logger.info(tmp_str)
tmp_str = f"good_labels_three={good_labels_three:3}"
logger.info(tmp_str)

logger.info("end real hardware probes classification...")
