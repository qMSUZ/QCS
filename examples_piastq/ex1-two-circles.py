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

import logging

import colorlog
import matplotlib.pyplot as plt
import numpy as np

from pcss_qapi import AuthorizationService
from pcss_qapi.aqt.provider import PCSS_AQTProvider

from qiskit import QuantumCircuit, transpile
from qiskit import version as qversion
from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_bloch_multivector
from qiskit_aqt_provider import AQTProvider
from qiskit_aqt_provider.primitives import AQTSampler

IN_RANGE = 0x1000
IN_RANGE_MARGIN = 0x1001

p0_0_min = 1
p0_0_max = 0
p1_0_min = 1
p1_0_max = 0

p0_0_avg = 0
p1_0_avg = 0

p0_1_min = 1
p0_1_max = 0
p1_1_min = 1
p1_1_max = 0

p0_1_avg = 0
p1_1_avg = 0

margin = 0.05

n_shots = 200


def trivial_classification(p0, p1):

    if (p0 >= p0_0_min) and (p0 <= p0_0_max) and (p1 >= p1_0_min) and (p1 <= p1_0_max):
        return (0, IN_RANGE)

    if (p0 >= p0_1_min) and (p0 <= p0_1_max) and (p1 >= p1_1_min) and (p1 <= p1_1_max):
        return (1, IN_RANGE)

    if (
        (p0 >= p0_0_min - margin)
        and (p0 <= p0_0_max + margin)
        and (p1 >= p1_0_min - margin)
        and (p1 <= p1_0_max + margin)
    ):
        return (0, IN_RANGE_MARGIN)

    if (
        (p0 >= p0_1_min - margin)
        and (p0 <= p0_1_max + margin)
        and (p1 >= p1_1_min - margin)
        and (p1 <= p1_1_max + margin)
    ):
        return (1, IN_RANGE_MARGIN)

    return None


def log_write_empty_line():
    for handler in logger.handlers:
        handler.stream.write("\n")
        handler.flush()


def log_close(logger_object):
    handlers = logger_object.handlers[:]
    for handler in handlers:
        handler.close()
        logger_object.removeHandler(handler)


#
# logger creation
#

logger = logging.getLogger("Ex1-TwoCircuits-Single-Qubit-ColorLogger")
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

file_handler = logging.FileHandler("ex1-two-circles.log", encoding="utf-8")
file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(file_formatter)


logger.addHandler(console_output_handler)
logger.addHandler(file_handler)

#
# data load
#


data = np.load("two-circles.npz")
d = data["d"]
l = data["l"]


logger.warning("begin ... ")
logger.info("-" * 60)
logger.info("ex1-two-circles example")
logger.info(f"data set loaded {l.shape[0]} probes")


testdata = d

d0 = d[l == 0]
d1 = d[l == 1]


plot_drawing = False

#
# circles drawing
#
if plot_drawing == True:
    plt.figure()
    plt.plot(d0[:, 0], d0[:, 1], ".b")
    plt.title("First circle")
    plt.axis("equal")

    plt.figure()
    plt.plot(d1[:, 0], d1[:, 1], ".g")
    plt.title("Second circle")
    plt.axis("equal")

    plt.figure()
    plt.plot(d0[:, 0], d0[:, 1], ".b")
    plt.plot(d1[:, 0], d1[:, 1], ".g")
    plt.title("Two circles")
    plt.axis("equal")

    # draw two circles
    plt.figure()
    plt.plot(d[:, 0], d[:, 1], ".r")
    plt.title("Two circles")
    plt.axis("equal")

    # draw two circles
    plt.figure()
    plt.plot(d0[:, 0], d0[:, 1], ".b", label="Circle A")
    plt.plot(d1[:, 0], d1[:, 1], "+r", label="Circle B")
    plt.title("Two circles")
    plt.legend()
    plt.axis("equal")
    plt.savefig("two-circles-dataset.eps")
    plt.show()
    
    # draw two circles
    plt.figure()
    plt.plot(testdata[0:399, 0], testdata[0:399, 1], ".", color="lightgray", label="Circle A")
    plt.plot(testdata[0:399, 0], testdata[0:399, 1], "+", color="lightgray", label="Circle B")  
    plt.plot(testdata[402, 0], testdata[402, 1], "o", color="r", label="")
    plt.plot(testdata[414, 0], testdata[414, 1], "o", color="r", label="")
    plt.plot(testdata[416, 0], testdata[416, 1], "o", color="r", label="")
    plt.plot(testdata[437, 0], testdata[436, 1], "o", color="r", label="")
    plt.plot(testdata[467, 0], testdata[467, 1], "o", color="r", label="")
    plt.plot(testdata[482, 0], testdata[482, 1], "o", color="r", label="")
    plt.plot(testdata[487, 0], testdata[487, 1], "o", color="r", label="")
    plt.plot(testdata[492, 0], testdata[492, 1], "o", color="r", label="")
    plt.title("Two circles")
    plt.legend()
    plt.axis("equal")
    plt.savefig("two-circles-dataset-class.eps")
    plt.show()    


#
# calculate probabilites
# on the learning set
#
#
for idx in range(400):
    qc = QuantumCircuit(1)
    qc.h(0)
    probe_arg = d[idx, 0] ** 2 + d[idx, 1] ** 2
    qc.ry(probe_arg, 0)
    state = Statevector(qc)
    # plot_bloch_multivector(state)

    if l[idx] == 0:
        p0_0_min = min(p0_0_min, state.probabilities()[0])
        p0_0_max = max(p0_0_max, state.probabilities()[0])
        p1_0_min = min(p1_0_min, state.probabilities()[1])
        p1_0_max = max(p1_0_max, state.probabilities()[1])

    if l[idx] == 1:
        p0_1_min = min(p0_1_min, state.probabilities()[0])
        p0_1_max = max(p0_1_max, state.probabilities()[0])
        p1_1_min = min(p1_1_min, state.probabilities()[1])
        p1_1_max = max(p1_1_max, state.probabilities()[1])

p0_0_avg = (p0_0_min + p0_0_max) / 2
p1_0_avg = (p1_0_min + p1_0_max) / 2

p0_1_avg = (p0_1_min + p0_1_max) / 2
p1_1_avg = (p1_1_min + p1_1_max) / 2


# p0_0_min p0_0_max p1_0_min p1_0_max p0_0_avg p1_0_avg

# log_write_empty_line()
logger.info("-" * 60)
logger.info("probabilites range for class zero and one")
tmp_str = (
    f"0: <{p0_0_min:.4f} {p0_0_max:.4f}>  1: <{p1_0_min:.4f} {p1_0_max:.4f}>",
    "for label 0",
)
logger.info(tmp_str)
tmp_str = (
    f"0: <{p0_1_min:.4f} {p0_1_max:.4f}>  1: <{p1_1_min:.4f} {p1_1_max:.4f}>",
    "for label 1",
)
logger.info(tmp_str)

# log_write_empty_line()
logger.info("-" * 60)
logger.info("executing on AQT")
logger.info(f"qiskit version: {qversion.get_version_info()}")

# simulation
logger.warning("simulation begin ...")

success_class = 0
failure_class = 0
no_in_range = 0
no_in_range_margin = 0

for idx in range(400, 500):
    rqc = QuantumCircuit(1)
    rqc.h(0)
    probe_arg = testdata[idx, 0] ** 2 + testdata[idx, 1] ** 2
    rqc.ry(probe_arg, 0)

    result = Statevector(rqc)

    quasi_dist = result.probabilities()[0], result.probabilities()[1]

    pl = trivial_classification(quasi_dist[0], quasi_dist[1])
    tmp_str = f"qd: {quasi_dist} probe: {idx} label {l[idx]} propose label {pl}"

    if pl != None:
        if pl[0] == l[idx]:
            success_class = success_class + 1
            if pl[1] == IN_RANGE:
                no_in_range = no_in_range + 1
            if pl[1] == IN_RANGE_MARGIN:
                no_in_range_margin = no_in_range_margin + 1
            logger.info(tmp_str)

        else:
            failure_class = failure_class + 1
            logger.error(tmp_str)
    else:
        failure_class = failure_class + 1
        logger.critical(tmp_str)


tmp_str = f"success: {success_class} failure {failure_class}"
logger.info(tmp_str)

tmp_str = f"success in range: {no_in_range} in range margin {no_in_range_margin}"
logger.info(tmp_str)


logger.warning("... simulation end")

# real hardware

logger.warning("real hardware begin ...")

provider = AQTProvider("ACCESS_TOKEN")
backend = provider.get_backend("offline_simulator_no_noise")

# PIAST-Q hardware access with token
#AuthorizationService.login("000-0-0-0-0-0-0-0-0")
#provider = PCSS_AQTProvider()
#backend = provider.get_direct_access_backend()


# transpiled_rqc = transpile(rqc, backend)
# transpiled_rqc.draw()

# Instantiate a sampler on the execution backend.
sampler = AQTSampler(backend)

# Optional: set the transpiler's optimization level.
# Optimization level 3 typically provides the best results.
sampler.set_transpile_options(optimization_level=3)

# Sample the circuit on the execution backend.

success_class = 0
failure_class = 0
no_in_range = 0
no_in_range_margin = 0

for idx in range(400, 500):
    rqc = QuantumCircuit(1)
    rqc.h(0)
    probe_arg = testdata[idx, 0] ** 2 + testdata[idx, 1] ** 2
    rqc.ry(probe_arg, 0)
    rqc.measure_all()

    result = sampler.run(rqc, shots=n_shots).result()

    quasi_dist = result.quasi_dists[0]

    pl = trivial_classification(quasi_dist[0], quasi_dist[1])
    tmp_str = f"qd: {quasi_dist} probe: {idx} label {l[idx]} propose label {pl}"

    if pl != None:
        if pl[0] == l[idx]:
            success_class = success_class + 1
            if pl[1] == IN_RANGE:
                no_in_range = no_in_range + 1
            if pl[1] == IN_RANGE_MARGIN:
                no_in_range_margin = no_in_range_margin + 1
            logger.info(tmp_str)

        else:
            failure_class = failure_class + 1
            logger.error(tmp_str)
    else:
        failure_class = failure_class + 1
        logger.critical(tmp_str)


tmp_str = f"success: {success_class} failure {failure_class}"
logger.info(tmp_str)

tmp_str = f"success in range: {no_in_range} in range margin {no_in_range_margin}"
logger.info(tmp_str)

logger.warning("... hardware end")
log_close(logger)
