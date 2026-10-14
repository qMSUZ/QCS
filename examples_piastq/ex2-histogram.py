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

import matplotlib.pyplot as plt
import numpy as np

# Ścieżka do pliku z danymi
#filename = "ex2-expr-1-extracted_values.txt"
filename = "ex2-expr-5-extracted_values.txt"

# Wczytanie danych
with open(filename, "r") as f:
    values = [float(line.strip()) for line in f if line.strip()]

# Konwersja do tablicy numpy
values = np.array(values)

# Utworzenie histogramu
plt.figure(figsize=(10, 6), layout='constrained')
plt.hist(values, bins=15)

# Opisy wykresu
plt.title("Histogram of error value for quantum switch", fontsize=14)
plt.xlabel("Value", fontsize=14)
plt.ylabel("Number of occurence", fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Siatka pomocnicza
plt.grid(True)

# Wyświetlenie wykresu
plt.savefig("ex2-hist-p1.eps",  bbox_inches='tight')
#plt.show()