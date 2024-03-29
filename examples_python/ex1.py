#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
/***************************************************************************
 *   Copyright (C) 2018, 2019, 2022 by Marek Sawerwain                     *
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
 """

import qcs

print()

print("Register after Had on qubit 0")
q = qcs.QuantumRegister()
q.Reset()
q.HadN(0)
q.Pr()
print()
q.PrSqr()
print()
q.PrFull()
print()
q.PrFullSqr()

print()

print("Register after SquareRootOfNot on qubit 0")
q.Reset()
q.SquareRootOfNotN(0)
q.Pr()
q.PrAsMatlab()
print()
q.PrFull()
print()
q.PrSqr()
print()
q.PrFullSqr()

print()

print("Register after set GHZ state for four qubits")
q = qcs.QuantumRegister(4)
q.Reset()
q.SetGHZState()
q.Pr()
q.PrAsMatlab()
