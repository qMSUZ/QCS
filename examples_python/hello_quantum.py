#! /usr/bin/env python
# -*- coding: utf-8 -*-

#***************************************************************************
#*   Copyright (C) 2005 -- 2022 by Marek Sawerwain                         *
#*                                         <M.Sawerwain@gmail.com>         *
#*                                         <M.Sawerwain@issi.uz.zgora.pl   *
#*                                                                         *
#*   Part of the Quantum Computing Simulator:                              *
#*    https://github.com/qMSUZ/QCS                                         *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 3 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************/

import qcs

print("Hello Quantum!")
print("Creation EPR pair aka Bell Pair")

q = qcs.QuantumRegister( 2 ) 
q.Reset()

q.HadN(0)
q.CNot(0, 1)

print("Quantum EPR state before measurement")
q.Pr()
