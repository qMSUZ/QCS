#***************************************************************************
#*   Copyright (C) 2005 -- 2011 by Marek Sawerwain                         *
#*                                         <M.Sawerwain@gmail.com>         *
#*                                                                         *
#*   Part of the Quantum Computing Simulator:                              *
#*   http://code.google.com/p/qcs                                          *
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

#! /usr/bin/python

import qcs

_D_LEVEL = 2

_TELEPORT_QUDIT = 0
_ALICE_QUDIT = 1
_BOB_QUDIT = 2

q=qcs.QuantumReg(3, _D_LEVEL)
q.Reset()

#q.SetKet("000")

q.NotN(_TELEPORT_QUDIT)			 # make qubit for teleportation
q.HadN(_TELEPORT_QUDIT)
print("begin print entry state")
q.Pr()
print("end print entry state")

				 #
q.HadN(_ALICE_QUDIT) 		 # make EPR pair
q.CNot(_ALICE_QUDIT, _BOB_QUDIT) #
				 #

q.CNot(_TELEPORT_QUDIT, _ALICE_QUDIT)
q.HadN(_TELEPORT_QUDIT)

# print "begin before the measuremment"
# q.Pr()
# print "end before the measuremment"

v1=q.MeasureOneQudit( _TELEPORT_QUDIT )
v2=q.MeasureOneQudit( _ALICE_QUDIT )
print("v1=", v1, "v2=", v2)
print("begin after the measuremment")
q.Pr()
print("end after the measuremment")

	
mat = q.CorrMatForTeleport(v1,v2)
q.ArbitrarySingleGate(_BOB_QUDIT, mat); 
		
print("begin after correction")
q.Pr()
print("end after correction")

del q
