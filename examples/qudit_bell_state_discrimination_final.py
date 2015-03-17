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

import qcs

_D_LEVEL = 3;

#
# |psi_{p,q}> = X^q_d(1) H(0)_d Z^p_d(0) CNOT_d(0,1) |00>
#
def make_psi(_r,p,q):
	_r.Reset()
	for i in range(0,q):
		_r.NotN(1)
	_r.HadN(0)
	for i in range(0,p):
		_r.PauliZ(0) 
	_r.CNot(0,1)	
	
def detect_bell_state( _r ):
	_r.HadN(2)
	_r.CNot(2,0)
	_r.CNot(2,1)
	_r.HadN_Conj(2) 
	#print "a1 measure"
	#_r.Pr()
	#print "a1 measure cleanup"
	#_r.Cleanup()
	#_r.Pr()
	a1=_r.MeasureOneQudit( 2 )
	#print "a1=", a1
	_r.CNot_Conj(0,3)  
	_r.CNot(1,3)
	#print "a2 measure
	#_r.Pr()
	#print "a2 measure cleanup"
	#_r.Cleanup()
	#_r.Pr()
	a2=_r.MeasureOneQudit(3)
	#print "a2=", a2
	return a1,a2

r=qcs.QuantumReg( 4, _D_LEVEL )
r.Reset()

make_psi(r,2,0)
print("d-level psi state")
r.Pr()

a1,a2=detect_bell_state(r)
print()
print("a1="+str(a1)+ " a2="+str(a2))
print()
r.Pr()


