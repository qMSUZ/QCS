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
import random
random.seed()

#q0 = qcs.Qubit() ; q0.SetThetaPsi(random.randint(0,180), random.randint(0,180)) #losowy stan
#q1 = qcs.Qubit() ; q1.Set0()
#q2 = qcs.Qubit() ; q2.Set0()

#print "teleported qubit:" ; print q0

#q_ar = qcs.QubitArray(3)
#q_ar.SetQubitN(0, q0)
#q_ar.SetQubitN(1, q1)
#q_ar.SetQubitN(2, q2)

_X = 0
_Y = 1
_Z = 2

q=qcs.QuantumReg(3)
q.Reset()
#q.SetKet("000");


q.NotN(0)
q.HadN(0)
#q.RotThetaN(0, 0.2)

#q.SetFromQubitArray(q_ar)
q.Pr() ; print "1st step"

q.HadN(_Y) 			 # make skew EPR pair
#q.RotThetaN(_Y, 0.2)
q.CNot(_Y, _Z)

q.Pr() ; print "d-epr making"

q.CNot(_X, _Y)
q.CNot(_Y, _Z)
q.HadN(_X)
#q.NotN(_X)

#[p0, p1] = q.ProbeQubitStdBase( _Z )
#print "before measure probe for Bob p0 = ", p0, " probe p1 =", p1
#q.PrSqr() ; print "end state print"
#q.CNot2_any(0, 1, 2, '0', '1');
#q.CNot2_any(0, 1, 2, '1', '1');
#q.PrSqr() ; print "end state after ampl corr"

q.Pr() ; print "before measurement"

v=q.MeasureN(_X, _Y)
#[p0, p1] = q.ProbeQubitStdBase( _Z )
#print "after measure probe for Bob p0 = ", p0, " probe p1 =", p1

print "v=",v

q.Pr()
	
if v==2:
	q.PauliZ(_Z)
if v==3:
	q.PauliZ(_Z)	



del q
