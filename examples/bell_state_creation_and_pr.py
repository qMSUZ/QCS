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

def make_psi_plus(_r):
    _r.SetKet("00")
    _r.HadN(0)
    _r.CNot(0,1)
    
def make_psi_minus(_r):
    _r.SetKet("10")
    _r.HadN(0)
    _r.CNot(0,1)

def make_phi_plus(_r):
    _r.SetKet("10")
    _r.HadN(1)
    _r.CNot(1,0)
    
def make_phi_minus(_r):
    _r.SetKet("11")
    _r.HadN(0)
    _r.CNot(0,1)

r=qcs.QuantumReg(2)
r.Reset()

make_psi_plus(r)
print "psi+"
r.Pr()

make_psi_minus(r)
print "psi-"
r.Pr()

make_phi_plus(r)
print "phi+"
r.Pr()

make_phi_minus(r)
print "phi-"
r.Pr()


