/***************************************************************************
 *   Copyright (C) 2021 by Marek Sawerwain                                 *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                         <M.Sawerwain@issi.uz.zgora.pl   *
 *                                                                         *
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
 
#include <stdio.h>

#include "qcs.h"

void qcs_rand_test1()
{
    int i;

    printf("Random number test:\n");

    for(i=0;i<10;i++)
        printf("[%02d] %f \n", i+1, qcs_randu());

}



int main( int argc, char *argv[] )
{
	qcs_core_library_initialization ( );

	printf(": qcs_rand_test\n");
	qcs_rand_test1();

	qcs_core_library_deinitialization();

	return 0;
}