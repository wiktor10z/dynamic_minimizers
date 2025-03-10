/**
 *    ius: indexing uncertain strings
 *    Copyright (C) 2023 E. Gabory, C. Liu, G. Loukides, S. P. Pissis, and W. Zuba.
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <cstring>
#include <ctime>
#include <chrono>
#include <malloc.h>

#include "input.h"
#include "minimizer_index.h"
#include "krfp.h"

using namespace std;
using namespace std::chrono;

using get_time = std::chrono::steady_clock;

int main (int argc, char ** argv ) 
{
   	Settings st = decode_switches(argc, argv);
   	istream& text = st.text.is_open()?st.text:cin;

	int ell = st.ell;
	MinimizerIndex M;
	text >> M;
   	M.build_index(st.z, ell);
    return 0;
}

