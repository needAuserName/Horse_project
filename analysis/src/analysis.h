/*    
 *    analysis.h		
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
#ifndef ANALYSIS
#define ANALYSIS

//#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <queue>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <sstream>
#include <vector>
#include <math.h>
#include <assert.h>
#include <queue>
#endif

using namespace std;



#endif

