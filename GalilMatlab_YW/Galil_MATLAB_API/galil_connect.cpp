
/*=================================================================
 *
 * This file implements Galil communication block for MATLAB
 *
 *=================================================================*/

#include "gclibo.h"

#include "mex.h"
#include <iostream> //std::cout
#include <string> //to_string, string, etc.
#include <cstdio> //sprintf, etc.
#include <cstring> //strlen, etc.

//galill connection function
int galil_connect()
{
    
    void* g ; //Galil connection //open connection to Galil controller
	GOpen("192.168.0.43 -d", &g);

	//GCommand(g, "\x12\x16", buf, sizeof(buf), &read_bytes);


}
//entry point to Mex function
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )  
{ 
   galil_connect();
}

