
#include "prelim.h"			// definition of all constants, global variables and functions

///////////////***** M  A  I  N *****///////////////
int main ()
{
	double t = 0.;				// time variable

    FILE *resf;							// file to save error
  	resf  = fopen ("L_2 error.txt", "wt"); // open the file to record results

  	double f[3][nx+2][ny+2] = {0.}; // f[0][][] pressure field variable (P); for nx cells + 2 ghost cells
									// f[1][][] x-direction velocity (u); for nx cells + 2 ghost cells
						 			// f[2][][] y-direction velocity (v); for nx cells + 2 ghost cells
									// for a better distinction, numbering of the variables matrix is  
									// different than other 3d arrays; [3][nx+2][ny+2] instead of [nx+1][ny+1][3]

	for(int j=1; j < ny+1; j++)		// Set initial values for domain variables						
		for(int i=1; i < nx+1; i++){  

			f[0][i][j] = p_0*cos(pi*(((i-1)*dx)+(dx/2.)))*cos(pi*(((j-1)*dy)+(dy/2.)));		// initialize (P)
			f[1][i][j] = u_0*sin(pi*(((i-1)*dx)+(dx/2.)))*sin(2.*pi*(((j-1)*dy)+(dy/2.)));	// initialize (u)
			f[2][i][j] = v_0*sin(2.*pi*(((i-1)*dx)+(dx/2.)))*sin(pi*(((j-1)*dy)+(dy/2.)));	// initialize (v)

			/*if (i > 2 && i < 6)		// add local pressure for testing
				if(j > 2 && j < 6)
					f[0][i][j] += 1.;*/
		}

										/////**** IMPLICIT EULER ****/////

	while ( t < tlimit){	// initiate time loop; continue up to a pre-set time
							// the loop will also end (break;) if the convergence limit is reached)

		double fi[nx+1][ny+1][3] = {0.}; 	// flux integral

		set_gcells(f);						// update ghost cells

		update_fi (f,fi);					// update flux integral

		appx_fac_dx (f,fi);					// update first stage approx. factorization
		appx_fac_dy (f,fi);					// update second stage approx. factorization ///* seperated for better tracking*/

    	double L2 [3] = {0.};				// records L2 error for each variable at each iter 
    	update_f (f, fi, L2);				// update f with calculated flux integral (fi) and return error (L2)

		fprintf (resf, "%6e %6e %6e %6e\n", t, L2[0], L2[1], L2[2]); // print L2 square for each variable at each iter 

	//	double ave_val[3] = {0.};		//record and print average valuee os the field variables
	//	ave_value(f, ave_val);			// returns average of P, u and v over the domain
	//	fprintf (resf, "%11e, %11e, %11e, %11e\n", t, ave_val[0], ave_val[1], ave_val[2]);

		t += dt ;    					// update time step

		if (L2[2] < Climit && 			// if l2 error norms are less than convergence limit --> exit loop 
			L2[1] < Climit &&
			L2[0] < Climit) break;
	}// end of the time loop

										/////**** POST PROCESSING ****/////

	generate_bfile_3var (t,f);			// output binary data for tecplot
	//symm_values (f[1]);				// return [*] varaibles' values on vertical symmetry line
	vortex_core(f); 					// find vortexies cores and print in file 
	//generate_vtkfile(t,f,1); 			// prints only 2d map of one variable (, ,*) in vtk format; 0:P, 1:U, 2:V
										// note that this function records data in 1-X,Y if (u_top < 0)

	cout << "Convergence reached at iter = " << t/dt << endl;

	double ave_val[3];
	ave_value(f, ave_val);		// print average value of each variable in the file
	cout << "Average results: " << " P = " << ave_val[0] 
								<< " U = " << ave_val[1]	
								<< " V = " << ave_val[2] << endl;

	fclose (resf);

	return 0;
}
