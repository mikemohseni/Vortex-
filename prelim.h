#ifndef _prelim_h_included_
#define _prelim_h_included_

   #include <fstream>
   #include <iostream>
   #include <stdio.h>
   #include <time.h>
   #include <cmath>
   #include <vector>
   #include <cstdlib>

   using namespace std;

   const double pi = 3.1415926535897;

   const int nx = 20;					// mesh size along x axis (j iter)
   const int ny = 60;					// mesh size along y axis (i iter)
   										//#### note that thomas solver assumes MAX(nx) = 200
   										//#### chnage if neccessary

   const double re   = 125.;			// Reynolds number
   const double p_0  = 1.0;				// pressure constant
   const double u_0  = 1.0;				// horizontal velocity constant
   const double v_0  = 1.0;				// vertical velocity constant
   const double beta = 1.0;				// artificial compressibility canstant 
   const double omega= 1.0;				// over-relaxation factor; set to 1.0 for no over-relaxation
   const double a_fac= -1.0;			// A-factor in filter term added to pressure equation rhs 
   										// filter off if a_fac = 0.;

   const double u_top = 1.; 			// velocity at top of the domain

   const double x_max = 1.0;			// length of the domain, x-direction
   const double y_max = 3.0;			// width of the domain, y-direction
   const double dx = x_max/nx;			// mesh size x-direction
   const double dy = y_max/ny;			// mesh size y-direction

   const double dt = 0.1;				// time step

   const double Climit = 0.000000000000001; 	// convergence limit
   const double tlimit = 2000.; 		// iteration time limit

   const double rcore = 0.3 ;			// vortices in this radius belong to a single core 

   void update_fi (const double aux[3][nx+2][ny+2], double aux1[nx+1][ny+1][3]); // update flux integral	 
   void set_gcells(double aux[3][nx+2][ny+2]);									// update ghost cells
   void SolveThomas(double LHS[][3], double RHS[],const int iSize);				// solve tri-diagonal matrix 
   void appx_fac_dx(const double aux_f[3][nx+2][ny+2], 
   					double aux[nx+1][ny+1][3]); 				// Approximate factorization - through rows (first step)
   void appx_fac_dy(const double aux_f[3][nx+2][ny+2],
   					double aux[nx+1][ny+1][3]); 				// Approximate factorization - through columns (second step)
   void update_f(double aux_f[3][nx+2][ny+2], 		
   				 double aux_fi[nx+1][ny+1][3], double err[3]); 	// update domain variables in time & return error 
   void Jacob_x(int i, int j, double aux1[3][3], 				// retun Jacobean matrix, x discretization
   				double aux2[3][3], double aux3[3][3],
				const double auxf[3][nx+2][ny+2]);		
   void Jacob_y(int i, int j, double aux1[3][3], 				// retun Jacobean matrix, y discretization
   				double aux2[3][3], double aux3[3][3],
				const double auxf[3][nx+2][ny+2]);
   double p (int q, int w);										// update u_x or velocity in x direction   
   double u (int q, int w);										// update u_x or velocity in x direction
   double v (int q, int w);										// update v_y or velocity in y direction

	////###POST PROCESSING FUNCTIONS###///////
   void error (double aux[nx+1][ny+1][3], double aux1[3]); 					// calculate error; differnce between current solution and fully-developed temp. field
   void exact_flux (int q, int w, double aux[3]);							// return analytical flux at cell q,w
   void ave_value (const double aux[3][nx+2][ny+2], double aux1[3]); 		// returns average of each field variable over the domain
   void generate_vtkfile(const double s, const double aux[3][nx+2][ny+2], int k); 	// print a vtk file of the solution for visualization
   void generate_bfile_3var(const double s, const double aux[3][nx+2][ny+2]); // print a vtk file of the solution for visualization (3 variables)
   void generate_bfile_1var(const double aux[nx+1][ny+1]); 					// print a vtk file of the solution for visualization (1 variables)
   void vortex_core (const double aux[3][nx+2][ny+2]); 						// returns the maximum vorticity strength in domain
   void flag_vortex (const double aux[3][nx+2][ny+2], int faux[nx+1][ny+1]);// find and flag all vortecies
   void symm_values (const double aux[nx+2][ny+2]);							// returns the values on vertical symmetry line for a specific variable 
   double vortex_strg (const double aux[3][nx+2][ny+2], const int q, const int w);// returns vortex strength

#include "analysis.h"
#include "post_processing.h"	

#endif
