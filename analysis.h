// in this header main analysis functiones are defined
#ifndef _analysis_h_included_ /// just to make sure nothing was defined erlier or somewhere else//*not neccessary here*//
#define _analysis_h_included_

#include "prelim.h"
#include "thomas_solver.h"	// Thomas system solver

////////////////****** UPDATE FLUX *****////////////////
void update_fi (const double aux[3][nx+2][ny+2], double aux1[nx+1][ny+1][3])
{
	for(int j=1; j < ny+1; j++){ 
		for(int i=1; i < nx+1; i++){ 

			double f_i [2][3]; // f_i[0][] --> 3 components of F_(i+0.5,j) ; f_i[1][] --> 3 components of F_(i-0.5,j)
			double g_i [2][3]; // g_i[0][] --> 3 components of G_(i,j+0.5) ; g_i[1][] --> 3 components ofG_(i,j-0.5)

			f_i[0][0] =	  ((aux[1][i+1][j] + aux[1][i][j])/(2.*beta))		// first component of F_(i+0.5,j)
						+ (((aux[0][i+1][j] - aux[0][i][j])*a_fac*dy));     // a_fact for pressure laplacian added for filter	

			f_i[0][1] =   pow((0.5*(aux[1][i+1][j] + aux[1][i][j])),2.)		// second component of F_(i+0.5,j)
						+ (0.5*(aux[0][i+1][j] + aux[0][i][j]))
						- ((aux[1][i+1][j] - aux[1][i][j])/(re*dx));

			f_i[0][2] =   ((0.5*(aux[1][i+1][j] + aux[1][i][j]))*			// third component of F_(i+0.5,j)
							(0.5*(aux[2][i+1][j] + aux[2][i][j])))
						- ((aux[2][i+1][j] - aux[2][i][j])/(re*dx));			

			f_i[1][0] =   ((aux[1][i][j] + aux[1][i-1][j])/(2.*beta))		// first component of F_(i-0.5,j)
						+ (((aux[0][i][j] - aux[0][i-1][j])*a_fac*dy));			

			f_i[1][1] =   pow((0.5*(aux[1][i][j] + aux[1][i-1][j])),2.)		// second component of F_(i-0.5,j)
						+ (0.5*(aux[0][i][j] + aux[0][i-1][j]))
						- ((aux[1][i][j] - aux[1][i-1][j])/(re*dx));		

			f_i[1][2] =   ((0.5*(aux[1][i][j] + aux[1][i-1][j]))			// third component of F_(i-0.5,j)
						*  (0.5*(aux[2][i][j] + aux[2][i-1][j])))
						- ((aux[2][i][j] - aux[2][i-1][j])/(re*dx));			

			g_i[0][0] =   ((aux[2][i][j+1] + aux[2][i][j])/(2.*beta))		// first component of G_(i,j+0.5)
						+ (((aux[0][i][j+1] - aux[0][i][j])*a_fac*dx));
	
			g_i[0][1] =   ((0.5*(aux[1][i][j] + aux[1][i][j+1])) 			// second component of G_(i,j+0.5)
						*  (0.5*(aux[2][i][j] + aux[2][i][j+1]))) 
						- ((aux[1][i][j+1] - aux[1][i][j])/(re*dy));			
	
			g_i[0][2] =   pow((0.5*(aux[2][i][j] + aux[2][i][j+1])),2.)		// third component of G_(i,j+0.5)
						+ (0.5*(aux[0][i][j] + aux[0][i][j+1]))
						- ((aux[2][i][j+1] - aux[2][i][j])/(re*dy));			
	
			g_i[1][0] =   ((aux[2][i][j] + aux[2][i][j-1])/(2.*beta))		// first component of G_(i,j-0.5)
						+ (((aux[0][i][j] - aux[0][i][j-1])*a_fac*dx*dy)/dy);     				
	
			g_i[1][1] =   ((0.5*(aux[1][i][j-1] + aux[1][i][j])) 			// second component of G_(i,j-0.5)
						*  (0.5*(aux[2][i][j-1] + aux[2][i][j]))) 
						- ((aux[1][i][j] - aux[1][i][j-1])/(re*dy));			

			g_i[1][2] =   pow((0.5*(aux[2][i][j-1] + aux[2][i][j])),2.)		// third component of G_(i,j-0.5)
						+ (0.5*(aux[0][i][j-1] + aux[0][i][j]))
						- ((aux[2][i][j] - aux[2][i][j-1])/(re*dy));			

			aux1[i][j][0] = - ((f_i[0][0] - f_i[1][0])/dx + (g_i[0][0] - g_i[1][0])/dy) ;	// first component of flux integral
			aux1[i][j][1] = - ((f_i[0][1] - f_i[1][1])/dx + (g_i[0][1] - g_i[1][1])/dy) ;	// second component of flux integral
			aux1[i][j][2] = - ((f_i[0][2] - f_i[1][2])/dx + (g_i[0][2] - g_i[1][2])/dy) ;	// third component of flux integral	

			for (int q = 0; q < 3; q++) if (aux1[i][j][q] == 0.)  aux1[i][j][q] = 0.0; // just to get rid of minus zero
		}
	}
}

///////////////*****SET GHOST CELL******///////////////
void set_gcells(double aux[3][nx+2][ny+2])
{
   for(int j=0; j < ny+2; j++){ 					// loop in y diraction
		for(int i=0; i < nx+2; i++){ 				// loop in x diraction or in a row

			if (i==0){								/// check left boundary
	
				aux[0][i][j] = aux[0][i+1][j];		// no-slip condition for P
				aux[1][i][j] = - aux[1][i+1][j];			// no-slip condition for u
				aux[2][i][j] = - aux[2][i+1][j];	// no-slip condition for v

				if (aux[0][i+1][j] == 0.) aux[0][i][j] = 0.; // to make sure no variable set to -0
				if (aux[1][i+1][j] == 0.) aux[1][i][j] = 0.;
				if (aux[2][i+1][j] == 0.) aux[2][i][j] = 0.;
			}

			if (i==nx+1 ){							//check right boundary

				aux[0][i][j] = aux[0][i-1][j];		// fully developed field, neumann bdc for P
				aux[1][i][j] = - aux[1][i-1][j];	// no-slip condition for u
				aux[2][i][j] = - aux[2][i-1][j];	// no-slip condition for v

				if (aux[0][i-1][j] == 0.) aux[0][i][j] = 0.; // to make sure no variable set to -0
				if (aux[1][i-1][j] == 0.) aux[1][i][j] = 0.;
				if (aux[2][i-1][j] == 0.) aux[2][i][j] = 0.;
			}

	  		if (j==ny+1 ){							//check bottom boundary

				aux[0][i][j] = aux[0][i][j-1];		// fully developed field, neumann bdc for P
				aux[1][i][j] = - aux[1][i][j-1];	// no-slip bdc for u 
				aux[2][i][j] = - aux[2][i][j-1];	// no-slip bdc for v

				if (aux[0][i][j-1] == 0.) aux[0][i][j] = 0.; // to make sure no variable set to -0
				if (aux[1][i][j-1] == 0.) aux[1][i][j] = 0.;
				if (aux[2][i][j-1] == 0.) aux[2][i][j] = 0.;
			}

	 		if (j==0){  							//check top boundary

				aux[0][i][j] = aux[0][i][j+1];						// fully developed field, neumann bdc for P
				aux[1][i][j] = (2.* u_top) - aux[1][i][j+1];		// no-slip moving wall bdc with u = u_top
				aux[2][i][j] = - aux[2][i][j+1];					// no-slip bdc for v

				if (aux[0][i][j+1] == 0.) aux[0][i][j] = 0.; // to make sure no variable set to -0
				if (aux[1][i][j+1] == 0.) aux[1][i][j] = 0.;
				if (aux[2][i][j+1] == 0.) aux[2][i][j] = 0.;
			}
		}
	}
}

////////////////****** CALCULATE 3 JACOBEAN MATRIX OF THE DISCRETIZATION IN X DIRECTION *****////////////////
void Jacob_x(int i, int j, double aux1[3][3], 
				double aux2[3][3], double aux3[3][3], const double auxf[3][nx+2][ny+2])
{
	aux1[2][0] = aux1[1][2] = aux1[0][2] = 0.;   		// set A_x materix components
	aux1[0][0] = a_fac*(dy/dx); // non-zero term if filter is on; zero if filter off (a_fac = 0)
	aux1[2][1] = -(0.25*(auxf[2][i-1][j] + auxf[2][i][j]))/dx;
	aux1[2][2] = -(0.25*(auxf[1][i-1][j] + auxf[1][i][j]) + 1./(dx*re))/dx;
	aux1[1][0] = -1./2./dx;
	aux1[1][1] = -(0.5*(auxf[1][i-1][j] + auxf[1][i][j]) + 1./(dx*re))/dx;
	aux1[0][1] = -1./(2.*beta*dx);

	aux2[2][0] = aux2[1][2] = aux2[0][2] = aux2[1][0] = aux2[0][1] = 0.;		// set B_x materix components
	aux2[0][0] = 2. * a_fac*(dy/dx);// non-zero term if filter is on
	aux2[2][1] = (auxf[2][i+1][j] - auxf[2][i-1][j])/(4.*dx);
	aux2[2][2] = ((auxf[1][i+1][j] - auxf[1][i-1][j])/4./dx) + (2./re/dx/dx);
	aux2[1][1] = ((auxf[1][i+1][j] - auxf[1][i-1][j])/2./dx) + (2./re/dx/dx);

	aux3[2][0] = aux3[1][2] = aux3[0][2] = 0.;   		// set C_x materix components
	aux3[0][0] = a_fac*(dy/dx);// non-zero term if filter is on
	aux3[2][1] = (1./4./dx)*(auxf[2][i][j] + auxf[2][i+1][j]);
	aux3[2][2] = ((1./4./dx)*(auxf[1][i][j] + auxf[1][i+1][j])) - (1./(dx*dx*re));
	aux3[1][0] = 1./2./dx;
	aux3[1][1] = ((1./2./dx)*(auxf[1][i][j] + auxf[1][i+1][j])) - (1./(dx*dx*re));
	aux3[0][1] = 1./(2.*beta*dx);

	for (unsigned int q = 0; q < 3; q++){			// to get rid of minus zeros
		for (unsigned int w = 0; w < 3; w++){
			if (aux1[q][w] == 0.) aux1[q][w] = 0.;
			if (aux2[q][w] == 0.) aux2[q][w] = 0.;
			if (aux3[q][w] == 0.) aux3[q][w] = 0.;
		}
	}
}

////////////////****** CALCULATE 3 JACOBEAN MATRIX OF THE DISCRETIZATION IN X DIRECTION *****////////////////
void Jacob_y(int i, int j, double aux1[3][3], double aux2[3][3], double aux3[3][3], const double auxf[3][nx+2][ny+2])
{
	aux1[2][1] = aux1[1][0] = aux1[0][1] = 0.;			// set A_y materix components
	aux1[0][0] = a_fac*(dx/dy); // non-zero term if filter is on; zero if filter off (a_fac = 0)	
	aux1[2][0] = -1./dy/2.; 
	aux1[2][2] = ((-1./2./dy)*(auxf[2][i][j-1] + auxf[2][i][j])) - (1./re/dy/dy);
	aux1[1][1] = ((-1./4./dy)*(auxf[2][i][j-1] + auxf[2][i][j])) - (1./re/dy/dy);
	aux1[1][2] = (-1./4./dy)*((auxf[1][i][j-1] + auxf[1][i][j]));
	aux1[0][2] = -1./2./dy/beta ;  		

	aux2[2][0] = aux2[2][1] = aux2[1][0] = aux2[0][1] = aux2[0][2] = 0.; // set B_y materix components
	aux2[0][0] = 2. * a_fac*(dx/dy);// non-zero term if filter is on
	aux2[2][2] = (1./2./dy)*(auxf[2][i][j+1] - auxf[2][i][j-1]) + (2./re/dy/dy);
	aux2[1][1] = (1./4./dy)*(auxf[2][i][j+1] - auxf[2][i][j-1]) + (2./re/dy/dy);
	aux2[1][2] = (1./4./dy)*(auxf[1][i][j+1] - auxf[1][i][j-1]);

	aux3[2][1] = aux3[1][0] = aux3[0][1] = 0.; 			// set C_y materix components
	aux3[0][0] = a_fac*(dx/dy);// non-zero term if filter is on
	aux3[2][0] = 1./2./dy;
	aux3[2][2] = ((1./2./dy)*(auxf[2][i][j] + auxf[2][i][j+1])) - (1./re/dy/dy);
	aux3[1][1] = ((1./4./dy)*(auxf[2][i][j] + auxf[2][i][j+1])) - (1./re/dy/dy);
	aux3[1][2] = (1./4./dy)*(auxf[1][i][j] + auxf[1][i][j+1]);
	aux3[0][2] = 1./2./beta/dy;

	for (unsigned int q = 0; q < 3; q++){			// to get rid of minus zeros
		for (unsigned int w = 0; w < 3; w++){
			if (aux1[q][w] == 0.) aux1[q][w] = 0.;
			if (aux2[q][w] == 0.) aux2[q][w] = 0.;
			if (aux3[q][w] == 0.) aux3[q][w] = 0.;
		}
	}	
}

///////////////***** APPROXIMATE FACTORIZATION - FIRST STEP - EVERY ROW******///////////////
void appx_fac_dx (const double aux_f[3][nx+2][ny+2], double aux_fi[nx+1][ny+1][3])
{
	
    double I[3][3] = {0.};   		// define identity matrix
    I[0][0] = I[1][1] = I[2][2] = 1.;
	
			for(int j=1; j < ny+1; j++){				//*******// going through rows as first part of appx. factorization

			double rhs [nx+2][3];			//  right-hand-side; b in Ax =b
			rhs[0][0] = 
			rhs[0][1] = 
			rhs[0][2] =	0.;					// set values for bdc;  ghost cells
			rhs[nx+1][0] = 
			rhs[nx+1][1] = 
			rhs[nx+1][2] = 0.;				// set values for bdc;  ghost cells

			double lhs [nx+2][3][3][3];		// left-hand-side; A in Ax =b

			// assign identity matrix to wall bdc
			Copy3x3(I, lhs[0][1]);			// bdc on first block on main diagonal
			Copy3x3(I, lhs [nx+1][1]);		// bdc on last block on main diagonal
			Copy3x3(I, lhs[0][2]);			// bdc on first block above main diagonal
			lhs[0][2][0][0] = -1.;			// set neumaan bdc for pressure							
			Copy3x3(I, lhs [nx+1][0]);		// bdc on last block below main diagonal
			lhs[nx+1][0][0][0] = -1.;		// set neumaan bdc for pressure

			for(int i=1; i < nx+1; i++){						

				rhs[i][0] = dt * aux_fi[i][j][0];				// set RHS and LHS  for each cell in the row
				rhs[i][1] = dt * aux_fi[i][j][1];
				rhs[i][2] = dt * aux_fi[i][j][2];

				double a_x [3][3];								// set jacobean matrix for cell i,j
				double b_x [3][3];
				double c_x [3][3];	
				Jacob_x(i,j,a_x,b_x,c_x,aux_f);

				MultSca(a_x,dt);			// set blocks in vector below main diagonal
				Copy3x3(a_x, lhs[i][0]);
				Add3x3(I, b_x, dt, b_x);	// set blocks in main diagonal vector
				Copy3x3(b_x, lhs[i][1]);
				MultSca(c_x,dt);			// set blocks in vector above main diagonal
				Copy3x3(c_x, lhs[i][2]);													
			}

			SolveBlockTri( lhs, rhs , nx+2); // slove Ax=b with block Thomas algorithm

			for(int i=1; i < nx+1; i++) CopyVec(rhs[i], aux_fi[i][j]);	// fi (initially flux integral vec) is overwritten with new data (delta_f)
		}
}

///////////////***** APPROXIMATE FACTORIZATION - SECOND STEP_EVERY COLUMN******///////////////
void appx_fac_dy(const double aux_f[3][nx+2][ny+2], double aux_fi[nx+1][ny+1][3])
{ 
  
    double I[3][3] = {0.};   		// define identity matrix
    I[0][0] = I[1][1] = I[2][2] = 1.;

   	for(int i=1; i < nx+1; i++){  				//*******// going through columns as second part of appx. factorization			

		double rhs [ny+2][3];								//  right-hand-side; b in Ax =b
		rhs[0][0] = rhs[0][2] =	rhs[0][1] = 0.;				// set values for bdc;  ghost cells
		rhs[ny+1][0] = rhs[ny+1][1] = rhs[ny+1][2] = 0.;	// set values for bdc;  ghost cells

		double lhs [ny+2][3][3][3];							// left-hand-side; A in Ax =b

		// assign identity matrix to wall bdc
		Copy3x3(I, lhs[0][1]);								// bdc on first block on main diagonal
		Copy3x3(I, lhs [ny+1][1]);							// bdc on last block on main diagonal
		Copy3x3(I, lhs[0][2]);								// bdc on first block above main diagonal
		lhs[0][2][0][0] = -1.;								// set neumaan bdc for pressure
		Copy3x3(I, lhs [ny+1][0]);							// bdc on last block below main diagonal
		lhs[ny+1][0][0][0] = -1.;							// set neumaan bdc for pressure

		for(int j=1; j < ny+1; j++){	// set RHS and LHS  for each cell in the row					

			CopyVec(aux_fi[i][j], rhs[j]);  // copy the results from previous step (delta_f) into the rhs vector 

			double a_y [3][3];								// set jacobean matrix for cell i,j
			double b_y [3][3];
			double c_y [3][3];	
			Jacob_y(i,j,a_y,b_y,c_y,aux_f);

			MultSca(a_y,dt);			// set blocks for below main diagonal
			Copy3x3(a_y, lhs[j][0]);
			Add3x3(I, b_y, dt, b_y);	// set blocks for main diagonal
			Copy3x3(b_y, lhs[j][1]);
			MultSca(c_y,dt);			// set blocks for above main diagonal
			Copy3x3(c_y, lhs[j][2]);													
		}

		SolveBlockTri( lhs, rhs, ny+2);

		for(int j=1; j < ny+1; j++) CopyVec(rhs[j], aux_fi[i][j]);
	}
}

///////////////*****UPDATE VARIABLES IN TIME & RETURN ERROR******///////////////
void update_f (double aux_f[3][nx+2][ny+2], double aux_fi[nx+1][ny+1][3], double err[3])
{
	
	for(int j=1; j < ny+1; j++)		// update final solution of appx. factor. or update f as well as L2
    	for(int i=1; i < nx+1; i++){
			double aux [3];
			aux [0] = aux_f[0][i][j] ;
			aux [1] = aux_f[1][i][j] ;
			aux [2] = aux_f[2][i][j] ;

		if (omega == 1.){		// if omega = 1 --> no over-relaxation

    		aux_f[0][i][j] += aux_fi[i][j][0];	// update solutions without over-relaxation
    		aux_f[1][i][j] += aux_fi[i][j][1];
    		aux_f[2][i][j] += aux_fi[i][j][2];
		}
		else {
				aux_f[0][i][j] += omega * (aux_fi[i][j][0]-aux_f[0][i][j]);	// update solutions with over-relaxation
				aux_f[1][i][j] += omega * (aux_fi[i][j][1]-aux_f[1][i][j]);
				aux_f[2][i][j] += omega * (aux_fi[i][j][2]-aux_f[2][i][j]);
   		}

			err [0] += ((aux [0]-aux_f[0][i][j])*(aux [0]-aux_f[0][i][j]))/(nx*ny);  // update L2 norm for each cell
    		err [1] += ((aux [1]-aux_f[1][i][j])*(aux [1]-aux_f[1][i][j]))/(nx*ny);
    		err [2] += ((aux [2]-aux_f[2][i][j])*(aux [2]-aux_f[2][i][j]))/(nx*ny);
		}

	err [0] = sqrt (err[0]);
	err [1] = sqrt (err[1]);
	err [2] = sqrt (err[2]);
}

#endif
