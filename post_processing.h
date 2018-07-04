#ifndef _post_processing_h_included_ /// just to make sure nothing was defined erlier or somewhere else
#define _post_processing_h_included_

#include "prelim.h"

////////////////****** CALC AVERAGE VALUES *****////////////////
void ave_value (const double aux[3][nx+2][ny+2], double aux1[3])
{
	for(int j=1; j < ny+1; j++)
    	for(int i=1; i < nx+1; i++){
    		aux1[0] += aux[0][i][j]/(nx*ny);
			aux1[1] += aux[1][i][j]/(nx*ny);
			aux1[2] += aux[2][i][j]/(nx*ny);    			
		}	
}

////////////////****** PRINT VALUES ON VERTICAL SYMMETRY LINE *****////////////////
void symm_values (const double aux[nx+2][ny+2])
{	
	FILE *middval;
	middval  = fopen ("midd_values.txt", "wt"); // open the file to record results

	double aux1[ny+1];

	if ( nx%2 == 0){ // if nx is even, symmetry value is average of the 2 middle cells
		int midd = nx/2;
		for(int j=1; j < ny+1; j++)
			aux1[j] = 0.5*(aux[midd][j]+aux[midd+1][j]);
	}

	else{ // if nx is odd, symmetry value is vlaue at (nx-1)/2
		int midd = (nx-1)/2;
		for(int j=1; j < ny+1; j++)
			aux1[j] = aux[midd][j];
	}
	
	for(int j=1; j < ny+1; j++) fprintf (middval, "%e\t %e\n",(((j-1)*dy)+(dy/2)), aux1[j]);

	fclose (middval);
}

///////////////****** FIND & FLAG VORTICES *****////////////////
void flag_vortex (const double aux[3][nx+2][ny+2], int faux[nx+1][ny+1]){

	for(int j=2; j < ny; j++){								// interior elements are checked for vortex
    	for(int i=2; i < nx; i++){							//assuming no vortex core on the wall
			if (aux[1][i][j-1]*aux[1][i][j+1] < 0)			// check if horizontal velocities at top and bottom of the cell are in opposite direction  
				if (aux[2][i-1][j]*aux[2][i+1][j] < 0) 		// checks if vertical velocities at right and left sides of the cell are in opposite direction
					if (aux[1][i][j-1]*aux[2][i-1][j] < 0)	// check if vortex has actual spiral circulation
					faux[i][j] = 1;							// if so, flag the element a vortex
    	}
    }
}

////////////////****** CALC VORTEX STRENGTH *****////////////////
double vortex_strg (const double aux[3][nx+2][ny+2], const int q, const int w)
{				// vortex strength is calcualted as second invarient Q
				// 2ed order centered scheme used for gradient estimation
	double u_x = (aux[1][q+1][w] - aux[1][q][w])/dx;
	double u_y = (aux[1][q][w+1] - aux[1][q][w])/dy;
	double v_x = (aux[2][q+1][w] - aux[2][q][w])/dx;
	double v_y = (aux[2][q][w+1] - aux[2][q][w])/dy;
	
	double invar = -0.5*((v_y*v_y)+(2.*u_y*v_x)+(u_x*u_x));
	
	return invar;	
}

////////////////****** RETURN VORTEX CORE POSITIONS & STRENGTH*****////////////////
void vortex_core (const double aux_f[3][nx+2][ny+2])
{
	FILE *vort;
	vort = fopen ("vortex.txt", "wt");

	int fcore [2][nx+1][ny+1] = {0};	// this variable caries 2 flags for each element of the domain; initialized as zeros
										// varible to flag vortex cells; if a cell is a vortex:fcore [0][nx+1][ny+1]=1 else:fcore [0][nx+1][ny+1]=0
										// if a vortex, ie fcore [0][nx+1][ny+1]=1, has been already considered as part of
										// a single core:fcore [1][nx+1][ny+1]=i where
										// "i" is # of the core the vortex belongs to

	int m = 0; 							// vortices counter

	flag_vortex(aux_f,fcore[0]); 		// find and flag all vortices

	for(int j=2; j < ny; ++j){			// no vortex on boundary so loop is restericted to interior cells
    	for(int i=2; i < nx; ++i){
    		
    		if (fcore[0][i][j] == 1 && fcore[1][i][j] == 0){		//if cell is flaged as vortex (fcore[0]=1) AND not assigned to a core (fcore[1]=0)
				
    			m ++;	
       			
    			fcore[1][i][j] = m ; 				// asign a number to the core
    			
    			double pos_x = ((i-1)*dx)+(dx/2);	// core x-axis position
    			double pos_y = ((j-1)*dy)+(dy/2);	// core y-axis position
    			double core_pres = aux_f[0][i][j]; 	// core presure
    			double ave_strg = vortex_strg (aux_f, i,j);			// core strength (average of vortices strength)						

    			for(int w=2; w < ny; w++){			// check the domain to check adjacent cells; 
    				for(int q=2; q < nx; q++){
    					if (fcore[0][q][w] == 1){	// check if the comparison cell is flaged as vortex
    						if (i != q || j != w){	// check if the primary cell is not checked by itself
							
    							double dist = sqrt ((((w-j)*dy)*((w-j)*dy)) + (((q-i)*dx)*((q-i)*dx))); // distance between cell (i,j) and (q,w)
    							
    							if (dist < rcore){				// check if the flaged cell fits in pre-set core radius limit 
    								fcore[1][q][w] = m;			// cell belongs to the same core as cell i,j
    								pos_x = 0.5*(pos_x + (((q-1)*dx)+(dx/2))); // with every new cell update core position	
    								pos_y = 0.5*(pos_y + (((w-1)*dy)+(dy/2)));
    								core_pres = 0.5*(core_pres + aux_f[0][q][w]); // with every new cell update average pressure
    								ave_strg = 	0.5*(ave_strg + vortex_strg (aux_f, q,w));				// update average core vortices' strength
								}
							}
						}
    				}
    			}   
			fprintf (vort, "Core# %d : located @ (%e, %e), average pressure: %e, average strength: %e\n" , m, pos_x, pos_y, core_pres, ave_strg); 			
    		}
    	}
    }

	fclose(vort);	

}

///////////////*****UPDATE L_2 ERROR (COMPARISON WITH EXACT SOLUTION)******///////////////
////// error is defined as the differnce between numerical and analytical flux
void error (double aux[nx+1][ny+1][3], double aux1[3])
{
	aux1[0] = aux1[1] = aux1[2] = 0.;

	for(int j=1; j < ny+1; j++){
    	for(int i=1; i < nx+1; i++){

    		double y = (((j-1)*dy)+(dy/2)) ; // y coordinate defined at center
    		double x = (((i-1)*dx)+(dx/2)) ; // x coordinate defined at center

			double comp_1 = (-pi/beta)*((u_0*cos(pi*x)*sin(2.*pi*y)) + (v_0*sin(2.*pi*x)*cos(pi*y)));	// first component of the analytical flux integral

			double comp_2 =   (p_0*pi*sin(pi*x)*cos(pi*y))												// second component of the analytical flux integral
							- (u_0*u_0*pi*sin(2.*pi*x)*pow(sin(2.*pi*y),2.))
							- ((u_0*v_0*pi*sin(pi*x)*sin(2.*pi*x))
							*   (cos(pi*y)*sin(2.*pi*y) + 2.*cos(2.*pi*y)*sin(pi*y)))
							- ((u_0/re)*(5.*pi*pi*sin(pi*x)*sin(2.*pi*y)));

			double comp_3 =   (p_0*pi*cos(pi*x)*sin(pi*y))												// third component of the analytical flux integral
							- (v_0*v_0*pi*sin(2.*pi*y)*pow(sin(2.*pi*x),2.))
							- ((u_0*v_0*pi*sin(pi*y)*sin(2.*pi*y))
							*   ((cos(pi*x)*sin(2.*pi*x)) + (2.*cos(2.*pi*x)*sin(pi*x))))
							- ((v_0/re)*(5.*pi*pi*sin(2.*pi*x)*sin(pi*y)));

			aux1[0] += ((aux[i][j][0] - comp_1)*(aux[i][j][0] - comp_1)) / (nx*ny);
			aux1[1] += ((aux[i][j][1] - comp_2)*(aux[i][j][1] - comp_2)) / (nx*ny);
			aux1[2] += ((aux[i][j][2] - comp_3)*(aux[i][j][2] - comp_3)) / (nx*ny);
		} 
	}
	aux1[0] = sqrt(aux1[0]);
	aux1[1] = sqrt(aux1[1]);
	aux1[2] = sqrt(aux1[2]);
}

///////////////****** GENERATE VTK FILE *****////////////////
void generate_vtkfile(const double s, const double aux[3][nx+2][ny+2], int k)
{
    FILE *out1;

    char name[50];
	
	double aux2d[nx+2][ny+2]; // record resluts in a 2d array vtk generator input
		for(int j=1; j < ny+1; j++)
    		for(int i=1; i < nx+1; i++){
    		if (u_top < 0.)	aux2d[i][j] = aux[k][nx-i+1][j]; // function records data in 1-X,Y if (u_top < 0)
			else   aux2d[i][j] = aux[k][i][j];  			
			}
		
	int m = 0;
	
	m = (int) (s * 100.);
	
  //  if (m%10 == 0 || s == 0.){
    
    sprintf (name, "Results_%d.vtk", m);

    out1 = fopen (name, "wt");

    fprintf(out1, "# vtk DataFile Version 2.0\n");
    fprintf(out1, "vtk output\n");
    fprintf(out1, "ASCII\n");
    fprintf(out1, "DATASET STRUCTURED_POINTS\n");
    fprintf(out1, "DIMENSIONS %d %d 1\n",nx,ny); 
    fprintf(out1, "SPACING 1 1 1\n");
    fprintf(out1, "ORIGIN 0 0 0\n");
    fprintf(out1, "POINT_DATA %d\n",(nx)*(ny)); 
    fprintf(out1, "SCALARS C float 1\n");
    fprintf(out1, "LOOKUP_TABLE default\n");

	for (int j = 1; j < ny+1; j++) {
	   for (int i = 1; i < nx+1; i++)
	   {
	      fprintf (out1, "%2.10e\n ", aux2d[i][j]);
	    }
	fprintf (out1, "\n");
	}
    fclose (out1);
    //}
}

///////////////****** GENERATE dat FILE - 3 variables *****////////////////
void generate_bfile_3var(const double s, const double aux[3][nx+2][ny+2])
{	
    FILE *out1;

    char name[50];
	
	int m = (int) (s * 100.);
    
    sprintf (name, "Results_%d.dat", m);

    out1 = fopen (name, "wt");

    fprintf(out1, "TITLE = \" Example: Multi-Zone 2D Plot \" \n" );
    fprintf(out1, "VARIABLES = \"X\", \"Y\", \"Pressure\", \"U\", \"V\" \n");
    fprintf(out1, "ZONE T=\"BIG ZONE\", I= %d, J= %d, DATAPACKING=POINT \n",nx,ny); 

	for (int j = 1; j < ny+1; j++) {
	   for (int i = 1; i < nx+1; i++)
	   {
	      fprintf (out1, "%d %d %e %e %e\n", i, j, aux[0][i][j], aux[1][i][j], aux[2][i][j]);
		}
	}
    fclose (out1);
}

///////////////****** GENERATE dat FILE - 1 variables*****////////////////
void generate_bfile_1var(const double aux[nx+1][ny+1]){
	
  	FILE *out1;

    char name[50];
	
//	int m = (int) (s * 100.);
    
    sprintf (name, "Results_angles.dat");

    out1 = fopen (name, "wt");

    fprintf(out1, "TITLE = \" Example: Multi-Zone 2D Plot \" \n" );
    fprintf(out1, "VARIABLES = \"X\", \"Y\", \"Angles\" \n");
    fprintf(out1, "ZONE T=\"BIG ZONE\", I= %d, J= %d, DATAPACKING=POINT \n",nx,ny); 

	for (int j = 1; j < ny+1; j++) {
	   for (int i = 1; i < nx+1; i++)
	   {
	      fprintf (out1, "%d %d %1.10e\n", i, j, aux[i][j]);
		}
	}
    fclose (out1);	
}

#endif
