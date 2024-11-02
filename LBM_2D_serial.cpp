#include <bits/stdc++.h>
using namespace std;

// The domain is a square domain with side of 0.2
// u = 0.1 LU, rho_in = 5, alpha = 0.01 LU are fixed in LU according to example in AA Moh. Book and others are calculated accordingly...
double rho_in = 5.0, u_in = 0.1, tolerance = 1e-8, alpha = 0.01;
double dx = 1.0, dy = 1.0, dt = 1.0; // All in LU
int nx = 101, ny = 101;

double Re = u_in*nx/alpha;
double omega = 1.0/(3.*alpha + 0.5);

vector<double> w = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
vector<pair<double,double>> C = {{0,0}, {0,-1}, {1,0}, {0,1}, {-1,0}, {1,-1}, {1,1}, {-1,1}, {-1,-1}};

typedef struct F
{
	double f0;
	double f1; 
	double f2;
	double f3;
	double f4;
	double f5;
	double f6;
	double f7;
	double f8;
}F;

void f_eq (vector<vector<F>>& f_equil, vector<vector<double>>& rho, vector<vector<double>>& u, vector<vector<double>>& v)
{
	for(int i=0; i<ny; i++){
		for(int j=0; j<nx; j++){
			f_equil[i][j].f0 = w[0]*rho[i][j]*(1 + 3*(C[0].first*u[i][j]+C[0].second*v[i][j]) + 4.5*pow((C[0].first*u[i][j]+C[0].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) );
			f_equil[i][j].f1 = w[1]*rho[i][j]*(1 + 3*(C[1].first*u[i][j]+C[1].second*v[i][j]) + 4.5*pow((C[1].first*u[i][j]+C[1].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f2 = w[2]*rho[i][j]*(1 + 3*(C[2].first*u[i][j]+C[2].second*v[i][j]) + 4.5*pow((C[2].first*u[i][j]+C[2].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f3 = w[3]*rho[i][j]*(1 + 3*(C[3].first*u[i][j]+C[3].second*v[i][j]) + 4.5*pow((C[3].first*u[i][j]+C[3].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f4 = w[4]*rho[i][j]*(1 + 3*(C[4].first*u[i][j]+C[4].second*v[i][j]) + 4.5*pow((C[4].first*u[i][j]+C[4].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f5 = w[5]*rho[i][j]*(1 + 3*(C[5].first*u[i][j]+C[5].second*v[i][j]) + 4.5*pow((C[5].first*u[i][j]+C[5].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f6 = w[6]*rho[i][j]*(1 + 3*(C[6].first*u[i][j]+C[6].second*v[i][j]) + 4.5*pow((C[6].first*u[i][j]+C[6].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f7 = w[7]*rho[i][j]*(1 + 3*(C[7].first*u[i][j]+C[7].second*v[i][j]) + 4.5*pow((C[7].first*u[i][j]+C[7].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) ); 
			f_equil[i][j].f8 = w[8]*rho[i][j]*(1 + 3*(C[8].first*u[i][j]+C[8].second*v[i][j]) + 4.5*pow((C[8].first*u[i][j]+C[8].second*v[i][j]),2) - 1.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]) );  
		}
	}
}

void apply_BC (vector<vector<F>>& f_new)
{
	for(int i=0; i<ny; i++){
		f_new[i][0].f2 = f_new[i][0].f4;
		f_new[i][0].f6 = f_new[i][0].f8;
		f_new[i][0].f5 = f_new[i][0].f7;

		f_new[i][nx-1].f4 = f_new[i][nx-1].f2;
		f_new[i][nx-1].f8 = f_new[i][nx-1].f6;
		f_new[i][nx-1].f7 = f_new[i][nx-1].f5;
	}

	for(int j=0; j<nx; j++){
		f_new[0][j].f3 = f_new[0][j].f1;
		f_new[0][j].f6 = f_new[0][j].f8;
		f_new[0][j].f7 = f_new[0][j].f5;
	}

	for(int j=1; j<nx-1; j++){
		double rhon = f_new[ny-1][j].f0 + f_new[ny-1][j].f2 + f_new[ny-1][j].f4 + 2*(f_new[ny-1][j].f3 + f_new[ny-1][0].f6 + f_new[ny-1][j].f7);
		f_new[ny-1][j].f1 = f_new[ny-1][j].f3;
		f_new[ny-1][j].f5 = f_new[ny-1][j].f7 + (rhon*u_in)/6.0;
		f_new[ny-1][j].f8 = f_new[ny-1][j].f6 - (rhon*u_in)/6.0;
	}
}

int main()
{
	auto start = std::chrono::steady_clock::now();

	vector<vector<F>> f_old(nx, vector<F>(ny)), f_new(nx, vector<F>(ny));
	vector<vector<F>> f_equil(nx, vector<F>(ny));
	vector<vector<double>> u_old(nx, vector<double>(ny)), u_new(nx, vector<double>(ny));
	vector<vector<double>> v_old(nx, vector<double>(ny)), v_new(nx, vector<double>(ny));
	vector<vector<double>> rho_old(nx, vector<double>(ny)), rho_new(nx, vector<double>(ny));
	vector<vector<F>> f_inter(nx, vector<F>(ny));

	FILE *Field_2D, *Error, *Field_x_01, *Field_y_01;

	// Error = fopen("Error.dat", "w");
	// Field_2D = fopen("Field_2D.dat", "w");
	Field_x_01 = fopen("Field_x_01.dat", "w");
	Field_y_01 = fopen("Field_y_01.dat", "w");

	// fprintf(Error, "Iteration  uerror  verror\n");
	fprintf(Field_x_01, "y\tvalue\n");
	fprintf(Field_y_01, "x\tvalue\n");

	double uerror = 1.0, verror = 1.0;

// Initialising the u, v, rho vectors
	for(int i=0; i<ny; i++){
		for(int j=0; j<nx; j++){
			u_old[i][j] = 0.0;
			u_new[i][j] = 0.0;
			v_old[i][j] = 0.0;
			v_new[i][j] = 0.0;
			rho_old[i][j] = rho_in;
			rho_new[i][j] = rho_in;			 
		}
	}

	for(int j=1; j<nx-1; j++){
		u_old[ny-1][j] = u_in;
		u_new[ny-1][j] = u_in;
		v_old[ny-1][j] = 0.0;
		v_new[ny-1][j] = 0.0;
	}

	// Initialising the f vectors equal to the f at equilibrium
	f_eq(f_equil, rho_old, u_old, v_old);
	
	f_old = f_equil;
	f_new = f_equil;

	// CPU printing to console the message when init. done and time marching loop starts.
	printf("\nStarting the main calculation loop\n\n");

	int iter = 0;
	while((uerror > tolerance) || (verror > tolerance))
	{
		if (iter != 0)
		{
			// Calculate f_equil for the current iteration rho, u, v values.
			f_eq (f_equil, rho_old, u_old, v_old);
		}

		// Collision step
		for(int i=0; i<ny; i++){
			for(int j=0; j<nx; j++){
				f_inter[i][j].f0 = omega*f_equil[i][j].f0 + (1.0-omega)*f_old[i][j].f0;
				f_inter[i][j].f1 = omega*f_equil[i][j].f1 + (1.0-omega)*f_old[i][j].f1;
				f_inter[i][j].f2 = omega*f_equil[i][j].f2 + (1.0-omega)*f_old[i][j].f2;
				f_inter[i][j].f3 = omega*f_equil[i][j].f3 + (1.0-omega)*f_old[i][j].f3;
				f_inter[i][j].f4 = omega*f_equil[i][j].f4 + (1.0-omega)*f_old[i][j].f4;
				f_inter[i][j].f5 = omega*f_equil[i][j].f5 + (1.0-omega)*f_old[i][j].f5;
				f_inter[i][j].f6 = omega*f_equil[i][j].f6 + (1.0-omega)*f_old[i][j].f6;
				f_inter[i][j].f7 = omega*f_equil[i][j].f7 + (1.0-omega)*f_old[i][j].f7;
				f_inter[i][j].f8 = omega*f_equil[i][j].f8 + (1.0-omega)*f_old[i][j].f8;
			}
		}

		// Streaming step
		for(int i=0; i<ny; i++){
			for(int j=0; j<nx; j++){
				f_new[i][j].f0 = f_inter[i][j].f0;
			}
		}
		
		for(int i=0; i<ny; i++){
			// RIGHT TO LEFT!
			for(int j=nx-1; j>0; j--){
				f_new[i][j].f2 = f_inter[i][j-1].f2;
			}

			// LEFT TO RIGHT!
			for(int j=0; j<nx-1; j++){
				f_new[i][j].f4 = f_inter[i][j+1].f4;
			}
		}

		for(int i=ny-1; i>0; i--){

			for(int j=0; j<nx; j++){
				f_new[i][j].f3 = f_inter[i-1][j].f3;
			}

			for(int j=nx; j>0; j--){
				f_new[i][j].f6 = f_inter[i-1][j-1].f6;
			}

			for(int j=0; j<nx-1; j++){
				f_new[i][j].f7 = f_inter[i-1][j+1].f7;
			}
		}

		for(int i=0; i<ny-1; i++){

			for(int j=0; j<nx; j++){
				f_new[i][j].f1 = f_inter[i+1][j].f1;
			}

			for(int j=0; j<nx-1; j++){
				f_new[i][j].f8 = f_inter[i+1][j+1].f8;
			}

			for(int j=nx-1; j>0; j--){
				f_new[i][j].f5 = f_inter[i+1][j-1].f5;
			}
		}

		// apply boundary conditions after streaming
		apply_BC(f_new);

		// calculate u_new, v_new, rho_new, uerror & verror
		double sum_u = 0.0, sum_v = 0.0;

		for(int i=0; i<ny; i++){
			for(int j=0; j<nx; j++){
				rho_new[i][j] = f_new[i][j].f0 + f_new[i][j].f1 + f_new[i][j].f2 + f_new[i][j].f3 + f_new[i][j].f4 + f_new[i][j].f5 + f_new[i][j].f6 + f_new[i][j].f7 + f_new[i][j].f8;
			}
		}

		for(int j=1; j<nx-1; j++)
			rho_new[ny-1][j] = f_new[ny-1][j].f0 + f_new[ny-1][j].f2 + f_new[ny-1][j].f4 + 2*(f_new[ny-1][j].f3 + f_new[ny-1][j].f6 + f_new[ny-1][j].f7);

		for(int i=1; i<ny-1; i++){
			for(int j=1; j<nx-1; j++){
				u_new[i][j] = (f_new[i][j].f2 - f_new[i][j].f4 + f_new[i][j].f5 + f_new[i][j].f6 - f_new[i][j].f7 - f_new[i][j].f8)/rho_new[i][j];
				v_new[i][j] = (f_new[i][j].f3 - f_new[i][j].f1 - f_new[i][j].f5 + f_new[i][j].f6 + f_new[i][j].f7 - f_new[i][j].f8)/rho_new[i][j];

				sum_u += pow(u_new[i][j] - u_old[i][j], 2);
				sum_v += pow(v_new[i][j] - v_old[i][j], 2);
 			}
		}

		uerror = sqrt(sum_u) / (nx*ny);
		verror = sqrt(sum_v) / (nx*ny);

		// updating the f_old to f_new calcuated
		f_old = f_new;

		// update rho, u, v for current iteration from the rho, u, v we got from previous iteration for this iteration
		rho_old = rho_new;
		u_old = u_new;
		v_old = v_new;

		iter++;

		// fprintf(Error, "%d %.9lf %.9lf\n", iter, uerror, verror);

		// Printing the error at each 100th iteration which is performed..
		if (iter % 100 == 0)
			printf("uerror: %.9lf & verror: %.9lf for iteration - %d\n", uerror, verror, iter);
	}

	printf("\nProgram is converged with the Final error of: \n");
	printf("uerror: %.9lf & verror: %.9lf for iteration - %d\n", uerror, verror, iter);

	auto end = std::chrono::steady_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

    // Writing the output files.
	// for(int i=0; i<ny; i++){
	// 	for(int j=0; j<nx; j++){
	// 		double x = 0.0 + j*(0.2/(nx-1));
	// 		double y = 0.0 + i*(0.2/(ny-1));
	// 		double p = rho_new[i][j]*(1.0/3);

	// 		fprintf(Field_2D, "%.8lf %.8f %.8lf %.8lf %.8lf\n", x, y, u_new[i][j], v_new[i][j], p);
	// 	}
	// }

	for(int i=0; i<ny; i++){
		double x = 0.1;
		double y = 0.0 + i*(0.2/(ny-1));
		double p = rho_new[i][50]*(1.0/3);

		fprintf(Field_x_01, "%.8lf %.8f %.8lf %.8lf %.8lf\n", x, y, u_new[i][50], v_new[i][50], p);
	}

	for(int j=0; j<nx; j++){
		double x = 0.0 + j*(0.2/(nx-1));
		double y = 0.1;
		double p = rho_new[50][j]*(1.0/3);

		fprintf(Field_y_01, "%.8lf %.8f %.8lf %.8lf %.8lf\n", x, y, u_new[50][j], v_new[50][j], p);
	}


	// fclose(Field_2D);
	// fclose(Error);
	fclose(Field_x_01);
	fclose(Field_y_01);

	return 0;
}