#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
//#define N 10
// int N = 50;
#define random() ((double) rand() / (RAND_MAX))

void copy_3d(double ***src, double ***dest, int n);
//void *memcpy(void *, const void *, unsigned long);
void write_out_vtk(int nx, int ny, int nz, double dx, double dy, double dz, int istep, double ***data1);
void write_into_file(double ***a, int timestep, FILE *out);

int N, thread_num; 
int main(int argc, char **argv)
{
	N = atoi(argv[2]);
    	thread_num = atoi(argv[1]);
    	omp_set_num_threads(thread_num);
	srand(time(NULL));

	//FILE *out = fopen("out.txt", "w");
	double ***array_pred = (double ***)malloc(N*sizeof(double **));
	double ***array_post = (double ***)malloc(N*sizeof(double **));
	int i, j, k, n, timesteps = 10;
	int i_down, i_up, j_down, j_up, k_down, k_up;
	int H = 1000;
	double K = 1;
	float t = 1e-4, h = 8*1e-2, f;


        for (i = 0; i< N; i++) {

   		array_post[i] = (double **) malloc(N*sizeof(double *));
		array_pred[i] = (double **) malloc(N*sizeof(double *));

          	for (j = 0; j < N; j++) {

              		array_post[i][j] = (double *)malloc(N*sizeof(double));
			array_pred[i][j] = (double *)malloc(N*sizeof(double));
         	 }	

        }
	

	// initial conditions
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < N; k++)
			{
				array_pred[i][j][k] = 0;
				array_post[i][j][k] = 0;
			}
		}
	}
	
	n = 0;
//	array_pred[4][4][4] = 1;
//	array_pred[3][4][4] = 1;
//	array_pred[5][4][4] = 1;
//	array_pred[4][4][3] = 1;
//	array_pred[4][4][5] = 1;
	
	for (i = 4; i < 6; i++)
	{

		for (j = 4; j < 6; j++)
		{
			for (k = 4; k < 6; k++)
			{
				array_pred[i][j][k] = 1;
			}
		}
	}


	//write_into_file(array_pred, n, out);
	struct timeval start, end;
    	gettimeofday(&start, NULL);
	#pragma omp parallel for shared(array_pred, array_post) private(n, i, j, k)
	for (n = 0; n < timesteps; n++){
		for (i = 0; i < N; i++){
			i_down = i == 0 ? (N - 1) : i;
			i_up = i == (N - 1) ? 0 : i;
			for (j = 0; j < N; j++){
				j_down = j == 0 ? (N - 1) : j;
                        	j_up = j == (N - 1) ? 0 : j;
				for (k = 0; k < N; k++){
					k_down = k == 0 ? (N - 1) : k;
					k_up = k == (N - 1) ? 0 : k;
						
			 		f = t*H*(pow((array_pred[i][j][k]/1.1 - 1/2 + 0.17), 3) \
						- 1/4*(array_pred[i][j][k]/1.1 - 1/2 + 0.17) - 0.05);
		//			f = t*H*(pow((array_pred[i][j][k] - 1/2), 3) - 1/4*(array_pred[i][j][k] - 1/2));

		                        array_post[i][j][k] = K*(t/pow(h, 2) * (array_pred[i_down][j][k] - 6 * array_pred[i][j][k] \
                                         + array_pred[i_up][j][k] + array_pred[i][j_down][k] \
                                         + array_pred[i][j_up][k] + array_pred[i][j][k_down] \
                                         + array_pred[i][j][k_up]) - f) + array_pred[i][j][k]; 
				}
			}
		}
	//	write_into_file(array_post, n + 1, out);
	//	write_out_vtk(N, N, N, h, h, h, n, array_pred);
	//	memcpy(array_pred, array_post, sizeof(double **)*N);
	//	array_pred = array_post;
		copy_3d(array_post, array_pred, N);
	}

        gettimeofday(&end, NULL);
        double delta = ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
        //printf("%d %d %lf\n", thread_num, N, delta);
        printf("%lf\n", delta);
	//fclose(out);
	return 0;
}

void write_out_vtk(int nx, int ny, int nz, double dx, double dy, double dz, int istep, double ***data1)
{

   FILE *out;
   char fname[200];

   sprintf(fname, "equation_time_%d.vtu", istep);
   out = fopen(fname,"w");

   int npoin = nx*ny*nz;
   int i, j, k;
   double x, y, z;

 // ???????????? ASCII ?????????? ?? ???????????????? VTK:

 // ??????????????????  VTK ??????????

   fprintf(out,"# vtk DataFile Version 2.0\n");
   fprintf(out,"time_10.vtk\n");
   fprintf(out,"ASCII\n");
   fprintf(out,"DATASET STRUCTURED_GRID\n");

 // ???????????????????? ?????????? ??????????:

   fprintf(out,"DIMENSIONS %5d %5d %5d\n",nx,ny,nz);
   fprintf(out,"POINTS %7d float\n", npoin);

   for(i = 0; i < nx; i++)
   {
      for(j = 0; j < ny; j++)
      {
         for(k = 0; k < nz; k++)
         {

            x = i*dx;
            y = j*dy;
            z = k*dz;

            fprintf(out, "%14.6e %14.6e %14.6e\n",x,y,z);
         }
      }
   }

 // ???????????? ???????????? ???? ??????????:

   fprintf(out,"POINT_DATA %5d\n", npoin);

   fprintf(out,"SCALARS CON float 1\n");

   fprintf(out,"LOOKUP_TABLE default\n");

   for(i = 0; i < nx; i++)
   {
      for(j = 0; j < ny; j++)
      {
         for(k = 0; k < nz; k++)
         {


            if(fabs(data1[i][j][k]) < 1.e-15)
               fprintf(out,"%14.6e\n",0.0);
            else   
               fprintf(out,"%14.6e\n",data1[i][j][k]);
         }
      }
   }

   fclose(out);
}

void copy_3d(double ***src, double ***dest, int n){
	int i, j, k;
	for(i = 0; i< n; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < n; k++){
                		dest[i][j][k] = src[i][j][k];
			}
		}
	}
}


void write_into_file(double ***a, int timestep, FILE *out)
{	
	fprintf(out, "%d\n", timestep);
	int i, j, k;
	for (i = 0; i < N; i++)
        {
		fprintf(out, "\n");
                for (j = 0; j < N; j++)
                {
			fprintf(out, "\n");
                        for (k = 0; k < N; k++)
                        {
                                fprintf(out, "%lf' '", a[i][j][k]);
                        }
                }
        }
}
