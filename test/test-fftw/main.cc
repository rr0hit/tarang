

#include <iostream>
#include <cmath>
#include <complex>
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#include <fftw3-mpi.h>

using namespace std;
     
int main(int argc, char **argv)
{
         const ptrdiff_t N0 = 8, N1 = 8;
         fftw_plan plan;
         fftw_complex *data;
         ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
     
         MPI_Init(&argc, &argv);
         fftw_mpi_init();
     
         /* get local data size and allocate */
         alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                              &local_n0, &local_0_start);
         data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
     
         /* create plan for forward DFT */
         plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
     
         /* initialize data to some function my_function(x,y) */
       //  for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j) {
          //  data[i*N1 + j] = my_function(local_0_start + i, j);
		//	real(data[i*N1 + j]) =1.0 ; imag(data[i*N1 + j]) =1.0 ;
	//	}
     
         /* compute transforms, in-place, as many times as desired */
	
	real(data[1,1]) = 1.0;
         fftw_execute(plan);
     
         fftw_destroy_plan(plan);
     
		cout << " DONE" << endl;
         MPI_Finalize();
     }




