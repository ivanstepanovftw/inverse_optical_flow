#include <iostream>
#include <ctime>

using namespace std;

#include "backward_flow.h"

extern "C" {
#include "iio.h"
}

bool read_image(const char *fname, float **Ir, float **Ig, float **Ib, int &nx, int &ny, int &nz)
{
	
	float *f;
	
	*Ir = *Ig = *Ib = NULL;
	
	f = iio_read_image_float_vec(fname, &nx, &ny, &nz);
	
	if(nz > 0)
	{
	    *Ir = new float[nx * ny];
	    *Ig = new float[nx * ny];
	    *Ib = new float[nx * ny];
	    
	    int dg=1, db=2;
	    
	    if(nz == 1) dg=db=0;
	    
	    for(int i = 0; i < ny; i++)
		for(int j = 0; j < nx; j++)
		{
		      (*Ir)[i*nx+j] = f[(i*nx+j)*nz];
		      (*Ig)[i*nx+j] = f[(i*nx+j)*nz+dg];
		      (*Ib)[i*nx+j] = f[(i*nx+j)*nz+db];
		}
	}
	
	return *Ir ? true : false;
}


bool read_flow(const char *fname, float **u, float **v, int &nx, int &ny)
{
	
	float *f;
	int nz;
	
	*u = *v = NULL;
	
	f = iio_read_image_float_vec(fname, &nx, &ny, &nz);
	
	if(nz > 0)
	{
	    *u = new float[nx * ny];
	    *v = new float[nx * ny];
    
	    int dg=1, db=2;
	    
	    if(nz == 1) dg=db=0;
	    
	    for(int i = 0; i < ny; i++)
		for(int j = 0; j < nx; j++)
		{
		      (*u)[i*nx+j] = f[(i*nx+j)*2];
		      (*v)[i*nx+j] = f[(i*nx+j)*2+1];
		}
	}
	
	return *u ? true : false;
}

void save_flow(const char *fname, float *u, float *v, int nx, int ny)
{
	//save the flow 
	float *f = new float[nx * ny * 2];
	for (int i = 0; i < nx * ny; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	}
	iio_save_image_float_vec(fname, f, nx, ny, 2);
	delete []f;
}


int main(int argc, char *argv[])
{
	if(argc < 4)
		cout << "Usage: " << argv[0] << " I1 I2 flow_in [flow_out mask fill strategy]" << endl;
	else
	{
		int nx, ny, nz;

		int i = 1;

		const char *image1   = argv[i]; i++;
		const char *image2   = argv[i]; i++;
		const char *flow_in  = argv[i]; i++;
		const char *flow_out = (argc > i)? argv[i]: "inverse_flow.uv"; i++;
		const char *mask_out = (argc > i)? argv[i]: NULL; i++;

		const int   strategy = (argc > i)? atoi(argv[i]): AVG_IMAGE_METHOD; i++;
		const int   fill     = (argc > i)? atoi(argv[i]): ORIENTED_FILL; i++;
		const int   verbose  = (argc > i)? atoi(argv[i]): 0; i++;
		
		float *I1r, *I1g, *I1b, *I2r, *I2g, *I2b;
		
		if(read_image(image1, &I1r, &I1g, &I1b, nx, ny, nz) &&
		   read_image(image2, &I2r, &I2g, &I2b, nx, ny, nz))
		{
		    float *u, *v, *u_, *v_, *m;

		    read_flow(flow_in, &u, &v, nx, ny);
		    
		    u_ = new float[nx * ny];
		    v_ = new float[nx * ny];
		    m  = new float[nx * ny];
		    
		    if(verbose)
		      cout << "strategy = " << strategy << " fill = " << fill << endl;

		    float diff;
		    struct timespec start, end;
		    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		    backward_flow(I1r, I1g, I1b, I2r, I2g, I2b, u, v, u_, v_, m, nx, ny, strategy, fill);
		    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		    diff = (float) (end.tv_nsec - start.tv_nsec)/1000000;
		    cout.precision(8);
		    cout << "Time: " << diff << endl;
		    
		    save_flow(flow_out, u_, v_, nx, ny);
		    
		    if(mask_out) save_flow(mask_out, m, m, nx, ny);

		    delete []I1r;
		    delete []I1g;
		    delete []I1b;
		    delete []I2r;
		    delete []I2g;
		    delete []I2b;
		    delete []u;
		    delete []v;
		    delete []u_;
		    delete []v_;
		    delete []m;
		}
	}
}
