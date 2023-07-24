#ifndef FILL_DISOCCLUSIONS
#define FILL_DISOCCLUSIONS

#include <cmath>
#include <vector>
#include <algorithm>

#define MIN_FILL  1
#define AVERAGE_FILL 2
#define ORIENTED_FILL 3

#define NO_DISOCCLUSION 1
#define DISOCCLUSION 0
#define STEREO_DISOCCLUSION -1
#define STREETLAMP_DISOCCLUSION -2


//Two pass algorithm for labelling regions - connected-component labelling algorithm
//returns the number of layers found
int connected_component_labeling(
    float *mask,
    int   *regions,
    int    nx, 
    int    ny 
)
{
    int label = 1;
    
    std::vector<int> corresponding_label;
    
    corresponding_label.push_back(0);
    
    for(int i = 0; i < nx * ny; i++) regions[i] = 0;
        
    //first pass for initial labelling
    for(int i = 0; i < ny; i++)

	for(int j = 0; j < nx; j++)
	{
	    int p = i * nx + j;

	    if(mask[p] == DISOCCLUSION)
	    {
		int l1 = 0;
		int l2 = 0;

		//take neighbors labels
		if(i-1 >= 0) l1 = regions[p - nx];
		if(j-1 >= 0) l2 = regions[p - 1];

		if(l1 > 0 || l2 > 0)
		{
		    if(l1 > 0)
		    {
			if(l2 > 0)
			{
			    int c = std::min(corresponding_label[l1], corresponding_label[l2]); 
			 
			    corresponding_label[l1] = c;
			    corresponding_label[l2] = c;
			 
			    regions[p] = c;
			}
			else
			    regions[p] = l1;
		    }
		    else 
			regions[p] = l2;
		}
		else
		{
		    regions[p] = label;    
		    corresponding_label.push_back(label);
		    label++;
		}
	    }
	}
	
    //unify labels
    for(int i = 2; i < (int) corresponding_label.size(); i++)
    {
	int c = corresponding_label[i];
	
	while(c != corresponding_label[c]) 
	    c = corresponding_label[c];
	
	corresponding_label[i] = c;
    }

    //second pass for simplifying labels
    for(int i = 0; i < ny * nx; i++)
	regions[i] = corresponding_label[regions[i]];
    
    return *(std::max_element(corresponding_label.begin(), corresponding_label.end()));
}

void save_regions(int *regions, int nx, int ny, int labels)
{
    float *r1 = new float[nx * ny];
    float *r2 = new float[nx * ny];

    std::vector<float> ang;
    ang.push_back(0);
    
    for(int i = 1; i < labels; i++) ang.push_back((float)rand()/RAND_MAX);

	for(int i = 0; i < nx * ny; i++)
	{
	    r1[i] = (regions[i]==0)?0:sin(2*3.141592*regions[i]/labels);
	    r2[i] = (regions[i]==0)?0:cos(2*3.141592*ang[regions[i]]);
	}

    delete []r1;
    delete []r2;
}


inline void min_test(float *mask, float *u, float *v, float *mu, float *mv, float *d, int i, int l)
{
    if(mask[i] != DISOCCLUSION)
    {
	const float uu = u[i];
	const float vv = v[i];
	const float dd = uu * uu + vv * vv;
	
	if(dd < d[l])
	{
	    d[l]  = dd;
	    mu[l] = uu;
	    mv[l] = vv;
	}
    }   
}

/**
 * 
 *   Function to fill the empty regions with the min neighbor value
 * 
 */
void minfill(
    float *u, 
    float *v, 
    float *mask,
    int    nx, 
    int    ny 
)
{
    int *regions = new int[nx * ny];

    //classify the regions and assign labels
    const int labels = connected_component_labeling(mask, regions, nx, ny);

    save_regions(regions, nx, ny, labels);
    
    //search for the minimum around every region
    float *min_u = new float[labels];
    float *min_v = new float[labels];
    float *min_d = new float[labels];
    
    for(int i = 0; i < labels; i++)
	min_u[i] = min_v[i] = min_d[i] = 99999.9;
    
    for(int i = 0; i < ny; i++)
	
	for(int j = 0; j < nx; j++)
	{
	      int p = i * nx + j;
	      
	      if(mask[p] == DISOCCLUSION)
	      {
		  const int l = regions[p]-1;
		  if(i >    0) min_test(mask, u, v, min_u, min_v, min_d, p-nx, l);
		  if(i < ny-1) min_test(mask, u, v, min_u, min_v, min_d, p+nx, l);
		  if(j >    0) min_test(mask, u, v, min_u, min_v, min_d, p-1,  l);
		  if(j < nx-1) min_test(mask, u, v, min_u, min_v, min_d, p+1,  l);
	      }
	}
	
    //fill the holes with the minimum values around the region
    for(int i = 0; i < nx * ny; i++)

	if(mask[i] == DISOCCLUSION)
	{
	    const int l = regions[i]-1;
	    
	    u[i] = min_u[l];
	    v[i] = min_v[l];
	}
	    
    delete []min_u;
    delete []min_v;
    delete []min_d;
    delete []regions;
}




/**
 * 
 *   Function to fill the empty regions with the closest min value
 * 
 */
void restricted_minfill(
    float *u, 
    float *v, 
    float *mask,
    int    nx, 
    int    ny,
    int    radius = 5
)
{
    float *mask2 = new float[nx * ny];
    float *mask3 = new float[nx * ny];

    for(int i = 0; i < nx * ny; i++) mask2[i] = mask3[i] = mask[i];

    bool end;

    do 
    {
	end = true;

	for(int i = 0; i < ny; i++)
	    for(int j = 0; j < nx; j++)
	    {
		const int p = i * nx + j;
		//if there is a disocclusion, then try to fill it up
		if(mask2[p] == DISOCCLUSION)
		{
		    float min_d = 99999.9;
		    float min_u = 99999.9;
		    float min_v = 99999.9;

		    bool min_found = false;

		    for(int k = ((i-radius>0)? i-radius: 0); k < ((i+radius<ny)? i+radius: ny-1); k++)
			for(int l = ((j-radius>0)? j-radius: 0); l < ((j+radius<nx)? j+radius: nx-1); l++)
			{
			    const int q = k * nx + l;
			    
			    if(mask2[q] != DISOCCLUSION)
			    {
				const int d = u[q] * u[q] + v[q] * v[q];
				
				if(d < min_d)
				{
				    min_d = d;
				    min_u = u[q];
				    min_v = v[q];
				    
				    min_found = true;
				}
			    }
			}

		    if(min_found)
		    {
			u[p] = min_u;
			v[p] = min_v;
			mask3[p] = NO_DISOCCLUSION;
		    }
		    else end = false;
		}
	    }

	std::swap(mask2, mask3);

    }while(!end);

    delete []mask2;
    delete []mask3;
}



/**
 * 
 *   Function to fill the empty regions with the closest min value
 * 
 */
void average_fill(
    float *u, 
    float *v, 
    float *mask,
    int    nx, 
    int    ny,
    int    radius = 5
)
{
    float *mask2 = new float[nx * ny];
    float *mask3 = new float[nx * ny];

    for(int i = 0; i < nx * ny; i++) mask2[i] = mask3[i] = mask[i];

    bool end;

    do 
    {
	end = true;

	for(int i = 0; i < ny; i++)
	    for(int j = 0; j < nx; j++)
	    {
		const int p = i * nx + j;
		//if there is a disocclusion, then try to fill it up
		if(mask2[p] != NO_DISOCCLUSION)
		{
		    float avg_u = 0.0;
		    float avg_v = 0.0;
		    int n = 0;

		    for(int k = ((i-radius>0)? i-radius: 0); k < ((i+radius<ny)? i+radius: ny-1); k++)
			for(int l = ((j-radius>0)? j-radius: 0); l < ((j+radius<nx)? j+radius: nx-1); l++)
			{
			    const int q = k * nx + l;

			    if(mask2[q] != DISOCCLUSION)
			    {
				avg_u += u[q]; 
				avg_v += v[q]; 
				n++;
			    }
			}

		    if(n >= radius)
		    {
			u[p] = avg_u/n;
			v[p] = avg_v/n;
			mask3[p] = NO_DISOCCLUSION;
		    }
		    else end = false;
		}
	    }

	std::swap(mask2, mask3);

    }while(!end);

    delete []mask2;
    delete []mask3;
}



/**
 *
 *   Function to fill the empty regions with the closest min value
 *
 */
void oriented_fill(
    const float *u,
    const float *v,
    float *u_,
    float *v_,
    float *mask,
    int    nx,
    int    ny
)
{
    const int size = nx * ny;

    int   *index = new int[size];
    
    float *du  = new float[size];
    float *dv  = new float[size];
    float *duu = new float[size];
    float *dvv = new float[size];

    //compute filling direction
    for(int i = 0; i < size; i++)
    {
	//in case of disocclusion
	if(mask[i] != NO_DISOCCLUSION)
	{
	    //normalize direction
	    const float d = sqrt(u[i] * u[i] + v[i] * v[i]);
	    duu[i] = du[i] = -u[i] / d;
	    dvv[i] = dv[i] = -v[i] / d;
	    index[i] = -2;
	}
    }

    bool end;
    do
    {	
	end = true;

	for(int i = 0; i < ny; i++)

	    for(int j = 0; j < nx; j++)
	    {
		const int p = i * nx + j;

		//if there is a disocclusion and no value assigned, then try to fill it 
		if(mask[p] != NO_DISOCCLUSION && index[p] < 0)
		{
		    int k = (int)((float)i + dv[p] + 0.5);
		    int l = (int)((float)j + du[p] + 0.5);

		    if(k < 0 || k >= ny) {
		      dvv[p] = dv[p] = -dvv[p];
		      du[p]  = duu[p];
		      k = i + dv[p] + 0.5;
		    }
		    if(l < 0 || l >= nx) {
		      duu[p] = du[p] = -duu[p];
		      dv[p]  = dvv[p];
		      l = j + du[p] + 0.5;
		    }

		    const int p1 = k * nx + l;
		    if(mask[p1] == NO_DISOCCLUSION)
			index[p] = p1;
		    else 
		    {
			end = false;

			//test the direction of both disocclusions
			const float d  = sqrt(u[p]  * u[p]  + v[p]  * v[p]);
			const float d1 = sqrt(u[p1] * u[p1] + v[p1] * v[p1]);
			const float uv = u[p] * u[p1] + v[p] * v[p1];
			if(uv / (d * d1) < 0.9 && index[p] != -1) {
			    if(d1 > d) {
				duu[p] = du[p] = -u[p1] / d1;
				dvv[p] = dv[p] = -v[p1] / d1;
				index[p] = -1;
			    }
			}
			du[p] = du[p] + duu[p];
			dv[p] = dv[p] + dvv[p];
		    }
		}
	    }
    } while(!end);
    
    //we fill the inverse flow
    for(int i = 0; i < size; i++)

      if(mask[i] != NO_DISOCCLUSION)
      {
	  u_[i] = u_[index[i]];
	  v_[i] = v_[index[i]];
      }    

    delete []index;
    delete []du;
    delete []dv;
    delete []duu;
    delete []dvv;
}


#endif
