#ifndef _INVERSE_OPTIC_FLOW_H
#define _INVERSE_OPTIC_FLOW_H

#include "fill_disocclusions.h"
#include <cfloat>

//constants definition for inverse optical flow algorithms
#define MAX_FLOW_METHOD  1   
#define MAX_IMAGE_METHOD 2
#define AVG_FLOW_METHOD  3
#define AVG_IMAGE_METHOD 4  

#define OCCLUSION 99999.0  //FLT_MAX
#define OCCLUSION_MAX_FLOW 0
#define WEIGHT_TH 0.25
#define MOTION_TH 0.25

/**
 * 
 *   Function to compute the backward flow from the forward flow
 * 
 */
void inverse_flow(
    const float *u, 
    const float *v, 
    float       *u_, 
    float       *v_,
    float       *mask,
    const int    nx, 
    const int    ny 
)
{
      for(int y = 0; y < ny; y++)

	for(int x = 0; x < nx; x++)
	{
	    //warpeamos el flujo
	    const int   pos = x + nx * y;
	    const float xw = (float) (x + u[pos]);
	    const float yw = (float) (y + v[pos]);

	    const int sx = (xw < 0)? -1: 1;
	    const int sy = (yw < 0)? -1: 1;

	    int xi  = (int) xw;
	    int yi  = (int) yw;

	    int dx = xi + sx;
	    int dy = yi + sy;

	    if(xi < 0) xi = 0;
	    else if (xi >= nx) xi = nx - 1;
	    if(yi < 0) yi = 0;
	    else if (yi >= ny) yi = ny - 1;

	    if(dx < 0) dx = 0;
	    else if (dx >= nx) dx = nx - 1;
	    if(dy < 0) dy = 0;
	    else if (dy >= ny) dy = ny - 1;

	    const int pos1 = xi + nx * yi;
	    const int pos2 = dx + nx * yi;
	    const int pos3 = xi + nx * dy;
	    const int pos4 = dx + nx * dy;

	    const float e1 = ((float) sx * (xw - xi));
	    const float E1 = ((float) 1.0 - e1);
	    const float e2 = ((float) sy * (yw - yi));
	    const float E2 = ((float) 1.0 - e2);

	    //metemos en los cuatro puntos la proporcion correspondiente
	    const float w1 = E1 * E2;
	    const float w2 = e1 * E2;
	    const float w3 = E1 * e2;
	    const float w4 = e1 * e2;

	    const float d  = u[pos] * u[pos] + v[pos] * v[pos];
	    const float d1 = u_[pos1] * u_[pos1] + v_[pos1] * v_[pos1];
	    const float d2 = u_[pos2] * u_[pos2] + v_[pos2] * v_[pos2];
	    const float d3 = u_[pos3] * u_[pos3] + v_[pos3] * v_[pos3];
	    const float d4 = u_[pos4] * u_[pos4] + v_[pos4] * v_[pos4];
	    
	    if(w1 >= WEIGHT_TH && d >= d1)
	    {
		u_[pos1] = -u[pos];
		v_[pos1] = -v[pos];
		mask[pos1] = NO_DISOCCLUSION;
	    }

	    if(w2 >= WEIGHT_TH && d >= d2)
	    {
		u_[pos2] = -u[pos];
		v_[pos2] = -v[pos];
		mask[pos2] = NO_DISOCCLUSION;
	    }

	    if(w3 >= WEIGHT_TH && d >= d3)
	    {
		u_[pos3] = -u[pos];
		v_[pos3] = -v[pos];
		mask[pos3] = NO_DISOCCLUSION;
	    }

	    if(w4 >= WEIGHT_TH && d >= d4)
	    {
		u_[pos4] = -u[pos];
		v_[pos4] = -v[pos];
		mask[pos4] = NO_DISOCCLUSION;
	    }
	} 
}


/**
 * 
 *   Function to compute the backward flow from the forward flow
 * 
 */
void inverse_image_max_flow(
    const float *I1r,
    const float *I1g,
    const float *I1b,
    const float *I2r,
    const float *I2g,
    const float *I2b,
    const float *u, 
    const float *v, 
    float       *u_, 
    float       *v_,
    float       *mask,
    const int    nx, 
    const int    ny 
)
{
      int size = nx * ny;
      
      float *DI = new float[size];

      for(int i = 0; i < size; i++) DI[i] = FLT_MAX;
      
      for(int y = 0; y < ny; y++)

	for(int x = 0; x < nx; x++)
	{
	    //warpeamos el flujo
	    const int   pos = x + nx * y;
	    const float xw  = (float) (x + u[pos]);
	    const float yw  = (float) (y + v[pos]);

	    const int sx = (xw < 0)? -1: 1;
	    const int sy = (yw < 0)? -1: 1;

	    int xi = (int) xw;
	    int yi = (int) yw;

	    int dx = xi + sx;
	    int dy = yi + sy;

	    if(xi < 0) xi = 0;
	    else if (xi >= nx) xi = nx - 1;
	    if(yi < 0) yi = 0;
	    else if (yi >= ny) yi = ny - 1;

	    if(dx < 0) dx = 0;
	    else if (dx >= nx) dx = nx - 1;
	    if(dy < 0) dy = 0;
	    else if (dy >= ny) dy = ny - 1;

	    const int pos1 = xi + nx * yi;
	    const int pos2 = dx + nx * yi;
	    const int pos3 = xi + nx * dy;
	    const int pos4 = dx + nx * dy;

	    const float e1 = ((float) sx * (xw - xi));
	    const float E1 = ((float) 1.0 - e1);
	    const float e2 = ((float) sy * (yw - yi));
	    const float E2 = ((float) 1.0 - e2);

	    //metemos en los cuatro puntos la proporcion correspondiente
	    const float w1 = E1 * E2;
	    const float w2 = e1 * E2;
	    const float w3 = E1 * e2;
	    const float w4 = e1 * e2;

	    const float d1 = (I1r[pos] - I2r[pos1]) * (I1r[pos] - I2r[pos1]) + 
			     (I1g[pos] - I2g[pos1]) * (I1g[pos] - I2g[pos1]) + 
			     (I1b[pos] - I2b[pos1]) * (I1b[pos] - I2b[pos1]);
	    const float d2 = (I1r[pos] - I2r[pos2]) * (I1r[pos] - I2r[pos2]) + 
			     (I1g[pos] - I2g[pos2]) * (I1g[pos] - I2g[pos2]) + 
			     (I1b[pos] - I2b[pos2]) * (I1b[pos] - I2b[pos2]);
	    const float d3 = (I1r[pos] - I2r[pos3]) * (I1r[pos] - I2r[pos3]) + 
			     (I1g[pos] - I2g[pos3]) * (I1g[pos] - I2g[pos3]) + 
			     (I1b[pos] - I2b[pos3]) * (I1b[pos] - I2b[pos3]);
	    const float d4 = (I1r[pos] - I2r[pos4]) * (I1r[pos] - I2r[pos4]) + 
			     (I1g[pos] - I2g[pos4]) * (I1g[pos] - I2g[pos4]) + 
			     (I1b[pos] - I2b[pos4]) * (I1b[pos] - I2b[pos4]);

	    if(w1 >= WEIGHT_TH && DI[pos1] >= d1)
	    {
		u_[pos1] = -u[pos];
		v_[pos1] = -v[pos];
		DI[pos1] = d1;
		mask[pos1] = NO_DISOCCLUSION;
	    }

	    if(w2 >= WEIGHT_TH && DI[pos2] >= d2)
	    {
		u_[pos2] = -u[pos];
		v_[pos2] = -v[pos];
		DI[pos2] = d2;
		mask[pos2] = NO_DISOCCLUSION;
	    }

	    if(w3 >= WEIGHT_TH && DI[pos3] >= d3)
	    {
		u_[pos3] = -u[pos];
		v_[pos3] = -v[pos];
		DI[pos3] = d3;
		mask[pos3] = NO_DISOCCLUSION;
	    }

	    if(w4 >= WEIGHT_TH && DI[pos4] >= d4)
	    {
		u_[pos4] = -u[pos];
		v_[pos4] = -v[pos];
		DI[pos4] = d4;
		mask[pos4] = NO_DISOCCLUSION;
	    }
	} 
}


inline void select_motion(
    const float d,
    const float u, 
    const float v,
    const float wght,
    float &d_,
    float &u_,
    float &v_,
    float &wght_,
    float &mask
)
{
    if(wght >= WEIGHT_TH)
    {
	if(fabs(d-d_) <= MOTION_TH)
	{
	    u_    += u * wght;
	    v_    += v * wght;
	    wght_ += wght;
	    mask   = NO_DISOCCLUSION;
	}
	else if(d >= d_)   
	{
	    //if it is an oclussion we retain the highest value
	    d_    = d;
	    u_    = u * wght;
	    v_    = v * wght;
	    wght_ = wght;
	    mask  = NO_DISOCCLUSION;
	}
    }
}


/**
 *
 *   Function to compute the backward flow from the forward flow
 *
 */
void inverse_average_flow(
    const float *u, 
    const float *v, 
    float       *u_, 
    float       *v_,
    float       *mask,
    const int    nx, 
    const int    ny 
)
{
      int size = nx * ny;
      
      float *avg_u = new float[size];
      float *avg_v = new float[size];

      float *wgt_ = new float[size];

      float *d_ = new float[size];

      for(int i = 0; i < size; i++)
      {
	avg_u[i] = avg_v[i] = wgt_[i] = 0;
	d_[i] = -999;
      }
      
      for(int y = 0; y < ny; y++)

	for(int x = 0; x < nx; x++)
	{
	    //warpeamos el flujo
	    const int   pos = x + nx * y;
	    const float xw = (float) (x + u[pos]);
	    const float yw = (float) (y + v[pos]);

	    const int sx = (xw < 0)? -1: 1;
	    const int sy = (yw < 0)? -1: 1;

	    int xi = (int) xw;
	    int yi = (int) yw;

	    int dx = xi + sx;
	    int dy = yi + sy;

	    if(xi < 0) xi = 0;
	    else if (xi >= nx) xi = nx - 1;
	    if(yi < 0) yi = 0;
	    else if (yi >= ny) yi = ny - 1;

	    if(dx < 0) dx = 0;
	    else if (dx >= nx) dx = nx - 1;
	    if(dy < 0) dy = 0;
	    else if (dy >= ny) dy = ny - 1;

	    const int pos1 = xi + nx * yi;
	    const int pos2 = dx + nx * yi;
	    const int pos3 = xi + nx * dy;
	    const int pos4 = dx + nx * dy;

	    const float e1 = ((float) sx * (xw - xi));
	    const float E1 = ((float) 1.0 - e1);
	    const float e2 = ((float) sy * (yw - yi));
	    const float E2 = ((float) 1.0 - e2);

	    //metemos en los cuatro puntos la proporcion correspondiente
	    const float w1 = E1 * E2;
	    const float w2 = e1 * E2;
	    const float w3 = E1 * e2;
	    const float w4 = e1 * e2;

	    const float d = u[pos] * u[pos] + v[pos] * v[pos];
	    
	    select_motion(
	      d, u[pos], v[pos], w1, d_[pos1], avg_u[pos1], 
	      avg_v[pos1], wgt_[pos1], mask[pos1]
	    );

	    select_motion(
	      d, u[pos], v[pos], w2, d_[pos2], avg_u[pos2], 
	      avg_v[pos2], wgt_[pos2], mask[pos2]
	    );

	    select_motion(
	      d, u[pos], v[pos], w3, d_[pos3], avg_u[pos3], 
	      avg_v[pos3], wgt_[pos3], mask[pos3]
	    );

	    select_motion(
	      d, u[pos], v[pos], w4, d_[pos4], avg_u[pos4], 
	      avg_v[pos4], wgt_[pos4], mask[pos4]
	    );   
	}

	for(int i = 0; i < size; i++) 
	{	    
	  if(mask[i] == NO_DISOCCLUSION)
	  {
	    u_[i] = -avg_u[i] / wgt_[i];
	    v_[i] = -avg_v[i] / wgt_[i];
	  }
	}

	delete []avg_u;
	delete []avg_v;
	delete []wgt_;
	delete []d_;
}

    
inline void select_image_motion(
    const float dI,
    const float u, 
    const float v,
    const float wght,
    float &d_,
    float &dI_,
    float &u_,
    float &v_,
    float &wght_,
    float &mask
)
{
    if(wght >= WEIGHT_TH)
    {
	const float d = u * u  + v * v ;
	
	if(fabs(d-d_) <= MOTION_TH)
	{
	    u_    += u * wght;
	    v_    += v * wght;
	    wght_ += wght;
	    mask   = NO_DISOCCLUSION;
	}
	else if(dI_ >= dI) 
	//if it is an oclussion we retain the value with the most similar image colors
	{
	    d_    = d;
	    dI_   = dI;
	    u_    = u * wght;
	    v_    = v * wght;
	    wght_ = wght;
	    mask  = NO_DISOCCLUSION;
	}
    }
}

/**
 * 
 *   Function to compute the backward flow from the forward flow
 * 
 */
void inverse_image_average_flow(
    const float *I1r,
    const float *I1g,
    const float *I1b,
    const float *I2r,
    const float *I2g,
    const float *I2b,
    const float *u, 
    const float *v, 
    float       *u_, 
    float       *v_,
    float       *mask,
    const int    nx, 
    const int    ny 
)
{
      int size = nx * ny;
      
      float *avg_u = new float[size];
      float *avg_v = new float[size];
      float *wgt_  = new float[size];
      float *d_    = new float[size];
      float *dI    = new float[size];

      for(int i = 0; i < size; i++)
      {
	avg_u[i] = avg_v[i] = wgt_[i] = 0;
	d_[i] = -999;
	dI[i] = FLT_MAX;
      }
      
      for(int y = 0; y < ny; y++)

	for(int x = 0; x < nx; x++)
	{
	    //warpeamos el flujo
	    const int   pos = x + nx * y;
	    const float xw  = (float) (x + u[pos]);
	    const float yw  = (float) (y + v[pos]);

	    const int sx = (xw < 0)? -1: 1;
	    const int sy = (yw < 0)? -1: 1;

	    int xi = (int) xw;
	    int yi = (int) yw;

	    int dx = xi + sx;
	    int dy = yi + sy;

	    if(xi < 0) xi = 0;
	    else if (xi >= nx) xi = nx - 1;
	    if(yi < 0) yi = 0;
	    else if (yi >= ny) yi = ny - 1;

	    if(dx < 0) dx = 0;
	    else if (dx >= nx) dx = nx - 1;
	    if(dy < 0) dy = 0;
	    else if (dy >= ny) dy = ny - 1;

	    const int pos1 = xi + nx * yi;
	    const int pos2 = dx + nx * yi;
	    const int pos3 = xi + nx * dy;
	    const int pos4 = dx + nx * dy;

	    const float e1 = ((float) sx * (xw - xi));
	    const float E1 = ((float) 1.0 - e1);
	    const float e2 = ((float) sy * (yw - yi));
	    const float E2 = ((float) 1.0 - e2);

	    //metemos en los cuatro puntos la proporcion correspondiente
	    const float w1 = E1 * E2;
	    const float w2 = e1 * E2;
	    const float w3 = E1 * e2;
	    const float w4 = e1 * e2;

	    const float dI1 = (I1r[pos] - I2r[pos1]) * (I1r[pos] - I2r[pos1]) + 
			     (I1g[pos] - I2g[pos1]) * (I1g[pos] - I2g[pos1]) + 
			     (I1b[pos] - I2b[pos1]) * (I1b[pos] - I2b[pos1]);
	    const float dI2 = (I1r[pos] - I2r[pos2]) * (I1r[pos] - I2r[pos2]) + 
			     (I1g[pos] - I2g[pos2]) * (I1g[pos] - I2g[pos2]) + 
			     (I1b[pos] - I2b[pos2]) * (I1b[pos] - I2b[pos2]);
	    const float dI3 = (I1r[pos] - I2r[pos3]) * (I1r[pos] - I2r[pos3]) + 
			     (I1g[pos] - I2g[pos3]) * (I1g[pos] - I2g[pos3]) + 
			     (I1b[pos] - I2b[pos3]) * (I1b[pos] - I2b[pos3]);
	    const float dI4 = (I1r[pos] - I2r[pos4]) * (I1r[pos] - I2r[pos4]) + 
			     (I1g[pos] - I2g[pos4]) * (I1g[pos] - I2g[pos4]) + 
			     (I1b[pos] - I2b[pos4]) * (I1b[pos] - I2b[pos4]);

	    select_image_motion(
	      dI1, u[pos], v[pos], w1, d_[pos1], dI[pos1], avg_u[pos1], 
	      avg_v[pos1], wgt_[pos1], mask[pos1]
	    );

	    select_image_motion(
	      dI2, u[pos], v[pos], w2, d_[pos2], dI[pos2], avg_u[pos2], 
	      avg_v[pos2], wgt_[pos2], mask[pos2]
	    );

	    select_image_motion(
	      dI3, u[pos], v[pos], w3, d_[pos3], dI[pos3], avg_u[pos3], 
	      avg_v[pos3], wgt_[pos3], mask[pos3]
	    );

	    select_image_motion(
	      dI4, u[pos], v[pos], w4, d_[pos4], dI[pos4], avg_u[pos4], 
	      avg_v[pos4], wgt_[pos4], mask[pos4]
	    );   
	}

	for(int i = 0; i < size; i++) 
	{	    
	  if(mask[i] == NO_DISOCCLUSION)
	  {
	    u_[i] = -avg_u[i] / wgt_[i];
	    v_[i] = -avg_v[i] / wgt_[i];
	  }
	}

	delete []avg_u;
	delete []avg_v;
	delete []wgt_;
	delete []d_;
	delete []dI;
}



/**
 * 
 *   Function to compute the backward flow from the forward flow
 * 
 */
int backward_flow(
    const float *I1r,
    const float *I1g,
    const float *I1b,
    const float *I2r,
    const float *I2g,
    const float *I2b,
    const float *u, 
    const float *v, 
    float       *u_, 
    float       *v_,
    float       *mask,
    int 	 	nx, 
    int 	 	ny,
    int		strategy,
    int		fill
)
{
    int size = nx * ny;

    float occlusion = (strategy == MAX_FLOW_METHOD)? OCCLUSION_MAX_FLOW: OCCLUSION;

    for(int i = 0; i < size; i++)
    {
	u_[i] = v_[i] = occlusion;
	mask[i] = DISOCCLUSION;
    }

    switch(strategy) {
      case MAX_FLOW_METHOD:  
	      inverse_flow(u, v, u_, v_, mask, nx, ny);
	      break;
	      
      case MAX_IMAGE_METHOD: 
	      inverse_image_max_flow(I1r, I1g, I1b, I2r, I2g, I2b, u, v, u_, v_, mask, nx, ny);
	      break;
	      
      case AVG_FLOW_METHOD: 
	      inverse_average_flow(u, v, u_, v_, mask, nx, ny);
	      break;
	      
      case AVG_IMAGE_METHOD: default:
	      inverse_image_average_flow(I1r, I1g, I1b, I2r, I2g, I2b, u, v, u_, v_, mask, nx, ny);
	      break;
    }
    
    if(fill == MIN_FILL) 
      restricted_minfill(u_, v_, mask, nx, ny);
    else if(fill == AVERAGE_FILL)
      average_fill(u_, v_, mask, nx, ny);
    else if(fill == ORIENTED_FILL)
      oriented_fill(u, v, u_, v_, mask, nx, ny);
    else if(strategy==MAX_FLOW_METHOD)
      for(int i = 0; i < size; i++)
	if(mask[i] == DISOCCLUSION)
	  u_[i] = v_[i] = OCCLUSION;

    return 0;
}

#endif
