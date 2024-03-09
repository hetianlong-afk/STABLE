#include "mex.h"
#include <math.h>
// input (TbAng_coef_real,TbAng_coef_imag,V_load_0_real,V_load_0_imag, V_load_real,V_load_imag) 
// output (V_load_real,V_load_imag,V_load_0_real,V_load_0_imag)
// 

double complex_multiply_Real(double a1, double b1, double a2, double b2)
{
	return a1*a2-b1*b2;
}

double complex_multiply_Imag(double a1, double b1, double a2, double b2)
{
	return a1*b2+a2*b1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *V_load_real_in, *V_load_imag_in, *pattern;
	double V_load_0_real, V_load_0_imag, TbAng_coef_real, TbAng_coef_imag;
	size_t V_num, P_num;
	double *V_load_real_out, *V_load_imag_out;
	double *V_load_0_real_out, *V_load_0_imag_out;
	int i;
	int j;
	// input data
	TbAng_coef_real = mxGetScalar(prhs[0]);
	TbAng_coef_imag = mxGetScalar(prhs[1]);
	V_load_0_real = mxGetScalar(prhs[2]);
	V_load_0_imag = mxGetScalar(prhs[3]);
	
	V_load_real_in = mxGetPr(prhs[4]);
	V_load_imag_in = mxGetPr(prhs[5]);
	
	V_num = (int)(mxGetM(prhs[4])*mxGetN(prhs[4]));
	
	pattern        = mxGetPr(prhs[6]);
	
	P_num = (int)(mxGetM(prhs[6])*mxGetN(prhs[6]));
	
	// output data
	plhs[0] = mxCreateDoubleMatrix(1,V_num,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,V_num,mxREAL);	
	V_load_real_out = mxGetPr(plhs[0]);
	V_load_imag_out = mxGetPr(plhs[1]);
	
	plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);	
	V_load_0_real_out = mxGetPr(plhs[2]);
	V_load_0_imag_out = mxGetPr(plhs[3]);
	j=0;
	for (i=0;i<P_num;i++)
	{
		
		if (pattern[i] ==1)
		{
			
			V_load_real_out[j] = complex_multiply_Real(TbAng_coef_real,TbAng_coef_imag,V_load_0_real,V_load_0_imag);
			V_load_imag_out[j] = complex_multiply_Imag(TbAng_coef_real,TbAng_coef_imag,V_load_0_real,V_load_0_imag);	
		
			V_load_0_real = V_load_real_in[j] + V_load_real_out[j];
			V_load_0_imag = V_load_imag_in[j] + V_load_imag_out[j];	
			j++;
		}
		else 
		{
			V_load_0_real = complex_multiply_Real(TbAng_coef_real,TbAng_coef_imag,V_load_0_real,V_load_0_imag);
			V_load_0_imag = complex_multiply_Imag(TbAng_coef_real,TbAng_coef_imag,V_load_0_real,V_load_0_imag);
		}
		
	}
	
	V_load_0_real_out[0] = V_load_0_real;
	V_load_0_imag_out[0] = V_load_0_imag;
}