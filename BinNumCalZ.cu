/*
 * Bins fast calculation for later convolution.
 * Takes a gpuArray input and returns a gpuArray output
 * e.g.       bin_num_q           = BinNumCalZ(q);
 * or        [bin_num_q,bin_num_x]= BinNumCalZ(q,x);
 * [bin_num_q,bin_num_x,bin_num_y]= BinNumCalZ(q,x,y);
 * Author: Tianlong He. 
 * Time  : 20200519
 */
 
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "cuda_runtime.h"
/*
 * Device code
 */
void __global__ BinNumCal(float const * const Z,
                         float * const binZ,
                         int const N,
						 int const rowsZ,
						 int const binnum)
{
    /* Calculate the global linear index, assuming a 1-d grid. */
	int i;
	int index,indexz;
	int const iL = threadIdx.x * N;
	int const iR = iL + N;
	extern __shared__ float binZshare[];
	
	for (i=0;i<binnum;i++) {
		binZshare[threadIdx.x * binnum + i] = 0;
	}
	/*__syncthreads();*/
	for (i=iL;i<iR;i+=5) {
		
		if (i < rowsZ) {
			indexz = blockIdx.x * rowsZ + i;
			index = (int)(Z[indexz]);
			binZshare[threadIdx.x*binnum+index] += 1;
			index = (int)(Z[indexz + 1]);
			binZshare[threadIdx.x*binnum+index] += 1;
			index = (int)(Z[indexz + 2]);
			binZshare[threadIdx.x*binnum+index] += 1;
			index = (int)(Z[indexz + 3]);
			binZshare[threadIdx.x*binnum+index] += 1;
			index = (int)(Z[indexz + 4]);
			binZshare[threadIdx.x*binnum+index] += 1;
		}
	}
	/*__syncthreads();*/
	for (i=0;i<binnum;i++) {
		index= (i + blockIdx.x * binnum)*blockDim.x + threadIdx.x;
		binZ[index] = binZshare[threadIdx.x * binnum + i];
	}
	
}

/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare all variables.*/
    mxGPUArray const *Z;
    mxGPUArray *binZ;
	float const *d_Z;
    float *d_binZ;
	size_t binnum;
    int N;
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 5;
    int blocksPerGrid;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
    if ((nrhs!=2) || !(mxIsGPUArray(prhs[1]))) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }
	
	binnum = (int)(mxGetScalar(prhs[0]));
    Z = mxGPUCreateFromMxArray(prhs[1]);

    /*
     * Verify that Z really is a double array before extracting the pointer.
     */
    if (mxGPUGetClassID(Z) != mxSINGLE_CLASS) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    /* 
     * Now that we have verified the data type, extract a pointer to the input
     * data on the device.
     */
    d_Z = (float const *)(mxGPUGetDataReadOnly(Z));
	
	const mwSize *dimsZ = mxGPUGetDimensions(Z);
	size_t nrowsZ = dimsZ[0];
	size_t ncolsZ = dimsZ[1];	
			
	mwSize dims[2] = {threadsPerBlock, ncolsZ * binnum}; // note here	
	
    /* Create a GPUArray to hold the result and get its underlying pointer. */
    binZ = mxGPUCreateGPUArray(2,
                            dims,
                            mxGPUGetClassID(Z),
                            mxGPUGetComplexity(Z),
                            MX_GPU_INITIALIZE_VALUES);
							
    d_binZ = (float *)(mxGPUGetData(binZ));

    /*
     * Call the kernel using the CUDA runtime API. We are using a 1-d grid here,
     * and it would be possible for the number of elements to be too large for
     * the grid. For this example we are not guarding against this possibility.
     */
    N = ((int)(nrowsZ)+threadsPerBlock-1)/threadsPerBlock;
	
    blocksPerGrid = (int)(ncolsZ);
    BinNumCal<<<blocksPerGrid, threadsPerBlock, threadsPerBlock*binnum*sizeof(float)>>>(d_Z, d_binZ, N, nrowsZ, binnum);

    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(binZ);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(Z);
    mxGPUDestroyGPUArray(binZ);
}
 