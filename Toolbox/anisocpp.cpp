#include <matrix.h>
#include <mex.h>
#include <cmath>
#include <string.h>


/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif


void fillzero(double *output, int length){
    for (int i = 0; i < length; i++){
        output[i] = 0;
    }
    
}// end of function


void setdiffl(const double *input, double *output, int length)
{
    int k = 513;// actually should be 515 but for mod 0 is 0 so decrease -2
    
    for (int i = 0; i < length; i++){
        
        
        if( (i%512)  == 0 )
        {
            k +=2;
        }
        
        
        output[i+k] = input[i];
    }
    
    
}// end of function


void setdelta1(const double *dff, const double *dffl, double *output, int nswe, int length)
{
    int k = 0;
    
    switch(nswe)
    {
        case 1 :
            k=512;
            break;
        case 2 :
            k=514;
            break;
        case 3 :
            k=1027;
            break;
        case 4 :
            k=-1;
            break;
    }
    
    
    for (int i = 0; i < length; i++){
        
        if( (i%512)  == 0 )
        {
            k +=2;
        }
        
        output[i] =  dffl[i+k] - dff[i];
    }
    
    
    
}// end of function



void setdelta2(const double *dff, const double *dffl, double *deltaN, double *deltaS,
        double *deltaE, double *deltaW, int nswe, int dimx, int dimy)
{
    int k=0;
    
    for (int i = 0; i < dimx; i++)
    {
        
        for (int j = 0; j < dimy; j++)
        {
            
            deltaN[i*dimy+j] =  dffl[(i+1)*dimy+(j+2+k)] - dff[i*dimy+j];
            deltaS[i*dimy+j] =  dffl[(i+1)*dimy+(j+4+k)] - dff[i*dimy+j];
            deltaE[i*dimy+j] =  dffl[(i+2)*dimy+(j+5+k)] - dff[i*dimy+j];
            deltaW[i*dimy+j] =  dffl[(i)*dimy+(j+1+k)] - dff[i*dimy+j];
            
        }
        
        k +=2;
    }
    
}// end of function

void setdelta3(const double *dff, const double *dffl, double *deltaN, double *deltaS,
        double *deltaE, double *deltaW, int nswe, int length)
{
    int k=0;
    
    for (int i = 0; i < length; i++){
        
        if( (i%512)  == 0 )
        {
            k +=2;
        }
        
        deltaN[i] =  dffl[i+k+512] - dff[i];
        deltaS[i] =  dffl[i+k+514] - dff[i];
        deltaE[i] =  dffl[i+k+1027] - dff[i];
        deltaW[i] =  dffl[i+k-1] - dff[i];
    }
    
    
}// end of function


void doOption1(const double *delta, double *output, double kappa, int length)
{
    for (int i = 0; i < length; i++){
        
        output[i] =  exp( -(pow(delta[i]/kappa, 2)) );
    }
    
}// end of function

void doOption12(const double *deltaN, const double *deltaS,const double *deltaE,
        const double *deltaW, double *cN, double *cS, double *cE, double *cW,
        double kappa, int length)
{
    for (int i = 0; i < length; i++){
        
        cN[i] =  exp( -(pow(deltaN[i]/kappa, 2)) );
        cS[i] =  exp( -(pow(deltaS[i]/kappa, 2)) );
        cE[i] =  exp( -(pow(deltaE[i]/kappa, 2)) );
        cW[i] =  exp( -(pow(deltaW[i]/kappa, 2)) );
        
    }
    
}// end of function



void doOption2(const double *delta, double *output, double kappa, int length)
{
    for (int i = 0; i < length; i++){
        
        output[i] =  1/( 1 + (pow(delta[i]/kappa, 2)) );
    }
    
}// end of function


void doOption22(const double *deltaN, const double *deltaS,const double *deltaE,
        const double *deltaW, double *cN, double *cS, double *cE, double *cW,
        double kappa, int length)
{
    for (int i = 0; i < length; i++){
        
        cN[i] =  1/( 1 + (pow(deltaN[i]/kappa, 2)) );
        cS[i] =  1/( 1 + (pow(deltaS[i]/kappa, 2)) );
        cE[i] =  1/( 1 + (pow(deltaE[i]/kappa, 2)) );
        cW[i] =  1/( 1 + (pow(deltaW[i]/kappa, 2)) );
    }
    
}// end of function



void doCalculation(const double *deltaN,const double *deltaS,const double *deltaE,const double *deltaW,
        const double *cN,const double *cS,const double *cE,const double *cW,
        double *output, double lambda, int length)
{
    for (int i = 0; i < length; i++)
    {
        output[i] = output[i] + lambda*(cN[i]*deltaN[i] + cS[i]*deltaS[i] + cW[i]*deltaW[i] + cE[i]*deltaE[i]);
        
    }
    
}// end of function




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //declare variables
    mxArray *img_in, *img_out, *diffl;
    mxArray *deltaN, *deltaS, *deltaW, *deltaE;
    mxArray *cN, *cS, *cW, *cE;
    const mwSize *dims;
    double  kappa, lambda, *imgP, *diffP;
    int  num_iter, option;
    double  *deltaNP, *deltaSP, *deltaWP, *deltaEP;
    double  *cNP, *cSP, *cWP, *cEP;
    size_t dimx, dimy, numdims;
    int i,j;
    
    
//associate inputs
    img_in =  mxDuplicateArray(prhs[0]);
    kappa =   (double) mxGetScalar(prhs[1]);
    option =    mxGetScalar(prhs[2]);
    num_iter =  mxGetScalar(prhs[3]);
    lambda =  (double) mxGetScalar(prhs[4]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    
//associate outputs
    img_out = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    diffl  =  mxCreateDoubleMatrix(dimy+2,dimx+2,mxREAL);
    deltaN  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    deltaS  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    deltaE  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    deltaW  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    cN  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    cS  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    cE  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    cW  =  mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    
//associate pointers
    imgP =  mxGetPr(img_in);
    //outP =  mxGetPr(img_out);
    diffP = mxGetPr(diffl);
    deltaNP = mxGetPr(deltaN);
    deltaSP = mxGetPr(deltaS);
    deltaEP = mxGetPr(deltaE);
    deltaWP = mxGetPr(deltaW);
    cNP = mxGetPr(cN);
    cSP = mxGetPr(cS);
    cEP = mxGetPr(cE);
    cWP = mxGetPr(cW);
    
    
    int lImSize = (dimx+2)*(dimy+2);
    int ImSize = dimx * dimy;
    
    
    /* Test
     * fillzero(diffP, lImSize);
     *
     * setdiffl(imgP, diffP, dimx, dimy);
     *
     * setdelta(imgP, diffP, deltaNP, 1, dimx, dimy);
     *
     * setdelta(imgP, diffP, deltaSP, 2, dimx, dimy);
     *
     * setdelta(imgP, diffP, deltaEP, 3, dimx, dimy);
     *
     * setdelta(imgP, diffP, deltaWP, 4, dimx, dimy);
     *
     * doOption2(deltaNP, cNP, kappa, dimx, dimy);
     * doOption2(deltaSP, cSP, kappa, dimx, dimy);
     * doOption2(deltaEP, cEP, kappa, dimx, dimy);
     * doOption2(deltaWP, cWP, kappa, dimx, dimy);
     *
     * doCalculation(deltaNP,deltaSP,deltaEP,deltaWP,cNP,cSP,cEP,cWP,imgP,lambda,dimx,dimy);
     *
     * mexPrintf("\n%f\n",imgP[102655]);
     *
     */
        
    
    for(i=0;i<num_iter;i++)
    {
        memset (diffP,0,lImSize-1); // alternative to fillzero a bit faster
        
        //fillzero(diffP, lImSize);      
        
        setdiffl(imgP, diffP, ImSize);
        
        setdelta1(imgP, diffP, deltaNP, 1, ImSize);
        
        setdelta1(imgP, diffP, deltaSP, 2, ImSize);
        
        setdelta1(imgP, diffP, deltaEP, 3, ImSize);
        
        setdelta1(imgP, diffP, deltaWP, 4, ImSize);
        
        
        //setdelta2(imgP, diffP, deltaNP, deltaSP, deltaEP, deltaWP, 1, dimx, dimy);
        
        //setdelta3(imgP, diffP, deltaNP, deltaSP, deltaEP, deltaWP, 1, ImSize);
        
        
        if(option == 1)
        {
            
            doOption1(deltaNP, cNP, kappa, ImSize);
            doOption1(deltaSP, cSP, kappa, ImSize);
            doOption1(deltaEP, cEP, kappa, ImSize);
            doOption1(deltaWP, cWP, kappa, ImSize);
            
            //doOption12(deltaNP,deltaSP,deltaEP,deltaWP, cNP, cSP, cEP, cWP, kappa, ImSize);
            
        }// end of if
        else
        {
            
            doOption2(deltaNP, cNP, kappa, ImSize);
            doOption2(deltaSP, cSP, kappa, ImSize);
            doOption2(deltaEP, cEP, kappa, ImSize);
            doOption2(deltaWP, cWP, kappa, ImSize);
            
            //doOption22(deltaNP,deltaSP,deltaEP,deltaWP, cNP, cSP, cEP, cWP, kappa, ImSize);
            
            
        }//end of else
        
        
        doCalculation(deltaNP,deltaSP,deltaEP,deltaWP,cNP,cSP,cEP,cWP,imgP,lambda,ImSize);
        
    } // end of for
    
    
    // set output
    plhs[0] = img_in;
    
     //mexPrintf("\n lImSize %f\n",lImSize);
    /**
     * mexPrintf("\n lImSize %d\n",lImSize);
     * mexPrintf("\n Kappa %f\n",kappa);
     * mexPrintf("\n Option %d\n",option);
     * mexPrintf("\n Num Iter %d\n",num_iter);
     * mexPrintf("\n Lambda %f\n",lambda);
     **/
    
    return;
}// end of function