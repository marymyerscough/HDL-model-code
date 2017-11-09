#include "auto_f2c.h"

int func(integer ndim, const doublereal *u, const integer *icp, const doublereal *par, integer ijac, doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
	//p
	f[0] = u[1];
	f[1] = -par[1]/par[0]*u[4]*u[6]/(1.0+u[4])+par[2]/par[0]*u[0];
	//q
	f[2] = u[3];
	f[3] = -par[4]/par[3]*u[4]*u[6]/(1.0+u[4])+par[5]/par[3]*u[2];
	//l
	f[4] = u[5];
	f[5] = par[7]/par[6]*u[4]*u[6]/(1.0+u[4])+par[8]/par[6]*u[4];
	//m
	f[6] = u[7];
	f[7] = par[10]/par[9]*(u[5]*u[7]+u[6]*(par[7]/par[6]*u[4]*u[6]/(1.0+u[4])+par[8]/par[6]*u[4]))+par[11]/par[9]*u[4]*u[6]/(1.0+u[4])-par[12]/par[9]*u[8]*u[10]/(par[17]+u[8])+par[13]/par[9]*u[6];
	//h
	f[8] = u[9];
	f[9] = par[15]/par[14]*u[8]*u[10]/(par[17]+u[8])+par[16]/par[14]*u[8];
	//N
	f[10] = u[11];
	f[11] = -par[11]/par[18]*u[4]*u[6]/(1.0+u[4])+par[12]/par[18]*u[8]*u[10]/(par[17]+u[8]);
	
	return 0;
}

int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)
{
	par[0] = 1.0e6;			//Dp
	par[1] = 1.0e6;			//Mup
	par[2] = 1.0e3;			//d_p
    
	par[3] = 1.0e6;			//Dq
	par[4] = 1.0e6;			//Muq
	par[5] = 1.0e3;			//d_q
    
	par[6] = 1.0e4;			//Dl
	par[7] = 1.0e5;			//Mul
	par[8] = 1.0e2;			//d_l
	
	par[9] = 1.0e2;			//Dm
	par[10] = 1.0e4;		//Chim
	par[11] = 1.0e2;		//Mum
	par[12] = 1.0e1;		//NuN - - - - - 1.0e1
	par[13] = 1.0e0;		//d_m
	
	par[14] = 1.0e4;		//Dh
	par[15] = 1.0e4;		//Nuh - - - - - 1.0e4
	par[16]	= 1.0e2;		//d_h
	par[17]	= 1.0e0;		//Kappa
	
	par[18] = 1.0e-2;		//DN
		
	par[19] = 1.0e5;		//Sigmap1
	par[20] = 1.0e4;		//Sigmap2
	par[21] = 1.0e0;		//Betap
	
	par[22] = 1.0e0;		//Sigmaq
    
	par[23] = 1.0e3;		//Sigmal
	par[24] = 1.0e1;		//Alphal - - - - - 1.0e1
	par[25] = 1.0e-1;		//Gammal
    
	par[26] = 4.0e-3;		//Sigmam
	par[27] = 1.0e1;		//Alpham - - - - - 1.0e1
	par[28] = 1.0e-1;		//Gammam	
	par[29] = 1.0e-1;		//P_0
	par[30] = 1.0e-1;		//A_0
	
	par[31] = 8.0e2;		//Sigmah
	
	par[32] = 1.0e1;		//NuN2
	
	u[0] = 0.0;
	u[1] = 0.0;
	u[2] = 0.0;
	u[3] = 0.0;
	u[4] = 0.0;
	u[5] = 0.0;
	u[6] = 0.0;
	u[7] = 0.0;
	u[8] = 0.0;
	u[9] = 0.0;
	u[10] = 0.0;
	u[11] = 0.0;
	
	return 0;
}

int pvls(integer ndim, const doublereal *u, doublereal *par)
{
	par[33]=getp("BV0",1,u);	//p(0)
	par[34]=getp("BV0",3,u);	//q(0)
	par[35]=getp("BV0",5,u);	//l(0)
	par[36]=getp("BV0",7,u);	//m(0)
	par[37]=getp("BV0",9,u);	//h(0)
	par[38]=getp("BV0",11,u);	//N(0)
	
	par[39]=getp("BV0",12,u);	//N'(0)
	par[40]=getp("BV1",12,u);	//N'(1)
	
	par[41]=getp("INT",11,u);	//int(N)
	
	//The first argument of GETP may be one of the following:
	//'NRM' (L2-norm),     'MAX' (maximum),
	//'INT' (integral),    'BV0 (left boundary value),
	//'MIN' (minimum),     'BV1' (right boundary value).
	//'MNT' (t value for minimum)
	//'MXT' (t value for maximum)
	//'NDIM', 'NDX' (effective (active) number of dimensions)
	//'NTST' (NTST from constant file)
	//'NCOL' (NCOL from constant file)
	//'NBC'  (active NBC)
	//'NINT' (active NINT)
	//'DTM'  (delta t for all t values, I=1...NTST)
	//'WINT' (integration weights used for interpolation, I=0...NCOL)
	
	//Also available are
	//'STP' (Pseudo-arclength step size used).
	//'FLD' (`Fold function', which vanishes at folds).
	//'BIF' (`Bifurcation function', which vanishes at singular points).
	//'HBF' (`Hopf function'; which vanishes at Hopf points).
	//'SPB' ( Function which vanishes at secondary periodic bifurcations).
	//'EIG' ( Eigenvalues/multipliers, I=1...2*NDIM, alternates real/imag parts).
	//'STA' ( Number of stable eigenvalues/multipliers).
	
    return 0;
}

int bcnd(integer ndim, const doublereal *par, const integer *icp, integer nbc, const doublereal *u0, const doublereal *u1, integer ijac, doublereal *fb, doublereal *dbc)
{
	fb[0] = u0[1]+par[19]/par[0]*u0[4]/(par[21]+u0[4])+par[20]/par[0]*u0[2];
    fb[1] = u1[1];
	fb[2] = u0[3]-par[22]/par[3]*u0[2];
	fb[3] = u1[3];
	fb[4] = u0[5]+par[23]/par[6]*(1.0+u0[8]/par[24])/(1.0+u0[8]/par[25]);
	fb[5] = u1[5];
	fb[6] = u0[7]+par[26]/par[9]*(1+u0[8]/par[27])/(1+u0[8]/par[28])*(u0[0]-par[29])*(1.0+par[30]*u0[2])+par[10]*par[23]/(par[9]*par[6])*u0[6]*(1.0+u0[8]/par[24])/(1.0+u0[8]/par[25]);
	fb[7] = u1[7];
	fb[8] = u0[9]+par[31]/par[14];
	fb[9] = u1[9];
	fb[10] = u0[11];
	fb[11] = u1[11];
	
	return 0;
}

int icnd(integer ndim, const doublereal *par, const integer *icp, integer nint, const doublereal *u, const doublereal *uold, const doublereal *udot, const doublereal *upold, integer ijac, doublereal *fi, doublereal *dint)
{
    return 0;
}

int fopt(integer ndim, const doublereal *u, const integer *icp, const doublereal *par, integer ijac, doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
    return 0;
}