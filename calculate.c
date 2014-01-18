// calculate the jitter parameter
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf.h>
#include "ptime_jitter.h"

#define ITMAX 10000  // Maximum allowed number of iterations.
#define EPS 3.0e-16 // Machine double floating-point precision.
//#define EPS 3.0e-8 // Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}

	return temp[n-1];
}

int remove_baseline (double *in, double frac_off, int n, double *out)
{
	// define the off_pulse range of std, frac_off is the fraction of the phase
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse(n, in, off_0, frac_off);

	// normalize std
	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
		//printf ("%lf \n", off_0[i]);
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int off_pulse (int nphase, double *in, double *out, double frac)
// define the on_pulse range, 10% of the phase
{
	int n = nphase;
	int num = (int)(n*frac);
	int i,j;
	double small;
	double ave;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num; j++)
			{
				small += in[j];
			}
			small = small/num;
		}
			
		ave = 0.0;
		for(j = 0; j < num; j++)
		{
			if ((i+j) > n-1)
			{
				ave += in[(i+j)-(n-1)];
			}
			else 
			{
				ave += in[i+j];
			}
		}
		ave = ave/num;

		small = (ave <= small ? ave : small);
		index = (ave <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	for (i = 0; i < num; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

int pre (double *s, int nphase, double *real_s, double *ima_s)
// dft std, output real_s, ima_s  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i = 0; i < nphase; i++)
	{
		test[i] = s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

    fftw_complex *out_s;
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	

	{
	    dft_profiles(nphase,s,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];
		/*
		for (i = 0; i < nphase; i++)                                
		{                                                        
			printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		}
		*/

		for (i = 0; i < nphase/2+1; i++)                                
		{                                                        
			real_s[i]=out_s[i][0];                                       
			ima_s[i]=out_s[i][1];                                        
		}
										
	}

	fftw_free(out_s); 
	fftw_free(out_t); 

	return 0;
}

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new)
// inverse dft, and add pi/2 phase
{
	double *dp;
    fftw_plan plan;
	fftw_complex *cp;

    dp = (double *)malloc(sizeof (double) * ncount);
	cp = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp, 0, sizeof (double) * ncount);
	memset(cp, 0, sizeof (fftw_complex) * ncount);

	// initialize the dft...
	double *dp_t;
    fftw_plan plan_t;
	fftw_complex *cp_t;

    dp_t = (double *)malloc(sizeof (double) * ncount);
	cp_t = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp_t, 0, sizeof (double) * ncount);
	memset(cp_t, 0, sizeof (fftw_complex) * ncount);

	int i;
    double real,ima,amp,cosina,sina;

	for (i = 0; i < ncount; i++)
	{
		if (i < ncount/2+1)
		//if (i < ncount/2)
		{
            real = real_p[i];
            ima = ima_p[i];
			amp = sqrt(real*real+ima*ima);

			if ( amp != 0.0 )
			{
				cosina = real/amp;
				sina = ima/amp;

				cp[i][0] = amp*(cosina*cos(-i*3.1415926)-sina*sin(-i*3.1415926));
				cp[i][1] = amp*(sina*cos(-i*3.1415926)+cosina*sin(-i*3.1415926));

				cp_t[i][0] = real_p[i];
				cp_t[i][1] = ima_p[i];
			}
			else
			{
				cp[i][0]=0.0;
				cp[i][1]=0.0;
				cp_t[i][0]=0.0;
				cp_t[i][1]=0.0;
			}
		}
		else
		{
			cp[i][0]=0.0;
			cp[i][1]=0.0;
			cp_t[i][0]=0.0;
			cp_t[i][1]=0.0;
		}
	}

    plan_t = fftw_plan_dft_c2r_1d(ncount, cp_t, dp_t, FFTW_MEASURE);

    fftw_execute(plan_t);

    fftw_destroy_plan(plan_t);

	/////////////////////////////////////////////////////////////////

    plan = fftw_plan_dft_c2r_1d(ncount, cp, dp, FFTW_MEASURE);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

	for (i = 0; i < ncount; i++)
	{
		p_new[i] = dp[i]/ncount;  // normalized by the ncount
		//printf ("%lf\n", p_new[i]);
	}

	return 0;
}

double u_0 (double *y, int nphase, double psrfreq)
{
	// move the std to the center
	double real_s[nphase/2+1],ima_s[nphase/2+1];
	pre(y, nphase, real_s, ima_s);

	/*
	int j;
	for (j = 0; j < nphase/2+1; j++)
	{
		printf ("%d %lf\n", j, real_s[j]);
	}
	*/

	double std_0[nphase];
	inverse_dft (real_s, ima_s, nphase, std_0);
	
	// remove the baseline
	double std_1[nphase];
	remove_baseline (std_0, 0.2, nphase, std_1);

	// normalize the std
	double peak;
	peak = find_peak_value(nphase, std_1);
	printf ("peak is: %lf\n", peak);

    double u_0[nphase];
    double t[nphase];

    double t_step=(1.0/psrfreq)/nphase;
    //double t_step=(1.0/173.687946184761)/1024;
    
    int i;
    for (i = 0; i < nphase; i++)
    {
        u_0[i]=std_1[i]/peak;
		t[i]=i*t_step-t_step*nphase/2;
    }

	/*
	int j;
	for (j = 0; j < nphase; j++)
	{
		printf ("%d %lf\n", j, std_0[j]);
	}
	*/

    double u_1=0.0;
    double u_2=0.0;
   
    for (i = 0; i < nphase-1; i++)
    //for (i = 0; i < nphase; i++)
    {
		u_1 += 0.5*(u_0[i]+u_0[i+1])*t_step;
		u_2 += 0.5*(u_0[i]+u_0[i+1])*pow(0.5*(t[i]+t[i+1]),2)*t_step;
		/*
		{
			if (i==0 || i==1023)
			{
				u_1 += 0.5*u_0[i]*t_step;
				//u_2+=0.5*u_0[i]*t[i]*t_step;
				u_2 += 0.5*u_0[i]*t[i]*t[i]*t_step;
			}
			else
			{
				u_1 += u_0[i]*t_step;
				//u_2+=u_0[i]*t[i]*t_step;
				u_2 += u_0[i]*t[i]*t[i]*t_step;
			}
		}
		*/
	}

    double U0;
    double N = floor(60*psrfreq); // number of pulse of each subint
    //double N=floor(64*60*psrfreq);

    //U0=u_2/(u_1);
    U0 = u_2/(u_1*N);

	printf ("u_1: %.10lf \n", u_1);
	printf ("u_2: %.10lf \n", u_2);
	printf ("sqrt((u_2/u_1)) (microsecond): %lf \n", sqrt((u_2/u_1))*1.0e+6);
	printf ("U0: %lf \n", U0*1e+6);

    return U0;
}

double chi (double sigma_J, double k, int nint, double *resi, double *e, double reduced_chi)
{
    //printf("Number of subint is: %d \n", nint);
    double n=nint;

    /////////////////////////////////////////////////////////////////////////////////
    double chisqr=0.0;
    //double sigma_J=1.896*1.35e-07;
    //double sigma_J=1.0e-07;
    int i;
    
    for (i = 0; i < n; i++)
    {
        //chisqr+=(3.0*resi[i]/(sigma_J+e[i]))*(3.0*resi[i]/(sigma_J+e[i]));
        chisqr+=(resi[i]/(sigma_J+e[i]))*(resi[i]/(sigma_J+e[i]));
        //chisqr+=(resi[i]/e[i])*(resi[i]/e[i]);
        //printf ("p_bar is: %f\n", p_bar);
    }

	//double reduced_chi = 1.0;
	double F;

    F=(chisqr/k)-reduced_chi;
    //F=(chisqr/63.0)-1.075;
    //F=(chisqr/63.0)-1.0;
    //F=(chisqr/63.0)-1.0;
    //printf ("chisqr/d.o.f is: %f\n", chisqr);

    ////////////////////////////////////////////////////////////////////////////////

    return F;
}

double zbrent(double (*func)(double sigma_J, double k, int nint, double *resi, double *error, double reduced_chi), double x1, double x2, double tol, double k, int nint, double *resi, double *error, double reduced_chi)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, k, nint, resi, error, reduced_chi),fb=(*func)(b, k, nint, resi, error, reduced_chi),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, k, nint, resi, error, reduced_chi);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

double chi_distribution (double *k, size_t dim, void *p)
// calculate the error of chi, then use it to derive the error of jitter parameter
{
	double f;
	//double k0;
	//double Gamma=8.222839e+33;

	struct my_params * fp=(struct my_params *)p;

	//printf ("chi_distribution: %lf %lf\n", fp->a, fp->b);
	//f=1.0;
	f=(1.0/((fp->b)*pow(2.0,(fp->a)/2.0)))*pow(k[0],(fp->a)/2.0-1.0)*exp(-k[0]/2.0);

	return f;
}

double error (double x, void *p, size_t N)
// calculate the error of chi, then use it to derive the error of jitter parameter
{
	//printf ("Integrating...\n");
	struct my_params *fp=(struct my_params *)p;
	//printf ("error: %lf %lf\n", fp->a, fp->b);

	double F;
    double res,err;

    double xl[1]={0.0};
    double xu[1]={x};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&chi_distribution,1,fp};

    size_t calls = N;
    gsl_rng_env_setup ();

    T=gsl_rng_default;

    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G, xl, xu, 1, 10000, r, s, &res, &err);
	//printf ("result0: %lf %lf\n", res, err);

    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 1, calls/5, r, s, &res, &err);
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

	//printf ("result of integration is: %e\n", res);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

	F=res-0.682; // 1 sigma
	//F=res-0.9;

    return F;
}

int Gamma_cal (double k, void *p)
// calculate the value of \Gamma(k/2), which will be used by error.c
{
	struct my_params *fp=(struct my_params *)p;
	//double k=63.0;
	//double k=64.0;
    //printf("k is: %lf\n", k);
    //printf("Gamma(k/2) is: %e\n", gsl_sf_gamma (k/2.0));
	fp->a=k;
	fp->b=gsl_sf_gamma (k/2.0);
	//printf ("Gamma_cal: %lf %lf\n", fp->a, fp->b);
	//(*Gamma)=gsl_sf_gamma (k/2.0);

    return 0;
}

double zbrent_err(double (*func)(double x, void *para, size_t N), double x1, double x2, double tol, void *para, size_t N)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	//struct my_params *fp=(struct my_params *)para;
	//printf ("zbrent_err: %lf %lf\n", fp->a, fp->b);
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, para, N),fb=(*func)(b, para, N),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, para, N);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

int jitter (double *resi, double *e, double *stdy, double psrfreq, int nint, int nphase, double degree, FILE *fp)
// calculate the jitter parameter and its error
{
    //puts (argv[1]);
    //puts (argv[2]);

    //int nint=0; 
    //int nsample=0; 
    //double resi[64],e[64];
    //double x[1024],y[1024];

    //read (argv[1], &nint, resi, e);
    //read (argv[2], &nsample, x, y);

    //printf("Number of subint is: %d \n", nint);

    /////////////////////////////////////////////////////////////////////////////////
	// calculate jitter parameter
	double k;
	//k=nint-10.0;
	k = degree;

	double reduced_chi;

	double sigma_J;

	reduced_chi = 1.0;

	// initial guess of the result
	/////////////////////////////////////////////////////////////////////
	double up = 1.0e-5;
	double step = up/100.0;
	double low = 1.0e-5;

	while( chi(up, k, nint, resi, e, reduced_chi)*chi(low, k, nint, resi, e, reduced_chi) > 0 )
	{
		low = low - step;
	}
	/////////////////////////////////////////////////////////////////////

	sigma_J = zbrent(chi, low, up, 1.0e-16, k, nint, resi, e, reduced_chi);

    printf("low: %e \n", low);
    printf("sigma_J: %e \n", sigma_J);
    printf("reduced_chisqr: %e \n", chi(sigma_J, k, nint, resi, e, reduced_chi));

    ////////////////////////////////////////////////////////////////////////////////

    double U;
    U = u_0(stdy, nphase, psrfreq);

    ///////////////////////////////////////////////////////////////////////////////

    double jitter;
    //double N=floor(60*173.687946184761);

    //jitter=sqrt(N)*sigma_J/U;
    jitter = sigma_J/sqrt(U);
    printf ("jitter parameter is: %f\n", jitter);
    //printf ("jitter parameter is: %f\n", sqrt(N)*sigma_J/0.14e-3);
    //printf ("jitter parameter is: %f\n", sqrt(N)*sigma_J/1.02e-3);
	
	//////////////////////////////////////////////////////////////////////////////////
	// calculte the error of jitter parameter
	struct my_params params;
    double chisqr_new;

	Gamma_cal(k, &params);

	// initial guess of the result
	/////////////////////////////////////////////////////////////////////
	double up_e = k+10;
	double step_e = up_e/100.0;
	double low_e = up_e;
    //printf ("%lf %lf\n", low_e, up_e);

	while( error(up_e, &params, 500000)*error(low_e, &params, 500000) > 0 )
	{
		low_e = low_e - step_e;
	}
    //printf ("%lf %lf\n", low_e, up_e);
	/////////////////////////////////////////////////////////////////////

	//chisqr_new = zbrent_err(error, low_e-10.0, up_e+10.0, 1.0e-16, &params); // get the 1 sigma chisqr
	chisqr_new = zbrent_err(error, low_e, up_e, 1.0e-16, &params, 500000); // get the 1 sigma chisqr

	reduced_chi = chisqr_new/k;
    printf ("chisqr_new: %lf\n", chisqr_new);

	/////////////////////////////////////////////////////////////////////

	// initial guess of the result
	/////////////////////////////////////////////////////////////////////
	up = 1.0e-6;
	step = up/100.0;
	low = 1.0e-6;

	while( chi(up, k, nint, resi, e, reduced_chi)*chi(low, k, nint, resi, e, reduced_chi) > 0 )
	{
		low = low - step;
	}
	/////////////////////////////////////////////////////////////////////

	double sigma_J_error, jitter_error;
	//sigma_J_error = zbrent(chi, 1e-7, 1e-5, 1.0e-16, k, nint, resi, e, reduced_chi);   // for 1 sigma reduced_chi, get the sigma_J
	sigma_J_error = zbrent(chi, low, up, 1.0e-16, k, nint, resi, e, reduced_chi);   // for 1 sigma reduced_chi, get the sigma_J

	jitter_error = fabs(sigma_J_error/sqrt(U)-jitter);
    printf ("The error of jitter parameter is: %f\n", jitter_error);
	
	// print out the results
    fprintf (fp, "%lf %lf\n", jitter, jitter_error);

    return 0;
}
