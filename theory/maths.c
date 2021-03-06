#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "fftw3.h"
#include "maths.h"


con constants = {
     3.14159265358979323846, //pi
     9.86960440108935861883, 
     6.28318530717958647693, //twopi
     0.69314718,
     2.90888208665721580e-4, //arcmin
     299792.458 //speed of light km/s
};
 
pre precision= {
	1e-5, //eps1
	1e-15, //eps2
	1e-10 //eps3
      	
};

lim limits = {
	0.01, //a_min
	1.e-3, //Delta_L_klin_min
	1.e5, //Delta_L_klin_max
	1.0e-3, //Pdelta_halo_k_min
	1.0e7,  //Pdelta_halo_k_max
	1.0,//P_2_s_min
	1.0e6,//P_2_s_max
	3.0e-5,//xi_via_hankel_theta_min
	1.0, //xi_via_hankel_theta_max
	1.e+4, //M_min
	1.e+16,//M_max
	20.0, //halo_s_min
	10000.0 //halo_s_max
};


 /* ============================================================ *
 * GSL_zeugs.					*
 * ============================================================ */
 

/* ============================================================ *
 * Complex gamma function.					*
 * ============================================================ */

void cdgamma(fftw_complex x, fftw_complex *res)
{
	double          xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

	xr = (double) x[0];
	xi = (double) x[1];

	if (xr<0) {
		wr = 1 - xr;
		wi = -xi;
	} else {
		wr = xr;
		wi = xi;
	}

	ur = wr + 6.00009857740312429;
	vr = ur * (wr + 4.99999857982434025) - wi * wi;
	vi = wi * (wr + 4.99999857982434025) + ur * wi;
	yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
			0.293729529320536228;
	yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
	ur = vr * (wr + 4.00000003016801681) - vi * wi;
	ui = vi * (wr + 4.00000003016801681) + vr * wi;
	vr = ur * (wr + 2.99999999944915534) - ui * wi;
	vi = ui * (wr + 2.99999999944915534) + ur * wi;
	yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
	yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
	ur = vr * (wr + 2.00000000000603851) - vi * wi;
	ui = vi * (wr + 2.00000000000603851) + vr * wi;
	vr = ur * (wr + 0.999999999999975753) - ui * wi;
	vi = ui * (wr + 0.999999999999975753) + ur * wi;
	yr += ur * 10.5400280458730808 + vr;
	yi += ui * 10.5400280458730808 + vi;
	ur = vr * wr - vi * wi;
	ui = vi * wr + vr * wi;
	t = ur * ur + ui * ui;
	vr = yr * ur + yi * ui + t * 0.0327673720261526849;
	vi = yi * ur - yr * ui;
	yr = wr + 7.31790632447016203;
	ur = log(yr * yr + wi * wi) * 0.5 - 1;
	ui = atan2(wi, yr);
	yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
	yi = ui * (wr - 0.5) + ur * wi;
	ur = yr * cos(yi);
	ui = yr * sin(yi);
	yr = ur * vr - ui * vi;
	yi = ui * vr + ur * vi;
	if (xr<0) {
		wr = xr * 3.14159265358979324;
		wi = exp(xi * 3.14159265358979324);
		vi = 1 / wi;
		ur = (vi + wi) * sin(wr);
		ui = (vi - wi) * cos(wr);
		vr = ur * yr + ui * yi;
		vi = ui * yr - ur * yi;
		ur = 6.2831853071795862 / (vr * vr + vi * vi);
		yr = ur * vr;
		yi = ur * vi;
	}

	(*res)[0]=yr; (*res)[1]=yi;
}




/*=============Besselfunctions==========================*/

double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
				+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
				+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
				+y*(-0.6911147651e-5+y*(0.7621095161e-6
				-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}


double bessj1(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
				+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
				+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
				+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
				+y*(0.8449199096e-5+y*(-0.88228987e-6
				+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
double bessj(int n, double x)
{
	int j,jsum,m;
	double ax,bj,bjm,bjp,sum,tox,ans;

	if (n < 2) sm2_error("Index n less than 2 in bessj");
	ax=fabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > (double) n) {
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+(int) sqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans : ans;
}
#undef ACC
#undef BIGNO
#undef BIGNI 


void sm2_error(char *s)
{
	printf("error:%s\n ",s);
	exit(1);
}

double **sm2_matrix(long nrl, long nrh, long ncl, long nch)
			/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) sm2_error("allocation failure 1 in sm2_matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) sm2_error("allocation failure 2 in sm2_matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void sm2_free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
		/* free a double matrix allocated by sm2_matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* ============================================================ *
 * dfridr.c							*
 * NR page 188. Returns derivate of func at x, initial step is  *
 * h. Error estimate in err.					*
 * Modified! func depends on two double! (like P_L)		*
 * ============================================================ */

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

double *sm2_vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) sm2_error("allocation failure in double vector()");
	return v-nl+NR_END;
}

int *int_vector(long nl, long nh)
			/* allocate a int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) sm2_error("allocation failure in int vector()");
	return v-nl+NR_END;
}

long *long_vector(long nl, long nh)
			/* allocate a int vector with subscript range v[nl..nh] */
{
	long *v;

	v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) sm2_error("allocation failure in int vector()");
	return v-nl+NR_END;
}

void sm2_free_vector(double *v, long nl, long nh)
		/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double sm2_dfridr(double (*func)(double,double), double x, double h, 
		  double *err, double aa)
{
	int i,j;
	double errt,fac,hh,**a,ans;

	ans = 1e30; /* dummy initialization */
	if (h == 0.0) sm2_error("h must be nonzero in sm2_dfridr.");
	a=sm2_matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	sm2_free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=sm2_vector(1,n);
	d=sm2_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
				sm2_error("Error in routine sm2_polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	sm2_free_vector(d,1,n);
	sm2_free_vector(c,1,n);
}


void sm2_spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=sm2_vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	sm2_free_vector(u,1,n-1);
}

void sm2_splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
		sm2_error("Bad xa input to routine sm2_splint");
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+
			((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#define FUNC(x) ((*func)(x))

double sm2_trapzd(double (*func)(double), double a, double b, int n, double *s)
{
	double x,tnm,sum,del;
	int it,j;

	if (n == 1) {
		return (*s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
			sum += FUNC(x);
		}
		*s=0.5*(*s+(b-a)*sum/tnm);
		return *s;
	}
}
#undef FUNC

/* ============================================================ *
 * Interpolates f at the value x, where f is a double[n] array,	*
 * representing a function between a and b, stepwidth dx.	*
 * 'lower' and 'upper' are powers of a logarithmic power law	*
 * extrapolation. If no	extrapolation desired, set these to 0	*
 * ============================================================ */
double sm2_interpol(double *f, int n, double a, double b, double dx, double x, double lower, double upper)
{
	double r;
	int  i;
	if (x < a) {
		if (lower==0.) {
			sm2_error("value too small in sm2_interpol");
			return 0.0;
		}
		return f[0] + lower*(x - a);
	}
	r = (x - a)/dx;
	i = (int)(floor(r));
	if (i+1 >= n) {
		if (upper==0.0) {
			if (i+1==n) {
				return f[i];  /* constant extrapolation */
			} else {
				sm2_error("value too big in sm2_interpol");
				return 0.0;
			}
		} else {
			return f[n-1] + upper*(x-b); /* linear extrapolation */
		}
	} else {
		return (r - i)*(f[i+1] - f[i]) + f[i]; /* interpolation */
	}
}


/* ============================================================ *
 * like interpol, but f beeing a 2d-function			*
 * 'lower' and 'upper' are the powers of a power law extra-	*
 * polation in the first argument				*
 * ============================================================ */
double sm2_interpol2d(double **f,
		      int nx, double ax, double bx, double dx, double x,
		      int ny, double ay, double by, double dy, double y,
		      double lower, double upper)
{
	double t, dt, s, ds;
	int i, j;
	if (x < ax-0.00001) {
		sm2_error("value too small in sm2_interpol2d");
	}
	if (x > bx+0.00001) {
	printf("%le %le\n",x,bx);	
	  sm2_error("value too big in sm2_interpol2d");
	}
	t = (x - ax)/dx; 
	i = (int)(floor(t));
	dt = t - i; 
	if (y < ay) {
		return ((1.-dt)*f[i][0] + dt*f[i+1][0]) + (y-ay)*lower;
	} else if (y > by) {
		return ((1.-dt)*f[i][ny-1] + dt*f[i+1][ny-1]) + (y-by)*upper;
	}
	s = (y - ay)/dy; 
	j = (int)(floor(s)); 
	ds = s - j; 
	if ((i+1==nx)&&(j+1==ny)) {
	  //printf("%d %d\n",i+1,j+1);
	  return (1.-dt)*(1.-ds)*f[i][j];
	  
	}
	if (i+1==nx){
	  //printf("%d %d\n",i+1,j+1);
	  return (1.-dt)*(1.-ds)*f[i][j]+ (1.-dt)*ds*f[i][j+1];
	}
	if (j+1==ny){
	  //printf("%d %d\n",i+1,j+1);
	  return (1.-dt)*(1.-ds)*f[i][j]+ dt*(1.-ds)*f[i+1][j];
	}
	return (1.-dt)*(1.-ds)*f[i][j] +(1.-dt)*ds*f[i][j+1] + dt*(1.-ds)*f[i+1][j] + dt*ds*f[i+1][j+1];
}


/* ============================================================ *
 * new interpol, but f beeing a 2d-function no extrapolation			*
 * ============================================================ */
double sm2_interpol2d_noextra(double **f,
		      int nx, double ax, double bx, double dx, double x,
		      int ny, double ay, double by, double dy, double y,
		      double lower, double upper)
{
	double t, dt, s, ds;
	int i, j;
	if (x < ax-0.00001) {
		sm2_error("value too small in sm2_interpol2d_noextra");
	}
	if (x > bx+0.00001) {
	printf("%le %le\n",x,bx);	
	  sm2_error("value too big in sm2_interpol2d_noextra");
	}
	if (y < ay-0.00001) {
		sm2_error("value too small in sm2_interpol2d_noextra");
	}
	if (y > by+0.00001) {
	printf("%le %le\n",y,by);	
	  sm2_error("value too big in sm2_interpol2d_noextra");
	}
	// below is only for roundoff uncertainties below the 0.00001 tolerance criteria above
	if ((x < ax) && (y < ay)) return f[0][0];
	if ((x > bx) && (y > by)) return f[nx-1][ny-1];
	
	t = (x - ax)/dx; 
	i = (int)(floor(t));
	dt = t - i; 
	
	s = (y - ay)/dy; 
	j = (int)(floor(s)); 
	ds = s - j; 

	if ((i+1==nx)&&(j+1==ny)) return (1.-dt)*(1.-ds)*f[i][j];
	
	if (i+1==nx) return (1.-dt)*(1.-ds)*f[i][j]+ (1.-dt)*ds*f[i][j+1];
	
	if (j+1==ny) return (1.-dt)*(1.-ds)*f[i][j]+ dt*(1.-ds)*f[i+1][j];	

	return (1.-dt)*(1.-ds)*f[i][j] +(1.-dt)*ds*f[i][j+1] + dt*(1.-ds)*f[i+1][j] + dt*ds*f[i+1][j+1];
}
double sm2_interpol_coyote(double **f,
		      int nx, double ax, double bx, double dx, double x,
		      int ny, double ay, double by, double dy, double y)
{
	double t, dt, s, ds;
	int i, j;
	if (x < ax) {
		sm2_error("value too small in sm2_interpol2d");
	}
	if (x > bx) {
		sm2_error("value too big in sm2_interpol2d");
	}
	t = (x - ax)/dx;

	i = (int)(floor(t));
//	printf("%d %d\n",i,nx);
	if (i+1 > nx || i < 0) {
		printf("%d %d\n",i,nx);
		sm2_error("index out of range in sm2_interpol_coyote");
		}
	dt = t - i;
	s = (y - ay)/dy;
	j = (int)(floor(s));
	ds = s - j;
	if ((i+1==nx)&&(j+1==ny)) return (1.-dt)*(1.-ds)*f[i][j];
	if (i+1==nx) return (1.-dt)*(1.-ds)*f[i][j]+ (1.-dt)*ds*f[i][j+1];
	if (j+1==ny) return (1.-dt)*(1.-ds)*f[i][j]+ dt*(1.-ds)*f[i+1][j];
	return (1.-dt)*(1.-ds)*f[i][j] +
			(1.-dt)*ds*f[i][j+1] +
			dt*(1.-ds)*f[i+1][j] +
			dt*ds*f[i+1][j+1];
}



double sm2_gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
		int j;

		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*log(tmp);
		ser=1.000000000190015;
		for (j=0;j<=5;j++) ser += cof[j]/++y;
		return (double)(-tmp+log(2.5066282746310005*ser/x));
}



/* ============================================================ *
 * qromb.c							*
 * Romberg Integration. Uses trapzd. NR p. 140			*
 * ============================================================ */

#define EPS 1.0e-6
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

double sm2_qromb(double (*func)(double), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1], strap;
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=sm2_trapzd(func,a,b,j,&strap);
		if (j >= K) {
			sm2_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	sm2_error("Too many steps in routine sm2_qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define EPS 1.0e-7
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

/* ============================================================ *
 * Romberg integration of an open intervall [a,b]. choose is a  *
 * pointer to a routine using an open quadrature formula.	*
 * Uses polint (Neville-Aitken) to extrapolate. NR p. 143	*
 * ============================================================ */

double sm2_qromo(double (*func)(double), double a, double b,
		 double (*choose)(double(*)(double), double, double, int))
{
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			sm2_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	sm2_error("Too many steps in routing sm2_qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define FUNC(x) ((*func)(x))

double sm2_midpnt(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC



void delta_theta_bin_mid(int Ndata,double *theta_binmid,double *deltatheta, int Lin)
{
	int i;
	double bin_min, bin_max, binbreite,logbinbreite;
	if (Lin ==1){
	    binbreite=(theta_binmid[Ndata-1] - theta_binmid[0])/(Ndata-1);
		for (i=0;i<Ndata;i++){
			deltatheta[i]=binbreite;
		}
		
	}
	// Note: deltatheta[i] has different extensions to both sides of theta_binmid[i], more precisely theta_binmid[i]+0.5deltatheta[i] != exp(bin_max)=exp(log(theta_binmid[i])+0.5*logbinbreite;)
	if (Lin ==0){
		logbinbreite=(log(theta_binmid[Ndata-1])-log(theta_binmid[0]))/(Ndata-1);
		for (i=0;i<Ndata;i++){
		    bin_min=log(theta_binmid[i])-0.5*logbinbreite;
		    bin_max=log(theta_binmid[i])+0.5*logbinbreite;
		    deltatheta[i]=exp(bin_max)-exp(bin_min);
		}
	}
}




