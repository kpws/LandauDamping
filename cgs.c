#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <fftw3.h>
#include <sys/time.h>

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KRESET "\x1B[0m"

#define MAX_PG 10

//0=whole phi range, small r range
//1=small phi range, large r range
double complex *polarDatas[MAX_PG];
fftw_complex *Gs, *Gso;
int nphis[MAX_PG],nrs[MAX_PG],n,cn,nThreads,nPolarGrid;
double *phiss[MAX_PG], *rss[MAX_PG],Nf,sqrtNf,L,cL,alpha,alpha3,csigma;
double cs=1.;
int acn;
double complex *acs;

int nFactor, cn2;
double delta, cL2, *ker, cDelta, Delta;

int cubic=1,useAsymCoef=0, debug=0;

static inline double complex cubicInterp(double complex p0,double complex p1,double complex p2,double complex p3, double x) {
	if(cubic)
		return p1 + 0.5 * x*(p2 - p0 + x*(2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + x*(3.0*(p1 - p2) + p3 - p0)));
	else
		return p1*(1.-x)+x*p2;
}

double complex asymExponent(double tau, double x, double lsqrtNf){
	double c1=0.0243615448342451073975190249505015144289430415363799596395;
	complex double res=(  3/(36*M_PI*M_PI)  -  c1*pow(fabs(x)*lsqrtNf,1./3)*cpow(1+I*tau/x,-2./3)  ) / lsqrtNf;
	if(useAsymCoef){
		double r=sqrt(tau*tau + x*x);
		double phi=atan(fabs(x/tau));
		double s=(acn-1)*phi/(M_PI/2.);
		int i=s;
		s-=i;
		double complex coef=acs[i]*(1.-s)+acs[i+1]*s;
		if(x*tau<0.)coef=conj(coef);
		
		res=r*pow(lsqrtNf*r,-4./3)*coef/phi;
		res=r*pow(lsqrtNf,-28./3)*pow(lsqrtNf*r,-8./3)*coef/phi;
		res=pow(lsqrtNf*r,-1./3)*coef/phi / lsqrtNf;
	}
	return res;
}

double complex GEnum(double tau, double x, double lsqrtNf){
	//tau=tau*(I+0.01);
	//return cexp((tau+I*x)*(tau+I*x)/(12*M_PI*csqrt(tau*tau+x*x)))/(2.*M_PI*(x*I - tau));
	
	//return sin(tau)*sin(tau)/(tau*tau)*exp(-x*x);// Benchmarking function
	if(tau==0. || x==0. || lsqrtNf==0.){
		fprintf(stderr, "Tried to evaluate GE for one of the arguments = 0.\n");
		exit(EXIT_FAILURE);
	}
	
	double sr=x*x+tau*tau;
    int g=nPolarGrid-1;
    
    double phi=atan(fabs(x/tau));
	
	while(phi>phiss[g][nphis[g]-2])
		g--;
    
    double complex *polarData;
    double rmax,phimin,phimax;
    int radii,angles;
    radii=nrs[g];
    angles=nphis[g];
    polarData=polarDatas[g];
    rmax=rss[g][radii-1] / lsqrtNf;
    phimin=phiss[g][0];
    phimax=phiss[g][angles-1];
    
    double complex e;
    if(rmax*rmax*cs*cs<sr)
        e=asymExponent(tau, x, sqrtNf);
    else{
        double r=sqrt(sr);
        
        double s=(radii-1)*r/rmax;
        double t=(angles-1)*(phi-(phimin))/(phimax-(phimin));
        
        int i=s; s-=i;
        int j=t; t-=j;
        
        //Bicubic interpolation
        double complex p;
        if(i-1<0 || i+2>=radii)
            p=(polarData[i+radii*(j)]*(1-s)+polarData[i+1+radii*(j)]*s)*(1-t) +
            (polarData[i+radii*(j+1)]*(1-s)+polarData[i+1+radii*(j+1)]*s)*t;
        else
            p=cubicInterp(
                          cubicInterp(polarData[i-1+radii*(j-1)],polarData[i+radii*(j-1)],polarData[i+1+radii*(j-1)],polarData[i+2+radii*(j-1)],s),
                          cubicInterp(polarData[i-1+radii*(j+0)],polarData[i+radii*(j+0)],polarData[i+1+radii*(j+0)],polarData[i+2+radii*(j+0)],s),
                          cubicInterp(polarData[i-1+radii*(j+1)],polarData[i+radii*(j+1)],polarData[i+1+radii*(j+1)],polarData[i+2+radii*(j+1)],s),
                          cubicInterp(polarData[i-1+radii*(j+2)],polarData[i+radii*(j+2)],polarData[i+1+radii*(j+2)],polarData[i+2+radii*(j+2)],s),t);
        
        e=r*p;
        if(x*tau<0.)e=conj(e);
    }
	//return /*cexp(e)*/ exp(-.005*(x*x+tau*tau)); //1./(2.*M_PI*(x*I - tau));
	return cexp(e)/(2.*M_PI*(x*I - tau));
	//return  1./(2.*M_PI*(x*I - tau));
}

void *calculateGs(void *rStartP)
{
	int iStart=*(int*)rStartP*n/nThreads;
	int iEnd=(*(int*)rStartP+1)*n/nThreads;
	
	int iSize=cn2+(iEnd-iStart-1)*nFactor;
	complex double *buffer=malloc(sizeof(complex double)*cn2*iSize);
	int bufferJ=0;	//next written line in buffer. If the kernel center 'rolls' in the buffer,
					//then we can avoid some read/writes
	int jj=0;//our place in Gs
	int operc=-1;
	for(int j=0;j<cn2+(n-1)*nFactor;j++){
		for(int i=0;i<iSize;i++){
			buffer[i+bufferJ*iSize]=GEnum( -L-cL2 + (iStart*nFactor + i)*cDelta, -L-cL2 + j*cDelta,sqrtNf);
		}
		bufferJ++;
		if(bufferJ>=cn2)bufferJ-=cn2;
		int fromStart=j+1-cn2;
		
		if(0<=fromStart  &&  fromStart%nFactor==0)
		{
			for(int ii=iStart;ii<iEnd;ii++){
				complex double acc=0;
				for(int cj=0;cj<cn2;cj++){
					int bj=(nFactor*jj+cj)%cn2;
					for(int ci=0;ci<cn2;ci++){
						int bi=nFactor*(ii-iStart)+ci;
						acc+=ker[ci+cn2*cj] * buffer[bi+bj*iSize];
					}
				}
				double sx=-1.+ii*2./(n-1);
				double sy=-1.+jj*2./(n-1);
				Gs[ii+jj*n]=acc*exp(-alpha3*(sx*sx+sy*sy));
			}
			jj++;
			if(*(int*)rStartP==0){
				int perc=100.*(jj)/(n-1);
				if(perc!=operc){
					printf("   %d%%   \r", perc );
					fflush(stdout);
					operc=perc;
				}
			}
		}
    }
    return NULL;
}

int main(int argc, char *argv[]){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int begin = tp.tv_sec * 1000 + tp.tv_usec / 1000;
	
    if(argc!=2){
        printf("Usage: cgs run_description\n");
        return 1;
    }
    FILE *fp;
	
	if(fftw_init_threads()==0)
	{
		fprintf(stderr, "Multithreaded fftw not available.\n");
		return 1;
	}

	printf(KGRN"Loading input file...\n"KRESET);
    if(! (fp = fopen(argv[1], "r")))
    {
        fprintf(stderr, "Unable to open input file \"%s\".\n",argv[1]);
        return 1;
    }
	fscanf(fp, "%d\n", &nPolarGrid);
	if(nPolarGrid>MAX_PG)
	{
		fprintf(stderr, "Can at most have 10 input files.\n");
		return 1;
	}
    char fileName[MAX_PG+1][1024];
    for(int g=0;g<nPolarGrid+1;g++)
        fscanf(fp, "%s\n", fileName[g]);
    fscanf(fp, "%d\n", &nThreads);
    fscanf(fp, "%lf\n", &Nf);
    sqrtNf=sqrt(Nf);
    fscanf(fp, "%lf\n", &L);
    fscanf(fp, "%d\n", &n);  n=1<<n;
	fclose(fp);
	
	alpha=3.7;//2.5; //UV aliasing suppressed by exp(-2.5^2) ~ .2%
	double alpha2=alpha; //Cut off kernel at exp(-2.5^2) ~ .2%
	alpha3=alpha; //Have window function be exp(-2.5^2) at spatial edges
	csigma=4*alpha*L/(n*M_PI);
	cL=alpha2*csigma;
	double pointsPerWidth=12.;
	cn=2*pointsPerWidth*cL/csigma;//51; //21 still seems to show some aliasing

	printf(KBLU" * L*sqrt(Nf) = %g\n"KRESET,L*sqrt(Nf));
	printf(KBLU" * 2*L*sqrt(Nf)/(n-1) = %g\n"KRESET,2*L*sqrt(Nf)/(n-1));
	printf(KBLU" * cn = %d\n"KRESET,cn);
	
	if(useAsymCoef){
		printf(KGRN"Loading asymptotic coefficents...\n"KRESET);
		if(! (fp = fopen("/Users/petter/asymCoef","r")))
		{
			fprintf(stderr, "Unable to open asymptotic coefficents file.\n");
			return 1;
		}
		fread(&acn,sizeof(int),1,fp);
		printf(KBLU" * Loading %d coefficients...\n"KRESET,acn);
		acs=malloc(sizeof(complex double)*acn);
		fread(acs,sizeof(complex double),acn,fp);
		fclose(fp);
		printf(KBLU" * Done\n"KRESET);
	}

    for(int g=0;g<nPolarGrid;g++){
        if(! (fp = fopen(fileName[g], "r")))
        {
            fprintf(stderr, "Unable to open polar grid file \"%s\"\n",fileName[g]);
            return 1;
        }
        printf(KGRN"Loading polar grid %d to memory...\n"KRESET,g);
        
        fscanf(fp, "%d\n", &nrs[g]);
        fscanf(fp, "%d\n", &nphis[g]);
        rss[g]=malloc(sizeof(double)*nrs[g]);
        phiss[g]=malloc(sizeof(double)*nphis[g]);
        for(int i=0;i<nrs[g];i++)
              fscanf(fp, "%lf\n", &rss[g][i]);
        for(int i=0;i<nphis[g];i++)
              fscanf(fp, "%lf\n", &phiss[g][i]);
        printf(KBLU" * %d radii between %g and %g\n"KRESET,nrs[g],rss[g][0],rss[g][nrs[g]-1]);
		printf(KBLU" * %d angles between %g and %g\n"KRESET,nphis[g],phiss[g][0],phiss[g][nphis[g]-1]);
        polarDatas[g]=malloc(sizeof(double complex)*nrs[g]*nphis[g]);
        
        for(int i=0;i<nrs[g];i++)
            for(int j=0;j<nphis[g];j++){
                double re,im;
                fscanf(fp, "%lf\n", &re);
                fscanf(fp, "%lf\n", &im);
                polarDatas[g][i+nrs[g]*j]=re+I*im;
            }
        fclose(fp);
		
		if(debug && useAsymCoef){
			double sqrtNf=sqrt(Nf);
			double r=rss[g][nrs[g]-1]/sqrtNf;
			printf(KCYN"polarGrid\tasym2\tasym3\n"KRESET);
			for(int j=0;j<nphis[g];j++){
				double phi=phiss[g][j];

				complex double e1=r*polarDatas[g][nrs[g]-1+nrs[g]*j];
				useAsymCoef=0;
				complex double e2=asymExponent(r*cos(phi), r*sin(phi), sqrtNf);
				useAsymCoef=1;
				complex double e3=asymExponent(r*cos(phi), r*sin(phi), sqrtNf);
				
				e1=e1-e2;
				e2=e3/e1;
				printf(KCYN"%8.8g + i%8.8g \t%8.8g + i%8.8g\n"KRESET, creal(e1),cimag(e1), creal(e2),cimag(e2));
				//printf(KCYN"%8.8g + i%8.8g \t%8.8g + i%8.8g \t%8.8g + i%8.8g\n"KRESET, creal(e1),cimag(e1), creal(e2),cimag(e2), creal(e3),cimag(e3));
			}
			
		}
    }
	
	printf(KGRN"Checking if polar grid sufficiently large...\n"KRESET);
	double maxErr=-1, maxVal=-1, phiWorst;
	for(int g=0;g<nPolarGrid;g++)
	for(int i=0;i<nphis[g];i++){
		double phi=phiss[g][i];
		if(g<nPolarGrid-1 && phi<phiss[g+1][nphis[g+1]-2])
			continue;
		double r=rss[g][nrs[g]-1]/sqrtNf;
		complex double e1=asymExponent(r*cos(phi), r*sin(phi), sqrtNf);
		complex double e2=r*polarDatas[g][nrs[g]-1+nrs[g]*i];
		double val=cabs(cexp(e2));
		double err=cabs(cexp(e1)-cexp(e2));
		if(maxErr<err){
			maxErr=err;
			phiWorst=phi;
		}
		if(maxVal<val)maxVal=val;
		
		//printf("%lg / %lg\n",err,val);
	}
	double maxRelErr=maxErr/maxVal;
	printf(KBLU" * Max relative difference of expansions is %lg%%, at phi=%g\n"KRESET,maxRelErr*100,phiWorst);
	if(0.01<maxRelErr){
		fprintf(stderr,"Error too large, increase polar grid size.\n");
		return 1;
	}
	
	printf(KGRN"Preparing sampling...\n"KRESET);
	nFactor=1+((int)(cn*L/cL/n))/2*2;  // (int) how much denser we sample for convolution than fft
	delta=2*L/((n-1)*nFactor); //about what delta would be
	cn2=cL/delta;// How many points we actually use for the kernel, supposed to be
	cn2=1+2*cn2;				// slightly larger than cn so points end up at points we can reuse.
	// We also want cn2 to be odd.
	cL2=(cn2-1)/2*delta; //How large our kernel actually is
	
	printf(KBLU" * Extra factor of samples needed for low-pass: %d*%d\n"KRESET,nFactor,nFactor);
	printf(KBLU" * Kernel size in samples: %d*%d\n"KRESET,cn2,cn2);
	
	int nTotal=cn2+(n-1)*nFactor; //How many samples we will actually have to calculate with GE per dimension
	if(nTotal/2*2!=nTotal)
	{
		fprintf(stderr, "Internal error, nTotal should be even.\n");
		return 1;
	}
	
	//precompute kernel
	double kerNorm=0;
	ker=malloc(sizeof(double)*cn2*cn2);
	cDelta=2*cL2/(cn2-1);
	
	for(int ci=0;ci<cn2;ci++)
		for(int cj=0;cj<cn2;cj++){
			double dx=-cL2 + ci*cDelta;
			double dy=-cL2 + cj*cDelta;
			ker[ci+cj*cn2]=exp(-(dx*dx+dy*dy)/(csigma*csigma));
			//if(ci!=cn2/2 && cj!=cn2/2)ker[ci+cj*cn2]=0;
			kerNorm+=ker[ci+cj*cn2];
		}
	for(int ci=0;ci<cn2*cn2;ci++) //normalize
		ker[ci]/=kerNorm;
	
	Delta=2*L/(n-1);
	

	printf(KGRN"Allocating memory to hold real space G...\n"KRESET);
	Gs=fftw_malloc(sizeof(fftw_complex)*n*n);
	Gso=fftw_malloc(sizeof(fftw_complex)*n*n);//maybe fftw can actually do it in-place for these sizes?
	if(!Gso){
		fprintf(stderr, "Not enough memory\n");
		return 1;
	}
	printf(KBLU" * Sucessfully allocated %d MB\n"KRESET,(int)(sizeof(double complex)*2*n*n/(1024*1024)));
	printf(KBLU" * Using this memory to find best fft scheme to use later\n"KRESET);
	fftw_plan_with_nthreads(nThreads);
	fftw_plan fftPlan = fftw_plan_dft_2d(n,n, Gs, Gso, FFTW_FORWARD, FFTW_ESTIMATE);
	
	printf(KGRN"Sampling G...\n"KRESET);
    printf(KBLU" * Starting %d threads...\n"KRESET, nThreads);
    int *startPoints=malloc(sizeof(int)*nThreads);;
    pthread_t *workers=malloc(sizeof(pthread_t)*nThreads);
    for(int i=0;i<nThreads;i++){
        startPoints[i]=i;
        if(pthread_create(&workers[i], NULL, calculateGs, &startPoints[i])) {
            fprintf(stderr, "Error creating thread\n");
            return 1;
        }
    }
    
    for(int i=0;i<nThreads;i++){
        printf(KYEL"              Waiting for thread %d to join...\r"KRESET, i);
        if(pthread_join(workers[i], NULL)) {
            fprintf(stderr, "Error joining thread\n");
            return 1;
        }
    }
	printf(KBLU" * Done                                             \n");
	
	printf(KGRN"Discrete Fourier transform...\n"KRESET);
	double omegaMax=M_PI*(n-1)*(n-1)/(2*L*n);
	double complex *phases=malloc(sizeof(double complex)*n);
	
	printf(KBLU" * Prefixing phases\n"KRESET);
	for(int i=0;i<n;i++)
		phases[i]=cexp(I*(omegaMax*L*2*i/(n-1)-L*omegaMax*.5));
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Gs[i+j*n]*=phases[i]*phases[j];
	
	printf(KBLU" * Running fftw3\n"KRESET);
	fftw_execute(fftPlan);
	//fftw_free(Gs);
	
	printf(KBLU" * Postfixing phases\n"KRESET);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			Gso[i+j*n]*=phases[i]*phases[j];
	
	printf(KBLU" * Downsampling, normalizing, and correcting for UV filter\n"KRESET);
	//fftw_destroy_plan(fftPlan);
	int stride=1<<2;
	int prunedN=n/2/stride;
	double deltaOmegaPruned=stride*2*omegaMax/(n-1);
	double omegaMaxPruned=deltaOmegaPruned*(prunedN-1)/2;
	double complex *GsPruned = malloc(sizeof(double complex)*prunedN*prunedN);
	if(!GsPruned){
		fprintf(stderr, "Not enough memory\n");
		return 1;
	}
	
	for(int i=0;i<prunedN;i++)
		for(int j=0;j<prunedN;j++){
			double omega=-omegaMaxPruned+i*deltaOmegaPruned;
			double kx=-omegaMaxPruned+j*deltaOmegaPruned;
			GsPruned[(prunedN-1-i)+j*prunedN]=( Gso[ (n/4 + stride/2 -1 +i*stride)  +  (n/4+  stride/2 -1 + j*stride)*n ]+
												Gso[ (n/4 + stride/2 -0 +i*stride)  +  (n/4+  stride/2 -1 + j*stride)*n ]+
												Gso[ (n/4 + stride/2 -1 +i*stride)  +  (n/4+  stride/2 -0 + j*stride)*n ]+
												Gso[ (n/4 + stride/2 -0 +i*stride)  +  (n/4+  stride/2 -0 + j*stride)*n ]  )/4*Delta*Delta * exp((omega*omega + kx*kx)*(csigma*csigma)/4);
			//GsPruned[i+j*prunedN]=Gs[ (n/4 + stride/2  +i*stride)  +  (n/4+ stride/2  + j*stride)*n ];// * exp((omega*omega + kx*kx)/(4*csigma*csigma));
		}
	//fftw_free(Gso);
	
    printf(KGRN"Saving result...\n"KRESET);
	
	char fn[200];
	snprintf(fn,200,"%s_%g_%d_%g",fileName[nPolarGrid],Nf,n,L);
    fp=fopen(fn,"w");
	fwrite(&Nf,sizeof(double),1,fp);
	fwrite(&L,sizeof(double),1,fp);
	fwrite(&n,sizeof(int),1,fp);
	fwrite(&csigma,sizeof(double),1,fp);
	fwrite(&cL,sizeof(double),1,fp);
	fwrite(&cn,sizeof(int),1,fp);

	
	fwrite(&prunedN,sizeof(int),1,fp);
	fwrite(&omegaMaxPruned,sizeof(double),1,fp);
	fwrite(GsPruned,sizeof(double complex),prunedN*prunedN,fp);

    fclose(fp);
    printf(KBLU" * Done\n"KRESET);
	gettimeofday(&tp, NULL);
	long int end = tp.tv_sec * 1000 + tp.tv_usec / 1000;
	printf(KYEL"Total wall time: %.1f seconds\n"KRESET,((double)(end - begin))/1000);

    
    return 0;
}
