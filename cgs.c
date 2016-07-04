#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <pthread.h>

//0=whole phi range, small r range
//1=small phi range, large r range
double complex *polarDatas[2],*Gs;
int nphis[2],nrs[2],n,cn,nThreads;
double *phiss[2], *rss[2],Nf,sqrtNf,L,cL,alpha,csigma;
double cs=1.;

double complex cubicInterp(double complex p0,double complex p1,double complex p2,double complex p3, double x) {
    return p1 + 0.5 * x*(p2 - p0 + x*(2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + x*(3.0*(p1 - p2) + p3 - p0)));
}

double complex asymExponent(double tau, double x, double sqrtNf){
	double c1=0.0243615448342451073975190249505015144289430415363799596395;
	return (  3/(36*M_PI*M_PI)  -  c1*pow(fabs(x)*sqrtNf,1./3)*cpow(1+I*tau/x,-2./3)  ) / sqrtNf;
}

double complex GEnum(double tau, double x, double sqrtNf){
    double sr=x*x+tau*tau;
    int g;
    
    double phi=atan(fabs(x/tau));
    
    if(phi<phiss[1][nphis[1]-2])
        g=1;
    else
        g=0;
    
    double complex *polarData;
    double rmax,phimin,phimax;
    int radii,angles;
    radii=nrs[g];
    angles=nphis[g];
    polarData=polarDatas[g];
    rmax=rss[g][radii-1];
    phimin=phiss[g][0];
    phimax=phiss[g][angles-1];
    
    double complex e;
    if(rmax*rmax*cs*cs<sr)
        e=asymExponent(tau, x, sqrtNf);
    else{
        double r=sqrt(sr);
        
        double s=(radii-1)*r*sqrtNf/rmax;
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
	return cexp (e)/(2.*M_PI*(x*I - tau));
}

void *calculateGs(void *rStartP)
{
	double cDelta=2*cL/(cn-1);
	
	double kernNorm=0;
	for(int ci=0;ci<cn;ci++)
	for(int cj=0;cj<cn;cj++){
		double dx=-cL + ci*cDelta;
		double dy=-cL + cj*cDelta;
		kernNorm+=exp(-(dx*dx+dy*dy)/(csigma*csigma));
	}
	kernNorm=1/kernNorm;
	
    double Delta=2*L/(n-1);
    int operc=-1;
    for(int i=*(int*)rStartP;i<n;i+=nThreads)
    for(int j=0;j<n;j++){
        complex double acc=0;
        for(int ci=0;ci<cn;ci++)
		for(int cj=0;cj<cn;cj++){
			double dx=-cL + ci*cDelta;
			double dy=-cL + cj*cDelta;
			acc+=exp(-(dx*dx+dy*dy)/(csigma*csigma))*GEnum( -L + i*Delta+dx, -L + j*Delta+dy,sqrtNf);
		}
        
        Gs[i+j*n]=acc*kernNorm;
        if(*(int*)rStartP==0){
            int perc=100.*(i*n+j)/(n*n-1);
            if(perc!=operc){
                printf("   %d%%   \r", perc );
                fflush(stdout);
                operc=perc;
            }
        }
    }
    return NULL;
}

int main(int argc, char *argv[]){
    if(argc!=2){
        printf("Usage: cgs run_description\n");
        return 1;
    }
    FILE *fp;
    
    if(! (fp = fopen(argv[1], "r")))
    {
        fprintf(stderr, "Unable to open input file \"%s\".\n",argv[1]);
        return 1;
    }
    char fileName[3][255];
    for(int g=0;g<3;g++)
        fscanf(fp, "%s\n", fileName[g]);
    fscanf(fp, "%d\n", &nThreads);
    fscanf(fp, "%lf\n", &Nf);
    sqrtNf=sqrt(Nf);
    fscanf(fp, "%lf\n", &L);
    fscanf(fp, "%d\n", &n);  n=1<<n;
	//fscanf(fp, "%lf\n", &cL);
	//fscanf(fp, "%d\n", &cn);
	//fscanf(fp, "%lf\n", &alpha);
	cn=71; //21 still seems to show some aliasing, odd seems better
	alpha=2.; //UV aliasing suppresed by exp(-2.5^2) ~ .2%
	double alpha2=alpha; //Cut off kernel at exp(-2.5^2) ~ .2%
	csigma=4*alpha*L/(n*M_PI);
	cL=alpha2*csigma;

	
    fclose(fp);

    for(int g=0;g<2;g++){
        if(! (fp = fopen(fileName[g], "r")))
        {
            fprintf(stderr, "Unable to open polar grid file \"%s\"\n",fileName[g]);
            return 1;
        }
        printf("Loading polar grid %d to memory...\n",g);
        
        fscanf(fp, "%d\n", &nrs[g]);
        fscanf(fp, "%d\n", &nphis[g]);
        rss[g]=malloc(sizeof(double)*nrs[g]);
        phiss[g]=malloc(sizeof(double)*nphis[g]);
        for(int i=0;i<nrs[g];i++)
              fscanf(fp, "%lf\n", &rss[g][i]);
        for(int i=0;i<nphis[g];i++)
              fscanf(fp, "%lf\n", &phiss[g][i]);
        printf(" * %d radii between %g and %g\n",nrs[g],rss[g][0],rss[g][nrs[g]-1]);
		printf(" * %d angles between %g and %g\n",nphis[g],phiss[g][0],phiss[g][nphis[g]-1]);
        polarDatas[g]=malloc(sizeof(double complex)*nrs[g]*nphis[g]);
        
        for(int i=0;i<nrs[g];i++)
            for(int j=0;j<nphis[g];j++){
                double re,im;
                fscanf(fp, "%lf\n", &re);
                fscanf(fp, "%lf\n", &im);
                polarDatas[g][i+nrs[g]*j]=re+I*im;
            }
        fclose(fp);
    }
	
	printf("Checking if polar grid sufficiently large...\n");
	double maxErr=-1, maxVal=-1, phiWorst;
	for(int g=0;g<2;g++)
	for(int i=0;i<nphis[g];i++){
		double phi=phiss[g][i];
		if(g<1 && phi<phiss[g+1][nphis[g+1]-2])
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
	printf(" * Max relative difference of expansions is %lg%%\n",maxRelErr*100);
	if(0.01<maxRelErr){
		fprintf(stderr,"Error too large at phi=%g, increase polar grid size.\n",phiWorst);
		return 1;
	}
	
    printf("Sampling G...\n");
    Gs=malloc(sizeof(double complex)*n*n);

    printf(" * Starting %d threads...\n", nThreads);
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
        printf("              Waiting for thread %d to join...\r", i);
        if(pthread_join(workers[i], NULL)) {
            fprintf(stderr, "Error joining thread\n");
            return 1;
        }
    }
    printf(" * Done\n");
    printf("Saving result...\n");
	
	char fn[200];
	snprintf(fn,200,"%s_%g_%d_%g",fileName[2],Nf,n,L);
    fp=fopen(fn,"w");
	fwrite(&Nf,sizeof(double),1,fp);
	fwrite(&L,sizeof(double),1,fp);
	fwrite(&n,sizeof(int),1,fp);
	fwrite(&csigma,sizeof(double),1,fp);
	fwrite(&cL,sizeof(double),1,fp);
	fwrite(&cn,sizeof(int),1,fp);
	fwrite(Gs,sizeof(double complex),n*n,fp);

    fclose(fp);
    printf(" * Done\n");

    
    return 0;
}