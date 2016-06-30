#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <pthread.h>

//0=whole phi range, small r range
//1=small phi range, large r range
double complex *polarDatas[2];
int nphis[2],nrs[2],n,nThreads;
double *phiss[2], *rss[2],sqrtNf,L;
double cs=1.;

double complex cubicInterp(double complex p0,double complex p1,double complex p2,double complex p3, double x) {
    return p1 + 0.5 * x*(p2 - p0 + x*(2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + x*(3.0*(p1 - p2) + p3 - p0)));
}

double complex GEnum(double tau, double x, double sqrtNf)
{
    double c1=0.0243615448342451073975190249505015144289430415363799596395;
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
        e=(  3/(36*M_PI*M_PI)  -  c1*pow(fabs(x)*sqrtNf,1./3)*cpow(1+I*tau/x,-2./3)  ) / sqrtNf;
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
    fscanf(fp, "%lf\n", &sqrtNf);
    sqrtNf=sqrt(sqrtNf);
    fscanf(fp, "%lf\n", &L);
    fscanf(fp, "%d\n", &n);
    fclose(fp);

    for(int g=0;g<2;g++){
        if(! (fp = fopen(fileName[g], "r")))
        {
            fprintf(stderr, "Unable to open polar grid file \"%s\"\n",fileName[g]);
            return 1;
        }
        printf("Loading polar grid %d to memory...\n",g);
        int nr,nphi;
        
        fscanf(fp, "%d\n", &nrs[g]);
        fscanf(fp, "%d\n", &nphis[g]);
        rss[g]=malloc(sizeof(double)*nrs[g]);
        phiss[g]=malloc(sizeof(double)*nphis[g]);
        for(int i=0;i<nrs[g];i++)
              fscanf(fp, "%lf\n", &rss[g][i]);
        for(int i=0;i<nphis[g];i++)
              fscanf(fp, "%lf\n", &phiss[g][i]);
        
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
    printf("Sampling G...\n");
    
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
        printf("              Waiting for thread %d to join...\r", i+1);
        if(pthread_join(workers[i], NULL)) {
            fprintf(stderr, "Error joining thread\n");
            return 1;
        }
    }
    printf(" * Done\n");
    printf("Saving result...\n");
    
    fp=fopen(fileName[3],"w");
    for(int i=0;i<n*n;i++)
        for(int j=0;j<n;j++){
            fprintf(fp,"%1.15f\n",creal(Gs[i]));
            fprintf(fp,"%1.15f\n",cimag(Gs[i]));
        }
    fclose(fp);
    printf(" * Done\n");

    
    return 0;
}