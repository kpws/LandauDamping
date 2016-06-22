#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include <complex.h>
#include <pthread.h>

const unsigned long inverseMaxError=1000000000000;

double rmaxD, *pg;
int nr, nk, nphi, nThreads;
mpf_t *coefR,*coefI;

void *calculateRadii(void *rStartP)
{
    int operc=-1;
    mpf_t r,prodR,prodI;
    mpf_init(r);
    mpf_init(prodR);
    mpf_init(prodI);
    mpf_t *rPowers=malloc(sizeof(mpf_t)*nk);
    mpf_init_set_ui(rPowers[0],1ul);
    for(int i=*(int*)rStartP;i<nr;i+=nThreads){
        mpf_set_d(r,rmaxD*i/(nr-1));
        for(int j=1;j<nk;j++){
            mpf_init(rPowers[j]);
            mpf_mul(rPowers[j],rPowers[j-1],r);
        } 
        mpf_t accR,accI;
        mpf_init(accR);
        mpf_init(accI);
        for(int j=0;j<nphi;j++){
            mpf_set_ui(accR,0ul);
            mpf_set_ui(accI,0ul);
            for(int k=0;k<nk;k++)
            {
                mpf_mul(prodR, rPowers[k], coefR[j+k*nphi]);
                mpf_mul(prodI, rPowers[k], coefI[j+k*nphi]);
                mpf_add(accR, accR, prodR);
                mpf_add(accI, accI, prodI);
                pg[i + j*nr] = mpf_get_d(accR) + I*mpf_get_d(accI);
            }
            if(*(int*)rStartP==0){
                int perc=100.*(i*nphi+j)/(nr*nphi-1);
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
    if(argc!=4){
        printf("Usage: pgs number_of_threads input_file output_file\n");
        return 1;
    }
    FILE *fp;
    printf("Loading coefficients...\n");
    if(! (fp = fopen(argv[2], "r")))
    {
        fprintf(stderr, "Unable to open input file.\n");
        return 1;
    }
    int accDec;

    fscanf (fp, "%d", &accDec);
    fscanf (fp, "%d", &nphi);
    fscanf (fp, "%d", &nk);

    printf(" * Accuracy: %d\n",accDec);
    printf(" * Number of angles: %d\n",nphi);
    printf(" * Number of coefficients: %d\n",nk);
    
    double *phis=malloc(sizeof(double)*nphi);
    for(int i=0;i<nphi;i++)
        fscanf(fp, "%lf", &phis[i]);
    
    int prec=accDec*log(10.)/log(2.)+30;
    mpf_set_default_prec(prec*3);
    coefR=malloc(sizeof(mpf_t)*nphi*nk);
    coefI=malloc(sizeof(mpf_t)*nphi*nk);

    if(!coefR || !coefI){
        fprintf(stderr, "Not enough memory.\n");
        return 1;
    }

    nk=100;
    int precDec=accDec+30, operc=-1;
    for(int j=0;j<nk;j++){
        for(int i=0;i<nphi;i++)
        {
            mpf_init(coefR[i + nphi*j]); 
            //gmp_fscanf(fp, "%Ff", coefR[i + nphi*j]);
            mpf_init(coefI[i + nphi*j]); 
            //gmp_fscanf(fp, "%Ff", coefI[i + nphi*j]);
            mpf_inp_str(coefR[i + nphi*j],fp,10);
            mpf_inp_str(coefI[i + nphi*j],fp,10);
        }
        int perc=100.*j/(nk-1);
        if(perc!=operc){
            printf("   %d%%     \r", perc );
            fflush(stdout);
            operc=perc;
        }
    }
    fclose(fp); 
    printf(" * Done\n");
    printf("Calculating max radius...\n");
    mpf_t rmax, r2, i2;
    mpf_init(rmax);
    mpf_init(r2);
    mpf_init(i2);
    for(int i=0;i<nphi;i++){
        mpf_mul(r2, coefR[i+nphi*(nk-1)], coefR[i+nphi*(nk-1)]);
        mpf_mul(i2, coefI[i+nphi*(nk-1)], coefI[i+nphi*(nk-1)]);
        mpf_add(r2,r2,i2);
        if(0<mpf_cmp(r2,rmax))
            mpf_set(rmax, r2);
    }
    mpf_sqrt(rmax,rmax);
    mpf_mul_ui(rmax,rmax,inverseMaxError);
    long exp;
    mpf_get_d_2exp(&exp,rmax);
    int nsqrt=log2(-exp);
    for(int i=0;i<nsqrt;i++)
        mpf_sqrt(rmax,rmax);
    
    rmaxD=pow(mpf_get_d(rmax),-pow(2.,nsqrt)/(nk-1.));
    printf(" * Maximum radius: %1.5f\n",rmaxD);
    
    nr=nk+10;
    
    printf("Series evaluation...\n");
    
    pg=malloc(sizeof(complex double)*nr*nphi);
    if(!pg){
        fprintf(stderr,"Not enough memory\n");
        return 1;
    }
    
    if(sscanf(argv[1],"%d",&nThreads)!=1 || !(0<nThreads)){
        fprintf(stderr,"Couldn't read number of threads\n");
        return 1;
    }

    printf(" * Starting %d threads...\n", nThreads);
    int *startPoints=malloc(sizeof(int)*nThreads);;
    pthread_t *workers=malloc(sizeof(pthread_t)*nThreads);
    for(int i=0;i<nThreads;i++){
        startPoints[i]=i;
        if(pthread_create(&workers[i], NULL, calculateRadii, &startPoints[i])) {
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
    
    fp=fopen(argv[3],"w");
    for(int i=0;i<nr;i++)
    for(int j=0;j<nphi;j++)
        fprintf(fp,"%1.15f\n",creal(pg[i+j*nr]));
    fclose(fp);
    printf(" * Done\n");
    //mpf_clear(r);  don't deallocate since we just quit anyways
    return 0;
}
