#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include <complex.h>

unsigned long inverseMaxError=1000000000000;

int main(int argc, char *argv[]){
    if(argc!=2){
        printf("Please provide an input file.\n");
        return 1;
    }
    FILE *fp;
    printf("Loading coefficients...\n");
    if(! (fp = fopen(argv[1], "r")))
    {
        printf("Unable to open input file.\n");
        return 1;
    }
    int accDec,nphi,nk;

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
    mpf_t *coefR=malloc(sizeof(mpf_t)*nphi*nk);
    mpf_t *coefI=malloc(sizeof(mpf_t)*nphi*nk);

    if(!coefR || !coefI)
    {
        printf("Not enough memory.\n");
        return 1;
    }
    //nk=200;
    int precDec=accDec+30, operc=-1;
    for(int j=0;j<nk;j++){
        for(int i=0;i<nphi;i++)
        {
            mpf_init(coefR[i + nphi*j]); 
            gmp_fscanf(fp, "%Ff", coefR[i + nphi*j]);
            mpf_init(coefI[i + nphi*j]); 
            gmp_fscanf(fp, "%Ff", coefI[i + nphi*j]);
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
    printf("Generating radius power matrix...\n");
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
    //gmp_printf ("%.*Ff \n", precDec, rmax);
    
    double rmaxD=pow(mpf_get_d(rmax),-pow(2.,nsqrt)/(nk-1.));
    printf(" * Maximum radius: %1.5f\n",rmaxD);
    //printf("%d,   %1.5f\n",nsqrt,-pow(2.,nsqrt)/(nk-1.));
    //printf("hm %1.50f",mpf_get_d(rmax));
    
    int nr=nk+10;
    mpf_t *rs=malloc(sizeof(mpf_t)*nr);
    mpf_t *rmat=malloc(sizeof(mpf_t)*nr*nk);
    for(int i=0;i<nr;i++){
        mpf_init_set_d(rs[i],rmaxD*i/(nr-1));
        mpf_init_set_ui(rmat[i+0*nr],1ul);
        for(int j=1;j<nk;j++){
            mpf_init(rmat[i+j*nr]);
            mpf_mul(rmat[i+j*nr],rmat[i+(j-1)*nr],rs[i]);
        }
    }
    
    printf(" * Done.\n");
    printf("Matrix multiply...\n");

    complex double *pg=malloc(sizeof(complex double)*nr*nphi);
    if(!pg){
        printf("Not enough memory.\n");
        return 1;
    }
    mpf_t accR,accI;
    mpf_init(accR);
    mpf_init(accI);
    for(int i=0;i<nr;i++)
    for(int j=0;j<nphi;j++){
        mpf_set_ui(accR,0ul);
        mpf_set_ui(accI,0ul);
        for(int k=0;k<nk;k++)
        {
            mpf_mul(r2, rmat[i+k*nr], coefR[j+k*nphi]);
            mpf_mul(i2, rmat[i+k*nr], coefI[j+k*nphi]);
            mpf_add(accR, accR, r2);
            mpf_add(accI, accI, i2);
            pg[i + j*nr] = mpf_get_d(accR) + I*mpf_get_d(accI);
        }
        int perc=100.*(i*nphi+j)/(nr*nphi-1);
        if(perc!=operc){
            printf("   %d%%     \r", perc );
            fflush(stdout);
            operc=perc;
        }
    }
    
    printf(" * Done\n");
    printf("Saving result...\n");
    
    fp=fopen("smallPhiGrid","w");
    for(int i=0;i<nr;i++)
    for(int j=0;j<nphi;j++)
        fprintf(fp,"%1.15f\n",creal(pg[i+j*nr]));
    fclose(fp);
    printf(" * Done\n");
    //mpf_clear(r);  don't deallocate since we just quit anyways
    return 0;
}
