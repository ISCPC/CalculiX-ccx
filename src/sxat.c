/*
 * sxat.c
 *
 *  Created on: Jun 17, 2021
 *      Author: uno
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "timelog.h"
#include "sxat.h"

/*
 * CalculiX solver for SXAT
 */

//typedef int32_t ITG;

static vesolver_handle_t hdl;

static inline void set_matrix(matrix_desc_t* desc, double *ad, double *au, double *adb, double *aub, const double sigma,
        ITG *icol, ITG *irow, const ITG nzs) {
    ITG i,j,k,l;

    ITG* pointers = desc->pointers;
    ITG* indice = desc->indice;
    double* value = desc->values;
    const ITG neq = desc->neq;

    k=desc->nnz;
    l=nzs;
    pointers[neq]=k;

    if (sigma==0.0) {
        for(i=neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                indice[--k]=irow[--l]-1;
                value[k]=au[l];
            }
            pointers[i]=--k;
            indice[k]=i;
            value[k]=ad[i];
        }
    } else {
        for(i=neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                indice[--k]=irow[--l]-1;
                value[k]=au[l]-(sigma)*aub[l];
            }
            pointers[i]=--k;
            indice[k]=i;
            value[k]=ad[i]-(sigma)*adb[i];
        }
    }
}

/*
 * symmetric matrix; the subdiagonal entries are stored
 *  column by column in au, the diagonal entries in ad;
 *  pardiso needs the entries row per row
 */
#define NNEW(a,b,c) a=(b *)calloc((c),sizeof(b))
#define RENEW(a,b,c) a=(b *)realloc((b *)(a),(c)*sizeof(b))
#define SFREE(a) free(a)

#define FORTRAN(func, parm)

static inline void set_matrix_symmetric(matrix_desc_t* desc, double *ad, double *au, double *adb, double *aub, const double sigma,
        ITG *icol, ITG *irow, const ITG nzs) {
    ITG* pointers = desc->pointers;
    ITG* icolpardiso = desc->indice;
    double* aupardiso = desc->values;
    const ITG neq = desc->neq;
    ITG i,j;

    ITG k = desc->nnz;
    ITG l=nzs;
    pointers[neq]=k+1;

    if(sigma==0.){
        for(i=neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                icolpardiso[--k]=irow[--l];
                aupardiso[k]=au[l];
            }
            pointers[i]=k--;
            icolpardiso[k]=i+1;
            aupardiso[k]=ad[i];
        }
    } else {
        for(i=neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                icolpardiso[--k]=irow[--l];
                aupardiso[k]=au[l]-sigma*aub[l];
            }
            pointers[i]=k--;
            icolpardiso[k]=i+1;
            aupardiso[k]=ad[i]-sigma*adb[i];
        }
    }
}

/*
 * off-diagonal terms  are stored column per
 * column from top to bottom in au;
 * diagonal terms are stored in ad
 */
static inline void set_matrix_asymmetric_3(matrix_desc_t* desc, double *ad, double *au, double *adb, double *aub, const double sigma,
        ITG *icol, ITG *irow, const ITG nzs) {
    const ITG neq = desc->neq;
    ITG* irowpardiso = NULL;
    ITG i,j;

    ITG ndim = desc->nnz;
    ITG* pointers = desc->pointers;
    ITG* icolpardiso = desc->indice;
    double* aupardiso = desc->values;

    NNEW(irowpardiso,ITG,ndim);	  

    ITG k=0;
    ITG k2=0;
    for(i=0;i<neq;i++){
        for(j=0;j<icol[i];j++){
            if(au[k]>1.e-12||au[k]<-1.e-12){
                icolpardiso[k2]=i+1;
                irowpardiso[k2]=irow[k];
                aupardiso[k2]=au[k];
                k2++;		  
            }
            k++;	      
        }	  
    }  

    /* diagonal terms */  
    for(i=0;i<neq;i++){
        icolpardiso[k2]=i+1;
        irowpardiso[k2]=i+1;
        aupardiso[k2]=ad[i];
        k2++;	  
    }
    ndim = k2;

    /* pardiso needs the entries row per row; so sorting is
       needed */ 
    ITG kflag=2;
    FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
                &ndim,&kflag));

    /* sorting each row */
    k=0;
    pointers[0]=1;
    for(i=0;i<neq;i++){
        ITG j=i+1;
        ITG kstart=k;
        do {
            ITG n;

            if(irowpardiso[k]!=j ){
                n=k-kstart;		  
                FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
                            &n,&kflag));
                pointers[i+1]=k+1;
                break;  
            }else{
                if(k+1==ndim){
                    n=k-kstart+1;	  
                    FORTRAN(isortid,(&icolpardiso[kstart],
                                &aupardiso[kstart],&n,&kflag));
                    break;	       
                }else{
                    k++;	       
                }  
            }
        } while(1);
    }
    pointers[neq]=ndim+1;
    SFREE(irowpardiso);
}

/*
 * lower triangular matrix is stored column by column in
 * au, followed by the upper triangular matrix row by row;
 * the diagonal terms are stored in ad
 */
static inline void set_matrix_asymmetric_1(matrix_desc_t* desc, double *ad, double *au, double *adb, double *aub, const double sigma,
        ITG *icol, ITG *irow, const ITG nzs, ITG *jq, const ITG nzs3) {
    const ITG neq = desc->neq;

    ITG* irowpardiso = NULL;
    ITG i,j;

    ITG ndim = desc->nnz;
    ITG* pointers = desc->pointers;
    ITG* icolpardiso = desc->indice;
    double* aupardiso = desc->values;

    /* reordering lower triangular matrix */
    NNEW(irowpardiso,ITG,ndim);

    ITG k=0;
    for(i=0;i<neq;i++){
        for(j=0;j<icol[i];j++){
            icolpardiso[k]=i+1;
            irowpardiso[k]=irow[k];
            aupardiso[k]=au[k];
            k++;
        }
    }

    /* pardiso needs the entries row per row; so sorting is
       needed */
    ITG kflag=2;
    FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
                &ndim,&kflag));

    /* sorting each row */
    k=0;
    pointers[0]=1;
    if(ndim>0){
        for(i=0;i<neq;i++){
            ITG j=i+1;
            ITG kstart=k;
            do{
                ITG n;
                /* end of row reached */
                if(irowpardiso[k]!=j){
                    n=k-kstart;
                    FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
                                &n,&kflag));
                    pointers[i+1]=k+1;
                    break;
                }else{
                    /* end of last row reached */
                    if(k+1==ndim){
                        n=k-kstart+1;
                        FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
                                    &n,&kflag));
                        break;
                    }else{
                        /* end of row not yet reached */
                        k++;
                    }
                }
            } while(1);
        }
    }
    pointers[neq]=ndim+1;
    SFREE(irowpardiso);

    /* composing the matrix: lower triangle + diagonal + upper triangle */
    ndim=neq+2*nzs;
    RENEW(icolpardiso,ITG,ndim);
    RENEW(aupardiso,double,ndim);
    k=ndim;
    for(i=neq-1;i>=0;i--){
        ITG l=k+1;
        for(j=jq[i+1]-1;j>=jq[i];j--){
            icolpardiso[--k]=irow[j-1];
            aupardiso[k]=au[j+nzs3-1];
        }
        icolpardiso[--k]=i+1;
        aupardiso[k]=ad[i];
        for(j=pointers[i+1]-1;j>=pointers[i];j--){
            icolpardiso[--k]=icolpardiso[j-1];
            aupardiso[k]=aupardiso[j-1];
        }
        pointers[i+1]=l;
    }
    pointers[0]=1;
}


int sxat_ve_init() {
    vesolver_init();
    hdl = vesolver_activate();
    return 0;
}

void sxat_ve_fini() {
    vesolver_deactivate(hdl);
    vesolver_finalize();
    return;
}

#if 0
void sxat_ve_factor(double *ad, double *au, double *adb, double *aub,
        double *sigma,ITG *icol, ITG *irow, ITG *neq, ITG *nzs) {
    vesolver_set_option();
    vesolver_set_matrix_csr();
}
static void sxat_ve_factor(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
        ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
        ITG *jq, ITG *nzs3){

    char *env;
    /*  char env1[32]; */
    ITG i,j,k,l;
    ITG ndim;
    ITG *irowpardiso=NULL,kflag,kstart,n,ifortran,
        lfortran,index,id,k2;
    double *b=NULL,*x=NULL;

    matrix_desc_t* desc;
    if(*symmetryflag==0){
        set_matrix_symmetric(desc, ad, au, adb, aub, *sigma, icol, irow, *nzs);
    }else{
        if(*inputformat==3){
            set_matrix_asymmetric_3(desc, ad, au, adb, aub, *sigma, icol, irow, *nzs);
        }else if(*inputformat==1){
            set_matrix_asymmetric_1(desc, ad, au, adb, aub, *sigma, icol, irow, *nzs);
        }
    }

    //FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
    //            pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
    //            b,x,&error));

    return;
}
#endif

void sxat_ve_factor(double *ad, double *au, double *adb, double *aub,
        const double sigma, ITG *icol, ITG *irow, const ITG neq, const ITG nzs,
	    const ITG symmetryflag, const ITG inputformat, ITG *jq, const ITG nzs3) {
    matrix_desc_t* desc;

    //vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_ITER_CG);
    vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_DIRECT_HS);

    if(symmetryflag==0) {
        printf(" Solving the system of equations using the symmetric VE solver\n\n");
        desc = vesolver_alloc_matrix(hdl, neq, neq+nzs,
            MATRIX_TYPE_CSC|MATRIX_TYPE_INDEX0|MATRIX_TYPE_SYMMETRIC);

        set_matrix_symmetric(desc, ad, au, adb, aub, sigma, icol, irow, nzs);
    } else {
        printf(" Solving the system of equations using the unsimmetric VE solver (format=%d)\n\n", inputformat);

        if(inputformat==3) {
            desc = vesolver_alloc_matrix(hdl, neq, neq+nzs,
                MATRIX_TYPE_CSR|MATRIX_TYPE_INDEX1|MATRIX_TYPE_ASYMMETRIC);
            set_matrix_asymmetric_3(desc, ad, au, adb, aub, sigma, icol, irow, nzs);
        } else if(inputformat==1) {
            desc = vesolver_alloc_matrix(hdl, neq, nzs,
                MATRIX_TYPE_CSR|MATRIX_TYPE_INDEX1|MATRIX_TYPE_ASYMMETRIC);
            set_matrix_asymmetric_1(desc, ad, au, adb, aub, sigma, icol, irow, nzs, jq, nzs3);
        }
    }
}

void sxat_ve_solve(double *b) {
    double res = 1.0e-10;

    vesolver_solve_sync(hdl, b, b, res);
}

void sxat_ve_cleanup() {
    vesolver_free_matrix(hdl);
}

#if 1
void sxat_ve_main(double *ad, double *au, double *adb, double *aub, const double sigma,
        double *b, ITG *icol, ITG *irow, const ITG neq, const ITG nzs,
	    const ITG symmetryflag, const ITG inputformat, ITG *jq, const ITG nzs3) {
    sxat_ve_factor(ad, au, adb, aub, sigma, icol, irow, neq, nzs,
        symmetryflag, inputformat, jq, nzs3);
    sxat_ve_solve(b);
    sxat_ve_cleanup();
}
#else
void sxat_ve_main(double *ad, double *au, double *adb, double *aub, const double sigma,
        double *b, ITG *icol, ITG *irow, const ITG neq, const ITG nzs,
	    const ITG symmetryflag, const ITG inputformat, ITG *jq, const ITG nzs3) {
    double res = 1.0e-10;
    matrix_desc_t* desc;

    //vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_ITER_CG);
    vesolver_set_option(hdl, VESOLVER_OPTION_SOLVER, VESOLVER_DIRECT_HS);

    if(symmetryflag==0) {
        printf(" Solving the system of equations using the symmetric VE solver\n\n");
        desc = vesolver_alloc_matrix(hdl, neq, neq+nzs,
            MATRIX_TYPE_CSC|MATRIX_TYPE_INDEX0|MATRIX_TYPE_SYMMETRIC);

        set_matrix_symmetric(desc, ad, au, adb, aub, sigma, icol, irow, nzs);
    } else {
        printf(" Solving the system of equations using the unsimmetric VE solver (format=%d)\n\n", inputformat);

        if(inputformat==3) {
            desc = vesolver_alloc_matrix(hdl, neq, neq+nzs,
                MATRIX_TYPE_CSR|MATRIX_TYPE_INDEX1|MATRIX_TYPE_ASYMMETRIC);
            set_matrix_asymmetric_3(desc, ad, au, adb, aub, sigma, icol, irow, nzs);
        } else if(inputformat==1) {
            desc = vesolver_alloc_matrix(hdl, neq, nzs,
                MATRIX_TYPE_CSR|MATRIX_TYPE_INDEX1|MATRIX_TYPE_ASYMMETRIC);
            set_matrix_asymmetric_1(desc, ad, au, adb, aub, sigma, icol, irow, nzs, jq, nzs3);
        }
    }
    //set_matrix(desc, ad, au, adb, aub, *sigma, icol, irow, *nzs);

    vesolver_set_matrix(hdl, desc);
    vesolver_solve_sync(hdl, b, b, res);
    vesolver_free_matrix(hdl);
}
#endif

#if 0
static void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3){

  char *env;
/*  char env1[32]; */
  ITG i,j,k,l,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
      msglvl=0,error=0,*irowpardiso=NULL,kflag,kstart,n,ifortran,
      lfortran,index,id,k2;
  ITG ndim,nthread,nthread_v;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0){
      printf(" Factoring the system of equations using the symmetric pardiso solver\n");
  }else{
      printf(" Factoring the system of equations using the unsymmetric pardiso solver\n");
  }

  iparm[0]=0;
  iparm[1]=3;
/* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
   must be done once  */
  if (nthread_mkl == 0) {
    nthread=1;
    env=getenv("MKL_NUM_THREADS");
    if(env) {
      nthread=atoi(env);}
    else {
      env=getenv("OMP_NUM_THREADS");
      if(env) {nthread=atoi(env);}
    }
    env=getenv("CCX_NPROC_EQUATION_SOLVER");
    if(env) {
      nthread_v=atoi(env);
      if (nthread_v <= nthread) {nthread=nthread_v;}
    }
    if (nthread < 1) {nthread=1;}
    sprintf(envMKL,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);  
    putenv(envMKL);
    nthread_mkl=nthread;
  }
    
  printf(" number of threads =% d\n\n",nthread_mkl);

  for(i=0;i<64;i++){pt[i]=0;}

  if(*symmetryflag==0){

      /* symmetric matrix; the subdiagonal entries are stored
         column by column in au, the diagonal entries in ad;
         pardiso needs the entries row per row */      

      mtype=-2;
      
      ndim=*neq+*nzs;
      
      NNEW(pointers,ITG,*neq+1);
      NNEW(icolpardiso,ITG,ndim);
      NNEW(aupardiso,double,ndim);
      
      k=ndim;
      l=*nzs;
      
      if(*sigma==0.){
	  pointers[*neq]=ndim+1;
	  for(i=*neq-1;i>=0;--i){
	      for(j=0;j<icol[i];++j){
		  icolpardiso[--k]=irow[--l];
		  aupardiso[k]=au[l];
	      }
	      pointers[i]=k--;
	      icolpardiso[k]=i+1;
	      aupardiso[k]=ad[i];
	  }
      }
      else{
	  pointers[*neq]=ndim+1;
	  for(i=*neq-1;i>=0;--i){
	      for(j=0;j<icol[i];++j){
		  icolpardiso[--k]=irow[--l];
		  aupardiso[k]=au[l]-*sigma*aub[l];
	      }
	      pointers[i]=k--;
	      icolpardiso[k]=i+1;
	      aupardiso[k]=ad[i]-*sigma*adb[i];
	  }
      }
  }else{
//      mtype=11;
      mtype=1;

      if(*inputformat==3){

          /* off-diagonal terms  are stored column per
             column from top to bottom in au;
             diagonal terms are stored in ad  */

	  ndim=*neq+*nzs;
	  NNEW(pointers,ITG,*neq+1);
	  NNEW(irowpardiso,ITG,ndim);	  
	  NNEW(icolpardiso,ITG,ndim);
	  NNEW(aupardiso,double,ndim);
	  
	  k=0;
	  k2=0;
	  for(i=0;i<*neq;i++){
	      for(j=0;j<icol[i];j++){
		  if(au[k]>1.e-12||au[k]<-1.e-12){
		      icolpardiso[k2]=i+1;
		      irowpardiso[k2]=irow[k];
		      aupardiso[k2]=au[k];
		      k2++;		  
		  }
		  k++;	      
	      }	  
	  }  
	  /* diagonal terms */  
	  for(i=0;i<*neq;i++){
	      icolpardiso[k2]=i+1;
	      irowpardiso[k2]=i+1;
	      aupardiso[k2]=ad[i];
	      k2++;	  
	  }
	  ndim=k2;
	  
	  /* pardiso needs the entries row per row; so sorting is
	     needed */ 
	  
	  kflag=2;
	  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			    &ndim,&kflag));
	  
	  /* sorting each row */
	  
	  k=0;
	  pointers[0]=1;
	  for(i=0;i<*neq;i++){
	      j=i+1;
	      kstart=k;
	      do{
		  if(irowpardiso[k]!=j ){
		      n=k-kstart;		  
		      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				       &n,&kflag));
		      pointers[i+1]=k+1;
		      break;  
		  }else{
		      if(k+1==ndim){
			  n=k-kstart+1;	  
			  FORTRAN(isortid,(&icolpardiso[kstart],
                                  &aupardiso[kstart],&n,&kflag));
			  break;	       
		      }else{
			  k++;	       
		      }  
		  }
	      }while(1);
	  }
	  pointers[*neq]=ndim+1;
	  SFREE(irowpardiso);

      }else if(*inputformat==1){
	  
          /* lower triangular matrix is stored column by column in
             au, followed by the upper triangular matrix row by row;
             the diagonal terms are stored in ad */

          /* reordering lower triangular matrix */

	  ndim=*nzs;
	  NNEW(pointers,ITG,*neq+1);
	  NNEW(irowpardiso,ITG,ndim);
	  NNEW(icolpardiso,ITG,ndim);
	  NNEW(aupardiso,double,ndim);
	  
	  k=0;
	  for(i=0;i<*neq;i++){
	      for(j=0;j<icol[i];j++){
		  icolpardiso[k]=i+1;
		  irowpardiso[k]=irow[k];
		  aupardiso[k]=au[k];
		  k++;
	      }
	  }
	  
	  /* pardiso needs the entries row per row; so sorting is
	     needed */
	  
	  kflag=2;
	  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
	  &ndim,&kflag));
	  
	  /* sorting each row */
	  
	  k=0;
	  pointers[0]=1;
	  if(ndim>0){
	      for(i=0;i<*neq;i++){
		  j=i+1;
		  kstart=k;
		  do{

	  	      /* end of row reached */

		      if(irowpardiso[k]!=j){
			  n=k-kstart;
			  FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					   &n,&kflag));
			  pointers[i+1]=k+1;
			  break;
		      }else{

		          /* end of last row reached */

			  if(k+1==ndim){
			      n=k-kstart+1;
			      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					       &n,&kflag));
			      break;
			  }else{

			      /* end of row not yet reached */

			      k++;
			  }
		      }
		  }while(1);
	      }
	  }
	  pointers[*neq]=ndim+1;
	  SFREE(irowpardiso);

	  /* composing the matrix: lower triangle + diagonal + upper triangle */

	  ndim=*neq+2**nzs;
	  RENEW(icolpardiso,ITG,ndim);
	  RENEW(aupardiso,double,ndim);
	  k=ndim;
	  for(i=*neq-1;i>=0;i--){
	      l=k+1;
	      for(j=jq[i+1]-1;j>=jq[i];j--){
		  icolpardiso[--k]=irow[j-1];
		  aupardiso[k]=au[j+*nzs3-1];
	      }
	      icolpardiso[--k]=i+1;
	      aupardiso[k]=ad[i];
	      for(j=pointers[i+1]-1;j>=pointers[i];j--){
		  icolpardiso[--k]=icolpardiso[j-1];
		  aupardiso[k]=aupardiso[j-1];
	      }
	      pointers[i+1]=l;
	  }
	  pointers[0]=1;
      }
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  return;
}
#endif

