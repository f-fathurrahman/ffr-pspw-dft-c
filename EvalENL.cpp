// eFeFeR (20910015), December 2011

#include "common_pspw_cuda.h"

void EvalENL()
{
  // Local
  int is, ia, ik, ig, i_lm, ist, lm_end;
  double complex sf;
  double complex *t_eigr_pkg=NULL;
  double enl_m;
  char trans;

  // Allocate memory
  t_eigr_pkg = (double complex*)malloc(NGWX*NLMAX*sizeof(double complex));
  for(int i=0; i<NGWX*NLMAX; i++) t_eigr_pkg[i] = Z_ZERO;

  // Reset ENL
  ENL = 0.0;

  for(is=1; is<=NSP; is++) {
    for(ia=1; ia<=NA[is-1]; ia++) {
      for(ik=1; ik<=NKPT; ik++) {
        lm_end = L_MAX[is-1]*L_MAX[is-1] - 2*L_LOC[is-1] + 1;
        for(i_lm=1; i_lm<=lm_end; i_lm++) {
          for(ig=1; ig<=NGW[ik-1]; ig++) {
            t_eigr_pkg[IDX2F(ig,i_lm,NGWX)] = PKG_A[IDX4F(i_lm,ig,is,ik,NLMAX,NGWX,NSP)]*
              conj(EIGR[IDX4F(ig,ia,is,ik,NGWX,NAX,NSP)]);
          }
          trans='T';
          zgemv_(&trans,&NGW[ik-1],&NX,&Z_ONE,&C0[IDX3F(1,1,ik,NGWX,NX)],&NGWX,
              &t_eigr_pkg[IDX2F(1,i_lm,NGWX)],&I_ONE, &Z_ZERO,
              &FNL[IDX5F(1,ik,is,ia,i_lm,NX,NKPT,NSP,NAX)],&I_ONE);
          enl_m = 0.0;
          for(ist=1; ist<=NX; ist++) {
            sf = FNL[IDX5F(ist,ik,is,ia,i_lm,NX,NKPT,NSP,NAX)];
            //enl_m = enl_m + FOCC[IDX2F(ist,ik,NX)]*(sf.x*sf.x + sf.y*sf.y);
            enl_m = enl_m + FOCC[IDX2F(ist,ik,NX)]*SQUARE(cabs(sf));
          }
          // No symmetry
          ENL = ENL + WKPT[ik-1]*enl_m*WNL[IDX2F(is,i_lm,NSP)];
        } // i_lm
      } // ik
    } // ia
  } // is

  // Free memory
  free(t_eigr_pkg); t_eigr_pkg = NULL;
}

