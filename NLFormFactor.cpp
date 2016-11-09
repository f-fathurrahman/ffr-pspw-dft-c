// eFeFeR (20910015), November 2011

#include "common_pspw.h"

void NLFormFactor()
{
  // Constants
  const double FPIBO = 4.0*M_PI/OMEGA;
  const double SQRT3 = sqrt(3.0);

  double ds;
  // Local variables
  int i, j, is, ir, i_l, m_l, ik, ig;
  double *vref=NULL, *fun=NULL, *fun1=NULL;
  double rr, rr2, dv, q1, arg, fac, ag, t, cosx, cosy, cosz;

  vref = (double*)malloc(MAXR*NSP*sizeof(double)); 
  for(i=0; i<MAXR*NSP; i++) vref[i] = 0.0;
  fun = (double*)malloc(MAXR*NSP*sizeof(double));
  for(i=0; i<MAXR*NSP; i++) fun[i] = 0.0;
  fun1 = (double*)malloc(MAXR*sizeof(double));
  for(i=0; i<MAXR; i++) fun1[i] = 0.0;

  for(is=1; is<=NSP; is++) {
    for(ir=1; ir<=MMAX[is-1]; ir++) {
      vref[IDX2F(ir,is,MAXR)] = VION[IDX3F(ir,is,L_LOC[is-1],MAXR,NSP)];
    }
  }

  for(is=1; is<=NSP; is++) {
    m_l = 0;
    for(i_l=1; i_l<=L_MAX[is-1]; i_l++) { // sum over AM
      if(i_l != L_LOC[is-1]) { // don't include local potential
        for(ir=1; ir<=MMAX[is-1]; ir++) {
          rr = R[IDX2F(ir,is,MAXR)];
          rr2 = rr*rr;
          dv = VION[IDX3F(ir,is,i_l,MAXR,NSP)] - vref[IDX2F(ir,is,MAXR)];
          fun[IDX2F(ir,is,MAXR)] = dv*PSI[IDX3F(ir,is,i_l,MAXR,NSP)]*rr2;
          fun1[ir-1] = dv*SQUARE( PSI[IDX3F(ir,is,i_l,MAXR,NSP)] ) * rr;
        }
        q1 = simpson(MMAX[is-1], fun1, CLOG[is-1]);
        if(i_l < 4) {
          for(j=1; j<=2*i_l-1; j++) {
            WNL[IDX2F(is,m_l+j,NSP)] = (2.0*i_l - 1.0)*FPIBO/q1;
          }
        } else {
          printf("ERROR: not implemented i_l = %d\n", i_l);
          abort();
        }
        //
        for(ik=1; ik<=NKPT; ik++) {
        for(ig=1; ig<=NGW[ik-1]; ig++) {
          t = sqrt( XKG[IDX2F(ig,ik,NGWX)] );
          arg = t*TPIBA;
          //
          if(i_l == 1) {
            if(t < 1.0e-4) {
              for(ir=1; ir<=MMAX[is-1]; ir++) {
                fun1[ir-1] = fun[IDX2F(ir,is,MAXR)];
              }
            } else {
              for(ir=1; ir<=MMAX[is-1]; ir++) {
                fac = arg*R[IDX2F(ir,is,MAXR)];
                fac = sin(fac)/fac;
                fun1[ir-1] = fac*fun[IDX2F(ir,is,MAXR)];
              }
            }
            q1 = simpson(MMAX[is-1], fun1, CLOG[is-1]);
            PKG[IDX4F(ig,is,ik,m_l+1,NGWX,NSP,NKPT)] = q1;
          }
          //
          else if(i_l == 2) {
            if(t < 1.e-4) {
              for(j=1; j<=3; j++) PKG[IDX4F(ig,is,ik,m_l+j,NGWX,NSP,NKPT)] = 0.0;
            } else {
              for(ir=1; ir<=MMAX[is-1]; ir++) {
                fac = arg*R[IDX2F(ir,is,MAXR)];
                fac = ( sin(fac)/fac - cos(fac) )/fac;
                fun1[ir-1] = fac*fun[IDX2F(ir,is,MAXR)];
              }
              q1 = simpson(MMAX[is-1], fun1, CLOG[is-1]);
              PKG[IDX4F(ig,is,ik,m_l+1,NGWX,NSP,NKPT)] = q1*( GGK[IDX3F(1,ig,ik,3,NGWX)] + XK[IDX2F(1,ik,3)] )/t;
              PKG[IDX4F(ig,is,ik,m_l+2,NGWX,NSP,NKPT)] = q1*( GGK[IDX3F(2,ig,ik,3,NGWX)] + XK[IDX2F(2,ik,3)] )/t;
              PKG[IDX4F(ig,is,ik,m_l+3,NGWX,NSP,NKPT)] = q1*( GGK[IDX3F(3,ig,ik,3,NGWX)] + XK[IDX2F(3,ik,3)] )/t;
            }
          }
          //
          else if(i_l == 3) {
            if(t < 1.e-4) {
              for(j=1; j<=5; j++) PKG[IDX4F(ig,is,ik,m_l+j,NGWX,NSP,NKPT)] = 0.0;
            } else {
              for(ir=1; ir<=MMAX[is-1]; ir++) {
                ag = arg*R[IDX2F(ir,is,MAXR)];
                fac = (3.0/SQUARE(ag)-1.0)*sin(ag) - 3.0/ag*cos(ag);
                fac = fac/ag;
                fun1[ir-1] = fac*fun[IDX2F(ir,is,MAXR)];
              }
              q1 = simpson(MMAX[is-1], fun1, CLOG[is-1]);
              cosx = ( GGK[IDX3F(1,ig,ik,3,NGWX)] + XK[IDX2F(1,ik,3)] )/t;
              cosy = ( GGK[IDX3F(2,ig,ik,3,NGWX)] + XK[IDX2F(2,ik,3)] )/t;
              cosz = ( GGK[IDX3F(3,ig,ik,3,NGWX)] + XK[IDX2F(3,ik,3)] )/t;
              PKG[IDX4F(ig,is,ik,m_l+1,NGWX,NSP,NKPT)] = q1*0.5*(3.0*SQUARE(cosz)-1.0);
              PKG[IDX4F(ig,is,ik,m_l+2,NGWX,NSP,NKPT)] = q1*SQRT3*cosz*cosx;
              PKG[IDX4F(ig,is,ik,m_l+3,NGWX,NSP,NKPT)] = q1*SQRT3*cosz*cosy;;
              PKG[IDX4F(ig,is,ik,m_l+4,NGWX,NSP,NKPT)] = q1*SQRT3*cosx*cosy;
              PKG[IDX4F(ig,is,ik,m_l+5,NGWX,NSP,NKPT)] = q1*SQRT3/2.0*(SQUARE(cosx) - SQUARE(cosy));
            }
          }
          //
          else {
            printf("ERROR: G-local potential for this angular momentum is not implemented\n");
            abort();
          }
          for(j=1; j<=2*i_l-1; j++) {
            PKG_A[IDX4F(m_l+j,ig,is,ik,NLMAX,NGWX,NSP)] = PKG[IDX4F(ig,is,ik,m_l+j,NGWX,NSP,NKPT)];
          }
        } // ig
        } // ik
        m_l = m_l + 2*i_l - 1;
      } // ENDIF l != L_LOC
    } // i_l
    // Display the prefactors
    /*printf("WNL = \n");
    for(i=1; i<=SQUARE(L_MAX[is-1])-2*L_LOC[is-1] + 1; i++) {
      printf("is = %d, %d %f\n", is, i, WNL[IDX2F(is,i,NSP)]);
    }*/
  } // is

  // Free memory
  free(vref); vref=NULL;
  free(fun); fun=NULL;
  free(fun1); fun1=NULL;
}

