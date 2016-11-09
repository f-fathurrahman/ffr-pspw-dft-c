// eFeFeR (20910015), October 2011

#include "common_pspw.h"

void GkVectorsGen()
{
  int ik, ig, igp;
  double gx, gy, gz, t;

  for(ik=1; ik<=NKPT; ik++) {
    NGW[ik-1] = 0;
  }

  for(ik=1; ik<=NKPT; ik++) {
    for(ig=1; ig<=NG; ig++) {
      gx = GG[IDX2F(ig,1,NGX+1)] + XK[IDX2F(1,ik,3)];
      gy = GG[IDX2F(ig,2,NGX+1)] + XK[IDX2F(2,ik,3)];
      gz = GG[IDX2F(ig,3,NGX+1)] + XK[IDX2F(3,ik,3)];
      t = gx*gx + gy*gy + gz*gz;
      if(t <= GCUTW) {
        NGW[ik-1] = NGW[ik-1] + 1;
        igp = NGW[ik-1];
        if(igp > NGWX) { 
          printf("ERROR in GkVectorsGen ngw(ik) > NGWX\n");
          printf("NGW(%d) = %d, NGWX=%d\n",ik,NGW[ik-1],NGWX);
          abort();
        }
        IGK[IDX2F(igp,ik,NGWX)] = ig;
        XKG[IDX2F(igp,ik,NGWX)] = t;
        N123[IDX2F(igp,ik,NGWX)] = N1[ig-1] + NR1*( (N2[ig-1]-1) + NR2*(N3[ig-1]-1) );
        GGK[IDX3F(1,igp,ik,3,NGWX)] = GG[IDX2F(ig,1,NGX+1)];
        GGK[IDX3F(2,igp,ik,3,NGWX)] = GG[IDX2F(ig,2,NGX+1)];
        GGK[IDX3F(3,igp,ik,3,NGWX)] = GG[IDX2F(ig,3,NGX+1)];
      }
    }
  }

}

