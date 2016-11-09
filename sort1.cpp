// eFeFeR (20910015), October 2011
//
// Converted naively from FHI98MD subroutine in order.f
//

void sort1(int N, double *arrin, int *index)
{
  int i, j, l, ir, indext;
  double q;

  for(i=1; i<=N; i++) {
    index[i-1] = i;
  }

  l = N/2 + 1;
  ir = N;

LABEL10:
  if(l > 1) {
    l = l - 1;
    indext = index[l-1];
    q = arrin[indext-1];
  } else {
    indext = index[ir-1];
    q = arrin[indext-1];
    index[ir-1] = index[0];
    ir = ir - 1;
    if(ir == 1) {
      index[0] = indext;
      goto LABEL150;
    }
  }
  i = l;
  j = l + l;

LABEL20:
  if(j <= ir) {
    if(j < ir) {
      if(arrin[index[j-1]-1] < arrin[index[j]-1]) j = j + 1;
    }
    if(q < arrin[index[j-1]-1]) {
      index[i-1] = index[j-1];
      i = j;
      j = j + j;
    } else {
      j = ir + 1;
    }
    goto LABEL20;
  }
  index[i-1] = indext;
  goto LABEL10;

LABEL150: return;

}

