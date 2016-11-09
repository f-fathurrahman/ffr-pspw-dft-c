// eFeFeR (20910015), October 2011

void SortGVectors(int N, double *Gvec, int *indexG)
{
  int i, j, flag=1;
  int temp;

  for(i=1; i<=N; i++) {
    indexG[i-1] = i;
  }

  for(i=1; (i <= N) && flag; i++) {
    flag = 0;
    for (j=1; j<=(N-1); j++) {
      if(Gvec[indexG[j+1]-1] < Gvec[indexG[j]-1]) {
        temp = indexG[j];
        indexG[j] = indexG[j+1];
        indexG[j+1] = temp;
        flag = 1;               // indicates that a swap occurred.
      }
    }
  }
}
