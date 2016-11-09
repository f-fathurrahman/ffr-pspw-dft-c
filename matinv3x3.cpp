#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

void matinv3x3(double *HM, double *HI)
{
  // Local variables
  double D11, D12, D13, D22, D23, D33, D21, D31, D32, DH;
  
  D11 = HM[IDX2F(2,2,3)]*HM[IDX2F(3,3,3)] - HM[IDX2F(2,3,3)]*HM[IDX2F(3,2,3)];
  D12 = HM[IDX2F(2,3,3)]*HM[IDX2F(3,1,3)] - HM[IDX2F(2,1,3)]*HM[IDX2F(3,3,3)];
  D13 = HM[IDX2F(2,1,3)]*HM[IDX2F(3,2,3)] - HM[IDX2F(3,1,3)]*HM[IDX2F(2,2,3)];
  D22 = HM[IDX2F(1,1,3)]*HM[IDX2F(3,3,3)] - HM[IDX2F(1,3,3)]*HM[IDX2F(3,1,3)];
  D23 = HM[IDX2F(1,1,3)]*HM[IDX2F(3,2,3)] - HM[IDX2F(3,1,3)]*HM[IDX2F(1,2,3)];
  D33 = HM[IDX2F(1,1,3)]*HM[IDX2F(2,2,3)] - HM[IDX2F(1,2,3)]*HM[IDX2F(2,1,3)];
  D21 = HM[IDX2F(3,2,3)]*HM[IDX2F(1,3,3)] - HM[IDX2F(1,2,3)]*HM[IDX2F(3,3,3)];
  D31 = HM[IDX2F(1,2,3)]*HM[IDX2F(2,3,3)] - HM[IDX2F(2,2,3)]*HM[IDX2F(1,3,3)];
  D32 = HM[IDX2F(1,3,3)]*HM[IDX2F(2,1,3)] - HM[IDX2F(1,1,3)]*HM[IDX2F(2,3,3)];
  DH  = HM[IDX2F(1,1,3)]*D11 + HM[IDX2F(1,2,3)]*D12 + HM[IDX2F(1,3,3)]*D13;
  D23 = -D23;
  HI[IDX2F(1,1,3)] = D11/DH;
  HI[IDX2F(2,2,3)] = D22/DH;
  HI[IDX2F(3,3,3)] = D33/DH;
  HI[IDX2F(1,2,3)] = D21/DH;
  HI[IDX2F(1,3,3)] = D31/DH;
  HI[IDX2F(2,3,3)] = D32/DH;
  HI[IDX2F(2,1,3)] = D12/DH;
  HI[IDX2F(3,1,3)] = D13/DH;
  HI[IDX2F(3,2,3)] = D23/DH;
}


