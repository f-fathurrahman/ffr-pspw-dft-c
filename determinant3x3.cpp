// eFeFeR (20910015), October 2011

double determinant3x3(double A1[3], double A2[3], double A3[3])
{ 
  return A1[0]*( A2[1]*A3[2] - A2[2]*A3[1] ) -
         A1[1]*( A2[0]*A3[2] - A2[2]*A3[0] ) +
         A1[2]*( A2[0]*A3[1] - A2[1]*A3[0] );
}

