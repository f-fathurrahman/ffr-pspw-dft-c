#include <stdio.h>
#include <stdlib.h>

/* Find out the prime factors
 * of a given number and print
 * them on the screen */
void factorize(int n)
{
  int d = 2;

  if(n < 2) return;
  
  printf("Prime factors of '%d': ", n);
  /* while the factor being tested
   * is lower than the number to factorize */
  while(d < n) {
    /* if valid prime factor */
    if(n % d == 0) {
      printf("%d x ", d);
      n /= d;
    }
    /* else: invalid prime factor */
    else {
      if(d == 2) d = 3;
      else d += 2;
    }
  }
     
  /* print last prime factor */
  printf("%d\n", d);
}

int main(int argc, char **argv)
{
  int number;

  if(argc < 2) {
    printf("ERROR: Need exactly one argument\n");
    return -1;
  }
  sscanf(argv[1],"%d",&number);
  factorize(number);

  return 0;
}


