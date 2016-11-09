// evaluation of three simple random number generators
// Written in Turbo C, modified for Pelles C

#include <stdio.h> // printf(), getchar(), puts()
     
// three different custom random number generators
int rand1(int lim);
int rand2(int lim);
int rand3(int lim);
     
int main(void)
{
	int k, rnd1, rnd2, rnd3, sum1, sum2, sum3;
	
	sum1 = sum2 = sum3 = 0;
	puts("Random numbers between 1 and 44:\n");
	printf("%15s%10s %10s %10s","Custom functions","rnd1()","rnd2()","rnd3()");
	for (k = 0; k < 40; k++) {
		rnd1 = rand1(44); // method one
		rnd2 = rand2(44); // method two
		rnd3 = rand3(44); // method three
		printf("\n%25d %10d %10d",rnd1,rnd2,rnd3);
		sum1 += rnd1;
		sum2 += rnd2;
		sum3 += rnd3;
	}
	
	printf("\n\n%15s%10d %10d %10d (ideal 22)\n"," Average:",sum1/k,sum2/k,sum3/k);
	
	return 0;
}
     
//
// returns random integer from 1 to lim
//
int rand1(int lim)
{
	static long a = 100001;
     
  a = (a * 125) % 2796203;
	return ((a % lim) + 1);
}
     
//
// returns random integer from 1 to lim (Gerhard's generator)
//
int rand2(int lim)
{
	static long a = 1; // could be made the seed value
     
  a = (a * 32719 + 3) % 32749;
	return ((a % lim) + 1);
}
     
//
// returns random integer from 1 to lim (Bill's generator)
//
int rand3(int lim)
{
	static long a = 3;
	
	a = (((a * 214013L + 2531011L) >> 16) & 32767);
	
	return ((a % lim) + 1);
}
     
