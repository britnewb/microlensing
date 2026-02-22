//compile with gcc -lm -o lastone lastone.c

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include <stdlib.h>

int main()
{
	double R_a;//Primary and secondary radii
	double A_a;//Area of the primary
	double r,theta;
	float midpoint;//midpoint of ring
	int j;
	double DENOMINATOR=0;
	double DENOMINATORA = 0;
	double DENOMINATORB = 0;
	double dr;
	double area=0;//area of segment
	float intensity;//intensity as a function of r,theta
	double dtheta = 2*pi/100000;
	int n =0;//counter
	double lastarea;//area of segment which isn't counted.
	double totalarea;//sum all area segments to verify it matches

	R_a = r_sun;
	A_a = pi*R_a*R_a;

	dr = R_a/100000;

	for(r=0;r<R_a;r += dr)
	{
	DENOMINATORA = 0;
		for (theta=0;theta<2*pi;theta +=dtheta)
		{
			DENOMINATORA += dtheta*1;
			if(r!=0)
			{
				area = (n*dr*n*dr - (n-1)*(n-1)*dr*dr)/(2/dtheta);
				midpoint = (dr*n + dr*(n+1))/2;
				//printf("area of section = %f\n",area);
			}
			if(r==0)
			{
				area = 0;
			}
			totalarea+=area;
		}	
		DENOMINATORB += r*DENOMINATORA*dr;
	n++;
	}

	lastarea = (n*dr*n*dr - (n-1)*(n-1)*dr*dr)/(2/dtheta);
	totalarea+=(lastarea*2*pi/dtheta);
	printf("lastarea = %f\n",lastarea);
	printf("Area of disc = %f\n",A_a);
	printf("ratio = %f\n",A_a/totalarea);
	printf("totalarea=%f\n",totalarea);
printf("DENOMINATORA= %f\n",DENOMINATORA);

	return (0);
}
