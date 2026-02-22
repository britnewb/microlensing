//AREA TEST
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"

float AREA(float R_a,float R_b, float d)
{
	float A,B,C,Area;
	A = pow(R_b,2.0)*acos((d*d + R_b*R_b - R_a*R_a)/(2*d*R_b));
	B = pow(R_a,2.0)*acos((d*d + R_a*R_a - R_b*R_b)/(2*d*R_a));
	C = 0.5*pow((R_a - d + R_b)*(R_a + d - R_b)*(R_a+d+R_b)*(R_b + d - R_a),0.5);
	Area= A+B-C;
	printf("Area= %lf\td= %lf\tA=%lf\tB=%lf\tC=%lf\n",Area,d,A,B,C);
	return (Area);
}

int main()
{

	float R_a, R_b,d,A;
	R_a = 2.0;
	R_b = 3.0;
	d = 4.0;
	printf("Area of a is %lf\n",pi*R_a*R_a);
	int i;
	for (i=0;i<17;i++)
	{
		A = AREA(R_a,R_b,d);
		d = d-0.25;
	}
return (0);
}
