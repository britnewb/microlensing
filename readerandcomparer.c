//Program to read in a data file and compare the values with the values generated herein.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include <stdlib.h>

int main(int argc, char **argv)
{
	FILE *infile1;//KOI74 data
	FILE *infile2;//generated data
	FILE *outfile;//chi2output.dat
	char *mode1 = "rt";
	char *mode2 = "wt";
	infile1 = fopen("KOI74.dat",mode1);
	infile2 = fopen("generated.dat",mode1);
	outfile = fopen("chi2output.dat",mode2);
	if (infile1 ==NULL)
		{
			fprintf(stderr,"Can't open KOI74.dat\n");
			exit (1);
		}
	if (infile2 ==NULL)
		{
			fprintf(stderr,"Can't open generated.dat\n");
			exit (1);
		}
	int i,j;
	double gendate;
	double genflux;
	double readdate;//date read from infile
	double readflux;//flux read from infile
	double currentreaddate;
	double currentreadflux;
	readdate = 55100.55998687;
	gendate = 55100.0195312500;
	char line1[17];
	char line2[17];
	char date[17];
	char *tmp;
	j=0;
	int num1,num2;
	num1=0;
	num2=0;
	int permanum;
	double difference;//between gendate and readdate
	while((fgets(line2,17,infile2)!=NULL)&&num2<=55)
	{
		sscanf (line2, "%lf", &gendate);
		//printf ("gendate = %.10lf\t", gendate);
		fgets(line2,17,infile2);
		sscanf (line2, "%lf", &genflux);
		//printf ("genflux = %.10lf\n", genflux);
		if(num2>=27)
		{
			while((fgets(line1, 17, infile1) != NULL))
			{
				sscanf (line1, "%lf", &readdate);
				difference=(gendate-readdate);
				fgets(line1, 17, infile1);//reads flux
				sscanf (line1, "%lf", &readflux);
				printf("readflux=%lf\n",readflux);
				if(fabs(difference)<=0.02)
				{
			 	fprintf(outfile,"gendate=%.10lf\tgenflux=%.10lf\treaddate = %.10lf\treadflux=%.10lf\n",gendate,genflux,readdate,readflux);
				}
				fgets(line1, 17, infile1);//for the error (3rd column)
				if(fabs(difference)<=0.02)
				{break;}
				num1++;			
			}
			num1=0;//reset num1 for the next iteration
			rewind(infile1);//starts at top of infile1
		}
   		num2++;
	}
   	fclose(infile1);  /* close the file prior to exiting the routine */
	fclose(infile2);
	fclose(outfile);
	return (0);
}
