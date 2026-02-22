//Program to generate a data file 

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include <stdlib.h>

int main(int argc, char **argv)
{
	FILE *outfile;
	char *mode1 = "wt";
	outfile = fopen("generated.dat",mode1);
	if (outfile ==NULL)
		{
			fprintf(stderr,"Can't open generated.dat\n");
			exit (1);
		}
	int i,j;
	float genval;//value generated for comparison with text value
	genval = 55100.0;
	double readdate;//date read from infile
	double readflux;//flux read from infile
	char line[17];
	char date[17];
	char *tmp;
	j=0;
	int num;
	num=0;
	for(i=0;i<1860;i++)
	{
		genval += 0.0214845;
		fprintf(outfile,"%.10lf\t0.9998456\n",genval);
	}
   	fclose(outfile);  /* close the file prior to exiting the routine */
	return (0);
}
