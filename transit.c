/**************NEW PROGRAM TO NUMERICALLY INTEGRATE FLUX***************/
//Just transit, no limb darkening, no microlensing
//store data in transit.dat
//compile with gcc -lm -o transit transit.c
//run as ./transit -h for the help information



#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h" //defines pi, G, au, m_sun, SQR(x)
#include <stdlib.h>

float AREA(float R_a,float R_b, float d)
{
	float A,B,C,Area;
	A = R_b*R_b*acos((d*d + R_b*R_b - R_a*R_a)/(2*d*R_b));
	B = R_a*R_a*acos((d*d + R_a*R_a - R_b*R_b)/(2*d*R_a));
	C = 0.5*pow((R_a - d + R_b)*(R_a + d - R_b)*(R_a+d+R_b)*(R_b + d - R_a),0.5);
	Area= A+B-C;
	return (Area);
}

int main(int argc, char *argv[])
{
	float A,Z,L_a,L_b; // A? Z? Luminosity of the primary and secondary
	float M_a,M_b;//mass of the primary and secondary
	float phi,theta;//zenith and polar angles respectively
	float psi;//angle of rotation in its own circular orbit.
	float mind,maxd;//minimum value of z == d =0
	float d;//projected distance between the two centres
	float R_a,R_b;//Primary and secondary radii
	float Flux; //Flux as a function of time
	float A_a,A_b;//Area of the first and second
	float v, w; //actual and observed rotational velocities
	float i;//inclination in radians
	float r;//radius of orbit
	int j,t,ttime;
	int Phase=1;
	int mint=0;
	float x,y,z;//new coordinates;
	int ti;//time at intersection
/**************************
Default setting:*/
v = 1; //actual rotational velocity in m/s?
w = 1; //observed rotational velocity
Flux = 1;
R_b = 1;
R_a = 1;
i = 0;
ttime = 10000;
psi=pi;
//should I set a default value for L_a and L_b?
/***********************/

	for (j = 0; j<argc; j++)
	{
   		switch(argv[j][1])
		{
			case 'a': // radius of a
   			sscanf(argv[j+1],"%f",&R_a);
   			break;

   			case 'b': // radius of b called r_s in SG
    			sscanf(argv[j+1],"%f",&R_b);
    			break;

			case 'c'://mass of a
			sscanf(argv[j+1], "%f",&M_a);
			break;
	
			case 'd'://mass of b
			sscanf(argv[j+1],"%f",&M_b);
			break;

  			case 't': // total time
  			sscanf(argv[j+1], "%d", &ttime);
  			break;
				  
			case 'r': //radius of orbit called D_LS in SG
			sscanf(argv[j+1], "%f", &r);
			break;
				   
			case 'l'://luminosity of the primary
			sscanf(argv[j+1], "%f", &L_a);
			break;
				   
			case 'm'://luminosity of the secondary
			sscanf(argv[j+1], "%f", &L_b);
			break;

			case 'v'://rotational velocity
			sscanf(argv[j+1], "%f", &v);
			break;

			case 'i'://inclination
			sscanf(argv[j+1], "%f", &i);
			break;
			
			case 'y'://min time
			sscanf(argv[j+1], "%d", &mint);
			break;
				   
			default:
			break;
   				
			case 'h': // HELP! 
			printf("Requires r,a,b,t,l,m,v,i\n"); //Can I make one large printf? 
			printf("Must be in solar radii, solar luminosities, solar masses, separations in AU, velocity in m/s\n");
  			fprintf(stderr, "Quick Helps: \n");
   			fprintf(stderr, "\n");
   			fprintf(stderr, " -r radius of orbit in AU\n");
   			fprintf(stderr, " -a radius_a in r_sun \n");
   			fprintf(stderr, " -b radius_b in r_sun \n");
			fprintf(stderr, " -c mass of a in m_sun\n");
			fprintf(stderr, " -d mass of c in m_sun\n");
   			fprintf(stderr, " -t total time in seconds \n");
			fprintf(stderr, " -l luminosity of the primary in l_sun\n");
			fprintf(stderr, " -m luminosity of the secondary in l_sun\n");
			fprintf(stderr, " -v actual rotational velocity in m/s\n");
			fprintf(stderr, " -i inclination in radians\n");
	//		fprintf(stderr, " -f distance from source to us in AU\n"); //Not implemented yet
   			fprintf(stderr, "\n");
   			break;
		}
	}

/******************UNIT CONVERSION TO SOLAR UNITS******************************************************/
R_a = R_a*r_sun;
R_b = R_b*r_sun;
M_a = M_a*m_sun;
M_b = M_b*m_sun;
r = au*r; //this is a recursive definition
L_a = L_a*l_sun;
L_b = L_b*l_sun;

FILE *outfile;
char *mode= "w";
outfile = fopen("transit.dat",mode); //Prepare the file transit.dat for writing output
if (outfile ==NULL)
{
	fprintf(stderr,"Can't open output.txt\n");
	exit (1);
}
	A_a = pi*R_a*R_a; //Calculate the area of a from the radius
	A_b = pi*R_b*R_b;
	const float firstz = z;//first value given of d.
	const float firstFlux = L_a/A_a;// + L_b/A_b;

psi = mint*v/r + pi; //Figure out the internal angle of the orbit given the start time & speed
printf("psi  = %f\n",psi);
if((psi>=pi)&&(psi<3*pi/2))//I think Phase 1 is when the transiting object is approaching the transit
	{Phase = 1;}
if((psi>=pi*3/2)&&(psi<2*pi)) //Phase 2 would be while it's leaving transit
	{Phase = 2;}
if((psi>=pi*2)&&(psi<5*pi/2))
	{Phase = 3;}
if((psi>=5*pi/2)&&(psi<3*pi))
	{Phase = 4;}
printf("Intersection at %f seconds\n",pi*r/(v*2)); //What's the difference between these two?
printf("Intersection is imminent when i = %f\n",asin((R_b + R_a)/r)); 
ti = pi*r/(2*v);
	for(t=mint;t<ttime;t++)
	{	
		psi = t*v/r + pi;
		x = -r*cos(i)*sin(psi);
		y = r*cos(psi); //projected coordinates
		z = tan(i)*x;//plane of the orbit
		d = sqrt(z*z + y*y);//R_S
		if(x>=0)
		{
			if (R_a + R_b <= d)//separated
			{
				Flux = 1;
				fprintf(outfile,"%d\t%f\t%d\tseparated\n",(t - ti),Flux,Phase);
			}
			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
				A = AREA(R_a,R_b,d); // Calculate the area of the current intersection
				Z = (A_a - A)/(A_a*A_a);
				Flux = L_a*Z/firstFlux;
				//Flux = (L_b/A_b + L_a*Z)/firstFlux;
				fprintf(outfile,"%d\t%.12f\t%d\tintersecting\n",(t - ti),Flux,Phase);
			}	
			if ((d+R_b)<R_a)//inside
			{
				Z = (A_a - A_b)/(A_a*A_a);
				Flux = L_a*Z/firstFlux;
				//Flux = (L_b/A_b + L_a*Z)/firstFlux;
				fprintf(outfile,"%d\t%.12f\t%d\tinside\n",(t - ti),Flux,Phase);
				if (((psi>=pi*3/2)&&(psi<2*pi))||(y>=0))
				{
					Phase = 2;
				}
				if ((y<-r*cos(2*v/r))||(y>r*cos(2*v/r)))
				{
					Phase = 3;
				}
			}
		}
		if(x<0)
		{
			if ((d+R_b)<R_a)//inside
			{
			Z = (A_a - A_b)/(A_a*A_a);
			Flux = ((L_b/A_b + L_a*Z)/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tinside\n",(t - ti),Flux,Phase);
			}
			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
			A = AREA(R_a,R_b,d);
			Z = (A_a - A)/(A_a*A_a);
			Flux = ((L_b/A_b + L_a*Z)/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tintersecting\t%f\n",(t - ti),Flux,Phase);
			}
			if (R_a + R_b <= d)//separated
			{
			Flux = 1;
			fprintf(outfile,"%d\t%f\t%d\tseparated\n",(t - ti),Flux,Phase);
			}	
			if(((psi>=5*pi/2)&&(psi<3*pi))||(y<=0))
			{
				Phase = 4;
			}
			if((psi>=pi)&&(psi<3*pi/2))
			{
				Phase = 1;
			}
		}
	}
	fclose(outfile);
	
	return (0);
}
