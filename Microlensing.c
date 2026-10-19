/**************NEW PROGRAM TO NUMERICALLY INTEGRATE FLUX***************/
//Function SG is based on the Sahu & Gillilan paper from 2003 called "Near-Field Microlensing and its Effects on Stellar Transit Observations by Kepler" in the Astrophysical Journal, 584:1042-1052

//Written by Britta Hansen in 2012 while working at the Canadian Institute for Theoretical Astrophysics under the supervision of Dr. Kelsey Hoffman
//compile with 
//gcc -lm -o Microlensing Microlensing.c
//
//Try running with test parameters:
//Microlensing -a 1 -b 0.00010051789 -c 1 -d 0.0009546 -v 5500000 -r 5 -l 1 -m 0.000000000000009 -y 213450 -t 213799 -i 0.0000001
//
//This parameter set gives a rough calculation of Jupiter transiting in front of the Sun



//D_LS = Distance from the lens to the star
//D_S = Distance to the source (the Sun)
//D_L = Distance from the lens to the observer (lens being a compact object causing microlensing such as a white dwarf)
















#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include <stdlib.h>

double SG(double R_E,double M_b,double R_a, double R_b, double d, double y,double z) //calculates the limb darkening from SG paper
{
	double r, 		//radius of the     in solar units
	theta, 			//angle of       in radians
	ydist,zdist, 		//
	yp, 			//
	ym, 			//
	NUMERATORA =0, 		//
	NUMERATORB =0, 		//
	Ap, 			//
	Am, 			//
	dr = R_a/10000, 	//
	dtheta = 2*pi/10000, 	//
	l,			//
	face,			//==rsin(theta)
	c1, c2, c3, c4,		//
	mu, I,rnorm,		//
	ratio;			//
				//
				//

	c1 = 0.5905;//0.701;	//G0V: 0.701;	//0.5118; 	//KOI81: 0.5905
	c2 = 0.0486;//0.149;	//0.149;	//0.0525;	//0.0486
	c3 = -0.0476;//0.277;	//0.277;	//0.4590;	//-0.0476
	c4 = 0.0074;//-0.297;//-0.297;	//-0.2727;	//0.0074

	for(r=0;r<R_a;r += dr)
	{
		NUMERATORA = 0;
		for (theta=0;theta<2*pi;theta +=dtheta)
		{
			Ap = 0;
			Am = 0;
			face = r*fabs(sin(theta));
			if(y>=0)
			{
				y = fabs(y);
				ydist = fabs(y - r*cos(theta));
			}
			if(y<0)
			{
				y = fabs(y);
				ydist = fabs(r*cos(theta) + y);
			}
			zdist = fabs(face - z);//determine zdist
			l = sqrt(zdist*zdist + ydist*ydist);//determine l
			yp = 0.5*(l+sqrt(l*l + 4*R_E*R_E));
			ym = 0.5*(l-sqrt(l*l + 4*R_E*R_E));
			ratio = l/R_E;
			if((yp>R_b)&&(l!=0))
			{
       				Ap = (pow(ratio,2) + 2.0)/((2.0*ratio)*sqrt(pow(ratio,2) + 4.0)) + 0.5;
			}
			if((fabs(ym)>R_b)&&(l!=0))
			{
				Am = (pow(ratio,2) + 2.0)/((2.0*ratio)*sqrt(pow(ratio,2) + 4.0)) - 0.5;
			}
			if(l==0)
			{
				printf("here\n");
			}	
			rnorm = r/R_a;
			mu = sqrt(1 - rnorm*rnorm);
			I = 1.0 - c1*(1 - pow(mu,0.5)) - c2*(1 - mu) - c3*(1 - pow(mu,1.5)) - c4*rnorm*rnorm;
			NUMERATORA += (Ap+Am)*dtheta*I;
		}	
		NUMERATORB += r*NUMERATORA*dr;
	}
	return(NUMERATORB);
}
float AREA(float R_a,float R_b, float d) //calculates the area of the primary light source which is unobstructed
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
	double A,Z,L_a,		//,L_b; Z is to compute intersection area the old-school way
	M_b,			//M_a,M_b;//mass of the primary and secondary
	theta,			//zenith and polar angles respectively
	psi,			//angle of rotation in it's own circular orbit.
	d,			//projected distance between the two centres
	R_a,R_b,		//Primary and secondary radii
	Flux, 			//Flux as a function of time
	A_a,A_b,		//Area of the first and second
	v,			//actual and observed rotational velocities
	i,			//inclination in radians
	r,			//radius of orbit
	D_S,			//D_LS + D_L
 	D_LS,			//lens-source distance
	D_L,			//lens-observer distance
	D, 			//ratio LS*L/S
	R_E,			//distance from primary to us
	micro,			//value of microlensing added according to SG
	DENOMINATORA =0,	//only used for numerical integration
	DENOMINATORB =0,
	R,			//counter for denominator integration
	c1, c2, c3, c4,
	mu, I,rnorm;		
	int ti;			//time to intersection
	int j,t,ttime;
	int Phase=1;		//t will be in minutes
	int mint=0;		//minimum time in system's orbit to begin the simulation
	float x,y,z;		//new coordinates;

	c1 = 0.5905;//0.701;//G0V: 0.701;//0.5118; //KOI81: 0.5905
	c2 = 0.0486;//0.149;//0.149;	//0.0525;	//0.0486
	c3 = -0.0476;//0.277;//0.277;	//0.4590;	//-0.0476
	c4 = 0.0074;//-0.297;//-0.297;	//-0.2727;	//0.0074

/**************************Default setting:*/
v = 0;
Flux = 0;
R_b = 0;
R_a = 0;
i = 0;
ttime = 10000;
psi=pi;
D_S = 600000;
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

			case 'v'://rotational velocity
			sscanf(argv[j+1], "%f", &v);
			break;

			case 'i'://inclination
			sscanf(argv[j+1], "%f", &i);
			break;

			case 'f'://distance from primary to us
			sscanf(argv[j+1], "%f", &D_S);
			break;
			
			case 'y'://min time
			sscanf(argv[j+1], "%d", &mint);
			break;
				   
			default:
			break;		
		}
	}

	R_a = R_a*r_sun;
	R_b = R_b*r_sun;
	M_b = M_b*m_sun;
	r = au*r;
	L_a = L_a*l_sun;
	D_S = D_S*au;

	FILE *outfile;
	char *mode= "w";
	outfile = fopen("output.dat",mode);
	if (outfile ==NULL)
	{
		fprintf(stderr,"Can't open output.txt\n");
		exit (1);
	}
	A_a = pi*R_a*R_a;
	A_b = pi*R_b*R_b;
	const float firstz = z;//first value given of d.
	const float firstFlux = L_a/A_a;
	micro = 1;
	psi = mint*v/r + pi;
	if((psi>=pi)&&(psi<3*pi/2))
		{Phase = 1;}
	if((psi>=pi*3/2)&&(psi<2*pi))
		{Phase = 2;}
	if((psi>=pi*2)&&(psi<5*pi/2))
		{Phase = 3;}
	if((psi>=5*pi/2)&&(psi<3*pi))
		{Phase = 4;}
	printf("Intersection at %f seconds\n",pi*r/(v*2));
	double dr = R_a/10000;
	double dtheta = 2*pi/10000;
	for(R=0;R<R_a;R += dr)
	{
		DENOMINATORA = 0;
		for (theta=0;theta<2*pi;theta +=dtheta)
		{
			rnorm = R/R_a;
			mu = sqrt(1 - rnorm*rnorm);
			I = 1.0 - c1*(1 - pow(mu,0.5)) - c2*(1 - mu) - c3*(1 - pow(mu,1.5)) - c4*rnorm*rnorm;
			DENOMINATORA += dtheta*I;
		}	
		DENOMINATORB += R*DENOMINATORA*dr;
	}
	ti = pi*r/(2*v);
	for(t=mint;t<ttime;t++)
	{	
		psi = t*v/r + pi;
		x = -r*cos(i)*sin(psi);
		D_LS = x;	
		D_L = D_S - D_LS;
		D = D_LS*D_L/D_S;
		R_E = 2*sqrt(M_b*G*D)/c;
		y = r*cos(psi); //projected coordinates
		z = r*sin(i)*sin(psi);//plane of the orbit
		d = sqrt(z*z + y*y);//R_S
		if(x>=0)
		{
			if (R_a + R_b <= d)//separated
			{
				micro = SG(R_E,M_b,R_a,R_b,d,y,z);
				Flux = micro/DENOMINATORB;
				fprintf(outfile,"%d\t%.10f\tseparated\n",(t - ti),Flux);
			}
			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
				micro = SG(R_E,M_b,R_a,R_b,d,y,z);
				Flux = micro/DENOMINATORB;
				fprintf(outfile,"%d\t%.10f\tintersecting\n",(t - ti),Flux);
			}
			if ((d+R_b)<R_a)//inside
			{
				micro = SG(R_E,M_b,R_a,R_b,d,y,z);
				Flux = micro/DENOMINATORB;
				fprintf(outfile,"%d\t%.10f\tinside\n",(t - ti),Flux);
			}
			if (((psi>=pi*3/2)&&(psi<2*pi))||(y>=0))
			{
				Phase = 2;
			}
			if ((y<-r*cos(2*v/r))||(y>r*cos(2*v/r)))
			{
				Phase = 3;
			}
		}
		if(x<0)
		{
			if (R_a + R_b <= d)//separated
			{
				Flux = 1;
				fprintf(outfile,"%d\t%f\t%d\tseparated\tr=%f\n",(t - ti),Flux,Phase);
			}
			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
				A = AREA(R_a,R_b,d);
				Z = (A_b - A)/(A_b*A_b);
				//Flux = ((L_a/A_a + L_b*Z)/firstFlux);
				fprintf(outfile,"%d\t%f\t%d\tintersecting\n",(t - ti),Flux,Phase);
			}
			if ((d+R_b)<R_a)//inside
			{
				Flux = (L_a/A_a/firstFlux);
				fprintf(outfile,"%d\t%f\t%d\tinside\n",(t - ti),Flux,Phase);
			}
			if ((y<r*(100000)*sin(v/(2*r)))&&(y>-r*(100000)*sin(v/(2*r))))
			{
				Phase = 4;
			}
			if ((y<-r*cos(2*v/r))||(y>r*cos(2*v/r)))
			{
				Phase = 1;
			}
		}
	}
	fclose(outfile);
	
	return (0);
}
