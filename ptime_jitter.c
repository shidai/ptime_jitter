//"Usage: ptime_jitter -f fname -std tname (-pt tname) -o oname \n"
//"Calculate the jitter parameter\n"
//fname: residuals.dat files produced by tempo2; tname: templates; oname: output .tim; -std: standard template format; -pt: ptime template;\n"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptime_jitter.h"

int main (int argc, char *argv[])
{
	//int h,i,j,k;

	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	char tname[128];   // name of template
	char oname[128];   // name of output .tim
	int mode; // to distinguish different type of templates

	int i;
	int index, n;
	for (i=0;i<argc;i++)
    {
		if (strcmp(argv[i],"-f") == 0)
		{
            index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 )
			{
				n++;
		    }
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 0; // standard template format
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
    }
	//printf ("%d\n", smode);

	// start to deal with different data file
	//
	// open file to write toa 
	FILE *fp;
	if ((fp = fopen(oname, "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
    //fprintf (fp, "S0    S    err\n");
	/////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////
	// read std
	char data[] = "[SUBINT]";
	char psrparam[] = "[PSRPARAM]";

	char std[50];
	char std_para[50];

	strcpy(std,tname);
	strcpy(std_para,tname);

	strcat(std_para, psrparam);

	if ( mode == 0)
	{
		strcat(std, data);
	}

	double psrfreq;
	psrfreq = read_psrfreq(std_para);
	printf ("PSR frequency: %.15lf\n", psrfreq);

	int nphase, nsub, nchn;
	nchn = get_nchan(std);	
	nsub = get_subint(std);	
	nphase = get_nphase(std);	

	double s[nphase];
	read_std(std,nsub,s,nphase,mode,nchn);
	/////////////////////////////////////////////////////////

	// read residual.dat and calculate jitter parameter
	double degree;
	int k;
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		strcpy(fname,argv[k]);

		// get the number of subint
		int nint;
		getNint (fname, &nint);
		degree = nint-1.0;

		// read residual.dat
		double resi[nint], error[nint];
		readResi (fname, resi, error);

		// calculate the jitter parameter and its error
		jitter (resi, error, s, psrfreq, nint, nphase, degree, fp);
	}

    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");

	return 0;
}
