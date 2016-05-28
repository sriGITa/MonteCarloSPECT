#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include <unistd.h>

#ifndef pi
#define pi 3.14159265358979323846
#endif  
 
int main(int argc, char *argv[])
{
	int   nx,ny,nz,i,j,k,m,n,l,sum;
	float sphere_X,sphere_Y,sphere_Z,sphere_XX,sphere_YY,max,X,Y;
	float cylinder_X,cylinder_Y,cylinder_Z,temp_XX,temp_YY,temp2_XX,temp2_YY;
	float *cylinder_XX,*cylinder_YY; 
	char  *s="standard";
	char  *s1="ultra_deluxe";
	char  *s2="deluxe";
	char  *s3="standard";
	char  *s4="benchmark";
	int   count[6]={0};
	float sphere_R[6]={0.95,1.27,1.59,1.91,2.54,3.18};
	float cylinder_R[6]={0.32,0.48,0.64,0.79,0.95,1.11};
	float cylinder_D[6]={8,6,5,5,4,4};
	int   cylinder_Count[6]={14,9,7,6,5,4};
	int   resolution,type;
    int   *tdata,*pdata;
    int   *ttdata,*ppdata;  
    FILE* slice[2000];
    FILE* fp;
    char  filename[128]={0};
    char  pathandname[256]={0};
    char  buf[80];
    char  bufup[80];
    getcwd(buf, sizeof(buf));
    cd("..");
    -
        
    if (argc != 3) {
       	printf("SYNTAX: %s <resolution> <choose phantom type>\n", argv[0]);
       	printf("resolution=100, 0.1mm per pixel; resolution=20, 0.5mm per pixel; \nresolution=10,1mm per pixel; resolution=1, 1cm per pixel.\n");
       //	printf("type=1 for Ultra Deluxe,type=2 for Deluxe,\ntype=3 for Standard and type=4 for Benchmark.\n");
      	printf("phantom name: ultra_deluxe, deluxe, standard, benchmark");
      	exit(1);
	    }

   sscanf(argv[1],"%d",&resolution);
   s=argv[2];
   for(i=0;i<strlen(s);i++)
   s[i]=tolower(s[i]);
  
   if (strcmp(s,s1)==0)
   	{type = 1;}
    if (strcmp(s,s2)==0)
   	{type = 2;}
   if (strcmp(s,s3)==0)
   	{type = 3;}
    if (strcmp(s,s4)==0)
   	{type = 4;}
      
   if(type<5)
    printf("type=%d\t\t%s\tphantom\n",type,s);   
    else if(type>=5||type<=0)
    {
	printf("error to type in the phantom name.\n");
	exit(1);
	}
    
  if(type==2)
   {
    cylinder_R[0]=0.48; cylinder_R[1]=0.64; cylinder_R[2]=0.79;
    cylinder_R[3]=0.95; cylinder_R[4]=1.11; cylinder_R[5]=1.27;
  
    cylinder_D[0]=6; cylinder_D[1]=5; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=4; cylinder_D[5]=4;  

    cylinder_Count[0]=10; cylinder_Count[1]=7; cylinder_Count[2]=6;
    cylinder_Count[3]=5; cylinder_Count[4]=4; cylinder_Count[5]=4;  

 	}
  if(type==3)
   {
	cylinder_R[0]=0.64; cylinder_R[1]=0.79; cylinder_R[2]=0.95;
    cylinder_R[3]=1.11; cylinder_R[4]=1.27; cylinder_R[5]=1.91;
    
    sphere_R[0]=1.27;sphere_R[1]=1.59;sphere_R[2]=1.91;
    sphere_R[3]=2.54;sphere_R[4]=3.18;sphere_R[5]=3.8;
    
    cylinder_D[0]=5; cylinder_D[1]=5; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=3; cylinder_D[5]=2.5;  

    cylinder_Count[0]=8; cylinder_Count[1]=6; cylinder_Count[2]=5;
    cylinder_Count[3]=4; cylinder_Count[4]=4; cylinder_Count[5]=3; 
   }
  if(type==4)
   {
	cylinder_R[0]=0.95; cylinder_R[1]=1.11; cylinder_R[2]=1.27;
    cylinder_R[3]=1.59; cylinder_R[4]=1.91; cylinder_R[5]=2.54;
    sphere_R[0]=1.27;sphere_R[1]=1.59;sphere_R[2]=1.91;
    sphere_R[3]=2.54;sphere_R[4]=3.18;sphere_R[5]=3.8;
    
    cylinder_D[0]=5; cylinder_D[1]=4; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=3; cylinder_D[5]=3;  

    cylinder_Count[0]=5; cylinder_Count[1]=4; cylinder_Count[2]=4;
    cylinder_Count[3]=3; cylinder_Count[4]=3; cylinder_Count[5]=2; 
   }

  // pixel size of the phantom
   nx=(21.6*resolution)+1;
   ny=(21.6*resolution)+1;
   nz=(14.6*resolution)+1;
  
 // center position of the spheres 
   sphere_X=21.6*resolution/2;
   sphere_Y=21.6*resolution/2;
   sphere_Z=(14.6-12.7)*resolution;
 
 // center position of the cylinders 
   cylinder_X=21.6*resolution/2;
   cylinder_Y=21.6*resolution/2;
   cylinder_Z=(14.6-8.8)*resolution;
  
 // check the input parameters
   if(resolution<0)
   {
	  printf("Resolution must be greater than 0.\n");
      exit(1);  
	}
 
 // create the output files
   for(i=0;i<nz;i++)
    {
	  if(i<10)
	  sprintf(filename,"slice0%d.txt",i);
	  else
	  sprintf(filename,"slice%d.txt",i);
	  
	  sprintf(pathandname,"%s/data/%s",buf,filename);        
	  if((slice[i]=fopen(pathandname,"w+"))==NULL)
	  return 0;
	}
	
	sprintf(filename,"phantom");
	sprintf(pathandname,"%s/%s",buf,filename);
	if((fp=fopen(pathandname,"w+"))==NULL)
	return 0;
 
 // malloc the memory 
   	if ((tdata = (int*)malloc(nx*ny*nz*sizeof(int)))==NULL)
	{
		fprintf(stderr,"Error in allocating memory\n");
		return 1;
	}
 //	malloc the memory
   if ((pdata = (int*)malloc(nx*ny*sizeof(int)))==NULL)
	{
		fprintf(stderr,"Error in allocating memory\n");
		return 1;
	}
 // malloc the memory for the center position for each cylindier
   if ((cylinder_XX = (float*)malloc(400*sizeof(int)))==NULL)
	{
		fprintf(stderr,"Error in allocating memory\n");
		return 1;
	}
   if ((cylinder_YY = (float*)malloc(400*sizeof(int)))==NULL)
	{
		fprintf(stderr,"Error in allocating memory\n");
		return 1;
	}
	
   ttdata=tdata;
   ppdata=pdata;
 // initialize  
   for(i=0;i<nx*ny*nz;i++)
      *tdata++=0;
   tdata=tdata-nx*ny*nz;
   
   for(i=0;i<nx*ny;i++)
      *pdata++=0;
   pdata=pdata-nx*ny;  
   
   max=0;
   for(m=0;m<6;m++)
   {
	  sphere_R[m]= sphere_R[m]*resolution/2; 
	  cylinder_R[m]=cylinder_R[m]*resolution/2;
	  if(sphere_R[m]>max)
	  max=sphere_R[m];
   }
 // generate the spheres
  for (m=0;m<6;m++)
  {   
	sphere_XX=sphere_X+(max*2+sphere_R[m])*cos(m*60*pi/180+30*pi/180);
    sphere_YY=sphere_Y+(max*2+sphere_R[m])*sin(m*60*pi/180+30*pi/180);
     for(i=0;i<nz;i++)
        {   
          for(j=0;j<ny;j++)
          {
			for(k=0;k<nx;k++)
			{
			  if(((k-sphere_XX)*(k-sphere_XX)+(j-sphere_YY)*(j-sphere_YY)+(i-sphere_Z)*(i-sphere_Z))<=(sphere_R[m]*sphere_R[m]))
			    { 
				 tdata=ttdata;
				 tdata=tdata+ny*nx*i+nx*j+k;
				 *tdata=1;	 		  	  
			    }			   	    
			}
		 }
	   }
 }
 
 // generate the cylinders
 for(m=0;m<6;m++)
 {
	 temp_XX=cylinder_X+cylinder_R[m]*cylinder_D[m]*cos(m*60*pi/180+30*pi/180);
	 temp_YY=cylinder_Y+cylinder_R[m]*cylinder_D[m]*sin(m*60*pi/180+30*pi/180); 
	 
     
	  for(n=0;n<cylinder_Count[m];n++)
 	  {
	   temp2_XX=temp_XX+4*cylinder_R[m]*n*cos(m*60*pi/180+60*pi/180);	 
	   temp2_YY=temp_YY+4*cylinder_R[m]*n*sin(m*60*pi/180+60*pi/180);
	  
	   *cylinder_XX=temp2_XX;
	   *cylinder_YY=temp2_YY;
	   //printf("%f\t%f\n",temp2_XX,temp2_YY);
	    count[m]++;
		for(l=0;l<cylinder_Count[m]-n-1;l++)
	      {  cylinder_XX++;
	         cylinder_YY++;
	         *cylinder_XX=temp2_XX+4*cylinder_R[m]*(l+1)*cos(m*60*pi/180);
	         *cylinder_YY=temp2_YY+4*cylinder_R[m]*(l+1)*sin(m*60*pi/180);
	        // printf("%f\t%f\n",*cylinder_XX,*cylinder_YY);
	         count[m]++;
          }
	   cylinder_XX++;
	   cylinder_YY++;  
	    
     }
}

//for(i=0;i<6;i++)
//printf("%d\t",count[i]);

sum=0;
for(i=0;i<6;i++)
sum=sum+count[i];
cylinder_XX=cylinder_XX-sum;
cylinder_YY=cylinder_YY-sum;


sum=0; 
 for(i=0;i<6;i++)
   { 
   for(m=sum;m<count[i]+sum;m++)
     { 
      for(j=0;j<ny;j++)
       {
	    for(k=0;k<nx;k++)
	 	  {
 			  if(((k-*cylinder_XX)*(k-*cylinder_XX)+(j-*cylinder_YY)*(j-*cylinder_YY))<=(cylinder_R[i]*cylinder_R[i]))
			  { 
				 pdata=ppdata;
				 pdata=pdata+nx*j+k;
				 *pdata=1;	 
				 //printf("%d\t%d\n",k,j);  	  
			   }
		   }
	   }
	   cylinder_XX++;
	   cylinder_YY++;
	}

	 sum=sum+count[i];
 }


for(i=(int)(cylinder_Z);i<nz;i++)
{
   tdata=ttdata;
   tdata=tdata+ny*nx*i;
     for(j=0;j<ny;j++)
     {
	  for(k=0;k<nx;k++)
  	   {	 
		  *tdata++=*ppdata++;
	   }
    }
  ppdata=ppdata-nx*ny;    
}

tdata=ttdata;
 // output the results
   for( i=0;i<nz;i++)
     { for( j=0;j<ny;j++)
        {
	     for(k=0;k<nx;k++)
	      { if(k<nx)
          {fprintf(slice[i],"%d\t",*ttdata++);}
          else
          fprintf(slice[i],"%d",*ttdata++);
	      
	      fprintf(fp,"%d\n",*tdata++);
	      }
         fprintf(slice[i],"\n");
	   }
     }    

    for(i=0;i<nz;i++)
    fclose(slice[i]);
    fclose(fp);
    printf("Generate a Jaszczak phantom!\n");
    printf("size:\t%dx%dx%d\n",nx,ny,nz);
    return 0;
    
}
