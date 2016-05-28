#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>


int main(int argc, char *argv[])
{
	int nx, ny, nz, type,i,j,k,t; 
	long  ii,jj,n;
    FILE* fp;
    FILE* file;
    FILE* totaldata;
    FILE* sino1;
    FILE* sino2;
    FILE* detector1;
    FILE* detector2;
    FILE* pixel_detector1;
    FILE* pixel_detector2;
    FILE* sino1_one;
    FILE* sino2_one;
    int *fdata,*pdata,*tdata,*sdata,*ddata,*odata,*sodata,*ssdata;
    int out;
    FILE* slice[400]; FILE* sino[400];
    char filename[128]={0};
    char pathandname[256]={0};
    char aa;
    char s[200];
    char buf[100];
    getcwd(buf, sizeof(buf));
    
    if (argc != 2) {
       	printf("SYNTAX: %s <input[Geant]> \n", argv[0]);
       	exit(1);
    }
  
    sprintf(pathandname,"%s/parameter.txt",buf);
    fp=fopen(pathandname,"r");//open the file in read mode
    if(fp==NULL)
    {
		printf("error to open parameter.txt\n.");
		return;
	}
	

	sprintf(filename,"%s",argv[1]);
	sprintf(pathandname,"%s/%s",buf,filename);
	file=fopen(pathandname,"r");//open the file in read mode
    if(file==NULL)
    {
		printf("error to open %s\n.",argv[1]);
		return;
	}
	
	sprintf(pathandname,"%s/results/sinogram1.txt",buf);
	sino1=fopen(pathandname,"w+");//open the file in read mode
    if(sino1==NULL)
    {
		printf("error to create sinogram1.txt\n.");
		return;
	}
	
	sprintf(pathandname,"%s/results/sinogram2.txt",buf);
	sino2=fopen(pathandname,"w+");//open the file in read mode
    if(sino2==NULL)
    {
		printf("error to create sinogram2.txt\n.");
		return;
	}
	
	sprintf(pathandname,"%s/results/sinogram1_one.txt",buf);
	sino1_one=fopen(pathandname,"w+");//open the file in read mode
    if(sino1_one==NULL)
    {
		printf("error to create sinogram1_one.txt\n.");
		return;
	}
	
	sprintf(pathandname,"%s/results/sinogram2_one.txt",buf);
	sino2_one=fopen(pathandname,"w+");//open the file in read mode
    if(sino2_one==NULL)
    {
		printf("error to create sinogram2_one.txt\n.");
		return;
	}
	
	sprintf(pathandname,"%s/results/pixel_detector1.txt",buf);
    pixel_detector1=fopen(pathandname,"w+");//open the file in read mode
    if(pixel_detector1==NULL)
    {
		printf("error to create pixel_detector1.txt\n.");
		return;
	}
	
	sprintf(pathandname,"%s/results/pixel_detector2.txt",buf);
	pixel_detector2=fopen(pathandname,"w+");//open the file in read mode
    if(pixel_detector2==NULL)
    {
		printf("error to create pixel_detector2.txt\n.");
		return;
	}
	
    sprintf(filename,"detector1.txt");
	sprintf(pathandname,"%s/results/%s",buf,filename);
	detector1=fopen(pathandname,"w+");
	if(detector1==NULL)
	{
	   printf("error to create %s\n",filename);
	   return;
	}
	
    sprintf(filename,"detector2.txt");
	sprintf(pathandname,"%s/results/%s",buf,filename);
	detector2=fopen(pathandname,"w+");
	if(detector2==NULL)
	{
	   printf("error to create %s\n",filename);
	   return;
	}


	if (ReadArgument(fp, "x=%s", s))
   	  {
		fclose(fp);
		return 2;
		}
	nx=atoi(s);
	
    if (ReadArgument(fp, "y=%s", s))
   	  {
		fclose(fp);
		return 2;
		}
	ny=atoi(s);

    if (ReadArgument(fp, "z=%s", s))
   	  {
		fclose(fp);
		return 2;
		}
	nz=atoi(s);	
  
   for(i=0;i<nz*2;i++)
    {
	  if(i<10)
	  sprintf(filename,"slice0%d.txt",i);
	  else
	  sprintf(filename,"slice%d.txt",i);
	  
	  sprintf(pathandname,"%s/data/%s",buf,filename);        
	  if((slice[i]=fopen(pathandname,"w+"))==NULL)
	  return 0;
	}
	
   for(i=0;i<nx*2;i++)
    {
	  if(i<10)
	  sprintf(filename,"sino0%d.txt",i);
	  else
	  sprintf(filename,"sino%d.txt",i);
	  
	  sprintf(pathandname,"%s/sino/%s",buf,filename);        
	  if((sino[i]=fopen(pathandname,"w+"))==NULL)
	  return 0;
	}
	

	sprintf(filename,"totaldata.txt");
	sprintf(pathandname,"%s/results/%s",buf,filename);
	if((totaldata=fopen(pathandname,"w+"))==NULL)
	    printf("error in creating totaldata.txt\n.");

  
    // Allocate memory for fdata 
	if ((fdata = (int*)malloc(nx*ny*nz*2*sizeof(int)))==NULL)
	{
		fprintf(stderr,"Error in allocating memory\n");
		return 1;
	}
    pdata=fdata;
    tdata=fdata;
    sdata=fdata;
    ddata=fdata;
    odata=fdata;
    sodata=fdata;
    ssdata=fdata;
    
    for(ii=0;ii<(2*(nx+1)*nz);ii++)
    {
		if(ii%258==0)
		{fscanf(file,"%s\n",&aa);
		//printf("%ld\t\t",ii);
		}
		if(ii!=0&&ii%129==0&&ii%258!=0)
		{fscanf(file,"%s\n",&aa);
		 printf("%ld\t\t",ii);
		}
		if(ii>0&&ii%129>0)
		{
		 for(jj=0;jj<ny;jj++)
		  {fscanf(file,"%d",&out);
		    *fdata++=out;
		  }
		  
		}
		
	}
   
	
	printf("\n%d\t%d\t%d\t",nx,ny,nz);
	n=nx*ny*nz;
	//fprintf(fp,"\n");
	for(i=0;i<nz*2;i++)
	{
     for(j=0;j<nx;j++)
	 {
      for(k=0;k<ny;k++)
	  {
		if(i%2==0)
	    {    
			t=i/2;
		    fprintf(slice[t],"%d\t",*pdata++);}
	    else
	    {
			 t=(i-1)/2+nz;
			 fprintf(slice[t],"%d\t",*pdata++);
		}
	    
	         fprintf(totaldata,"%d\t",*tdata++);
	    
	    if(i%2==0)
	    {
		     fprintf(detector1,"%d\t",*ddata++);
		     fprintf(pixel_detector1,"%d\n",*odata++);
		}
		else
		{
			fprintf(detector2,"%d\t",*ddata++);
			fprintf(pixel_detector2,"%d\n",*odata++);
		}
		
	    if(j==64)
	    {
			if(i%2==0)
	     {
			 fprintf(sino1,"%d\t",*sdata++);
			 fprintf(sino1_one,"%d\n",*sodata++);
		 }
	     else
	     {
	         fprintf(sino2,"%d\t",*sdata++);
	         fprintf(sino2_one,"%d\n",*sodata++);
	     }
	    }
	    else
	    {
			sdata++;
	        sodata++;
		} 
	 }
	    
	    fprintf(slice[t],"\n");  
	   
	    fprintf(totaldata,"\n");
	    if(i%2==0)
	    {
			fprintf(detector1,"\n");
			}
		else
		{
			fprintf(detector2,"\n");
			}
	    if(j==64)
	    {  if(i%2==0)
	          fprintf(sino1,"\n");
	       else
	          fprintf(sino2,"\n");
	    }
	 }
   }
   
 
   for(i=0;i<nz*2;i++)
	{
     for(j=0;j<nx;j++)
	 {
      for(k=0;k<ny;k++)
	  {

		  if(i%2==0)
	       {  fprintf(sino[j],"%d\t",*ssdata++); }
	      else
	       { 
	          fprintf(sino[j+nx],"%d\t",*ssdata++); }
	  }
	  if(i%2==0)
	       {  fprintf(sino[j],"\n"); }
	      else
	       { 
	          fprintf(sino[j+nx],"\n"); }
      
      }
    }

 
 
     
    for(i=0;i<nz*2;i++)
   {fclose(slice[i]);}
   for(i=0;i<nx*2;i++)
    {fclose(sino[i]);}
    printf("Convert Geant4 file to txt files!\n");
    fclose(totaldata);
    fclose(fp);
    fclose(file);
    fclose(sino1);
    fclose(sino2);
    fclose(detector1);
    fclose(detector2);
    fclose(pixel_detector1);
    fclose(pixel_detector2);
    fclose(sino1_one);
    fclose(sino2_one);
}


int ReadArgument(FILE *fp, char *Parameter, void *Data)
   {
	char s[200];
	do
   	{
		fgets(s,40,fp);
		} while(s[0] == '*');
	if (sscanf(s,Parameter,Data) != 1)
   	{
		printf("Apparent error in ReadArgument:\n");
		printf("\tInput Field: %s\n",s);
		printf("\tFormat Specifier: %s\n",Parameter);
		return 1;
		}
	return 0;
}


