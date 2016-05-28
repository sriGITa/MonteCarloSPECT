#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>   
#include <unistd.h>

#ifndef pi
#define pi 3.14159265358979323846
#endif 

int main(int argc, char *argv[])
{
    float angle,total_angle;
    int  number,i,j,counts,subnumber,type,gpsorgun,source_type,Jaszczak_type;
    int  sum,m,n,l,count[6]={0};
    FILE* slice;
    FILE* fp;
    FILE* test;
    char filename[128]={0};
    char pathandname[256]={0};
    char buf[80];
    char bufup[80];
    float sourcetocollimator,Length,offset,radius,halfZ;
	float cylinder_R[6]={0.32,0.48,0.64,0.79,0.95,1.11};
    float cylinder_D[6]={8,6,5,5,4,4};
	//int   cylinder_Count[6]={14,9,7,6,5,4};
	int   cylinder_Count[6]={6,4,4,3,3,2};
	float cylinder_X,cylinder_Y,cylinder_Z,temp_XX,temp_YY,temp2_XX,temp2_YY;
	float *cylinder_XX,*cylinder_YY; 
    
    // get the current dictionary
    getcwd(buf, sizeof(buf));
    chdir("..");
    getcwd(bufup,sizeof(bufup));
    chdir(buf);
   
    
    type=1;    
    if (argc != 4) {
       	printf("SYNTAX: %s <each_angle> <total_angle> <beamOn>\n", argv[0]);
       	exit(1);
    }
    
   sscanf(argv[1],"%f",&angle);
   sscanf(argv[2],"%f",&total_angle);
   sscanf(argv[3],"%d",&counts);
   
   sprintf(filename,"nuclide.txt");
   sprintf(pathandname,"%s/%s",buf,filename);
   	if (!(fp = fopen(pathandname,"r")))
   	{
		printf("Error opening nuclide file %s!\n", filename);
		fclose(fp);
		return 1;
	}
	fscanf(fp,"%d\n",&type);
	fscanf(fp,"%d\n",&source_type);
	fscanf(fp,"%f\n",&radius);
	fscanf(fp,"%f\n",&halfZ);
	fscanf(fp,"%d\n",&Jaszczak_type);
	printf("%d\n",type);   
	
   sprintf(filename,"collimator_detector.txt");
   sprintf(pathandname,"%s/%s",buf,filename);
   if(!(test=fopen(pathandname,"r")))
   {
	  printf("Error opening collimator_detector file %s!\n", filename);
      fclose(test);
	  return 1;
   }
   fscanf(test,"gpsorgun=%d\n",&gpsorgun);
      
   
   if(type>6||type<1)
   type=1;
   
   sprintf(filename,"vis.mac");
   number=(int)(total_angle/angle);
   subnumber=number/12;
   
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
   
  if(Jaszczak_type==2)
   {
    cylinder_R[0]=0.48; cylinder_R[1]=0.64; cylinder_R[2]=0.79;
    cylinder_R[3]=0.95; cylinder_R[4]=1.11; cylinder_R[5]=1.27;
  
    cylinder_D[0]=6; cylinder_D[1]=5; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=4; cylinder_D[5]=4;  

    cylinder_Count[0]=10; cylinder_Count[1]=7; cylinder_Count[2]=6;
    cylinder_Count[3]=5; cylinder_Count[4]=4; cylinder_Count[5]=4;  

 	}
  if(Jaszczak_type==3)
   {
	cylinder_R[0]=0.64; cylinder_R[1]=0.79; cylinder_R[2]=0.95;
    cylinder_R[3]=1.11; cylinder_R[4]=1.27; cylinder_R[5]=1.91;
    
    cylinder_D[0]=5; cylinder_D[1]=5; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=3; cylinder_D[5]=2.5;  

    cylinder_Count[0]=8; cylinder_Count[1]=6; cylinder_Count[2]=5;
    cylinder_Count[3]=4; cylinder_Count[4]=4; cylinder_Count[5]=3; 
   }
  if(Jaszczak_type==4)
   {
	cylinder_R[0]=0.95; cylinder_R[1]=1.11; cylinder_R[2]=1.27;
    cylinder_R[3]=1.59; cylinder_R[4]=1.91; cylinder_R[5]=2.54;
    
    cylinder_D[0]=5; cylinder_D[1]=4; cylinder_D[2]=4;
    cylinder_D[3]=4; cylinder_D[4]=3; cylinder_D[5]=3;  

    cylinder_Count[0]=5; cylinder_Count[1]=4; cylinder_Count[2]=4;
    cylinder_Count[3]=3; cylinder_Count[4]=3; cylinder_Count[5]=2; 
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

sum=0;
for(i=0;i<6;i++)
sum=sum+count[i];
cylinder_XX=cylinder_XX-sum;
cylinder_YY=cylinder_YY-sum;
     
 if(gpsorgun==1)
  {
   for(i=1;i<=12;i++)
   {  
   sprintf(pathandname,"%s/%d/%s",bufup,i,filename);
   
   if((slice=fopen(pathandname,"w+"))==NULL)
	  return 0;
	 
	fprintf(slice,"/run/initialize\n");
	fprintf(slice,"/gun/particle ion\n");
	if(type==1)
	{
	fprintf(slice,"/gun/ion 43 99 0 140.6833 keV\n");
	fprintf(slice,"/N01/fullChain false\n");
	fprintf(slice,"/grdm/nucleusLimits 100 100 43 43\n\n");
    }
    if(type==2)
    {
	fprintf(slice,"/gun/ion 53 123 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 123 123 53 53\n\n");
    }
    if(type==3)
    {
	fprintf(slice,"/gun/ion 27 57 0 0 keV\n\n");
    }
    if(type==4)
    {
	fprintf(slice,"/gun/ion 49 111 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 111 111 49 49\n\n");
    }
    if(type==5)
    {
     fprintf(slice,"/gun/ion 81 201 0 0 keV\n");
     fprintf(slice,"/grdm/nucleusLimits 201 201 81 81\n\n");
    }
    if(type==6)
    {
     fprintf(slice,"/gun/ion 53 131 0 0 keV\n");
     fprintf(slice,"/grdm/nucleusLimits 131 131 53 53\n\n");
    }                                                                                                                                                                                         

   for(j=0;j<subnumber;j++)
	   {
	     fprintf(slice,"/mydet/armAngle\t%f\tdeg\n",(i-1)*30+j*angle);
         fprintf(slice,"/run/beamOn\t%d\n",counts);
        }
       
    fclose(slice);
   }
 }
 if(gpsorgun==2)
 {
  for(i=1;i<=12;i++){ 
     sprintf(pathandname,"%s/%d/%s",bufup,i,filename);
     if((slice=fopen(pathandname,"w+"))==NULL)
	  return 0;	 
	fprintf(slice,"/run/initialize\n");
	fprintf(slice,"/gps/particle ion\n");
	
    fscanf(test,"sourcetocollimator=%f\n",&sourcetocollimator);
    fscanf(test,"halfLengofcollimator=%f\n",&Length); 
    offset=sourcetocollimator+Length;
    
  if(source_type==1){
    if(type==1){
    fprintf(slice,"/gps/ion 43 99 0 140.6833 keV\n");
    fprintf(slice,"/N01/fullChain false\n");
	fprintf(slice,"/grdm/nucleusLimits 100 100 43 43\n");
   }
   if(type==2){
    fprintf(slice,"/gps/ion 53 123 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 123 123 53 53\n\n");
   }
   if(type==3){
	fprintf(slice,"/gps/ion 27 57 0 0 keV\n\n"); 
   }
   if(type==4){
    fprintf(slice,"/gps/ion 49 111 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 111 111 49 49\n\n");   
   }
   if(type==5){
     fprintf(slice,"/gps/ion 81 201 0 0 keV\n");
     fprintf(slice,"/grdm/nucleusLimits 201 201 81 81\n\n");
   }
    if(type==6){
     fprintf(slice,"/gps/ion 53 131 0 0 keV\n");
     fprintf(slice,"/grdm/nucleusLimits 131 131 53 53\n\n");
   }   
	fprintf(slice,"/gps/pos/type Point\n");	
	fprintf(slice,"/gps/pos/center 0.0 0.0 %f cm\n",-offset);
    fprintf(slice,"/gps/ang/type iso\n\n");
}
    
  if(source_type==2){	
	if(type==1){
    fprintf(slice,"/gps/ion 43 99 0 140.6833 keV\n");
    fprintf(slice,"/N01/fullChain false\n");
	fprintf(slice,"/grdm/nucleusLimits 100 100 43 43\n");
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");    
   }
   	if(type==2){
	fprintf(slice,"/gps/ion 53 123 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 123 123 53 53\n");
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");    
   }
   if(type==3){
	fprintf(slice,"/gps/ion 27 57 0 0 keV\n");
    fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");  	   
   }
   if(type==4){
    fprintf(slice,"/gps/ion 49 111 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 111 111 49 49\n");
    fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");  	   
   }
   if(type==5){
	fprintf(slice,"/gps/ion 81 201 0 0 keV\n");
    fprintf(slice,"/grdm/nucleusLimits 201 201 81 81\n");
    fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	   
   }
   if(type==6){
	fprintf(slice,"/gun/ion 53 131 0 0 keV\n");
    fprintf(slice,"/grdm/nucleusLimits 131 131 53 53\n");
    fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f\n",radius);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position 0.0 0.0 %f cm\n",-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	 	   
   } 
  }
  
  if(source_type==3){
	if(type==1){
    sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
	fprintf(slice,"/gps/ion 43 99 0 140.6833 keV\n");
    fprintf(slice,"/N01/fullChain false\n");
	fprintf(slice,"/grdm/nucleusLimits 100 100 43 43\n");
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
	fprintf(slice,"/gps/ion 43 99 0 140.6833 keV\n");
    fprintf(slice,"/N01/fullChain false\n");
	fprintf(slice,"/grdm/nucleusLimits 100 100 43 43\n");		
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  }
}
   if(type==2){
    sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
    fprintf(slice,"/gps/ion 53 123 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 123 123 53 53\n");
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
    fprintf(slice,"/gps/ion 53 123 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 123 123 53 53\n");	
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  }
 }
   if(type==3){
	sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
	fprintf(slice,"/gps/ion 27 57 0 0 keV\n\n"); 
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
	fprintf(slice,"/gps/ion 27 57 0 0 keV\n\n"); 
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  }
}
   if(type==4){
	sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
	fprintf(slice,"/gps/ion 49 111 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 111 111 49 49\n\n"); 
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
	fprintf(slice,"/gps/ion 49 111 0 0 keV\n");
	fprintf(slice,"/grdm/nucleusLimits 111 111 49 49\n\n"); 
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  }  
   }
   if(type==5){
     sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
    fprintf(slice,"/gps/ion 81 201 0 0 keV\n");
    fprintf(slice,"/grdm/nucleusLimits 201 201 81 81\n\n");
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
    fprintf(slice,"/gps/ion 81 201 0 0 keV\n");
    fprintf(slice,"/grdm/nucleusLimits 201 201 81 81\n\n");
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  } 
}
    if(type==6){
     sum=0; 
  for(n=0;n<6;n++)
    { 
    for(m=sum;m<count[n]+sum;m++)
     { 
	 if(n==0&&m==0){
     fprintf(slice,"/gps/ion 53 131 0 0 keV\n");
     fprintf(slice,"/grdm/nucleusLimits 131 131 53 53\n\n");
    }
    else
    {
	fprintf(slice,"/gps/source/add 1.0\n");
    fprintf(slice,"/gps/ion 53 131 0 0 keV\n");
    fprintf(slice,"/grdm/nucleusLimits 131 131 53 53\n\n");
	}
	
	fprintf(slice,"/gps/pos/type Volume\n");
	fprintf(slice,"/gps/pos/shape Cylinder\n");
	fprintf(slice,"/gps/pos/radius %f cm\n",cylinder_R[n]);
	fprintf(slice,"/gps/pos/halfz %f cm\n",halfZ);
	fprintf(slice,"/gps/position %f %f %f cm\n",*cylinder_XX,*cylinder_YY,-offset);
	fprintf(slice,"/gps/ang/type iso\n");
	fprintf(slice,"/gps/pos/rot1 1 0 1\n");
	fprintf(slice,"/gps/pos/rot2 0 0 1\n\n");	
	cylinder_XX++;
	cylinder_YY++;
	}
	 sum=sum+count[i];
  }
   }   
   
}    
    for(j=0;j<subnumber;j++)
	 {
	  fprintf(slice,"/mydet/armAngle\t%f\tdeg\n",(i-1)*30+j*angle);
      fprintf(slice,"/run/beamOn\t%d\n",counts);
     }
       
    fclose(slice);   
 }	 
}
   
   fclose(fp);
   fclose(test);
   printf("Generate vis.mac files successfully!\n");
   return 0;
    
    
}

