#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

main(int argc, char *argv[]){
 FILE *inred, *innir, *incloud, *outfp;
 unsigned char red;
 unsigned char nir;
 unsigned char cloud;
 float ndvi; 
 int i,j,k,l;

 
 if(argc == 1){
   printf("\ncalc_ndvi.exe <red> <nir> <cloud> <output>\n\n");
   exit(1);
 }
 
 if((inred=fopen(argv[1],"rb"))==NULL){
   printf("\ninput RED file open error\n\n");
   exit(0);
 }

 if((innir=fopen(argv[2],"rb"))==NULL){
   printf("\ninput NIR file open error\n\n");
   exit(0);
 }

 if((incloud=fopen(argv[3],"rb"))==NULL){
   printf("\ninput CLOUD file open error\n\n");
   exit(0);
 }

 if((outfp=fopen(argv[4],"wb"))==NULL){
   printf("\noutput file open error\n\n");
   exit(0);
 }

 while(fread(&red,sizeof(unsigned char),1,inred)){
        fread(&nir,sizeof(unsigned char),1,innir);
	fread(&cloud,sizeof(unsigned char),1,incloud);

	if((red > 0) && (cloud < 2)){
		ndvi = ((float)nir - (float)red)/((float)nir + (float)red);
        }else{
		ndvi = 0.0;
  	}     
   	fwrite(&ndvi,sizeof(float),1,outfp);
 }
  
 fclose(inred);
 fclose(innir);
 fclose(incloud);
 fclose(outfp);
 
} /* end main */
