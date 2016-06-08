#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h> 

#define bands 3

main(int argc, char *argv[]){
  FILE *infp, *outfp;
  
  float input[bands]; /*dt probs, logit probs, class map stack*/
  float output;
  int i,j;
  
  if(argc == 1){
    printf("\n\nbayes.exe <input> <output>\n\n");
    exit(1);
  }
    
  if((infp=fopen(argv[1],"rb"))==NULL){
    printf("\n\ninput file open error\n\n");
    exit(1);
  }
    
  if((outfp=fopen(argv[2],"wb"))==NULL){
    printf("\n\noutput file open error\n\n");
    exit(1);
  }


  while(fread(&input, sizeof(float),bands,infp)){
        
        /* range check*/
        if(input[0] == 0){
           input[0] = 0.00001;
        }
	if(input[1] == 0){
	  input[1] = 0.00001;
	}

	/*if missing data in NBARS, use the probs from the logit model and disregard the dt probs layer. otherwise, use bayes to estimate urban probs
 */
	
	if(input[2] == 18){   /*assumes missing data class is 18 in class map*/

	  output = input[1];
	  }

	else {

	output = input[0]*input[1] / ((input[0] * input[1]) + ((1 - input[0]) * (1 - input[1]))); 
	
	}

        fwrite(&output,sizeof(float),1,outfp);
  }
 fclose(infp);
 fclose(outfp);
 
} // end main
