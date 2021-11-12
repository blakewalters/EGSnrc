/*
# DOSXYZnrc version: 2017.v2.rev1
#
# v1 includes the following additions to the basic 2017 version:
#     --ability to read binary .egsphant files
#     --ability to write binary .3dbindose files
#     --no. of parallel chunks as a user input (u_chunk)
#     --outputs time stamp to .lock file
# v2 includes above plus:
#     --implementation of an 2D intensity source based on the
#       analytical source model (AM) of Granton & Verhaegen,
#       PMB,58:3377-3395, 2013.
#
###############################################################################
#
#  EGSnrc functions to read and write partial dose files
#  Copyright (C) 2015 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:          Blake Walters, 1998
#
#  Contributors:    Iwan Kawrakow
#
###############################################################################
#
#  Blake Walters, 2004:
#
#  This file replaces files write_pardose.c and read_pardose.c. Since 2004 the
#  .pardose files are automatically recombined at the end of a parallel run
#  or when the user recombines them manually by running dosxyznrc with the
#  input flag IRESTART=4.
#
###############################################################################
*/


#include <stdio.h>
#include <stdlib.h>
#include "egs_config1.h"

#ifdef __cplusplus
extern "C" {
#endif

void F77_OBJ_(write_pardose,WRITE_PARDOSE)(
           double *temp2, int *imax, int *jmax, int *kmax,
           double *endep, double *endep2, const char *fname) {

    FILE *fp;
    int i,j,k;
    int ensize;
    float *sendep,*sendep2;

    if( !fname ) {
        printf("\n Null fname passed to write_pardose!\n");
        return;
    }
    ensize=(*imax)*(*jmax)*(*kmax);
    sendep = (float *) malloc(ensize*sizeof(float));
    sendep2 = (float *) malloc(ensize*sizeof(float));

    /*** open file for writing ***/

    if ((fp=fopen( fname,"wb"))==NULL) {
        printf("\n write_pardose: cannot open file %s\n",fname);
        return; /*exit(1);*/
    }

    if(fwrite((void *)temp2,
                sizeof(double),1,fp)!=1)
        printf("\n Error writing temp2. \n");

    if(fwrite((int *)imax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing imax. \n");

    if(fwrite((int *)jmax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing jmax. \n");

    if(fwrite((int *)kmax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing kmax. \n");

    for (i = 0; i < (*imax)*(*jmax)*(*kmax); i++) {
        sendep[i]=endep[i];
        
    }

    if(fwrite((float *)sendep,
                sizeof(float),ensize,fp)!=ensize)
        printf("\n Error writing endep. \n");

    if(fwrite((float *)sendep2,
                sizeof(float),ensize,fp)!=ensize)
        printf("\n Error writing endep2. \n");

    fclose(fp);
    free(sendep); free(sendep2);

    return;

}

void F77_OBJ_(read_pardose,READ_PARDOSE)(
           double *temp2, int *imax, int *jmax, int *kmax,
           float *endep, float *endep2, const char *fname) {

FILE *fp;
int i,j,k;
int ensize,xbsize,ybsize,zbsize;
float *sendep,*sendep2;

  /*** open file for reading ***/

  if ((fp=fopen(fname,"rb"))==NULL) {
    printf("\n Cannot open file \n");
    exit(1);
  }

  if(fread(temp2, sizeof(double),1,fp)!=1)
       printf("\n Error reading temp2. \n");

  if(fread((int *)imax,
      sizeof(int),1,fp)!=1)
       printf("\n Error reading imax. \n");

  if(fread((int *)jmax,
      sizeof(int),1,fp)!=1)
       printf("\n Error reading jmax. \n");

  if(fread((int *)kmax,
      sizeof(int),1,fp)!=1)
       printf("\n Error reading kmax. \n");

   /*** allocate memory for sendep,sendep2 ***/

  ensize=(*imax)*(*jmax)*(*kmax);
  sendep = (float *) malloc(ensize*sizeof(float));
  sendep2 = (float *) malloc(ensize*sizeof(float));

  if(fread((float *)sendep,
      sizeof(float),ensize,fp)!=ensize)
       printf("\n Error reading endep. \n");

  if(fread((float *)sendep2,
      sizeof(float),ensize,fp)!=ensize)
       printf("\n Error reading endep2. \n");

  for (i = 0; i < (*imax)*(*jmax)*(*kmax); i++) {
        endep[i]=sendep[i];
        endep2[i]=sendep2[i];
  }

  fclose(fp);
  free(sendep);free(sendep2);
  return;
}

void F77_OBJ_(read_binphant,READ_BINPHANT)(int *nmed,
              char *media, int *imax, int *jmax, int *kmax, float *xbound,
              float *ybound, float *zbound, short *med, double *rhor,
              const char* fname) {

FILE *fp;
int ensize;
float *srhor;
char *buf;

  /*** open file for reading ***/

  if ((fp=fopen(fname,"rb"))==NULL) {
    printf("\n Cannot open file %s\n",fname);
    exit(1);
  }

  if(fread((int *) nmed, sizeof(int),1,fp)!=1)
       printf("\n Error reading nmed. \n");

  buf = (char *) malloc(sizeof(char)*24);

  for (int i=0; i<*nmed; i++) {
    fread(buf,sizeof(char),24,fp);
    for (int j=0; j< 24; j++) {
      media[j*4+i*24*4]=buf[j];
      //fill the rest of the 4 chars with blank 
      for (int k=1; k<4; k++) media[j*4+i*24*4+k]=' ';
    }
    //read dummy char for some reason
    fread(buf,sizeof(char),1,fp);
  }

  //read no. of voxels in X,Y,Z
  if(fread( (int *) imax, sizeof(int),1,fp)!=1)
      printf("\n Error reading imax. \n");
  if(fread( (int *) jmax, sizeof(int),1,fp)!=1)
      printf("\n Error reading jmax. \n");
  if(fread( (int *) kmax, sizeof(int),1,fp)!=1)
      printf("\n Error reading kmax. \n");

  ensize = (*imax)*(*jmax)*(*kmax);

  //read voxel boundaries
  if(fread((float *) xbound, sizeof(float),*imax+1,fp)!=*imax+1)
      printf("\n Error reading xbound. \n");
  if(fread((float *) ybound, sizeof(float),*jmax+1,fp)!=*jmax+1)
     printf("\n Error reading ybound. \n");
  if(fread((float *) zbound, sizeof(float),*kmax+1,fp)!=*kmax+1)
      printf("\n Error reading zbound. \n");

  //read media numbers
  if(fread((short *) med, sizeof(short),ensize,fp)!=ensize)
     printf("\n Error reading med. \n");
  //read densities
  srhor = (float *) malloc(sizeof(float)*ensize);

  if(fread((float *) srhor, sizeof(float),ensize,fp)!=ensize)
     printf("\n Error reading srhor. \n");
  //rhor is offset by 1
  for (int i=0; i<ensize; i++) rhor[i+1]=srhor[i];

  free(srhor);free(buf);

  fclose(fp);

  return;
}

void F77_OBJ_(write_bindose,WRITE_BINDOSE)(int *imax, int *jmax, int *kmax, float *xbound,
              float *ybound, float *zbound, double *dd, double *ddun,
              const char* fname) {

FILE *fp;
int ensize;
float *sdd,*sddun;

  ensize=(*imax)*(*jmax)*(*kmax);
  sdd = (float *) malloc(ensize*sizeof(float));
  sddun = (float *) malloc(ensize*sizeof(float)); 

  /*** open file for writing ***/

  if ((fp=fopen(fname,"wb"))==NULL) {
    printf("\n Cannot open file %s\n",fname);
    exit(1);
  }

  //write dimensions

  if(fwrite((int *)imax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing imax. \n");

  if(fwrite((int *)jmax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing jmax. \n");

  if(fwrite((int *)kmax,
                sizeof(int),1,fp)!=1)
        printf("\n Error writing kmax. \n");

  //write voxel boundaries

  if(fwrite((float *)xbound,
                sizeof(float),*imax+1,fp)!=*imax+1)
        printf("\n Error writing xbound. \n"); 

  if(fwrite((float *)ybound,
                sizeof(float),*jmax+1,fp)!=*jmax+1)
        printf("\n Error writing ybound. \n");

  if(fwrite((float *)zbound,
                sizeof(float),*kmax+1,fp)!=*kmax+1)
        printf("\n Error writing zbound. \n");

  for (int i = 0; i < ensize; i++) {
        sdd[i]=dd[i];
        sddun[i]=ddun[i];

  }

  //write doses and uncertainties

  if(fwrite((float *)sdd,
                sizeof(float),ensize,fp)!=ensize)
        printf("\n Error writing doses. \n");

  if(fwrite((float *)sddun,
                sizeof(float),ensize,fp)!=ensize)
        printf("\n Error writing dose uncertainties. \n");

  free(sdd);free(sddun);
  fclose(fp);
  return;
}

void F77_OBJ_(read_i2d_file_bin,READ_I2D_FILE_BIN)(const char* fname,
              double *minx, double* maxx, int *nx, 
              double *miny, double* maxy, int *ny, 
              double *ivals) {

FILE *fp;
float *sivals,*sminx,*smaxx,*sminy,*smaxy;

  /*** open file for reading ***/

  printf("\n Reading binary file of 2D intensities, %s\n",fname);

  if ((fp=fopen(fname,"rb"))==NULL) {
    printf("\n Cannot open file %s\n",fname);
    exit(1);
  }

  sminx = (float *) malloc(sizeof(float));
  smaxx= (float *) malloc(sizeof(float));

  if(fread((float *) sminx, sizeof(float),1,fp)!=1)
       printf("\n Error reading minx. \n");
  if(fread((float *) smaxx, sizeof(float),1,fp)!=1)
       printf("\n Error reading maxx. \n");
  if(fread((int *) nx, sizeof(int),1,fp)!=1)
       printf("\n Error reading nx. \n");

  *minx=*sminx;
  *maxx=*smaxx;

  sminy = (float *) malloc(sizeof(float));
  smaxy= (float *) malloc(sizeof(float));

  if(fread((float *) sminy, sizeof(float),1,fp)!=1)
       printf("\n Error reading miny. \n");
  if(fread((float *) smaxy, sizeof(float),1,fp)!=1)
       printf("\n Error reading maxy. \n");
  if(fread((int *) ny, sizeof(int),1,fp)!=1)
       printf("\n Error reading ny. \n");

  *miny=*sminy;
  *maxy=*smaxy;

  sivals = (float *) malloc(sizeof(float)); 

  for (int i=0; i<(*nx)*(*ny); i++) {
      if(fread((float *) sivals, sizeof(float),1,fp)!=1)
       printf("\n Error reading intensity value %d.\n",i);
      ivals[i]=*sivals;
  }

  free(sivals);
  fclose(fp); 
}

#ifdef __cplusplus
}
#endif
