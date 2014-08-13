/*
 * io.cpp
 * Source Code to ildg Input/Output functions:
 *  - read in timeslices i until f from file
 * Created on: Aug 26, 2013
 * Author:
 */
#define _FILE_OFFSET_BITS 64

#include<lime.h>
#include<complex>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "io.h"
#include "io_utils.h"

using std::complex;

#define MAXBUF 1048576

int read_lime_gauge_field_doubleprec_timeslices(double *config,
    const char *filename, const int T,
    const int LX, const int LY,
    const int LZ, const int slice_i,
    const int slice_f) {
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  double tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n",
          status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp("ildg-binary-data",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)) {
    if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)/2) {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n",
          (n_uint64_t)bytes, filename,
          (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double));
      fprintf(stderr, "Aborting...!\n");
      fflush( stdout );
      exit(501);
    }
    else {
      fclose(ifs);
      fprintf(stderr, "single precision read!\n");

      fprintf(stderr, "Not implemented!\n");
      exit(EXIT_FAILURE);

      return( read_lime_gauge_field_singleprec(config, filename, T, LX, LY, LZ) );
    }
  }

  bytes = (n_uint64_t)72*sizeof(double);

  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
        for(x = 0; x < LX; x++) {

          // check for endianess and reading in data
          // the pointer limereader is internally increased by bytes
          // in the limeReaderReadData function
          if(!words_bigendian) {
            status = limeReaderReadData(tmp, &bytes, limereader);
            byte_swap_assign(tmp2, tmp, 72);
          }
          else
            status = limeReaderReadData(tmp2, &bytes, limereader);
          // check if reading was successfull
          if(status < 0 && status != LIME_EOR) {
            fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
                status, filename);
            exit(500);
          }

          // we just want to read in data in the specific range of timeslices
          // must be here because the limereader pointer must be increased correctly
          // could be done with much more performance but might be tricky to do correctly
          if(t<slice_i || t>slice_f)
            continue;

          // copy of link variables from tmp2 into config
          // ILDG has mu-order: x,y,z,t so it is changed here to: t,x,y,z !
          const size_t p = (size_t) ( ((t-slice_i)*LX*LY*LZ + x*LY*LZ + y*LZ + z) * 72); // position in config
          size_t k = 0;
          for(size_t mu = 1; mu <= 4; mu++) { // mu=4 is for the shift of U_t
            size_t index;
            if (mu != 4)
              index = p + mu*18; // for U_x, U_y and U_z
            else
              index = p; // U_t is copied into the beginning of
            // the (config+p) array

            for(size_t i = 0; i < 3; i++) {
              for(size_t j = 0; j < 3; j++) {
                config[index+6*i+2*j] = tmp2[2*k];
                config[index+6*i+2*j+1] = tmp2[2*k+1];
                k++;
              }
            }

          } // loop over mu ends here

        } // loop over position space ends here
      }
    }
  }
  if(status < 0 && status != LIME_EOR) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
        status, filename);
    exit(500);
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}



