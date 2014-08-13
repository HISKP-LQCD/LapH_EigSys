/*
 * io.h
 * Declares:
 * 	- read_lime_gague_field_doubleprec_timeslices
 *  Created on: Aug 26, 2013
 *      Author: helmes
 */

#ifndef IO_H_
#define IO_H_

int read_lime_gauge_field_doubleprec_timeslices(double *config,
    const char *filename,
    const int T, const int LX,
    const int LY, const int LZ,
    const int slice_i,
    const int slice_f);

int read_lime_gauge_field_singleprec(double *config, const char * filename,
    const int T, const int LX, const int LY,
    const int LZ);

#endif /* IO_H_ */
