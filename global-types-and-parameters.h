/*  Remove-Outlier-Sequences (version 0.9.6)
 *  Copyright 2023 by Christoph Mayer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BaitFisher.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For any enquiries send an Email to Christoph Mayer
 *  c.mayer.zfmk@uni-bonn.de
 *
 *  When publishing work that is based on the results please cite:
 *  ....
 *  
 */
#ifndef GLOBAL_TYPES_AND_PARAMETERS_H
#define GLOBAL_TYPES_AND_PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include "faststring2.h"
#include <climits>
#include <fstream>

#define DEBUG

#define PROGNAME "OliInSeq"
#define VERSION  "0.9.6"

struct outlier_parameters
{
  float     max_gap_proportion_in_window_in_non_gap_mode;
  int       windowsize;
  float     gap_open;
  float     gap_ext;
  float     factor_front_back_gaps;
  float     gap_match_score; //  = 5.63;   // mean blosum62 match score 
  int       model_and_gap_mode;
  float     outlier_IQR_factor;
  float     minimum_IQR;

  int       min_number_valid_seqs_in_window;

  float     threshold_remove_whole_seq_prop_outlier_windows;
};

struct general_parameters
{
  faststring input_fasta_file_name;
  faststring output_fasta_file_name;
  unsigned   VERBOSITY;
  bool       rm_gap_ambig_lower_case_sites;
  char       mask_mode; // 0: no masking, 1: lower case masking, 2: X masking
  bool       preserver_lower_case_for_input;

  faststring corresponding_nucleotide_alignment_file_name;
  int        rm_versus_detect_only_mode;   

};


#define macromax(x,y) ((x)<(y) ? (y) : (x))
#define macromin(x,y) ((x)<(y) ? (x) : (y))

#ifdef  DEBUG
#define DEBUGOUT1(x)        std::cerr << x                << '\n';
#define DEBUGOUT2(x,y)      std::cerr << x << y           << '\n';
#define DEBUGOUT3(x,y,z)    std::cerr << x << y << z      << '\n';
#define DEBUGOUT4(x,y,z,w)  std::cerr << x << y << z << w << '\n';
#else
#define DEBUGOUT1(x)
#define DEBUGOUT2(x,y)
#define DEBUGOUT3(x,y,z)
#define DEBUGOUT4(x,y,z,w)
#endif


void good_bye_and_exit(FILE *of, int);
void good_bye_and_exit(std::ostream&, int);
void init_param(general_parameters &, outlier_parameters &);
void read_and_init_parameters(int argc, char** argv, std::ostream &, general_parameters &, outlier_parameters &);
//void print_parameters(FILE*, const char *s,         const general_parameters &, const outlier_parameters &);
void print_parameters(std::ostream&, const char *s, const general_parameters &, const outlier_parameters &);
void print_calling_command_line(FILE*,unsigned, char**);
void print_calling_command_line(std::ostream&,unsigned, char**);

#endif
