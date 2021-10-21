/*  Remove-Outlier-Sequences (version 0.9.5)
 *  Copyright 2018 by Christoph Mayer
 *
 *  This source file is part of the Remove-Outlier-Sequences program.
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
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  For any enquiries send an Email to Christoph Mayer
 *  c.mayer.zfmk@uni-bonn.de
 *
 *  When publishing work that is based on the results please cite:
 *  ....
 *  
 */
#include "global-types-and-parameters.h"
#include "tclap/CmdLine.h"
#include <string>
#include "faststring2.h"

using namespace TCLAP;
using namespace std;


bool file_exists(const char* fn)
{
  ifstream is;
  is.open(fn);
  if (is.fail())
  {
    return false;
  }
  is.close();
  return true;
}

void good_bye_and_exit(FILE *of, int error)
{
  if (error == 0)
    fprintf(of, "#\n#\n## Finished successfully. ##\n");
  fprintf(of, "\n" PROGNAME " says goodbye.\n\n");
  exit(error);
}

void good_bye_and_exit(ostream &os, int error)
{
  if (error == 0)
    os << "#\n#\n## Finished successfully. ##\n";
  os << "\n" PROGNAME " says goodbye.\n\n";
  exit(error);
}

void init_param(general_parameters &general_param, outlier_parameters &out_param)
{
  out_param.max_gap_proportion_in_window_in_non_gap_mode     =  0.5;
  out_param.windowsize                                       =  30;
  /**/  out_param.gap_open                                         = -1.53;   // Only used in gap scoring mode
  /**/  out_param.gap_ext                                          = -1.00;   // Only used in gap scoring mode
  /**/  out_param.factor_front_back_gaps                           =  0.5;   
  /**/  out_param.gap_match_score                                  =  5.63;   // mean blosum62 match score 
  /**/  out_param.model_and_gap_mode                               =  1;      // 1: blosum62
  out_param.outlier_IQR_factor                               =  1.5;      // 3: mask only extreme outliers, 1.5: mask moderate outliers
  out_param.minimum_IQR                                      =  0.563;    //3*5.63/out_param.windowsize;
  out_param.min_number_valid_seqs_in_window                  =  5;
  out_param.threshold_remove_whole_seq_prop_outlier_windows  =  0.5;

  general_param.VERBOSITY                                    = 1;
  general_param.rm_gap_ambig_lower_case_sites                = true;

  general_param.mask_mode                                    = 0;
  general_param.preserver_lower_case_for_input               = false;
  general_param.rm_versus_detect_only_mode                   = 0;      // 0: rm, 1: detect
}

void read_and_init_parameters(int argc, char** argv, ostream &logerr, general_parameters &general_param, outlier_parameters &out_param)
{
  init_param(general_param, out_param);

  bool b_masklower = true;
  bool b_maskambig = false;
  bool b_nomasking = false;

  try
  {
    CmdLine cmd("**********************************************************\n"
		"The " PROGNAME  " program has been designed to identify, mask and/or remove outlier sequences in molecular sequence alignments."
		"Accepted input files are: Fasta files of aligned amino acid sequences. "
		"Sequences can either be analysed as a whole or using a sliding window. "
		"In most cases a sliding window is recommended, since it is more sensitive if sequences are outliers only in parts of the alignment. "
		"\n\n"
		"The algorithm works as follows: Pairwise sequence scores (blosum62 scores) are computed for all sequences in the input file. From this, a mean score is determined for each sequece to all other sequences. For the score the median and the quartiles (Q1, Q2=median, Q3) are determined. The distance IQR=Q3-Q1 is called the inter-quartile range (IQR). With a parameter, IQR can be set to a minimum value, which is necessary for very similar sequences, since for almost identical sequences one or a few substitutions are already identified as outliers. If this is the case, it is recommended to increase the minimum value of IQR. A sequence or partial sequence in the sliding window is identified as an outlier, if its blosum score is less than Q1-IQR*IQR-factor, where IQR-factor is a parameter with a default value of 1.5. If you want to mask more sequences use a smaller value, e.g. 1.0, if you want to mask fewer sequences use a larger value, e.g 2.0. The default value of 1.5 corresponds to the general recommendations for the outlier detection in statistics. Preferred values depend on the data set. For moderately dissimilar sequences, larger values such as 2.0 often lead to better results. For very dissimilar sequences, outliers will only be identified if the value is reduced to 1.0. Sequences for which a specified proportion of the windows are outliers are removed entirely. Gap/ambig/lower case sites can be removed automatically. A nucleotide alignment correpsonding to the amino acid alignment can be specified and resitues will be masked/removed corresponding to the amino acid alignment."
		""
		,
		' ', VERSION);

    ValueArg<string> correesp_nuc_fasta_file_name_Arg("", "corresponding-nuc-fasta-file-name",
					 "Name of file corresponding to the specified amino acid file in which outliers are detected. The corresponding file has to correspond to the ....",
					 false, "", "string");
    cmd.add( correesp_nuc_fasta_file_name_Arg );


    // This flag might need changes in the code which determines the score.
    //    SwitchArg respect_input_case_Arg("", "respect-input-case",
    //			     "By default, input sequences are converted to upper case symbols. If With this argument, outlier sequence segments are masked with lower case amino acid or nucleotide symbols. Default: Do not mask.",
    //			     false);
    //    cmd.add( mask_lower_Arg );


    // TODO: Add gap related paramters, after gap scoring has been made insistent.


    SwitchArg rm_gap_ambig_lowercase_sites_Arg("R", "remove-gap-ambig-lowerCase-sites",
					       "With this argument, alignment sites are removed which contain only gap, ambig, and lower case sequence symbols. ",
					       false);
    cmd.add( rm_gap_ambig_lowercase_sites_Arg );


    ValueArg<float> IQR_factor_Arg("f", "IQR-factor",
				   "For a given window, outliers are determined by computing Q1-IQT*IQR-factor. Scores lower than this are considered as outliers. "
				   "Default: 1.5",
				   false, out_param.outlier_IQR_factor, "floating point number");
    cmd.add( IQR_factor_Arg );

    ValueArg<float> min_IQR_Arg("I", "minimum-IQR",
				"In conserved or almost invariant regions, a sequence is identified as an outlier even if it differs by only one or two residues. "
				"To avoid recognizing such a sequence as an outlier, a minimum value for the IQR can be specified. "
				"If the true IQR is smaller than the minimum IQR, the minimum IRQ is used instead of the real IQR. "
				"Default: 0.563",
				false, out_param.minimum_IQR, "floating point number");
    cmd.add( min_IQR_Arg );

    

    ValueArg<float> max_gap_prop_Arg("g", "maximum-gap-proportion",
				     "If a sequence segment of window size contains only or mostly gaps, its score to all other sequences could be used to identify it as an outlier. "
				     "To avoid problems und meaningless scores, we do not consider sequences, which have a gap proportion of a value larger than this parameter. The maximum-gap-proportion "
				     "parameter must be in the range 0..1. With a value of 0.5, if a sequence has more than 50% gaps in a window, the score of this sequence will not be considered "
				     "in computing the inter quartile distance nor can this sequence be determined as an outlier in the particular window. This sequence segment is also not considered when computing "
				     "the proportion of windows in which a sequence is an outlier. Default: 0.5", 
				     false, out_param.max_gap_proportion_in_window_in_non_gap_mode, "floating point number");
    cmd.add( max_gap_prop_Arg );

    ValueArg<int> windows_size_Arg("w", "windowsSize",
				   "The window size for the sliding window outlier detection. For a value of 0 no sliding window will be used, but the sequence will be analyses as a whole. Default: 30",
				   false, out_param.windowsize, "unsigned integer");
    cmd.add( windows_size_Arg );
 
    ValueArg<unsigned> verbosity_Arg("", "verbosity",
				     "The verbosity option controls the amount of information " PROGNAME 
				     " writes to the console while running. 0: Print only welcome message and essential error messages that lead to exiting the program. "
				     ,
				     //				     "1: report also wanrings, 2: report also progress, 3: report more detailed progress, >10: debug output. Maximum 10000: write all possible diagnostic output."
				     //				     " A value of 2 is required if startup parameters should be reported.\n",
				     false, general_param.VERBOSITY, "unsigned integer");
    cmd.add( verbosity_Arg );


    // SwitchArg mask_lower_Arg("m", "mask-lower-case",
    // 			     "With this argument, outlier sequence segments are masked with lower case amino acid or nucleotide symbols. Default: true, i.e. do lower case masking.",
    // 			     true);
    // cmd.add( mask_lower_Arg );


    SwitchArg no_masking_Arg("N", "no-masking",
			     "The default is to mask outliers with lower case symbols. With the -M option, residues are masked with the ambiguity code characters N (nucleotide) or X (amino acids). "
			     "With the -N option, masking can be turned off completely.",
			     false);
    cmd.add( no_masking_Arg );


    SwitchArg mask_Ambig_Arg("M", "mask-ambig",
			     "With this argument, outlier sequence segments of nucleotide sequences are masked with Ns and amino acid sequences with X. Default: mask with lower case symbols.",
			     false);
    cmd.add( mask_Ambig_Arg );

    ValueArg<int> min_num_valid_sequences_Arg ("", "minNumValidSeqsInWindow",
					       "A meaningful outlier test can normally only be carried out, if the number of sequences is large. This is not always the case. "
					       "In order to be able to compute quartiles and a IQR, we decided to require 5 valid sequences in a sliding window to be able to detect outliers. "
					       "If you think a different value for the minimum number of sequences is appropriate, it can be specified with this parameter. Default: 5. ",
					       false, out_param.min_number_valid_seqs_in_window, "integer");
    cmd.add( min_num_valid_sequences_Arg );

    ValueArg<float> threshold_remove_whole_seq_prop_outlier_windows_Arg ("e", "thresholdPropOutlierWindows2removeSequences",
				     "A sequence can be completely removed from an MSA if it is identified as an outlier in at least the specified proportion of sequence windows. "
									 "The default is 0.5. So a sequence which is identified as an outlier in 50% or more of the sequence windows, "
									 "it is removed completely. A value >1 implies, that sequences are never removed. ",
									 false, out_param.threshold_remove_whole_seq_prop_outlier_windows, "floating point number");
    cmd.add( threshold_remove_whole_seq_prop_outlier_windows_Arg );
    
    //// Not yet implemented:
    //    ValueArg<string> reference_fasta_file_name_Arg("r", "reference-sequence-file",
    //					 "Name of fasta file from which standard distances are determined. NOT YET IMPLEMENTED!!",
    //					 true, "", "string");
    //    cmd.add( reference_fasta_file_name_Arg );

    ValueArg<string> output_fasta_file_name_Arg("o", "output-file",
						"Name of output file.",
						true, "", "string");
    cmd.add( output_fasta_file_name_Arg );


    ValueArg<string> input_fasta_file_name_Arg("i", "input-file",
					 "Name of input file in fasta format.",
					 true, "", "string");
    cmd.add( input_fasta_file_name_Arg );

    

//************************************************
    cmd.parse( argc, argv );
//************************************************

//    faststring tmp, tmp1;

    // Assigning parameters to variables:
    general_param.input_fasta_file_name                                 = input_fasta_file_name_Arg.getValue().c_str();
    general_param.output_fasta_file_name                                = output_fasta_file_name_Arg.getValue().c_str();
    general_param.VERBOSITY                                             = verbosity_Arg.getValue();
    general_param.rm_gap_ambig_lower_case_sites                         = rm_gap_ambig_lowercase_sites_Arg.getValue();

    general_param.corresponding_nucleotide_alignment_file_name          = correesp_nuc_fasta_file_name_Arg.getValue().c_str(); 

    //    b_masklower                                                         = mask_lower_Arg.getValue();
    b_maskambig                                                         = mask_Ambig_Arg.getValue();
    b_nomasking                                                         = no_masking_Arg.getValue();

    cerr << "Value of b_maskambig: " << (int) b_maskambig << endl;  
    cerr << "Value of b_nomasking: " << (int) b_nomasking << endl; 
    
    out_param.outlier_IQR_factor                                        = IQR_factor_Arg.getValue();
    out_param.minimum_IQR                                               = min_IQR_Arg.getValue();
    out_param.windowsize                                                = windows_size_Arg.getValue();
    out_param.threshold_remove_whole_seq_prop_outlier_windows           = threshold_remove_whole_seq_prop_outlier_windows_Arg.getValue();
    out_param.min_number_valid_seqs_in_window                           = min_num_valid_sequences_Arg.getValue(); 
  }

  catch (ArgException &e)
  {
    logerr << "ERROR: " << e.error() << " for arg " << e.argId() << "\n";
    good_bye_and_exit(stdout, -1);
  }

  // Do some checks here:
  if (b_maskambig && b_nomasking)
  {
    cerr << "ERROR: The parameters: --no-maksing (-N) and --mask-ambig (-M) cannot "
            "both be specified. Please remove one of the two parameters. Exiting.\n";
    exit(-3);
  }

  if (b_nomasking) // turn off default value
  {
    b_masklower = false;
  }

  if (b_maskambig) // overwrite default value
    b_masklower = false;
  
  if (out_param.min_number_valid_seqs_in_window <=2)
  {
    cerr << "ERROR: The parameter value of the option -minNumValidSeqsInWindow has to be >2. Please modify this. Exiting.\n";
    exit(-3);
  }

  // This should not be possible. Might be removed later.
  if ( (b_masklower && b_nomasking) || (b_nomasking && b_maskambig) || (b_maskambig && b_masklower) )
  {
    cerr << "INTERNAL ERROR: At least two of the three exclusive internal masking state variables are true.\n";
    cerr << "b_masklower: " << (b_masklower ? "true":"false") << '\n';
    cerr << "b_maskambig: " << (b_maskambig ? "true":"false") << '\n';
    cerr << "b_nomasking: " << (b_nomasking ? "true":"false") << '\n'; 
    cerr << "Please report this problem. Exiting.\n";
    exit(-3);
  }

  if (b_masklower)
    general_param.mask_mode = 1;

  if (b_maskambig)
    general_param.mask_mode = 2;

  if (out_param.windowsize == 0)
    out_param.windowsize = INT_MAX;
}



//void print_parameters(FILE *of, const char *s, const outlier_parameters &out_param, const general_parameters &general_param)
//{
//     if (!global_blast_evalue_commandline.empty())
//     {
//       fprintf(of, "%sBlast evalue cut-off:                             %s\n", s, global_blast_evalue_commandline.c_str());
//     }
//     else
//     {
//       fprintf(of, "%sBlast evalue cut-off:                             <not specified>\n", s);
//     }

//     if (!global_blast_result_file.empty())
//     {
//       fprintf(of, "%sDo not conduct a new blast search, but use the    \n"
// 	          "%sfollowing existing blast result file:             %s\n",  s,s, global_blast_result_file.c_str() );
//     }


//   } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
//}


void print_parameters(ostream &os, const char *s, const general_parameters &general_param, const outlier_parameters &out_param)
{
  os << s <<   PROGNAME " uses the following command line parameters/default values:\n";
  os << s <<   "Input fasta file name:                                       " <<  general_param.input_fasta_file_name  << '\n';
  os << s <<   "Output fasta file name:                                      " <<  general_param.output_fasta_file_name << '\n';
  os << s <<   "Window size:                                                 " <<  out_param.windowsize                 << '\n';
  os << s <<   "IQR factor:                                                  " <<  out_param.outlier_IQR_factor         << '\n';
  os << s <<   "Minimum IRQ:                                                 " <<  out_param.minimum_IQR                << '\n';
  os << s <<   "Maximum proportion of gaps in sequences in window            " <<  out_param.max_gap_proportion_in_window_in_non_gap_mode << '\n';
  os << s <<   "Minimum number of valid sequences in window                  " <<  out_param.min_number_valid_seqs_in_window << '\n';
  os << s <<   "Threshold for prop. of outlier windows to remove whole seq.: " <<  out_param.threshold_remove_whole_seq_prop_outlier_windows << '\n';

  if ( out_param.model_and_gap_mode == 2)
  {
    os << s <<   "Include gaps into computing the score:                       yes\n";
    //    os << s <<   "Parameters specific to modes including gap score:\n";
    //    os << s <<   "   gap opening score                                         " <<  out_param.gap_open        << '\n';
    //    os << s <<   "   gap extension score                                       " <<  out_param.gap_ext         << '\n';
    //    os << s <<   "   gap match score                                           " <<  out_param.gap_match_score << '\n';
  }
  else
  {
    os << s <<   "Include gaps into computing the score:                       no\n";
  }

  if ( general_param.mask_mode == 0)
    os << s <<   "Mask outlier sequence segments mode:                         no maksing\n";
  else if  ( general_param.mask_mode == 1)
    os << s <<   "Mask outlier sequence segments mode:                         mask with lower case residues\n";
  else if  ( general_param.mask_mode == 2)
    os << s <<   "Mask outlier sequence segments mode:                         mask with Ns (DNA) or Xs (protein)\n";

  if ( out_param.model_and_gap_mode == 1)
    os << s <<   "Scoring method:                                              blosum62 scores.\n";
  else 
    os << s <<   "Scoring method:                                              UNKNOWN\n";

  if ( general_param.rm_gap_ambig_lower_case_sites)
    os << s <<   "Removing gap, ambig, lower case only sites:                  yes\n";
  else
    os << s <<   "Removing gap, ambig, lower case only sites:                  no\n";

  os << s <<   "Verbosity:                                                   " <<  general_param.VERBOSITY              << '\n';




//   if (global_mode == 'a')
//   {
//     os << s << "Filter mode:                                      Retain only the best locus per alignment file; criterion: minimize number of baits.\n";
//   }
//   else if (global_mode == 'A')
//   {
//     os << s << "Filter mode:                                      Retain only the best locus per alignment file; criterion: maximize sequence coverage in bait region.\n";
//   }
//   else if (global_mode == 'f')
//   {
//     os << s << "Filter mode:                                      Retain only the best locus per feature; criterion: minimize number of baits in bait region.\n";
//   }
//   else if (global_mode == 'F')
//   {
//     os << s << "Filter mode:                                      Retain only the best locus per feature; criterion: maximize sequence coverage in bait region.x\n";
//   }
//   else if (global_mode == 'b')
//   {
//     os << s << "Filter mode:                                      Blast filter: Remove all bait regions from output that belong to the ALIGNMENT in which one bait has at least two good hits to the reference genome.\n";
//   }
//   else if (global_mode == 'B')
//   {
//     os << s << "Filter mode:                                      Blast filter: Remove all bait regions from output that belong to the FEATUE in which one bait has at least two good hits to the reference genome.\n";
//   }
//   else if (global_mode == 'x')
//   {
//     os << s << "Filter mode:                                      Blast filter: Remove all bait REGIONs in which one bait has at least two good hits to the reference genome.\n";
//   }
//   else if (global_mode == 'C')
//   {
//     os << s << "Filter mode:                                      Blast filter: Conduct converage filter without any other blast filter.\n";
//   }
//   else if (global_mode == 't')
//   {
//     os << s << "Filter mode:                                      Thinning mode: Retain only every " << global_thinning_step_width << " th bait region. Choose the starting offset such that the number of baits is minimized.\n";    
//   }
//   else if (global_mode == 'T')
//   {
//     os << s << "Filter mode:                                      Thinning mode: Retain only every " << global_thinning_step_width << " th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n";
//   }
//   else if (global_mode == 'h')
//   {
//     os << s << "Filter mode:                                      Thinning mode new algorithm: Retain only every " << global_thinning_step_width  << " th bait region. Choose the starting offset such that the number of baits is minimized.\n";
//   }
//   else if (global_mode == 'H')
//   {
//     os << s << "Filter mode:                                      Thinning mode new algorithm: Retain only every " << global_thinning_step_width  << " th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n";
//   }
//   else if (global_mode == 'c')
//   {
//     if (global_conversion_mode == 1)
//     {
//       os << s << "Conversion mode:                                Convert bait file to four-column-upload format.\n";
//     }
//     else
//     {
//       os << s << "ERROR: Invalid internal conversion mode: " << global_conversion_mode << " .\nPlease report this bug.\n";
//       good_bye_and_exit(cout, -99);
//     }
//   }
//   else if (global_mode == 'S')
//   {
//     os << s << "Mode:                                             Prints stats mode.\n";
//   }
//   else
//   {
//     os << s << "ERROR: Invalid internal mode:                     " << global_mode << '\n';
//     good_bye_and_exit(cout, -99);
//   }

//   if (global_mode == 'c' && !global_ProbeID_prefix.empty())
//   {
//     os << s << "ProbeID suffix string in conversion mode:         " << global_ProbeID_prefix.c_str() << '\n';
//   }

//   os << s << "Verbosity:                                        " << global_verbosity << '\n';


//   if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
//   {
//     faststring str;

//     os << s << "Blast search settings and parameters:\n";
//     os << s << "Blast executable:                                 " << global_blast_exe.c_str() << '\n';
//     os << s << "Blast data base file name:                        " << global_blast_db.c_str() << '\n';

//     os << s << "Blast minimum hit coverage of bait in first hit:  " << global_blast_min_hit_coverage_of_bait << '\n';

//     os << s << "Blast maximum evalue of first/best hit:           " << global_blast_max_first_hit_evalue << '\n';
//     os << s << "Blast maximum evalue of second best hit:          " << global_blast_max_second_hit_evalue << '\n';

//     if (!global_blast_extra_commandline.empty())
//     {
//       os << s << "Blast extra commandline argument:                 " << global_blast_extra_commandline.c_str() << '\n';
//     }
//     else
//     {
//       os << s << "Blast extra commandline argument:                 <not specified>\n";
//     }

//     if (!global_blast_evalue_commandline.empty())
//     {
//       os << s << "Blast evalue cut-off:                             " << global_blast_evalue_commandline.c_str() << '\n';
//     }
//     else
//     {
//       os << s << "Blast evalue cut-off:                             <not specified>\n";
//     }

//     if (!global_blast_result_file.empty())
//     {
//       os << s << "Do not conduct a new blast search, but use the    " << '\n'
// 	 << s << "following existing blast result file:             " << global_blast_result_file.c_str()  << '\n';
//     }
//   } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
}



void print_calling_command_line(FILE *of, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  fprintf(of, "\nThe BaitFilter program was called with the following command line:\n");
  fputs(cmd.c_str(), of);
  fputc('\n',of);  
}


void print_calling_command_line(ostream &os, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  os << "\nThe BaitFilter program was called with the following command line:\n";
  os << cmd.c_str() << '\n';  
}



// void print_output_creation_message(FILE *of, char *s)
// {
//   fprintf(of, "%sOutput created by " PROGNAME ", version " VERSION "\n", s);
//   fprintf(of, "%s\n", s);
//   print_parameters(of, s);
//   fprintf(of, "%s\n", s);
// }

// void print_output_creation_message(ostream &os, char *s)
// {
//   os << s << "Output created by " PROGNAME ", version " VERSION "\n";
//   os << s << '\n';
//   print_parameters(os, s);
//   os << s << '\n';
// }
