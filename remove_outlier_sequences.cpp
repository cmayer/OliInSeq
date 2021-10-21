#include "distances_scores.h"
#include <iostream>
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include "print_container.h"
#include <vector>
#include <map>
#include <iomanip>

#include "CFile/CFile2_1.h"
#include "faststring2.h"
#include "statistic_functions.h"
#include "global-types-and-parameters.h"

using namespace std;

unsigned VERBOSITY;


int one_round_of_outlier_removal(outlier_parameters out_param, CSequences2 **pseqs, vector<faststring> &seqnames_removed, vector<int> &below_or_eq_max_gap_proportion_in_mode_1)
{
  CSequences2 *seqs = *pseqs;

  CDists_collection *p_sc_all   = NULL;

  if (out_param.model_and_gap_mode == 1) // blosum62
  {
    p_sc_all = new CDists_collection (*seqs, out_param.gap_open, out_param.gap_ext, out_param.gap_match_score, out_param.factor_front_back_gaps, 0);
  }
  else if (out_param.model_and_gap_mode == 2) // blosum62+gaps
  {
    cerr << "This mode has to be revised." << endl;
    exit(-13);
    p_sc_all = new CDists_collection (*seqs, out_param.gap_open, out_param.gap_ext, out_param.gap_match_score, out_param.factor_front_back_gaps, 1);
  }
  else
  {
    cerr << "Internal error: Unknown mode in line 37." << endl;
    exit(-13);
  }

  std::vector<float> mean_dists;
  p_sc_all->get_vector_of_mean_dists_to_others(mean_dists);
  // p_sc_all->get_dists(mean_dists);

  std::vector<float> mean_scores_for_stats;
  mean_scores_for_stats.reserve(mean_dists.size());

  for (unsigned i=0; i<below_or_eq_max_gap_proportion_in_mode_1.size(); ++i)
  {
    if ( below_or_eq_max_gap_proportion_in_mode_1[i] )
    {
      mean_scores_for_stats.push_back(mean_dists[i]);
    }
  }

  if (VERBOSITY >= 10)
  {
    cout.setf(ios::fixed);
    cout.precision(2);
 
    cout << endl;
    cout << "Mean blosum score of sequences to all other sequences: " << endl;

    for (unsigned i=0; i<mean_dists.size(); ++i)
    {
      cout.width(20);
      cout  << std::left << seqs->get_seq_by_index(i)->getFullName() << " "
	    << mean_dists[i] << " " << below_or_eq_max_gap_proportion_in_mode_1[i] << endl;
    }
    cout << endl;

    if (VERBOSITY >= 50)
    {
      cout << "Scores for stats: " << endl;
      for (unsigned i=0; i<mean_scores_for_stats.size(); ++i)
      {
	cout  << mean_scores_for_stats[i] << ", ";
      }
      cout << endl;
    }
  }

  double Q1, Q2, Q3, lB, uB;

  vec_median_quartile_outlier_bounds_method3(mean_scores_for_stats, Q1, Q2, Q3, lB, uB, out_param.outlier_IQR_factor, out_param.minimum_IQR);

  if (VERBOSITY >= 10)
  {
    cerr << "Quartile 1: " << Q1 << endl;
    cerr << "Quartile 2: " << Q2 << endl;
    cerr << "Quartile 3: " << Q3 << endl;
    cerr << "Lower outlier fence for the blosum score: " << lB << endl;
    cerr << "Sequences with a blosum score lower than this will be considered to be outliers. " << lB << endl;
  //  cerr << "uB: " << uB << endl;
  }

  vector<short int> outlier_list; // Indices of outliers found in this round. Indices only valid in this round.

  int i1 = mean_dists.size()-1;
  while ( i1 >= 0 )
  {
    // Small Blosum mean scores correspond to a low similarity.
    // A score below lB is considered to be an outlier.    
    //    cerr << "Chekcing; " << i1 << " " << mean_dists[i1] << endl;

    if (mean_dists[i1] < lB && below_or_eq_max_gap_proportion_in_mode_1[i1])
    {
      //      cerr << "O: "<< i1 << endl;
      outlier_list.push_back(i1);              // Store index of outlier. Index only valid in this round.
      // Store the name of the outlier, since indices are not unique among different rounds.
      seqnames_removed.push_back(seqs->get_seq_by_index(i1)->getFullName());
    }
    --i1;
  }

  if (VERBOSITY >= 20)
  {
    cerr << "Removing outliers:" << endl;
    cerr << "List of outlier sequences with score below lower outlier fence:" << endl;
    i1 = outlier_list.size()-1;
    while (i1>=0)
    {
      cerr << ":" << outlier_list[i1]+1 << ": " << seqs->get_seq_by_index(outlier_list[i1])->getFullName() << endl;    
      --i1;
    }
  }

  for (unsigned i2=0; i2 < outlier_list.size(); ++i2)
  {
    // The two following erase commands keep the list of sequences and the
    // below_or_eq_max_gap_proportion_in_mode_1 vector in sync.
    seqs->remove_del(outlier_list[i2]);
    below_or_eq_max_gap_proportion_in_mode_1.erase(below_or_eq_max_gap_proportion_in_mode_1.begin()+outlier_list[i2]);
  }

  if (VERBOSITY >= 20)
  {
    cerr << endl;
    cerr << "Removing gap only positions." << endl;
  }

  CSequences2 *seqs_rmgo = new CSequences2(CSequence_Mol::protein);
  seqs->alignment_without_gap_only_positions(*seqs_rmgo);

  if (p_sc_all)
    delete p_sc_all;

  // Here we delete the original alignemnt and replace it with the modified alignment.
  delete *pseqs;

  *pseqs = seqs_rmgo;

  return outlier_list.size();
}

// This function should obtain a copy of the sequences, if outliers should not be remove from the **pseqs object.
void detect_and_remove_outliers_from_CSeqs(outlier_parameters out_param, CSequences2 **pseqs, vector<faststring> &seqnames_removed,  map<faststring, unsigned> &count_valid_seqs_in_window)
{
  CSequences2 *seqs;
  seqs = *pseqs;

  unsigned Nseqs = seqs->GetTaxaNum();

  vector<int> below_or_eq_max_gap_proportion_in_mode_1(Nseqs);
  int         num_valid_seqs_in_window = 0;

  // Determine gap proportions for all sequences
  for (unsigned i=0; i<Nseqs; ++i)
  {
    if (seqs->get_gap_proportion(i) > out_param.max_gap_proportion_in_window_in_non_gap_mode )
    {
      below_or_eq_max_gap_proportion_in_mode_1[i] = false;
    }
    else
    {
      ++num_valid_seqs_in_window;
      below_or_eq_max_gap_proportion_in_mode_1[i] = true;
    }
  }
  
//   // Code to investigate
//   bool Aligator_good_window = false;
//   if (below_or_eq_max_gap_proportion_in_mode_1[15])
//     Aligator_good_window = true;


  // Nothing needs to be done if we have an insufficiently mumber of valid sequences
  if (num_valid_seqs_in_window < out_param.min_number_valid_seqs_in_window)
  {
    return;
  }

  // Count total number of valid sequences in this window:
  for (unsigned i=0; i<Nseqs; ++i)
  {
    if (!below_or_eq_max_gap_proportion_in_mode_1[i] )
    {
      if (VERBOSITY >= 50)
	cout << seqs->get_seq_by_index(i)->getFullName() << " -" << endl; 
    }
    else
    {
      const faststring &s = seqs->get_Seq_FullName(i);
      add_or_count(count_valid_seqs_in_window, s );
      if (VERBOSITY >= 50)
	cout << seqs->get_seq_by_index(i)->getFullName() << " +" << endl; 
      if (VERBOSITY >= 100)
	cout << "Intermediate results: " << s << " " << count_valid_seqs_in_window[s] << endl; 
    }
  }

  // Detect the outliers: This modifies below_or_eq_max_gap_proportion_in_mode_1, so we cannot use this array for other purposes any more.
  int round = 1;
  int number_removed = one_round_of_outlier_removal(out_param, pseqs, seqnames_removed, below_or_eq_max_gap_proportion_in_mode_1);

  if (VERBOSITY >= 50)
    cerr << "Number of outliers removed in round " << round << ": " << number_removed << endl;

  while (number_removed > 0)
  {
    ++round;
    number_removed = one_round_of_outlier_removal(out_param, pseqs, seqnames_removed, below_or_eq_max_gap_proportion_in_mode_1);
    if (VERBOSITY >= 50)
      cerr << "Number of outliers removed in round " << round << ": " << number_removed << endl;
  }

//   // Code to investigate special case:
//   if (Aligator_good_window)
//   {
//     bool Alig_is_outlier=false;
//     for (unsigned i=0; i<seqnames_removed.size(); ++i)
//     {
//       if (seqnames_removed[i] == "Alligator_mi")
// 	Alig_is_outlier = true;
//     }
//     if (!Alig_is_outlier)
//     cerr << "FOUND window in which Alligator is no outlier. " << endl;
//   }

}

// Divides the alignment into windows and calls the detect_outliers_from_CSeqs for each windows.
// Reports results 
void sliding_window_check_outlier(outlier_parameters out_param, int rm_vs_detect_mode, CSequences2 *seqs_full_alignment, map<faststring, unsigned> &outlier_counter_map, map<faststring, vector<unsigned> > &outlier_coord_map, map<faststring, float> &outlier_windows_prop) 
{
  //  bool                      valid_seqgment;
  int                       alignment_length = seqs_full_alignment->GetPosNum();
  vector<faststring>        seqnames_removed;
  map<faststring, unsigned> count_valid_seqs_in_window;

  (void) rm_vs_detect_mode; // THIS MODE IS NOT IMPLEMENTED IN THIS FUNCTION.

  // For all starting points of windows:


  // If the sequence is shorter than the window size then we either remove the whole sequence or nothing.
  // We do not consider partial sequences at the beginning or the end.
  // In particular this is the case of the user specified a window size of 0 in order to deliberately for a whole sequence analysis.
  if (out_param.windowsize > alignment_length )
  {
   if (VERBOSITY >= 10)
      {
	cout << "******************************************************" << endl;
	cout << "Working on full siquence: "                             << endl;
	cout << "******************************************************" << endl;
      }

      CSequences2 *seqs_partial = new CSequences2(*seqs_full_alignment, 0, alignment_length);
      detect_and_remove_outliers_from_CSeqs(out_param, &seqs_partial, seqnames_removed, count_valid_seqs_in_window);

      vector<faststring>::iterator v_it=seqnames_removed.begin(), v_it_end=seqnames_removed.end();
      while (v_it != v_it_end)
      {
	add_or_count(outlier_counter_map, *v_it);
	outlier_coord_map[*v_it].push_back(0);
	++v_it;
      }
      seqnames_removed.clear();
      delete seqs_partial;
      if (VERBOSITY >= 10)
	cout << endl;
  }
  else
  {
    // For each starting position:
    for (int i=0; i<alignment_length-out_param.windowsize; ++i)
    {
      if (VERBOSITY >= 10)
      {
	cout << "******************************************************" << endl;
	cout << "Working on sliding window coordinate: " << i+1 << endl;
	cout << "******************************************************" << endl;
      }

      CSequences2 *seqs_partial = new CSequences2(*seqs_full_alignment, i, i+out_param.windowsize);
      detect_and_remove_outliers_from_CSeqs(out_param, &seqs_partial, seqnames_removed, count_valid_seqs_in_window);

      vector<faststring>::iterator v_it=seqnames_removed.begin(), v_it_end=seqnames_removed.end();
      while (v_it != v_it_end)
      {
	add_or_count(outlier_counter_map, *v_it);
	outlier_coord_map[*v_it].push_back(i);
	++v_it;
      }
      seqnames_removed.clear();
      delete seqs_partial;
      if (VERBOSITY >= 10)
	cout << endl;
    }
  }

  map<faststring, unsigned>::iterator m_it, m_it_end;

  if (VERBOSITY >= 1)
  {
    cout << "Outlier stats for all species:" << endl;
    cout << "Species\t valid windows\toutlier windows\tproportion" << endl;
  }

  m_it=outlier_counter_map.begin();
  m_it_end=outlier_counter_map.end();
  unsigned max_seq_name_length = 0;
  while (m_it != m_it_end)
  {
    if ( (m_it->first).length() > max_seq_name_length )
      max_seq_name_length = (m_it->first.length());
    ++m_it;
  }

  //  cerr << "Max length: " << max_seq_name_length << endl;

  m_it=outlier_counter_map.begin();
  m_it_end=outlier_counter_map.end();
  while (m_it != m_it_end)
  {
    outlier_windows_prop[m_it->first] = (double)m_it->second/count_valid_seqs_in_window[ m_it->first]; 
    // VALID windows OUTPUT and SCORES
    vector<unsigned> &local_v = outlier_coord_map[m_it->first];

    if (VERBOSITY >= 1)
    {
      cout.setf(ios::fixed);
      cout.precision(3);

      cout  << setw(max_seq_name_length) << m_it->first << "\t " << setw(8) << count_valid_seqs_in_window[ m_it->first] << "\t " << setw(8) << m_it->second << "\t "
	   << setw(8) << (double)m_it->second/count_valid_seqs_in_window[ m_it->first] << endl;

      //      cout << m_it->first << " valid windows " << setw(8) << count_valid_seqs_in_window[ m_it->first] << "   outlier windows " << setw(8) << m_it->second << " proportion "
      //	   << setw(8) << (double)m_it->second/count_valid_seqs_in_window[ m_it->first] << endl;

      if (VERBOSITY >= 2)
      {
	cout << "Coordinates: ";
	for (unsigned j=0; j < local_v.size(); ++j)
	{ 
	  cout <<  local_v[j]+1 << ", ";
	}
	cout << endl << endl;
      }
    }

    ++m_it;
  }
}


void remove_and_mask_outliers(general_parameters general_param, outlier_parameters out_param, CSequences2 *seqs, map<faststring, vector<unsigned> > &outlier_coord_map, map<faststring, float> &outlier_windows_prop)
{
  for (int i=seqs->GetTaxaNum()-1; i >= 0; --i)
  {
    CSequence_Mol  *seq    = seqs->get_seq_by_index(i);
    faststring     seqname = seq->getFullName();
    char           ambig   = 'X';

    if (outlier_windows_prop[seqname] > out_param.threshold_remove_whole_seq_prop_outlier_windows)
    {
       seqs->remove_del(i);
    }
    else
    {
      vector<unsigned> &coords = outlier_coord_map[seqname];
      for (unsigned j=0; j < coords.size(); ++j)
      {
	char *p1 = const_cast<char *>(seq->getSeqBegin()+coords[j]);
	char *p2 = p1+out_param.windowsize;

	//	cout << "Value of general_param.mask_mode: " << (int) general_param.mask_mode << endl;

	if (general_param.mask_mode == 1)
	{
	  //	  cout << "Calling lower masking." << endl;
	  seq->maskSeq_tolower(p1,p2);
	}
	else if (general_param.mask_mode == 2)
	{
	  //	  cout << "Calling X masking." << endl;
	  seq->maskSeq(p1,p2, ambig);
	}
      }
    }
  }
  /*
  ofstream os(general_param.output_fasta_file_name.c_str());
  if (os.fail())
  {
    cerr << "Could not open output file. Exiting." << endl;
    exit(-33);
  }

  CSequences2 *seqs_rmgo = new CSequences2(CSequence_Mol::protein);
  if (rm_gap_ambig_lower_case_sites)
    seqs->alignment_without_gap_N_lowerCase_only_positions(*seqs_rmgo);
  else
    seqs->alignment_without_gap_only_positions(*seqs_rmgo);

  seqs_rmgo->ExportSequences(os, 'f', 1000000000);
  //  seqs->ExportSequences(os, 'f', 1000000000);
  os.close();
  
  delete seqs;
  delete seqs_rmgo;
  */
}

void remove_and_mask_outliers_in_a_corresponding_nucleotide_alignment(general_parameters general_param, outlier_parameters out_param, CSequences2 *prot_seqs, CSequences2 *nuc_seqs, map<faststring, vector<unsigned> > &outlier_coord_map, map<faststring, float> &outlier_windows_prop)
{
  if (general_param.corresponding_nucleotide_alignment_file_name.empty())
    return;

  CFile         file;
  const char *  fn = general_param.corresponding_nucleotide_alignment_file_name.c_str();

  if (VERBOSITY >= 1)
    cerr << "Reading corresponding nucleotide fasta file: " << fn << endl;

  file.ffopen(fn);

  if (file.fail())
  {
    cerr << "Could not open specified file: " << fn << endl;
    exit(-1);
  }

  nuc_seqs = new CSequences2(CSequence_Mol::dna);

  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(CSequence_Mol::convert_toupper);
  nuc_seqs->read_from_Fasta_File(file, pflag, 0, -1, false);
  file.ffclose();

  if (prot_seqs->GetTaxaNum() != nuc_seqs->GetTaxaNum())
  {
    cerr << "ERROR: The number of sequences in the corresponding nucleotide file differs from the number of sequences in the amino acid file." << endl; 
    exit(-17);
  }

  if (VERBOSITY >= 1)
  {
    cerr << "Found this number of taxa in corresponding nucleotide file:      " << nuc_seqs->GetTaxaNum()  << endl;
  }

  for (int i = prot_seqs->GetTaxaNum()-1; i >= 0; --i)
  {
    CSequence_Mol  *prot_seq = prot_seqs->get_seq_by_index(i); // Use the names of the original Sequence file.
    faststring     seqname   = prot_seq->getFullName();

    CSequence_Mol  *nuc_seq  = nuc_seqs->get_seq_by_index(i);
    char           ambig     = 'N';
    faststring     seqname2  = nuc_seq->getFullName();

    // Check sequence names for being equal:
    if (seqname != seqname2)
    {
      cerr << "Warning: Sequence names do not correspond between protein and corresponding nucleotide sequences. "
	      "Sequence number: "
	   << i << endl;
    }
    if (prot_seq->length()*3 != nuc_seq->length() )
    {
      cerr << "ERROR: Sequence lengths do not correspond for protein sequences and corresponding nucleotide sequence file. "
	      "Sequence number: "
	   << i << endl;
      exit(-5);
    }
    
    if (outlier_windows_prop[seqname] > out_param.threshold_remove_whole_seq_prop_outlier_windows)
    {
       nuc_seqs->remove_del(i);
    }
    else
    {
      vector<unsigned> &coords = outlier_coord_map[seqname];
      for (unsigned j=0; j < coords.size(); ++j)
      {
	// Since coordinates are 0 based, we simply multiply be 3:
	const char *p1 = const_cast<char *>( nuc_seq->getSeqBegin() + coords[j]*3);
	const char *p2 = p1+out_param.windowsize*3;

	const char *seq_end = nuc_seq->getSeqEnd();
	if (p2 > seq_end)
	{
	  cerr << "WARNING: End coordinate when masking a partial sequence in the corresponding DNA sequence is "
	          "beyond the end of the sequence. Overhang: " << p2 - seq_end  << endl;
	  cerr << "A reason could be that wrong sequence files have been specified. Check whether the file specified by -i is a protein file "
	          "and the file specified by --corresponding-nuc-fasta-file-name is a nuceotide file." << endl; 
	  p2 = seq_end;
	}
	if (p1 > seq_end)
	{
	  cerr << "ERROR: Start coordinate when masking a partial sequence in the corresponding DNA sequence is "
	          "beyond the end of the sequence. Overhang: " << p2 - seq_end  << endl;
	  p2 = seq_end;
	}


	//	cout << "Value of general_param.mask_mode: " << (int) general_param.mask_mode << endl;

	if (general_param.mask_mode == 1)
	{
	  //	  cout << "Calling lower masking." << endl;
	  nuc_seq->maskSeq_tolower(p1,p2);
	}
	else if (general_param.mask_mode == 2)
	{
	  //	  cout << "Calling X masking." << endl;
	  nuc_seq->maskSeq(p1,p2, ambig);
	}
      }
    }
  }

  /*
  faststring output_name = general_param.corresponding_nucleotide_alignment_file_name.filename_dirname_and_basename_without_extension();

  output_name += "_orm.fas";

  ofstream os(output_name.c_str());
  if (os.fail())
  {
    cerr << "Could not open output file. Exiting." << endl;
    exit(-33);
  }

  CSequences2 *nuc_seqs_rmgo = new CSequences2(CSequence_Mol::dna);
  if (rm_gap_ambig_lower_case_sites)
    nuc_seqs->alignment_without_gap_N_lowerCase_only_positions(*nuc_seqs_rmgo);
  else
    nuc_seqs->alignment_without_gap_only_positions(*nuc_seqs_rmgo);

  nuc_seqs_rmgo->ExportSequences(os, 'f', 1000000000);
  //  nuc_seqs->ExportSequences(os, 'f', 1000000000);
  os.close();

  delete nuc_seqs;
  delete nuc_seqs_rmgo;
  */
}


// Reads the sequences and calls the sliding_window_check_outlier function.
void run_outlier_determination(general_parameters general_param, outlier_parameters out_param, int rm_vs_detect_mode)
{
  const char *fn      = general_param.input_fasta_file_name.c_str();
  //  const char *out_fn  = general_param.output_fasta_file_name.c_str();

  CFile         file;
  vector<float> score_all;
  vector<float> mean_score_to_other_sequences;

  if (VERBOSITY >= 1)
    cerr << "Reading fasta file: " << fn << endl;

  file.ffopen(fn);

  if (file.fail())
  {
    cerr << "Could not open specified file: " << fn << endl;
    exit(-1);
  }

  CSequences2 *prot_seqs = NULL;
  CSequences2 *nuc_seqs  = NULL;
  prot_seqs = new CSequences2(CSequence_Mol::protein);

  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(CSequence_Mol::convert_toupper);
  prot_seqs->read_from_Fasta_File(file, pflag, 0, -1, false);
  file.ffclose();

  if (VERBOSITY >= 1)
  {
    cerr << "Found this number of taxa:      " << prot_seqs->GetTaxaNum()  << endl;
    cerr << "Found this number of positions: " << prot_seqs->GetPosNum()  << endl;
  }

  map<faststring, unsigned>          outlier_counter_map;
  map<faststring, vector<unsigned> > outlier_coord_map;    // The coordinates we store are 0 based.
  map<faststring, float>             outlier_windows_prop; // Proportion of outlier windows among all valid windows for each sequence

  sliding_window_check_outlier(out_param, rm_vs_detect_mode, prot_seqs, outlier_counter_map,
			       outlier_coord_map, outlier_windows_prop);
  
  if (rm_vs_detect_mode == 0) // Remove mode
  {
    // Do we have to handle a corresponding nucleotide alignment?
    // This has to be done before handling the prot_seqs since sequences can be removed in this process
    // and the number of sequences in prot_seqs is not allowed to change before handling the nucleotides.
    
    if ( !general_param.corresponding_nucleotide_alignment_file_name.empty() )
    {
      CFile         file;
      const char   *fn = general_param.corresponding_nucleotide_alignment_file_name.c_str();

      if (VERBOSITY >= 1)
	cerr << "Reading corresponding nucleotide fasta file: " << fn << endl;

      file.ffopen(fn);

      if (file.fail())
      {
	cerr << "Could not open specified file: " << fn << endl;
	exit(-1);
      }


      nuc_seqs = new CSequences2(CSequence_Mol::dna);

      CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(CSequence_Mol::convert_toupper);
      nuc_seqs->read_from_Fasta_File(file, pflag, 0, -1, false);
      file.ffclose();

      if (prot_seqs->GetTaxaNum() != nuc_seqs->GetTaxaNum())
      {
	cerr << "ERROR: The number of sequences in the corresponding nucleotide file differs from the number "
	        "of sequences in the amino acid file." << endl; 
	exit(-17);
      }

      if (VERBOSITY >= 1)
      {
	cerr << "Found this number of taxa in corresponding nucleotide file:      " << nuc_seqs->GetTaxaNum()  << endl;
      }

      remove_and_mask_outliers_in_a_corresponding_nucleotide_alignment(general_param,
								       out_param, prot_seqs, nuc_seqs,
								       outlier_coord_map, outlier_windows_prop);

      // FINISHED: Removed and masked outlier in nuc_seqs
    } // END  if ( !general_param.corresponding_nucleotide_alignment_file_name.empty() )

    // In any case we remove and mask outliers in prot_seqs
    remove_and_mask_outliers(general_param, out_param, prot_seqs, outlier_coord_map, outlier_windows_prop);

    // Determine gap only positions in the amino acid alignment:
    vector<bool> vector_gap_ambig_lowerCase_only_positions;

    if (general_param.rm_gap_ambig_lower_case_sites)
      prot_seqs->get_vector_gap_ambig_lowerCase_only_positions(vector_gap_ambig_lowerCase_only_positions);
    else
      prot_seqs->get_vector_gap_only_positions(vector_gap_ambig_lowerCase_only_positions);
      
    CSequences2 *prot_seqs_rmgo = new CSequences2(CSequence_Mol::protein);

    prot_seqs->alignment_without_specified_positions(*prot_seqs_rmgo, vector_gap_ambig_lowerCase_only_positions);

    ofstream os(general_param.output_fasta_file_name.c_str());
    if (os.fail())
    {
      cerr << "Could not open output file. Exiting." << endl;
      exit(-33);
    }
    prot_seqs_rmgo->ExportSequences(os, 'f', 1000000000);
    os.close();
 
    delete prot_seqs;
    delete prot_seqs_rmgo;

    if (nuc_seqs)
    {
      CSequences2 *nuc_seqs_rmgo = new CSequences2(CSequence_Mol::dna);

      nuc_seqs->alignment_without_specified_positions_prot2nuc(*nuc_seqs_rmgo, vector_gap_ambig_lowerCase_only_positions);
      
      faststring output_name = general_param.corresponding_nucleotide_alignment_file_name.filename_dirname_and_basename_without_extension();
      output_name += "_orm.fas";

      ofstream os(output_name.c_str());
      if (os.fail())
      {
	cerr << "Could not open output file. Exiting." << endl;
	exit(-33);
      }

      /*
      if (rm_gap_ambig_lower_case_sites)
	nuc_seqs->alignment_without_gap_N_lowerCase_only_positions(*nuc_seqs_rmgo);
      else
	nuc_seqs->alignment_without_gap_only_positions(*nuc_seqs_rmgo);
      */

      nuc_seqs_rmgo->ExportSequences(os, 'f', 1000000000);
      os.close();

      delete nuc_seqs_rmgo;
      delete nuc_seqs;
    }
    
  }
  else // Detect mode
  {
    cerr << "Detect mode not fully implemented. Only remove and masking mode has been implemented." << endl;
    exit(-2);
  }
}

// Reads parameters and calls the run_outlier_determination function
int main(int argc, char ** argv)
{
  outlier_parameters  out_param;
  general_parameters  general_param;

  int        rm_vs_detect_mode=0;


  //  general_param.input_fasta_file_name  = "/Daten/Programmieren/C++/Programme/remove_outlier_sequences/remove_outlier_sequences-v0.9.2/example-data/extrem/EOG0907000P_pep_dialign.fas";
  //  general_param.output_fasta_file_name = "out.fas";
  //  init_param(general_param, out_param);

  read_and_init_parameters(argc, argv, cerr, general_param, out_param);
  VERBOSITY = general_param.VERBOSITY;

  if (VERBOSITY >=1)
    print_parameters(cerr, "", general_param, out_param);



//   faststring mode         = argv[1];
//   faststring infile       = argv[2];
//   faststring resultfile   = argv[3];
//   faststring special_mode = 0;

//   if (argc == 5)
//     special_mode = argv[4];

//   if (special_mode == "detect-only")
//   {
//     rm_vs_detect_mode = 7;
//   }
//   else
//   {
//     rm_vs_detect_mode = 0;
//   }
  
//   if (mode == "blosum62")
//   {
//     out_param.model_and_gap_mode = 1;

//   }
//   else if (mode == "blosum62+gaps")
//   {
//     out_param.model_and_gap_mode = 2;

//   }
//   else
//   {
//     cerr << "Unknown mode: " << mode << endl;
//     cerr << "Usage: " << argv[0] << " (dna|dna+gpas|blosum62|blosum62+gaps) fasta-input-filename fasta-output-filename [detect-only]." << endl;
//   }

  run_outlier_determination(general_param, out_param, rm_vs_detect_mode);

  if (VERBOSITY >= 1)
    cerr << "Finished detecting outlier sequences in the file " << general_param.input_fasta_file_name << endl << endl;
}
