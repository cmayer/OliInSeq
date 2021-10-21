#ifndef DISTANCES_SCORES_H
#define DISTANCES_SCORES_H

#include "faststring2.h"
#include "aa-matrices.h"

#include <iostream>
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include <vector>
#include "statistic_functions.h"
#include "basic-DNA-RNA-AA-routines.h"

// File in Klassen

int get_index_for_index_pair(int i, int j)
{
  // In all_aa_distances, all_aa_distances_gaps we push the distances in the order:
  //    (1,0), (2,0), (2,1), (3,0), (3,1), (3,2), (4,0), ...
  // They 0,    1,     2,     3,  
  return i*(i-1)/2+j;
}


// Fully implemented
inline void aaScore_BLOSUM62(const char* s1,
			      const char* s2,
			      unsigned N,
			      int      &blosum62_sum,
			      unsigned &len,
			      unsigned &unknown
			      )
{
  unsigned i;

  int      blosum62_sum_local  = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;

  char c1, c2;
  int  res;

  for (i=0; i<N; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( !(c1 == '-' || c2 == '-') )
    {
      res = get_BLOSUM62_score(c1, c2);

      //      std::cout << "c1 c2 dist " << c1 << " " << c2 << " " << res << std::endl;
 
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. " << c1 << " " << c2 << std::endl;
	++unknown_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << std::endl;
	blosum62_sum_local += res;
	++len_local;
      }
    }
  }
  blosum62_sum  = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}

inline void aaScore_BLOSUM62(const char* s1,
			      const char* s2,
			      unsigned range_beg,
			      unsigned range_end,
			      int      &blosum62_sum,
			      unsigned &len,
			      unsigned &unknown
			      )
{
  unsigned i;

  int      blosum62_sum_local  = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;

  char c1, c2;
  int  res;

  for (i=range_beg; i<range_end; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( !(c1 == '-' || c2 == '-') )
    {
      res = get_BLOSUM62_score(c1, c2);

      //      std::cout << "c1 c2 dist " << c1 << " " << c2 << " " << res << std::endl;
 
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. " << c1 << " " << c2 << std::endl;
	++unknown_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << std::endl;
	blosum62_sum_local += res;
	++len_local;
      }
    }
  }
  blosum62_sum  = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}


// Fully implemented
inline void aaScore_BLOSUM62(faststring &s1,
			     faststring &s2, 
			     int        &blosum62_sum,
			     unsigned   &len,
			     unsigned   &unknown
			     )
{
  if (s1.length() != s2.length())
  {
    blosum62_sum = 0;  // These values are indicative of this error or empty sequences.
    len           = 0;
    unknown       = 0;
    return;
  }
  aaScore_BLOSUM62(s1.c_str(), s2.c_str(), s1.length(), blosum62_sum, len, unknown);
}

// Fully implemented
inline void aaScore_BLOSUM62(faststring &s1,
			     faststring &s2,
			     unsigned   range_beg,
			     unsigned   range_end,
			     int        &blosum62_sum,
			     unsigned   &len,
			     unsigned   &unknown
			     )
{
  if (s1.length() < range_beg || s1.length() < range_end ||
      s2.length() < range_beg || s2.length() < range_end    )
  {
    blosum62_sum = 0;  // These values are indicative of this error or empty sequences.
    len           = 0;
    unknown       = 0;
    return;
  }
  aaScore_BLOSUM62(s1.c_str(), s2.c_str(), range_beg, range_end, blosum62_sum, len, unknown);
}

// not verified !!!!!!!!!!!!

// factor_front_back_gaps should be >0. Should not be 0.
// Scoring gaps:
// AA--AA:
// AA--AA: 0 gap-openings, 2 gap-matches
//
// ----AA--:
// --AAAA--: 1 gap-opeing, 4 gap-matches
//
// AA----AA:
// AAAA--AA: 1 gap-opening, 4 gap-matches


// Verify whether fully implemented.
inline void aaScore_BLOSUM62_gaps(const char* s1,
				   const char* s2,
				   unsigned N,
				   float    gap_open,
				   float    gap_ext,
		    		   float    gap_match_score,
				   float    factor_front_back_gaps,
				   float    &blosum62_sum,
				   unsigned &len,
				   unsigned &unknown
				   )
{
  float    blosum62_sum_local = 0;
  unsigned len_local           = 0;
  unsigned unknown_local       = 0;

  int      gap_s1 = 0;
  int      gap_s2 = 0;

  char c1, c2;
  int  res;

  unsigned start, end;

  bool gap_at_start_s1 = false;
  bool gap_at_start_s2 = false;

  bool in_gap;

  end   = N;

/*   // Determine start for comparison: */
/*   for (i=0; i< N; ++i) */
/*   { */
/*     if ( s1[i] != '-' || s2[i] != '-' ) */
/*       break; */
/*     else */
/*     { */
/*       blosum62_sum_local += gap_gap_match_score * factor_front_back_gaps; */
/*       ++len_local; */
/*     } */

/*   } */
/*   start = i; */
 
/*   for (i=N-1; i>=start; --i) */
/*   { */
/*     if ( s1[i] != '-' || s2[i] != '-' ) */
/*       break; */
/*     else */
/*     { */
/*       blosum62_sum_local += gap_gap_match_score * factor_front_back_gaps; */
/*       ++len_local; */
/*     } */
/*   } */
/*   end = i+1; // i is the position with */

  // In the range (start .. end-1) we start and end with a position that does
  // not have a gap in both sequences.

  start = 0;
  end   = N;

  if (start >= end)
  {
    blosum62_sum  = 0;
    len           = 0;
    unknown       = 0;
    return;
  }

  // Only one of the two ifs can be true:
  if (s1[start] == '-')
    gap_at_start_s1 = true;
  else if (s2[start] == '-')
    gap_at_start_s2 = true;

  for (unsigned i=start; i<end; ++i)
  {
    c1 = s1[i];
    c2 = s2[i];

    if ( c1 == '-' || c2 == '-' ) // => either gap-gap or gap-nongap or nongap-gap
    {
      if (c1 != c2)
      {
	++len_local; // We count this position
	if (c1 == '-')
	{
	  // This could be the end of a gap in s2:
	  if (gap_s2 > 0)
	  {
	    if (gap_at_start_s2)
	    {
	      blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
	      gap_s2 = 0;
	      gap_at_start_s2 = false;
	    }
	    else
	    {
	      blosum62_sum_local += (gap_open + gap_s2*gap_ext);
	      gap_s2 = 0;
	    }
	  }
	  if ( is_aa_blosum_code_index(c2) )
	  {
	    ++gap_s1;
	  }
	  else
	  {
	    ++unknown_local;
	    --len_local;
	  }
	}
	else // must be (c2 == '-')
	{
	  // This could be the end of a gap in s1:
	  if (gap_s1 > 0)
	  {
	    if (gap_at_start_s1)
	    {
	      blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
	      gap_s1 = 0;
	      gap_at_start_s1 = false;
	    }
	    else
	    {
	      blosum62_sum_local += (gap_open + gap_s1*gap_ext);
	      gap_s1 = 0;
	    }
	  }
	  if ( is_aa_blosum_code_index(c1) )
	  {
	    ++gap_s2;
	  }
	  else
	  {
	    ++unknown_local;
	    --len_local;
	  }
	}
      } // END if (c1 != c2)
      else // We have a gap in both sequences.
      {
	
      }
    }
    else // ! ( c1 == '-' || c2 == '-' ) => nongap-nongap
    {
      if (gap_s1 > 0)
      {
	if (gap_at_start_s1)
	{
	  blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
	  gap_s1 = 0;
	  gap_at_start_s1 = false;
	}
	else
	{
	  blosum62_sum_local += (gap_open + gap_s1*gap_ext);
	  gap_s1 = 0;
	}
      }
      else if (gap_s2 > 0) // DOUBLECHECK
      {
	if (gap_at_start_s2)
	{
	  blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
	  gap_s2 = 0;
	  gap_at_start_s2 = false;
	}
	else
	{
	  blosum62_sum_local += (gap_open + gap_s2*gap_ext);
	  gap_s2 = 0;
	}
      }

      ++len_local;
      res = get_BLOSUM62_score(c1, c2);
      if (res == 255)
      {
	std::cerr << "Found unexpected symbol(s). No score defined for these symbols. " << c1 << " " << c2 << std::endl;
	++unknown_local;
	--len_local;
      }
      else
      {
	//	std::cerr << c1 << " " << c2 << " " << (int)res << std::endl;
	blosum62_sum_local += res;
      }
    }
  } // END for (i=0; i<N; ++i)

  // Treat end gaps
  if (gap_s1 > 0)
  {
    blosum62_sum_local += (gap_open + gap_s1*gap_ext)*factor_front_back_gaps;
    gap_s1 = 0;
  }
  if (gap_s2 > 0)
  {
    blosum62_sum_local += (gap_open + gap_s2*gap_ext)*factor_front_back_gaps;
    gap_s2 = 0;
  }

  blosum62_sum = blosum62_sum_local;
  len           = len_local;
  unknown       = unknown_local;
}

// Fully implemented
inline float all_aa_distances(CSequences2 &seqs, std::vector<float> *dists)
{
  unsigned N      = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();
  
  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  unsigned  i,j;
  int       tmp;
  unsigned  l;
  unsigned  u;

  float sum=0;
  float d;

  if (dists != NULL)
  {
    dists->reserve( (N-1)*N/2 + 1 );

    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62(seqs.get_Seq_Data(i),
			     seqs.get_Seq_Data(j),
			     posNum,
			     tmp, l, u);
	if (l>0)
	{
	  d = (float) tmp/l;
	}
	else
	{
	  d = 0;
	}
	sum += d;

	//	std::cout << "(" << i << "," << j << "): " << d << std::endl; 
	//	std::cout << "Adding(1): " << d << std::endl;
	//	std::cout << seqs.get_Seq_Data(i) << std::endl;
	//	std::cout << seqs.get_Seq_Data(j) << std::endl;
	//	std::cout << "with tmp: " << tmp << " l: " << l << std::endl;
	//	std::cout << "New sum: " << sum << std::endl;
	//	std::cout << sum << std::endl;

	dists->push_back(d);
      }
    }
  }
  else
  {
    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62(seqs.get_Seq_Data(i),
			     seqs.get_Seq_Data(j),
			     posNum,
			     tmp, l, u);
	if (l>0)
	{
	  d = (float) tmp/l;
	}
	else
	{
	  d = 0;
	}
	//	std::cout << "(" << i << "," << j << "): " << d << std::endl; 
	//	std::cout << "Adding(2): " << (float)tmp/l << std::endl;
	//	std::cout << sum << std::endl;
	sum += d;
	//	std::cout << sum << std::endl;
      }
    }
  }
  return sum;
}

// Not verified ???
inline float all_aa_distances_gaps(CSequences2 &seqs, std::vector<float> *dists,
				   float gap_open, float gap_ext, float gap_match_score, float factor_front_back_gaps)
{
  unsigned N      = seqs.GetTaxaNum();
  unsigned posNum = seqs.GetPosNum();
  
  // CSequences2 is a class for aligned sequences, but we check whether all sequences
  // have equal length before we proceed.
  //  if (!seqs.equal_length_of_all_sequences() )
  //    return false;

  unsigned  i,j;
  float     tmp;
  unsigned  l;
  unsigned  u;

  float sum = 0;
  float d;

  if (dists != NULL)
  {
    dists->reserve( (N-1)*N/2 + 1 );

    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62_gaps(seqs.get_Seq_Data(i),
			      seqs.get_Seq_Data(j),
			      posNum,
			      gap_open, gap_ext, gap_match_score, factor_front_back_gaps,
			      tmp, l, u);
	if (l>0)
	{
	  d = (float) tmp/l;
	}
	else
	{
	  d = 0;
	}

	//	std::cout << "Adding: " << d << std::endl;
	sum += d;
	dists->push_back(d);
      }
    }
  }
  else // if (dists != NULL)
  {
    for (i=0; i<N; ++i)
    {
      for (j=0; j<i; ++j)
      {
	aaScore_BLOSUM62_gaps(seqs.get_Seq_Data(i),
			      seqs.get_Seq_Data(j),
			      posNum,
			      gap_open, gap_ext, gap_match_score, factor_front_back_gaps,
			      tmp, l, u);

	if (l>0)
	{
	  d = (float) tmp/l;
	}
	else
	{
	  d = 0;
	}
	//	std::cout << "Adding: " << d << std::endl;
	sum += d;
      }
    }
  }
  return sum;
}

// Count differences between two nucleotide (DNA) sequences:
// AA: identity is increased by 1. len is increased by 1.*
// AC: len is increased by 1.*
// RA: ambig is increased by 1. *
// RR: ambis is increased by 1.
// -A: gap is increased by 1.
// -R: gap is increased by 1. ambig is not increased by 1!
// --: gap is increased by 1.*
// OA: unknown is increased by 1.
// OR: unknown is increased. Ambig is not increases.
// O-: unknown is increased by 1. gap in not increased by one.
//
inline void nucStrictIdentity(const char *s1,
		       const char *s2,
		       unsigned N, 
		       int      &identity, 
		       unsigned &len,
		       unsigned &ambig,
		       unsigned &gap,
		       unsigned &unknown)
{
  unsigned i;

  int      identity_sum_local = 0;
  unsigned len_dna_local      = 0;
  unsigned ambig_local        = 0;
  unsigned gap_local          = 0;
  unsigned gap_or_ambig_local = 0;
  unsigned unknown_local      = 0;

  char c1, c2;
  int  res;

  for (i=0; i<N; ++i)
  {
    c1 = toupper(s1[i]);
    c2 = toupper(s2[i]);

    if ( is_DNA_base(c1) && is_DNA_base(c2) ) // This is the most probale case.
    {
      ++len_dna_local;
      if (c1 == c2)
	++identity_sum_local;
    }
    else if (is_DNA_or_DNA_ambig_or_GAP(c1) && is_DNA_or_DNA_ambig_or_GAP(c2)) // This is the next most probable case.
    {
      // c1 and c2 are valid symbols. They are not both DNA symbols, so one is a gap or an ambig symbol.
      ++gap_or_ambig_local;
      if (c1 == '-' || c2 == '-')
      {
	++gap_local;
	--len_dna_local;  // VERIFY
      }
      else
      {
	++ambig_local;
      }
    }
    else // c1 or c2 is unknown
    {
      ++unknown_local;
    }
  }
  identity      = identity_sum_local;
  len           = len_dna_local;
  unknown       = unknown_local;
}

// Code taken from distdna.c 

enum nt_changes{
    AA = 0,
    AC,
    AG,
    AT,
    CA,
    CC,
    CG,
    CT,
    GA,
    GC,
    GG,
    GT,
    TA,
    TC,
    TG,
    TT
};


/**
 * pairwise frequencies
 *
 * compute all pairwise nucleotide combination frequencies among two sequences.
 *
 * @param char *s1	The first sequence
 * @param char *s2	The second sequence
 * @param int slen	Sequences length
 * @param double f[16]	An array to store all 4x4 combination possibilities
 *
 */
inline void pairwise_frequencies(char *s1, char *s2, int slen, double *f)
{
    long int i, validpos;

    /* reset frequency counts */
    for (i = 0; i < 16; i++) f[i] = 0.;

    /* count pairwise frequencies and valid positions */
    validpos = slen;
    for (i = 0; i < slen; i++) {
        switch (s1[i]) {
	case 'A':
	    switch(s2[i]) {
	        case 'A': f[AA] += 1.0; break;
		case 'G': f[AG] += 1.0; break;
		case 'C': f[AC] += 1.0; break;
		case 'T': f[AT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'C':
	    switch(s2[i]) {
	        case 'A': f[CA] += 1.0; break;
		case 'G': f[CG] += 1.0; break;
		case 'C': f[CC] += 1.0; break;
		case 'T': f[CT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'G':
	    switch(s2[i]) {
	        case 'A': f[GA] += 1.0; break;
		case 'G': f[GG] += 1.0; break;
		case 'C': f[GC] += 1.0; break;
		case 'T': f[GT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	case 'T':
	    switch(s2[i]) {
	        case 'A': f[TA] += 1.0; break;
		case 'G': f[TG] += 1.0; break;
		case 'C': f[TC] += 1.0; break;
		case 'T': f[TT] += 1.0; break;
		default: validpos--; break;
	    }
	    break;
	default: validpos--; break;
	}
    }
    
    /* convert counts to frequencies */
    if (validpos != 0) {
        for (i = 0; i < 16; i++) f[i] /= (double) validpos;
    }
    /* else all f[i] == 0 */
}

/* double distance_p(char *s1, char *s2, char *s1_end) */
/* { */
/*   while (s1 < s1_end) */
/*   { */
/*     ++s1; */
/*     ++s2; */


/*   } */
/* } */

inline double dist_p_distdna(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int    gaps, k;
    double dist;

    /* compute unambiguous sum of matches and gaps */
    gaps = 0;
    matches = 0.;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero: return maximum distance */
    if (slen == 0) return FLT_MAX;

    /* compute uncorrected distance */
    dist = 1.0 - (matches / ((double) (slen - gaps) + ((double)gaps * gapwt)));

    return dist;
}

inline double dist_p_distPROTEIN(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int    gaps, k;
    double dist;

    /* compute unambiguous sum of matches and gaps */
    gaps = 0;
    matches = 0.;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero: return maximum distance */
    if (slen == 0) return FLT_MAX;

    /* compute uncorrected distance */
    dist = 1.0 - (matches / ((double) (slen - gaps) + ((double)gaps * gapwt)));

    return dist;
}


inline double dist_jukes_cantor_dnadist(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero */
    if (slen == 0) return FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / ((double) (slen - gaps) + ((double)gaps * gapwt)));

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     *
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     *
     */
    dist = -(3. / 4.) * log(1. - (dist / (3. / 4.)));

    return dist;
}

inline double dist_jukes_cantor_distPROTEIN(const char *s1, const char *s2, int slen, double gapwt=0, int ambig=0)
{
    double matches;
    int gaps, k;
    double dist;

    /* compute unambiguoug sum of matches and gaps */
    matches = gaps = 0;
    for (k = 0; k < slen; k++) {
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    gaps++;
	//	else if (ambig)
	//	    matches += nt_ambig_match(s1[k], s2[k]);
	else if (s1[k] == s2[k])
	    matches += 1.0;
    }

    /* avoid divide by zero */
    if (slen == 0) return FLT_MAX;

    /* compute uncorrected distance */
    /*
     * p = proportion of different nt
     * n = number of nt examined
     * k = number of nt pairs that differ
     * p ~= k/n = 1 - (matches/npos)
     */
    dist = 1.0 - (matches / ((double) (slen - gaps) + ((double)gaps * gapwt)));

    /* compute Jukes-Cantor distance */
    /*
     * d = -b log_e[1 - (p/b)]
     *
     * assuming equal frequencies for all nt, f(nt) = 1/4 and b=3/4
     *
     */
    dist = -(19. / 20.) * log(1. - (dist / (19. / 20.)));

    return dist;
}


inline double dist_tamura_dnadist(const char *s1, const char *s2, int slen)
{
    int k;
    int matches, transitions, transversions, npos;
    int gc1, at1, gc2, at2;
    double dist, P, Q, theta1, theta2, C;
    const char *nt = "AGCTU";
    const char *pur = "AGR";
    const char *pyr = "CTUY";

    /* compute transitions, transversions and GC content */
    matches = transitions = transversions = npos = 0;
    gc1 = at1 = gc2 = at2 = 0;
    for (k = 0; k < slen; k++) {
        /* this check is for compatibility with DISTMAT but 
	   shouldn't be done as it forces ignoring valid
	   nucleotides without a corresponding counterpart
	   in the other sequence *** JR *** */
        if (strchr("ATUGC", s1[k]) && strchr("ATUGC", s2[k])) { 

	/* check GC1 content */
	if (strchr("GCS", s1[k]))
	    gc1++;
	else if (strchr("ATUW", s1[k]))
	    at1++;
	/* check GC2 content */
	if (strchr("GCS", s2[k]))
	    gc2++;
	else if (strchr("ATUW", s2[k]))
	    at2++;

	}/* anything else cannot be identified as G+C and is ignored */

	/* gaps are ignored */
	if ((is_gap(s1[k])) || (is_gap(s2[k])))
	    continue;
	else if (strchr(nt, s1[k]) == strchr(nt, s2[k]))
	    /* both are the same unambiguous nucleotide */
	    matches++, npos++;
	else {
	    /* count transitions and transversions */
	    if (strchr(pur, s1[k]) && strchr(pur, s2[k]))
		/* R -> R */
		transitions++, npos++;
	    else if (strchr(pyr, s1[k]) && strchr(pyr, s2[k]))
		/* Y -> Y */
		transitions++, npos++;
	    else
		/* R -> Y */
		transversions++, npos++;
	}
    }
    
    /*
     * avoid divide by zero: if no matches return max. distance;
     * plus,if npos != o ==> there is at least one non-ambiguous nt
     * on each sequence and therefore (gc?+ac?) cannot be zero
     */
    if (npos == 0) return FLT_MAX;
    
    /*
     * P=transitions/npos;
     * Q=transversions/npos;
     * theta1 = GC fraction in sequence 1
     * theta2 = GC fraction in sequence 2
     * C = theta1 + theta2 - 2*theta1*theta2
     * 
     * distance = -C ln(1-P/C-Q) - 0.5(1-C) ln(1-2Q)
     */
    P = (double) transitions / (double) npos;
    Q = (double) transversions / (double) npos;
    theta1 = (double) gc1 / (double) (gc1 + at1);
    theta2 = (double) gc2 / (double) (gc2 + at2);
    C = theta1 + theta2 - (2 * theta1 * theta2);

    if ((2. * Q) > 1.) return -1.;	/* cannot use Tamura's distance */

    /* compute Tamura distance */
    dist = -C * log(1 - P / C - Q) - (0.5 * (1 - C) * log(1 - (2 * Q)));

    return dist;
}



class CDists_collection
{
 private:
  std::vector<float> dists;   // Entries in this vector have the order (1,0),(2,0),(2,1),(3,0),(3,1),(3,2),(4,0),(4,1),....

  char               mode;    // 0: no gaps in score, 1: gaps in score.

  unsigned           numTaxa;

  float              gap_open;
  float              gap_ext;
  float              gap_match_score;
  float              factor_front_back_gaps;

  static int vector_index_of_pair(int i, int j)
  {
    // We want to compute pairs (1,0),(2,0),(2,1),(3,0),(3,1),(3,2),(4,0),(4,1),....
    // to index values: 0,1,2,3,4,5,......
    // The index can be computed using i*(i-1)/2+j if i>j
    // or the scores are symetric:     j*(j-1)/2+i if j>i

    if (i==j)
      return -1;

    if (i>j)
      return i*(i-1)/2+j;
    else
      return j*(j-1)/2+i;
  }


  public:
 CDists_collection(CSequences2 &seqs, float pgap_open, float pgap_ext, float pgap_match_score, float pfactor_front_back_gaps, int pmode):
  mode(pmode), numTaxa(seqs.GetTaxaNum()), gap_open(pgap_open), gap_ext(pgap_ext), gap_match_score(pgap_match_score), factor_front_back_gaps(pfactor_front_back_gaps)
  {
    if (mode == 0)
    {
      all_aa_distances(seqs, &dists);
    }
    else if (mode == 1)
    {
      all_aa_distances_gaps(seqs, &dists, gap_open, gap_ext, gap_match_score, factor_front_back_gaps);
    }
    else
    {
      dists.clear();
    }
  }

  void print_debug(std::ostream &os)
  {
    os.setf(std::ios::fixed);
    os.precision(2);

    os << "Satus of CDists_collection object:" << std::endl;
    os << "Gap mode: " << (int)mode << std::endl;
    os << "Gap mode: " << (mode == 0 ? "without gaps" : (mode == 1 ? "with gaps" : "Unknown value")) << std::endl;
    os << "gap_open: " << gap_open  << std::endl;
    os << "gap_ext:  " << gap_ext << std::endl;
    os << "gap_match_score:        " << gap_match_score << std::endl;
    os << "factor_front_back_gaps: " << factor_front_back_gaps << std::endl;
    os << "numTaxa:                " << numTaxa  << std::endl;
    os << "Number of dists in vector: " << dists.size() << std::endl;
    os << "Dists: ";

    for (unsigned i=0; i<dists.size(); ++i)
      os << dists[i] << ", ";
    os << std::endl;

    os << "Dists as lower matrix: " << std::endl;

    for (unsigned i=0, k=0; i<numTaxa; ++i)
    {
      for (unsigned j=0; j<i; ++j)
      {
	os << dists[k] << " ";
	++k;
      }
      os << std::endl;
    }
    os << std::endl;
    os << std::endl;

  }

  bool successfully_populated()
  {
    return dists.size()>0;
  }

  float get_dist(int i, int j)
  {
    if (i < j)
      return dists[vector_index_of_pair(i,j)];
    else if (i > j)
      return dists[vector_index_of_pair(j,i)];
    else
      return 0; // Not a good value to indicate an "error".
  }

  void get_min_max(float &min_score, float &max_score)
  {
    vec_min_max(dists, min_score, max_score);
  }

  void get_mean_sd(double &mean_score, double &sd_score)
  {
    vec_mean_sd(dists, mean_score, sd_score);
  }

  void get_mean_sd(double &mean_score, double &sd_score, float &sum, float &sum_squares)
  {
    vec_mean_sd(dists, mean_score, sd_score, sum, sum_squares);
  }



  //######################
  // These element functions are too specific. The user should do this.
  void get_median_quartile(double &Q1_score, double &median_score, double &Q3_score)
  {
    vec_median_quartile_method3(dists, Q1_score, median_score, Q3_score);
  }

  void get_median_quartile_outlier_bounds(double &Q1_score, double &median_score, double &Q3_score, double &O_lower, double &O_upper)
  {
    vec_median_quartile_outlier_bounds_method3(dists, Q1_score, median_score, Q3_score, O_lower, O_upper);
  }
  //######################





  void get_dists(std::vector<float> &v)
  {
    v = dists;
  }

  void get_dists(unsigned i, std::vector<float> &v)
  {
    v.clear();
    
    unsigned j;
    for (j=0; j<numTaxa; ++j)
    {
      if (i==j)
      {
	
      }
      else
      {
	v.push_back(get_dist(i,j));
      }
    }
  }

  float get_mean_dist_to_others(unsigned i)
  {
    unsigned j;
    float sum=0;
    float d;

    //    std::cout << "Averaging: ";

    if (numTaxa < 2)
      return 0;

    for (j=0; j<numTaxa; ++j)
    {
      if (i==j)
      {
	//+0	
      }
      else
      {
	d = get_dist(i,j);
	//	std::cout << "(" << i << "," << j << "):" << d << ",";
	sum += d;
      }
    }
    //    std::cout << ": " << sum/(numTaxa-1) << std::endl; 
    return sum/(numTaxa-1);
  }

  void get_vector_of_mean_dists_to_others(std::vector<float> &v)
  {
    unsigned j;

    v.clear();
    for (j=0; j<numTaxa; ++j)
    {
      v.push_back(get_mean_dist_to_others(j));
    }
  }


};




#endif

