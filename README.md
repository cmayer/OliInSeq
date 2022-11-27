# OliInSeq: Detecting, masking and removing of outlier sequences in multiple sequence alignments

Table of contents:

- [About the OliInSeq program](#about-the-OliInSeq-package)
- [Compiling and installing OliInSeq](#compiling-and-installing)
  * [System requirements:](#system-requirements)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Frequently aksed questions](#Frequently-aksed-questions)

## About the OliInSeq program <a id="about-the-OliInSeq-package"></a>
Outlier detection, outlier masking and removing in multiple sequence alignments.

Outliers in multiple sequence alignments are a common source of error in all alignment bases sequence analyses problems.
Outliers add misleading signal in pyhlogeny or population genetic analyses. They can also result in a model misspecification in model based analyses.

OliInSeq stands for OutLIer Detection IN SEQuence alingnments. The idea for this project was the joint idea of Helena Vizan-Rico and Christoph Mayer. The OliInSeq program has been developed by Dr. Christoph Mayer and was first used in the analyses for the publication Vizan-Rico et al. 2019.

**When using the OliInSeq program, please cite:**
Will be updated soon.

## System requirements:  <a id="system-requirements"></a>
OliInSeq can be compiled on all platforms, which provide a C++ compiler.
In particular this includes Windows, MacOS and Linux operating systems.
Here I will only explain how to compile it on Mac and Linux computers.

## Compiling and installing OliInSeq <a id="compiling-and-installing"></a>
- Download the project or clone the project locally.
- On the command line go to the project folder and type "make".
- Make sure that you copy the OliInSeq-v0.9.6 program to a folder that listed in your $PATH variable, or copy it to the folder you want to use it from.

## Documentation <a id="documentation"></a>
You can always get online help on OliInSeq by typing:
OliInSeq-v0.x.y --help

The OliInSeq program has been designed to identify, mask and/or remove
   outlier sequences in molecular sequence alignments.Accepted input files
   are: Fasta files of aligned amino acid sequences. Sequences can either
   be analysed as a whole or using a sliding window. In most cases a
   sliding window is recommended, since it is more sensitive if sequences
   are outliers only in parts of the alignment. 

   The algorithm works as follows: Pairwise sequence scores (blosum62
   scores) are computed for all sequences in the input file. From this, a
   mean score is determined for each sequence to all other sequences. For
   the score the median and the quartiles (Q1, Q2=median, Q3) are
   determined. The distance IQR=Q3-Q1 is called the inter-quartile range
   (IQR). With a parameter, IQR can be set to a minimum value, which is
   necessary for very similar sequences, since for almost identical
   sequences one or a few substitutions are already identified as outliers.
   If this is the case, it is recommended to increase the minimum value of
   IQR. A sequence or partial sequence in the sliding window is identified
   as an outlier, if its blosum score is less than Q1-IQR*IQR-factor, where
   IQR-factor is a parameter with a default value of 1.5. If you want to
   mask more sequences use a smaller value, e.g. 1.0, if you want to mask
   fewer sequences use a larger value, e.g 2.0. The default value of 1.5
   corresponds to the general recommendations for the outlier detection in
   statistics. Preferred values depend on the data set. For moderately
   dissimilar sequences, larger values such as 2.0 often lead to better
   results. For very dissimilar sequences, outliers will only be identified
   if the value is reduced to 1.0. Sequences for which a specified
   proportion of the windows are outliers are removed entirely.
   Gap/ambig/lower case sites can be removed automatically. A nucleotide
   alignment corresponding to the amino acid alignment can be specified and
   residues will be masked/removed corresponding to the amino acid
   alignment.


This displays the following usage:
 ./OliInSeq-v0.9.6  -i <string> -o <string> [-e <floating point number>]
                      [--minNumValidSeqsInWindow <integer>] [-M] [-N]
                      [--verbosity <unsigned integer>] [-w <unsigned
                      integer>] [-g <floating point number>] [-I <floating
                      point number>] [-f <floating point number>] [-R]
                      [--corresponding-nuc-fasta-file-name <string>] [--]
                      [--version] [-h]

Where: 
   -i <string>,  --input-file <string>
     (required)  Name of input file in fasta format.

   -o <string>,  --output-file <string>
     (required)  Name of output file.

   -e <floating point number>, 
      --thresholdPropOutlierWindows2removeSequences <floating point
      number>
     A sequence can be completely removed from an MSA if it is identified
     as an outlier in at least the specified proportion of sequence
     windows. The default is 0.5. So a sequence which is identified as an
     outlier in 50% or more of the sequence windows, it is removed
     completely. A value >1 implies, that sequences are never removed. 

   --minNumValidSeqsInWindow <integer>
     A meaningful outlier test can normally only be carried out, if the
     number of sequences is large. This is not always the case. In order to
     be able to compute quartiles and a IQR, we decided to require 5 valid
     sequences in a sliding window to be able to detect outliers. If you
     think a different value for the minimum number of sequences is
     appropriate, it can be specified with this parameter. Default: 5. 

   -M,  --mask-ambig
     With this argument, outlier sequence segments of nucleotide sequences
     are masked with Ns and amino acid sequences with X. Default: mask with
     lower case symbols.

   -N,  --no-masking
     The default is to mask outliers with lower case symbols. With the -M
     option, residues are masked with the ambiguity code characters N
     (nucleotide) or X (amino acids). With the -N option, masking can be
     turned off completely.

   --verbosity <unsigned integer>
     The verbosity option controls the amount of information OliInSeq
     writes to the console while running. 0: Print only welcome message and
     essential error messages that lead to exiting the program. 

   -w <unsigned integer>,  --windowsSize <unsigned integer>
     The window size for the sliding window outlier detection. For a value
     of 0 no sliding window will be used, but the sequence will be analyses
     as a whole. Default: 30

   -g <floating point number>,  --maximum-gap-proportion <floating point
      number>
     If a sequence segment of window size contains only or mostly gaps, its
     score to all other sequences could be used to identify it as an
     outlier. To avoid problems und meaningless scores, we do not consider
     sequences, which have a gap proportion of a value larger than this
     parameter. The maximum-gap-proportion parameter must be in the range
     0..1. With a value of 0.5, if a sequence has more than 50% gaps in a
     window, the score of this sequence will not be considered in computing
     the inter quartile distance nor can this sequence be determined as an
     outlier in the particular window. This sequence segment is also not
     considered when computing the proportion of windows in which a
     sequence is an outlier. Default: 0.5

   -I <floating point number>,  --minimum-IQR <floating point number>
     In conserved or almost invariant regions, a sequence is identified as
     an outlier even if it differs by only one or two residues. To avoid
     recognizing such a sequence as an outlier, a minimum value for the IQR
     can be specified. If the true IQR is smaller than the minimum IQR, the
     minimum IRQ is used instead of the real IQR. Default: 0.563

   -f <floating point number>,  --IQR-factor <floating point number>
     For a given window, outliers are determined by computing
     Q1-IQT*IQR-factor. Scores lower than this are considered as outliers.
     Default: 1.5

   -R,  --remove-gap-ambig-lowerCase-sites
     With this argument, alignment sites are removed which contain only gap
     , ambig, and lower case sequence symbols. 

   --corresponding-nuc-fasta-file-name <string>
     Name of file corresponding to the specified amino acid file in which
     outliers are detected. The corresponding file has to correspond to the
     ....

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.



### Compiling OliInSeq

Please unpack the OliInSeq archive you downloaded 


## Quickstart <a id="quickstart"></a>



## Frequently asked questions <a id="Frequently-aksed-questions"></a>
No questions have been asked so far.
