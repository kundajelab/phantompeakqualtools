## Summary

This package computes informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data. It can also be used to obtain robust estimates of the predominant fragment length or characteristic tag shift values in these assays.

## Introduction

This set of programs operate on mapped Illumina single-end read datasets in tagAlign or BAM format.
They can be used to 

1. Compute the predominant insert-size (fragment length) based on strand cross-correlation peak
2. Compute Data quality measures based on relative phantom peak
3. Call Peaks and regions for punctate binding datasets

## Citations

If you are using the code or results in any formal publication please cite

[1] Landt SG1, Marinov GK, Kundaje A et al. ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res. 2012 Sep;22(9):1813-31. doi: 10.1101/gr.136184.111.

[2] Kharchenko PK, Tolstorukov MY, Park PJ, Design and analysis of ChIP-seq experiments for DNA-binding proteins Nat Biotechnol. 2008 Dec;26(12):1351-9

## Dependencies

NOTE: The current package does not run on a MacOS or Windows.
* unix, bash, R (>=3.1), awk, samtools, boost C++ libraries
* R packages: spp (>=1.14), caTools, snow, snowfall, bitops, Rsamtools (in Bioconductor)

## Files

1. `spp_1.14.tar.gz` : modified SPP peak-caller package (The original SPP-peak caller package was written by Peter Kharchenko [2], https://github.com/hms-dbmi/spp)
2. `run_spp.R` : The script to compute the frag length, data quality characteristics based on cross-correlation analysis and/or peak calling

## Installation

1. First make sure that you have installed R (version 3.1 or higher)

2. Also, you must have the Boost C++ libraries installed. Most linux distributions have these preinstalled.
If not, you can easily get these from your standard package manager for your linux distribution.
e.g synaptic package manager (`apt-get`) for ubuntu or emerge for gentoo.

3. Clone the repo.
   ```
   $ git clone https://github.com/kundajelab/phantompeakqualtools
   ```
4. Install the following R packages
   * snow (if you want parallel processing)
   * snowfall
   * bitops
   * caTools
   * Rsamtools (in the Bioconductor package)
   ```
   $ cd phantompeakqualtools
   $ R
   (From within R)
   > install.packages("snow", repos="http://cran.us.r-project.org")
   > install.packages("snowfall", repos="http://cran.us.r-project.org")
   > install.packages("bitops", repos="http://cran.us.r-project.org")
   > install.packages("caTools", repos="http://cran.us.r-project.org")
   > source("http://bioconductor.org/biocLite.R")
   > biocLite("Rsamtools",suppressUpdates=TRUE)
   > install.packages("./spp_1.14.tar.gz")
   ```

5. If your alignment files are BAM, you must have the samtools executable in your path so that the R script `run_spp.R` can call it using the `system()` command
You can get samtools from (here)[http://samtools.sourceforge.net/]
You can add the following line to your `~/.bashrc` file

   ```
   export PATH="<path_to_samtools_executable>:${PATH}"
   ```

6. Run `run_spp.R`

   ```
   Rscript run_spp.R <options>
   ```
   
## General usage

Usage: `Rscript run_spp.R <options>`

* Mandatory arguments

   | argument                | description                                 |
   |-------------------------|---------------------------------------------|
   |-c=\<ChIP_alignFile\>    | full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz) |
   
* Mandatory arguments for peak calling

   | argument                | description                                 |
   |-------------------------|---------------------------------------------|
   |-i=\<Input_alignFile\>   | full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz) |

* Optional arguments

   | argument                  | description                                     |
   |---------------------------|-------------------------------------------------|
   |-s=\<min\>:\<step\>:\<max\>| strand shifts at which cross-correlation is evaluated, default=-500:5:1500 |
   |-speak=\<strPeak\>         | user-defined cross-correlation peak strandshift |
   |-x=\<min\>:\<max\>         | strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10) |
   |-p=\<nodes\>               | number of parallel processing nodes, default=0  |
   |-fdr=\<falseDisoveryRate\> | false discovery rate threshold for peak calling |
   |-npeak=\<numPeaks\>        | threshold on number of peaks to call            |
   |-tmpdir=\<tempdir\>        | temporary directory (if not specified R function `tempdir()` is used) |
   |-filtchr=\<chrnamePattern\>| pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name |

* Output arguments

   | argument                   | description                                    |
   |----------------------------|------------------------------------------------|
   |-odir=\<outputDirectory\>   | name of output directory (If not set same as ChIP file directory is used) |
   |-savn=\<narrowpeakfilename\>| NarrowPeak file name (fixed width peaks)       |
   |-savn                       |                                                |
   |-savr=\<regionpeakfilename\>| RegionPeak file name (variable width peaks with regions of enrichment around peak summits) |
   |-savr                       |                                                |
   |-savd=\<rdatafile\>         | save Rdata file                                |
   |-savd                       |                                                |
   |-savp=\<plotdatafile\>      | save cross-correlation plot                    |
   |-savp                       |                                                |
   |-out=\<resultfile\>         | append peakshift/phantomPeak results to a file |
   |-rf                         | if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists |
   |-clean                      | if used it will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling `run_spp.R` is creating temporary files |

## Typical usage

1. Determine strand cross-correlation peak / predominant fragment length OR print out quality measures. `-out=<outFile>` will create and/or append to a file named <outFile> several important characteristics of the dataset.
   ```
   Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
   ```

   The file contains 11 tab delimited columns.

   |col.| abbreviation    | description                                                                                          |
   |----|-----------------|------------------------------------------------------------------------------------------------------|
   |1   | Filename        | tagAlign/BAM filename                                                                                |
   |2   | numReads        | effective sequencing depth i.e. total number of mapped reads in input file                           |
   |3   | estFragLen      | comma separated strand cross-correlation peak(s) in decreasing order of correlation.                 |
   |4   | corr_estFragLen | comma separated strand cross-correlation value(s) in decreasing order (COL2 follows the same order)  |
   |5   | phantomPeak     | Read length/phantom peak strand shift                                                                |
   |6   | corr_phantomPeak| Correlation value at phantom peak                                                                    |
   |7   | argmin_corr     | strand shift at which cross-correlation is lowest                                                    |
   |8   | min_corr        | minimum value of cross-correlation                                                                   |
   |9   | NSC             | Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8                                  |
   |10  | RSC             | Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)                  |
   |11  | QualityTag      | Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)       |

   The top 3 local maxima locations that are within 90% of the maximum cross-correlation value are output.
   In almost all cases, the top (first) value in the list represents the predominant fragment length.
   If you want to keep only the top value simply run

   ```
   sed -r 's/,[^\t]+//g' <outFile> > <newOutFile>
   ```
   You can run the program on multiple datasets in parallel and append all the quality information to the same <outFile> for a summary analysis.

   NSC values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. 
   Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type OR it could be due to poor quality)

   RSC values range from 0 to larger positive values. 1 is the critical threshold.
   RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.

   Qtag is a thresholded version of RSC.

2. Peak calling
   ```
   Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -fdr=<fdr> -odir=<peak_call_output_dir> -savr -savp -savd -rf
   Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -npeak=<npeaks> -odir=<peak_call_output_dir> -savr -savp -savd -rf
   ```
3. For IDR analysis you want to call a large number of peaks (relaxed threshold) so that the IDR model has access to a sufficient noise component.

   ```
   Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -npeak=300000 -odir=<peak_call_output_dir> -savr -savp -rf -out=<resultFile>
   ```

## Notes

* It is EXTREMELY important to filter out multi-mapping reads from the BAM/tagAlign files. Large number of multimapping reads can severly affect the phantom peak coefficient and peak calling results.

* If a dataset seems to have high PCR bottlenecking, then you might want to actually clamp the number of unique mappping reads per position to 1 or upto 5. If not the phantom peak coefficient can be artificially good.

* For the IDR rescue strategy, one needs to pool reads from replicates and then shuffle and subsample the mapped reads to create two balanced pseudoReplicates. This is much easier to implement on tagAlign/BED read-mapping files using the unix 'shuf' command. So it is recommended to use the tagAlign format.

* In most cases, you can simply use the maximum reported strand correlation peak as the predominant fragment length.
However, it is useful to manually take a look at the cross-correlation plot to make sure the selected max peak is not an artifact.

* Also, if there are problems with library size-selection, a dataset's cross-correlation profile can have multiple strong cross-correlation peaks. This is currently not autodetected.

## Input file formats

1. **BAM format**: This is a binary alignment format specified in http://samtools.sourceforge.net/SAM-1.3.pdf
You MUST have samtools installed to use run_spp.R with BAM files

2. **TagAlign files**: This a text-based BED3+3 alignment format that is easier to manipulate. It contains 6 tab delimited columns.

   |col.| abbrv.     | type   | description                                                        |
   |----|------------|--------|--------------------------------------------------------------------|
   |  1 | chrom      | string | Name of the chromosome                                             |
   |  2 | chromStart | int    | The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0. |
   |  3 | chromEnd   | int    | The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. |
   |  4 | sequence   | string | Sequence of this read                                              |
   |  5 | score      | int    | Indicates uniqueness or quality (preferably 1000/alignmentCount).  |
   |  6 | strand     | char   | Orientation of this read (+ or -)                                  |

NOTE: You dont have to store the sequence of reads in the sequence field as the peak caller never really uses that field. You can just put the letter 'N' in that field. This saves space significantly.

For the IDR rescue strategy, one needs to use shuffled and subsampled version of the alignment files. This is much easier to implement on tagAlign text files using the unix `shuf` command.
So it is recommended to preferably use the tagAlign format.

## Converting BAM to TAGALIGN files

It is very quick to convert pre-filtered BAM files (after removing unmapped, low quality reads, multimapping reads and duplicates) to gzipped tagAlign files using
```
samtools view -F 0x0204 -o - <bamFile> | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > <gzip_TagAlignFileName>
```

## Output file formats

1. **NarrowPeak/RegionPeak format**: The output peak file is in BED6+4 format known as tagAlign. It consists of 10 tab-delimited columns
   
   |col.| abbrv.      | type   | description                                                     |
   |----|-------------|--------|-----------------------------------------------------------------|
   |  1 | chrom       | string | Name of the chromosome                                          |
   |  2 | chromStart  | int    | The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0. |
   |  3 | chromEnd    | int    | The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. |
   |  4 | name        | string | Name given to a region (preferably unique). Use '.' if no name is assigned. |
   |  5 | score       | int    | Indicates how dark the peak will be displayed in the browser (1-1000). If '0', the DCC will assign this based on signal value. Ideally average signalValue per base spread between 100-1000. |
   |  6 | strand      | char   | +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned. |
   |  7 | signalValue | float  | Measurement of overall (usually, average) enrichment for the region. | 
   |  8 | pValue      | float  | Measurement of statistical signficance (-log10). Use -1 if no pValue is assigned. |
   |  9 | qValue      | float  | Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned. |
   | 10 | peak        | int    | Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called. |

## Releted pipelines and tools

* [AQUAS Peakcalling pipeline for TF and histone ChIP-seq data](https://github.com/kundajelab/TF_chipseq_pipeline)

## References

ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Landt SG`*`, Marinov GK`*`, Kundaje A`*`, Kheradpour P`*`, Pauli F, Batzoglou S, Bernstein BE, Bickel P, Brown JB, Cayting P, Chen Y, DeSalvo G, Epstein C, Fisher-Aylor KI, Euskirchen G, Gerstein M, Gertz J, Hartemink AJ, Hoffman MM, Iyer VR, Jung YL, Karmakar S, Kellis M, Kharchenko PV, Li Q, Liu T, Liu XS, Ma L, Milosavljevic A, Myers RM, Park PJ, Pazin MJ, Perry MD, Raha D, Reddy TE, Rozowsky J, Shoresh N, Sidow A, Slattery M, Stamatoyannopoulos JA, Tolstorukov MY, White KP, Xi S, Farnham PJ, Lieb JD, Wold BJ, Snyder M.
Genome Res. 2012 Sep;22(9):1813-31. doi: 10.1101/gr.136184.111.

## Contributors
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University (email: anshul AT kundaje DOT net)
* Youngsook Lucy Jung - Postdoctoral Fellow, Park Lab, Harvard Medical School
* Peter Kharchenko - Assistant Professor, Dept. of Biomedical Informatics, Harvard Medical School
* Peter Park - Associate Professor, Dept. of Biomedical Informatics, Harvard Medical School
