// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>

#ifndef CMDLINEOPTS_H
#define CMDLINEOPTS_H

#define VERSION_NUMBER	"1.97"
#define RELEASE_DATE    "22 May 2024"

class CmdLineOpts {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static bool parseCmdLineOptions(int argc, char **argv);
    static void printUsage(FILE *out, char *programName);

    //////////////////////////////////////////////////////////////////
    // public static fields : variables set by command-line options
    //////////////////////////////////////////////////////////////////

    // Genotype filename
    static char *genoFile;

    // Individual filename
    static char *indFile;

    // SNP data filename
    static char *markerFile;

    // Is the input in VCF format?
    static int vcfInput;

    // Output filename prefix
    static char *outPrefix;

    // Print text format haplotypes?
    static int txtOutput;

    // Print ped format haplotypes?
    static int pedOutput;

    // Print vcf format haplotypes?
    static int vcfOutput;

    // Print JSON format haplotypes?
    static int jsonOutput;

    // Print JSON format parent haplotypes?
    static int jsonParOutput;

    // Print inheritance vector values to a csv file?
    static int ivOutput;

    // Detect crossovers? If non-zero, gives the number of informative markers
    // required to establish an event as a crossover (vs. a conversion or
    // erroneous apparent recombination)
    static int detectCO;

    // When <detectCO> is non-zero, how many markers at the beginning and
    // end of the chromosome are needed for a crossover to be called? Might
    // want to set this to less than <detectCO>.
    static int edgeCO;

    // Print in IMPUTE2 format?
    static int useImpute2Format;

    // When reading a PLINK format BED file, do not print the family ids
    // in the output file
    static int noFamilyId;

    // What is the name of the X chromosome?
    static char *XchrName;

    // Analyze only specified chromosome -- if non-zero
    static char *onlyChr;

    // Starting position for analysis of partial chromosome
    static int startPos;

    // Starting position for analysis of partial chromosome
    static int endPos;

    // Force write of output files?  By default the following is false, and
    // the program quits if any of the output files exist (except the log file)
    static int forceWrite;

    // Introduce an error state if doing so saves at least this many
    // recombinations. 0 disables.
    static int minErrRecombs;

    // Number of successive informative markers that can be treated as errors
    // (Really the number of markers back to look for a potential previous
    // state where the intervening informative markers will be flagged as
    // errors)
    static uint8_t errorLength;

    // For detecting when the children's IV values have switched from a parent
    // with data to a parent without data.
    static int bothParHetThreshold;

    // For detecting when a missing data parent transmitted the same haplotype
    // to all the children, meaning that there are long stretches with zero
    // informative markers for that parent. If at least <forceInformInit>
    // markers occur without any informative markers for the parent, HAPI will
    // force in an informative marker. This will lead to a state with an IV
    // where all the children have the same haplotype.
    static int forceInformInit;

    // Once a force inform interval has triggered, HAPI will introduce a forced
    // informative marker separated by this many markers so long as there are
    // fewer than <forceInformTolerance> markers detected as informative for
    // the one-hap-trans parent.
    static int forceInformSeparation;

    // See comment on previous field: allow this many informative markers
    // between every <forceInformSeparation> number of markers. If there are
    // more than this, HAPI will not introduce a forced informative marker, but
    // will not break the force inform state (see below).
    static int forceInformTolerance;

    // See comment in previous two fields: this many informative markers in
    // a <forceInformSeparation> number of markers stops the tracking of such
    // intervals. This is in place to allow for a burst of erroneous markers in
    // the midst of a large one hap trans region without breaking the interval.
    static int numInformToBreakForceInform;

    // Verbose log?
    static int verbose;

    // How many parents (0-2) in a nuclear family must have data in order to
    // phase the family
    static int minNumParentsData;

    // How many children in a nuclear family must have data in order to phase
    // the family
    static int minNumChildrenData;

    // Bit mask for forcing missing data on parents. Values:
    // 0: (default) no forced missingness
    // 1: dad missing
    // 2: mom missing
    // 3: both dad and mom missing
    static uint8_t forceMissingParBits;
};

#endif // CMDOPTIONS_H
