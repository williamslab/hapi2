// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>

#ifndef CMDLINEOPTS_H
#define CMDLINEOPTS_H

#define VERSION_NUMBER	"1.87.8b"
#define RELEASE_DATE    "17 Oct 2018"

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

    // Introduce an error state at one marker if doing so saves at least this
    // many recombinations. 0 disables.
    static int max1MarkerRecomb;

    // Verbose log?
    static int verbose;

    // How many parents (0-2) in a nuclear family must have data in order to
    // phase the family
    static int minNumParentsData;

    // How many children in a nuclear family must have data in order to phase
    // the family
    static int minNumChildrenData;
};

#endif // CMDOPTIONS_H
