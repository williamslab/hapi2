// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <genetio/marker.h>
#include "cmdlineopts.h"

#define DEFAULT_NO_ERROR_MAX	2

////////////////////////////////////////////////////////////////////////////////
// define/initialize static members
char *CmdLineOpts::genoFile = NULL;
char *CmdLineOpts::indFile = NULL;
char *CmdLineOpts::markerFile = NULL;
int   CmdLineOpts::vcfInput = 0;
char *CmdLineOpts::outPrefix = NULL;
int   CmdLineOpts::txtOutput = 0;
int   CmdLineOpts::pedOutput = 0;
int   CmdLineOpts::vcfOutput = 0;
int   CmdLineOpts::ivOutput = 0;
int   CmdLineOpts::useImpute2Format = 0;
int   CmdLineOpts::noFamilyId = 0;
char *CmdLineOpts::XchrName = NULL;
char *CmdLineOpts::onlyChr = NULL;
int   CmdLineOpts::startPos = 0;
int   CmdLineOpts::endPos = INT_MAX;
int   CmdLineOpts::forceWrite = 0;
int   CmdLineOpts::max1MarkerRecomb = DEFAULT_NO_ERROR_MAX;
int   CmdLineOpts::oneHapTransThreshold = 100;
int   CmdLineOpts::detectCO = 0;
int   CmdLineOpts::edgeCO = 0;
int   CmdLineOpts::verbose = 0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    START_POS = CHAR_MAX + 1,
    END_POS,
    NO_ERROR_MAX,
    DETECT_CO,
    EDGE_CO,
  };

  static struct option const longopts[] =
  {
    {"geno", required_argument, NULL, 'g'},
    {"snp",  required_argument, NULL, 's'},
    {"ind",  required_argument, NULL, 'i'},
//    {"base", required_argument, NULL, 'b'},
    {"plink", required_argument, NULL, 'p'},
    // Probably won't ever do vcf input, but if I do have it be the -v option,
    // with --vcf being for output
//    {"vcf", required_argument, NULL, 'v'},
    {"out",  required_argument, NULL, 'o'},
    {"txt", no_argument, &CmdLineOpts::txtOutput, 1},
    {"ped", no_argument, &CmdLineOpts::pedOutput, 1},
    {"vcf", no_argument, &CmdLineOpts::vcfOutput, 1},
    {"iv", no_argument, &CmdLineOpts::ivOutput, 1},
    {"detect_co", required_argument, NULL, DETECT_CO},
    {"edge_co", required_argument, NULL, EDGE_CO},
    {"chr", required_argument, NULL, 'c'},
    {"start", required_argument, NULL, START_POS},
    {"end", required_argument, NULL, END_POS},
    {"no_err_max", required_argument, NULL, NO_ERROR_MAX},
//    {"impute2", no_argument, &CmdLineOpts::useImpute2Format, 1},
    {"no_family_id", no_argument, &CmdLineOpts::noFamilyId, 1},
    {"force", no_argument, &CmdLineOpts::forceWrite, 1},
    {"verbose", no_argument, &CmdLineOpts::verbose, 1},
    {0, 0, 0, 0}
  };

  // option index for getopt_long()
  int optionIndex = 0;
  int c;
  bool endPosSet = false; // did the user set the end position?

  bool haveGoodArgs = true;

  int prefixLength = 0;

  // Note: not currently supporting the -b and -v options corresponding to
  //    Eigenstrat/packed Ancestry Map format input and VCF input
  // Code to process these options still exists here, however (see cases below),
  // but the optstring being set without these options will prevent the code
  // from being run
//  char optstring[80] = "g:i:s:b:p:v:o:c:";
  char optstring[80] = "g:i:s:p:o:c:";
  while ((c = getopt_long(argc, argv, optstring, longopts, &optionIndex))
									!= -1) {
    errno = 0;    // initially: may get errors from strtol
    char *endptr; // for strtol
    switch (c) {
      case 0:
	// flag set by getopt_long()
	break;

      case 'g':
	if (genoFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of genotype filename\n");
	  haveGoodArgs = false;
	}
	genoFile = optarg;
	break;
      case 'i':
	if (indFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of individual filename\n");
	  haveGoodArgs = false;
	}
	indFile = optarg;
	break;
      case 's':
	if (markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of SNP filename\n");
	  haveGoodArgs = false;
	}
	markerFile = optarg;
	break;
      case 'b':
	if (genoFile != NULL || indFile != NULL || markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input files\n");
	  haveGoodArgs = false;
	}
	prefixLength = strlen(optarg);
	genoFile = new char[prefixLength + 5 + 1]; // + ".geno" + '\0'
	markerFile = new char[prefixLength + 4 + 1]; // + ".snp" + '\0'
	indFile = new char[prefixLength + 4 + 1]; // + ".ind" + '\0'
	sprintf(genoFile, "%s.geno", optarg);
	sprintf(markerFile, "%s.snp", optarg);
	sprintf(indFile, "%s.ind", optarg);
	break;
      case 'p':
	if (genoFile != NULL || indFile != NULL || markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input files\n");
	  haveGoodArgs = false;
	}
	prefixLength = strlen(optarg);
	genoFile = new char[prefixLength + 4 + 1]; // + ".bed" + '\0'
	markerFile = new char[prefixLength + 4 + 1]; // + ".bim" + '\0'
	indFile = new char[prefixLength + 4 + 1]; // + ".fam" + '\0'
	sprintf(genoFile, "%s.bed", optarg);
	sprintf(markerFile, "%s.bim", optarg);
	sprintf(indFile, "%s.fam", optarg);
	break;
      case 'v':
	if (genoFile != NULL || indFile != NULL || markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input files\n");
	  haveGoodArgs = false;
	}
	genoFile = optarg;
	vcfInput = 1;
	break;
      case 'o':
	outPrefix = optarg;
	break;
      case 'c':
	onlyChr = optarg;
	break;
      case START_POS:
	startPos = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse starting position option as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case END_POS:
	endPosSet = true;
	endPos = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse starting position option as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case NO_ERROR_MAX:
	max1MarkerRecomb = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse no_err_max option as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case DETECT_CO:
	detectCO = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse detect_co option as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	if (detectCO <= 0) {
	  fprintf(stderr, "ERROR: must provide strictly positive value for detect_co option\n");
	  haveGoodArgs = false;
	}
	break;
      case EDGE_CO:
	edgeCO = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse edge_co option as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	if (edgeCO <= 0) {
	  fprintf(stderr, "ERROR: must provide strictly positive value for edge_co option\n");
	  haveGoodArgs = false;
	}
	break;


      case '?':
	// bad option; getopt_long already printed error message
        printUsage(stderr, argv[0]);
	exit(1);
	break;

      default:
	exit(1);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Check for errors in command line options

  // For VCF format input, only genoFile and outPrefx should be set; others NULL
  if (vcfInput) {
    if (genoFile == NULL || outPrefix == NULL || indFile != NULL ||
							  markerFile != NULL) {
      if (haveGoodArgs)
	fprintf(stderr, "\n");
      if (outPrefix == NULL)
	fprintf(stderr, "ERROR: require output filenames\n");
      if (indFile != NULL || markerFile != NULL)
	fprintf(stderr, "ERROR: for VCF input, ind and snp should not be defined\n");
      haveGoodArgs = false;
    }
  }
  else if (genoFile == NULL || indFile == NULL || markerFile == NULL ||
							    outPrefix == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: required genotype, ind, snp, and output filenames\n");
    haveGoodArgs = false;
  }

  if (vcfInput) {
    if (!useImpute2Format)
      vcfOutput = 1; // default to VCF output unless explicitly set to IMPUTE2
  }

  if (startPos && !onlyChr) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "WARNING: using starting position without specified chromosome\n");
  }
  if (endPosSet && !startPos) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr,
	    "WARNING: using ending position, with no starting position\n");
  }

  if (detectCO && !edgeCO)
    // set edgeCO if detectCO enabled but edgeCO not specified
    edgeCO = detectCO;
  else if (edgeCO > detectCO) {
    if (!detectCO)
      fprintf(stderr, "ERROR: must enable --detect_co to use --edge_co\n");
    else
      fprintf(stderr, "ERROR: --edge_co value cannot be greater than --detect_co\n");
    haveGoodArgs = false;
  }

  if (!(txtOutput || pedOutput || vcfOutput || ivOutput || detectCO)) {
    fprintf(stderr, "ERROR: must choose at least one type of results to print\n");
    haveGoodArgs = false;
  }

  // TODO: allow user to specify X chromosome name; if unspecified, do:
  XchrName = new char[2];
  strcpy(XchrName, "X");

  if (!haveGoodArgs) {
    printUsage(stderr, argv[0]);
  }

  return haveGoodArgs;
}

// Prints usage message to <out>.  <programName> should be argv[0]
void CmdLineOpts::printUsage(FILE *out, char *programName) {
  fprintf(out, "\n");
  fprintf(out, "HAPI v%s!    (Released %s)\n\n", VERSION_NUMBER,
	  RELEASE_DATE);
  fprintf(out, "Usage:\n");
  fprintf(out, "%s [ARGUMENTS]\n", programName);
  fprintf(out, "\n");
  fprintf(out, "REQUIRED ARGUMENTS:\n");
  fprintf(out, "1. INPUT FILES - EITHER:\n");
  fprintf(out, "  -p, --plink <prefix>\tLoads <prefix>.bed, <prefix>.bim, <prefix>.fam\n");
  fprintf(out, " OR:\n");
  fprintf(out, "  -g, --geno <filename>\tGenotype file in PLINK BED format\n");
  fprintf(out, "  -s, --snp <filename>\tPLINK BIM file\n");
  fprintf(out, "  -i, --ind <filename>\tPLINK FAM file\n");
  // Note: not currently supporting EIGENSTRAT or packed ancestry map format:
//  fprintf(out, "  -g, --geno <filename>\tgenotype file in Eigenstrat or Ancestrymap format or\n");
//  fprintf(out, "\t\t\tPLINK BED file\n");
//  fprintf(out, "  -s, --snp <filename>\tSNP file or PLINK BIM file\n");
//  fprintf(out, "  -i, --ind <filename>\tindividual file or PLINK FAM file\n");
//  fprintf(out, " OR:\n");
//  fprintf(out, "  -b, --base <prefix>\tloads <prefix>.geno, <prefix>.snp, <prefix>.ind\n");
  fprintf(out, "\n");
  fprintf(out, "2. OUTPUT DIRECTORY:\n");
  fprintf(out, "  -o, --out <dir>\tOutput directory (creates directory <dir> and <dir>.log)\n");
  fprintf(out, "\n");
  fprintf(out, "3. RESULTS TO PRINT - ONE OR MORE OF:\n");
  fprintf(out, "  --vcf\t\t\tVCF format haplotypes\n");
  fprintf(out, "  --ped\t\t\tPLINK ped format haplotypes\n");
  fprintf(out, "  --txt\t\t\tText format haplotypes\n");
  fprintf(out, "  --iv\t\t\tInheritance vector data in CSV format\n");
  fprintf(out, "  --detect_co <#>\tDetect crossover events\n");
  fprintf(out, "\t\t\tNumeric argument specifies how many informative markers\n");
  fprintf(out, "\t\t\tmust have a switched haplotype for an event to be called\n");
//  fprintf(out, "\t\t\tValues >=5 are recommended. Smaller values may detect\n");
//  fprintf(out, "\t\t\tnon-crossovers or events driven by assembly errors\n");

  // Note: not currently supporting VCF format:
//  fprintf(out, " OR:\n");
//  fprintf(out, "  -v, --vcf <filename>\tbgzipped and indexed VCF file; - for stdin\n");
  fprintf(out, "\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  -c, --chr <string>\tOnly analyze specified chromosome\n");
  fprintf(out, "  --start <#>\t\tstart position on given chromosome\n");
  fprintf(out, "  --end <#>\t\tend position on given chromosome\n");
  fprintf(out, "\n");
  fprintf(out, "  --no_err_max <#>\tMaximum number of recombinations attributable to a\n");
  fprintf(out, "\t\t\tsingle marker before it is called an error. Default: %d.\n",
	  DEFAULT_NO_ERROR_MAX);
  fprintf(out, "\t\t\t0 disables; 1 calls single marker non-crossovers errors\n");
  fprintf(out, "  --edge_co <#>\t\tAllow detection of the \"background haplotype\" with fewer\n");
  fprintf(out, "\t\t\tthan --detect_co markers at start/end of chromosomes.\n");
  fprintf(out, "\t\t\tA child must have a sequence of this many informative\n");
  fprintf(out, "\t\t\tmarkers from one haplotype before calling crossovers\n");
  fprintf(out, "\n");
  // Note: not currently supporting VCF or IMPUTE2 output
//  fprintf(out, "  --impute2\t\tprint phase results in IMPUTE2 format\n");
//  fprintf(out, "  --vcf_out\t\toutput phase in bgzip VCF format (default for -v input)\n");
  fprintf(out, "  --force\t\tForce writing to output files (overwrite if they exist)\n");
  fprintf(out, "\n");
  fprintf(out, "  --no_family_id\tIgnore family ids from PLINK .fam file\n");
  fprintf(out, "\n");
  fprintf(out, "  --verbose\t\tPrint verbose messsages to log during phasing\n");
  fprintf(out, "\n");
}
