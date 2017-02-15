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

////////////////////////////////////////////////////////////////////////////////
// define/initialize static members
char *CmdLineOpts::genoFile = NULL;
char *CmdLineOpts::indFile = NULL;
char *CmdLineOpts::markerFile = NULL;
int   CmdLineOpts::vcfInput = 0; // TODO
char *CmdLineOpts::outFile = NULL;
int   CmdLineOpts::useImpute2Format = 0;
int   CmdLineOpts::vcfOutput = 0; // TODO
int   CmdLineOpts::noFamilyId = 0;
char *CmdLineOpts::XchrName = NULL; // TODO
char *CmdLineOpts::onlyChr = NULL; // TODO
int   CmdLineOpts::startPos = 0;
int   CmdLineOpts::endPos = INT_MAX;
int   CmdLineOpts::forceWrite = 0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    START_POS = CHAR_MAX + 1,
    END_POS,
  };

  static struct option const longopts[] =
  {
    {"geno", required_argument, NULL, 'g'},
    {"snp",  required_argument, NULL, 's'},
    {"ind",  required_argument, NULL, 'i'},
//    {"base", required_argument, NULL, 'b'},
    {"plink", required_argument, NULL, 'p'},
//    {"vcf", required_argument, NULL, 'v'},
    {"out",  required_argument, NULL, 'o'},
    {"chr", required_argument, NULL, 'c'},
    {"start", required_argument, NULL, START_POS},
    {"end", required_argument, NULL, END_POS},
    {"impute2", no_argument, &CmdLineOpts::useImpute2Format, 1},
//    {"vcf_out", no_argument, &CmdLineOpts::vcfOutput, 1},
    {"no_family_id", no_argument, &CmdLineOpts::noFamilyId, 1},
    {"force", no_argument, &CmdLineOpts::forceWrite, 1},
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
	outFile = optarg;
	break;
      case 'c':
	onlyChr = optarg;
	break;
      case START_POS:
	startPos = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse starting position option\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case END_POS:
	endPosSet = true;
	endPos = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse starting position option\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
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

  // For VCF format input, only genoFile and outFile should be set; others NULL
  if (vcfInput) {
    if (genoFile == NULL || outFile == NULL || indFile != NULL ||
							  markerFile != NULL) {
      if (haveGoodArgs)
	fprintf(stderr, "\n");
      if (outFile == NULL)
	fprintf(stderr, "ERROR: require output filenames\n");
      if (indFile != NULL || markerFile != NULL)
	fprintf(stderr, "ERROR: for VCF input, ind and snp should not be defined\n");
      haveGoodArgs = false;
    }
  }
  else if (genoFile == NULL || indFile == NULL || markerFile == NULL ||
							     outFile == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: required genotype, ind, snp, and output filenames\n");
    haveGoodArgs = false;
  }

  if (vcfOutput && useImpute2Format) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: --vcf_out and --impute2 contradictory: choose one only\n");
    haveGoodArgs = false;
  }
  else if (vcfInput) {
    if (!useImpute2Format)
      vcfOutput = 1; // default to VCF output unless explicitly set to IMPUTE2
  }

  if (startPos && !onlyChr) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "WARNING: using starting position, without specified chromosome\n");
  }
  if (endPosSet && !startPos) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr,
	    "WARNING: using ending position, with no starting position\n");
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
  fprintf(out, "  -o, --out <prefix>\toutput file prefix (suffix .phgeno,.phsnp,.phind,.log)\n");
  fprintf(out, " AND EITHER:\n");
  fprintf(out, "  -g, --geno <filename>\tgenotype file in PLINK BED format\n");
  fprintf(out, "  -s, --snp <filename>\tPLINK BIM file\n");
  fprintf(out, "  -i, --ind <filename>\tPLINK FAM file\n");
  // Note: not currently supporting EIGENSTRAT or packed ancestry map format:
//  fprintf(out, "  -g, --geno <filename>\tgenotype file in Eigenstrat or Ancestrymap format or\n");
//  fprintf(out, "\t\t\tPLINK BED file\n");
//  fprintf(out, "  -s, --snp <filename>\tSNP file or PLINK BIM file\n");
//  fprintf(out, "  -i, --ind <filename>\tindividual file or PLINK FAM file\n");
//  fprintf(out, " OR:\n");
//  fprintf(out, "  -b, --base <prefix>\tloads <prefix>.geno, <prefix>.snp, <prefix>.ind\n");
  fprintf(out, " OR:\n");
  fprintf(out, "  -p, --plink <prefix>\tloads <prefix>.bed, <prefix>.bim, <prefix>.fam\n");
  // Note: not currently supporting VCF format:
//  fprintf(out, " OR:\n");
//  fprintf(out, "  -v, --vcf <filename>\tbgzipped and indexed VCF file; - for stdin\n");
  fprintf(out, "\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  -c, --chr <#>\t\tonly analyze specified chromosome number\n");
  fprintf(out, "\n");
  fprintf(out, "  --start <#>\t\tstart position on given chromosome\n");
  fprintf(out, "  --end <#>\t\tend position on given chromosome\n");
  fprintf(out, "\n");
  fprintf(out, "  --impute2\t\tprint phase results in IMPUTE2 format\n");
  // Note: not currently supporting VCF output
//  fprintf(out, "  --vcf_out\t\toutput phase in bgzip VCF format (default for -v input)\n");
  fprintf(out, "  --force\t\tforce writing to output files (overwrite if they exist)\n");
  fprintf(out, "\n");
  fprintf(out, "  --no_family_id\tignore family ids from PLINK .fam file --\n");
  fprintf(out, "\t\t\tdefault PLINK ids are of the form \"family_id:person_id\"\n");
  fprintf(out, "\n");
}
