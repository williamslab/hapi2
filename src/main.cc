// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>
#include <genetio/marker.h>
#include <genetio/personbulk.h>
#include <genetio/personio.h>
#include <genetio/util.h>
#include "cmdlineopts.h"
#include "phaser.h"
#include "analysis.h"

void createOutDir();
bool openFilesToWrite(char *&filename, FILE *resultFiles[3], int chrIdx,
		      const char *parentIds[2], int famIdLen, int totalFileLen,
		      int &allocFilenameLen, FILE *log);
bool checkIfFileExists(char *filename, bool printWarning);

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  // Allocate a long string so we can have filenames with the parent ids, the
  // family ids, the chromosome name, and the file extension.
  // If needed, we'll increase this below
  int prefixLen = strlen(CmdLineOpts::outPrefix);
  int allocFilenameLen = prefixLen + 1024;
  char *filename = new char[ allocFilenameLen ];

  // open the .log file for writing
  sprintf(filename, "%s.log", CmdLineOpts::outPrefix);
  FILE *log = fopen(filename, "w");
  if (!log) {
    fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
    perror("open");
    exit(1);
  }

  // Get started: print status to stdout and the log
  FILE *outs[2] = { stdout, log };
  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "\n");
    fprintf(out, "HAPI v%s!    (Released %s)\n\n", VERSION_NUMBER,
	    RELEASE_DATE);

    if (CmdLineOpts::vcfInput) {
      fprintf(out, "VCF input file:\t\t%s\n", CmdLineOpts::genoFile);
    }
    else {
      fprintf(out, "Genotype file:\t\t%s\n", CmdLineOpts::genoFile);
      fprintf(out, "SNP file:\t\t%s\n", CmdLineOpts::markerFile);
      fprintf(out, "Individual file:\t%s\n", CmdLineOpts::indFile);
    }
    fprintf(out, "Output directory:\t%s\n", CmdLineOpts::outPrefix);

    fprintf(out, "\n");

    if (CmdLineOpts::onlyChr != NULL) {
      fprintf(out, "Chromosome:\t\t%s\n\n", CmdLineOpts::onlyChr);
    }

    fprintf(out, "Output to be generated in output directory:\n");
    if (CmdLineOpts::txtOutput) {
      fprintf(out, "  Text format haplotypes\n");
    }
    if (CmdLineOpts::ivOutput) {
      fprintf(out, "  Inheritance vectors for all markers\n");
    }
    if (CmdLineOpts::detectCO) {
      fprintf(out, "  Detected COs, %d informative markers required to call CO\n",
	      CmdLineOpts::detectCO);
    }
    if (CmdLineOpts::edgeCO != CmdLineOpts::detectCO) {
      fprintf(out, "    Background haplotype at chromosome start/end established with %d markers\n",
	      CmdLineOpts::edgeCO);
    }

    fprintf(out, "\n");
  }

  //////////////////////////////////////////////////////////////////////////
  // Read the genotype data
  PersonIO<PersonBulk>::readData(CmdLineOpts::genoFile,
				 CmdLineOpts::markerFile,
				 CmdLineOpts::indFile, CmdLineOpts::onlyChr,
				 CmdLineOpts::startPos, CmdLineOpts::endPos,
				 CmdLineOpts::XchrName,
				 CmdLineOpts::noFamilyId, log,
				 /*allowEmptyParents=*/ true,
				 /*bulkData=*/ true);

  createOutDir();

  //////////////////////////////////////////////////////////////////////////
  // Have output directory and log, can go forward!

  dynarray<NuclearFamily *> toBePhased;
  for(NuclearFamily::fam_ht_iter iter = NuclearFamily::familyIter();
			       iter != NuclearFamily::familyIterEnd(); iter++) {
    NuclearFamily *theFam = iter->second;
    if (theFam->numChildren() > 1) {
      int childrenWithData = 0;
      for(int c = 0; c < theFam->numChildren(); c++)
	if (theFam->_children[c]->hasData())
	  childrenWithData++;
      if (childrenWithData >= 2)
	toBePhased.append(theFam);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Phase!
  int numChrs = Marker::getNumChroms();
  int numFamsToBePhased = toBePhased.length();

  Phaser::init();
  int numFinished = 0;
  for(int f = 0; f < numFamsToBePhased; f++) {

    for(int o = 0; o < 2; o++) {
      FILE *out = outs[o];
      fprintf(out, "Phasing families with two or more children... %d / %d",
	      numFinished, numFamsToBePhased);
    }
    printf("\r");
    fflush(stdout);
    fprintf(log, "\n");


    NuclearFamily *theFam = toBePhased[f];
    const char *parentIds[2] = { theFam->_parents->first->getId(),
				 theFam->_parents->second->getId() };
    int famIdLen = theFam->_parents->first->getFamilyIdLength();
    // add 4 for prefix (e.g., /hap), 3 for dashes between ids, and 4 for ".ext"
    int famSpecificFileLen = strlen(parentIds[0]) + strlen(parentIds[1]) +
			     famIdLen + 4 + 3 + 4;

    // Phase the current family on each chromosome successively and print
    // the results the user requested
    theFam->initFam();
    for(int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
      FILE *resultsFiles[3];
      bool shouldPhase = openFilesToWrite(filename, resultsFiles, chrIdx,
			  parentIds, famIdLen,
			  /*totalFileLen=*/ prefixLen + famSpecificFileLen,
			  allocFilenameLen, log);

      if (!shouldPhase) {
	for(int o = 0; o < 2; o++) {
	  FILE *out = outs[o];
	  fprintf(out, "Unable to write any output for chr %s: will not phase it\n",
		  Marker::getChromName(chrIdx));
	}
	continue;
      }

      Phaser::run(theFam, chrIdx);

      if (CmdLineOpts::txtOutput && resultsFiles[0]) {
	theFam->printHapTxt(resultsFiles[0], chrIdx);
	fclose(resultsFiles[0]);
      }
      if (CmdLineOpts::ivOutput && resultsFiles[1]) {
	theFam->printIvCSV(resultsFiles[1], chrIdx);
	fclose(resultsFiles[1]);
      }
      if (CmdLineOpts::detectCO && resultsFiles[2]) {
	Analysis::findCOs(theFam, resultsFiles[2], chrIdx);
	fclose(resultsFiles[2]);
      }
    }

    numFinished++;
  }

  // Finished phasing all families!
  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "Phasing families with two or more children... done");
  }
  printf("               \n");
  fprintf(log, "\n");

  fclose(log);
}

// Creates directory that is the output prefix specified by the user
void createOutDir() {
  // make the directory to store the output in
  if (mkdir(CmdLineOpts::outPrefix, 0777) != 0) {
    if (errno == EEXIST) {
      struct stat buffer;
      if (stat( CmdLineOpts::outPrefix, &buffer) == 0) {
	if (!S_ISDIR(buffer.st_mode)) {
	  fprintf(stderr, "ERROR: %s exists and is not a directory\n",
		  CmdLineOpts::outPrefix);
	  exit(5);
	}
	// else: success -- directory already created
      }
      else {
	fprintf(stderr, "ERROR: unable to access %s\n", CmdLineOpts::outPrefix);
	perror("stat");
      }
    }
    else {
      fprintf(stderr, "ERROR: couldn't create directory %s\n",
	      CmdLineOpts::outPrefix);
      perror("mkdir");
      exit(1);
    }
  }
}

// Attempts to open the files to be printed to, and alerts the user if any
// exist. Returns true if at least one output file exists
bool openFilesToWrite(char *&filename, FILE *resultsFiles[3], int chrIdx,
		      const char *parentIds[2], int famIdLen, int totalFileLen,
		      int &allocFilenameLen, FILE *log) {
  const char *chrName = Marker::getChromName(chrIdx);
  int chrNameLen = strlen(chrName);

  if (totalFileLen + chrNameLen + 1 > allocFilenameLen) {
    delete [] filename;
    allocFilenameLen = (totalFileLen + chrNameLen + 1) * 2;
    filename = new char[ allocFilenameLen ];
  }

  bool haveAnOutput = false; // only if one of the outputs is OK to write
  const char *types[3] = { "hap", "iv", "co" };
  const char *ext[3] = { ".txt", ".csv", "" };
  const int  wantType[3] = { CmdLineOpts::txtOutput, CmdLineOpts::ivOutput,
			     CmdLineOpts::detectCO };
  for(int o = 0; o < 3; o++) {
    if (wantType[o]) {
      if (famIdLen == 0)
	sprintf(filename, "%s/%s-%s-%s.%s%s", CmdLineOpts::outPrefix,
		types[o], parentIds[0], parentIds[1], chrName, ext[o]);
      else
	sprintf(filename, "%s/%s-%.*s-%s-%s.%s%s", CmdLineOpts::outPrefix,
		types[o], famIdLen, parentIds[0],
		&parentIds[0][famIdLen+1], &parentIds[1][famIdLen+1],
		chrName, ext[o]);

      bool okToWrite = checkIfFileExists(filename, CmdLineOpts::forceWrite);
      haveAnOutput = haveAnOutput || okToWrite;
      if (okToWrite) {
	resultsFiles[o] = fopen(filename, "w");
	if(!resultsFiles[o]) {
	  fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
	  perror("open");
	  fprintf(log, "ERROR: couldn't open %s for writing!\n", filename);
	  fprintf(log, "open: %s\n", strerror(errno));
	}
      }
    }
  }

  return haveAnOutput;
}

bool checkIfFileExists(char *filename, bool printWarning) {
  struct stat buffer;
  if (stat( filename, &buffer) == 0) {
    if (printWarning) {
      fprintf(stderr, "WARNING: file %s exists; --force set, so overwriting\n",
	      filename);
      return true;
    }
    else {
      fprintf(stderr, "ERROR: file %s exists; use --force to overwrite\n",
	      filename);
      return false;
    }
  }
  return true;
}
