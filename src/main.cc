// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>
#include <vector>
#include <genetio/marker.h>
#include <genetio/personbulk.h>
#include <genetio/personio.h>
#include <genetio/util.h>
#include "cmdlineopts.h"
#include "phaser.h"
#include "analysis.h"

void  createOutDir();
FILE *setupJsonOutput(const char *jsonOutName, char *filename, FILE **outs,
		      bool noOtherJson);
void  getFamiliesToBePhased(FILE *log, dynarray<NuclearFamily *> &toBePhased);
bool  openFilesToWrite(char *&filename, FILE *resultsFiles[6], int chrIdx,
		       const char *parentIds[2], int famIdLen, int totalFileLen,
		       int &allocFilenameLen, FILE *log);
bool  checkIfFileExists(char *filename, bool printWarning);

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

    if (CmdLineOpts::onlyChr != NULL || CmdLineOpts::startPos != 0 ||
	CmdLineOpts::endPos != INT_MAX) {
      if (CmdLineOpts::onlyChr != NULL)
	fprintf(out, "Chromosome:\t\t%s\n", CmdLineOpts::onlyChr);
      if (CmdLineOpts::startPos != 0)
	fprintf(out, "Start position:\t\t%d\n", CmdLineOpts::startPos);
      if (CmdLineOpts::endPos != INT_MAX)
	fprintf(out, "End position:\t\t%d\n", CmdLineOpts::endPos);
      fprintf(out, "\n");
    }

    fprintf(out, "Errors:\t\t\t");
    if (CmdLineOpts::maxNoErrRecombs == 0)
      fprintf(out, "disabled\n\n");
    else {
      fprintf(out, "enabled\n");
      fprintf(out, "Number of recombinations to trigger error:\t%d\n",
	      CmdLineOpts::maxNoErrRecombs);
      fprintf(out, "Maximum informative markers per error:\t\t%d\n",
	      CmdLineOpts::errorLength);
      fprintf(out, "\n");
    }

    if (CmdLineOpts::minNumParentsData > 0 ||
	CmdLineOpts::minNumChildrenData > 2) {
      if (CmdLineOpts::minNumParentsData > 0)
	fprintf(out, "Minimum number of parents per family:\t\t%d\n",
		CmdLineOpts::minNumParentsData);
      if (CmdLineOpts::minNumChildrenData > 2)
	fprintf(out, "Minimum number of children per family:\t\t%d\n",
		CmdLineOpts::minNumChildrenData);
      fprintf(out, "\n");
    }

    if (CmdLineOpts::forceMissingParBits) {
      fprintf(out, "Omitting parent data for:\t\t\t");
      switch (CmdLineOpts::forceMissingParBits) {
	case 1:
	  fprintf(out, "dad\n\n");
	  break;
	case 2:
	  fprintf(out, "mom\n\n");
	  break;
	case 3:
	  fprintf(out, "both parents\n\n");
	  break;
	default:
	  fprintf(out, "ERROR: value of %d\n\n",
		  CmdLineOpts::forceMissingParBits);
	  exit(1);
	  break;
      }
    }

    fprintf(out, "Output to be generated in output directory:\n");
    if (CmdLineOpts::vcfOutput) {
      fprintf(out, "  VCF format haplotypes\n");
    }
    if (CmdLineOpts::pedOutput) {
      fprintf(out, "  PLINK ped format haplotypes\n");
    }
    if (CmdLineOpts::txtOutput) {
      fprintf(out, "  Text format haplotypes\n");
    }
    if (CmdLineOpts::jsonOutput) {
      fprintf(out, "  JSON format haplotypes\n");
    }
    if (CmdLineOpts::jsonParOutput) {
      fprintf(out, "  JSON format parent haplotypes\n");
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

  FILE *jsonParFile = (CmdLineOpts::jsonParOutput) ?
			  setupJsonOutput("parent_phase.json", filename, outs,
			      /*noOtherJson=*/ !CmdLineOpts::jsonOutput) :
			  NULL;
  FILE *jsonFile = (CmdLineOpts::jsonOutput) ?
			  setupJsonOutput("all.json", filename, outs,
			      /*noOtherJson=*/ !jsonParFile) :
			  NULL;

  //////////////////////////////////////////////////////////////////////////
  // Have output directory and log, can go forward!

  // Get list of families to be phased (based on minimum number of children and
  // parents with data):
  dynarray<NuclearFamily *> toBePhased;
  getFamiliesToBePhased(log, toBePhased);

  int numFamsToBePhased = toBePhased.length();

  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "\n");
//    fprintf(out, "Have data for %d individuals.",
//	    PersonBulk::_allIndivs.length());
    fprintf(out, "Have data for %d nuclear families with >= %d children",
	    numFamsToBePhased, CmdLineOpts::minNumChildrenData);
    if (CmdLineOpts::minNumParentsData == 1)
      fprintf(out, " and >= 1 parent");
    if (CmdLineOpts::minNumParentsData == 2)
      fprintf(out, " and 2 parents");
    fprintf(out, ".\n\n");
  }
  if (numFamsToBePhased == 0) {
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "ERROR: no families with %d or more children to phase so exiting.\n",
	      CmdLineOpts::minNumChildrenData);
    exit(5);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Phase!
  int numChrs = Marker::getNumChroms();

  Phaser::init();
  int numFinished = 0;
  for(int f = 0; f < numFamsToBePhased; f++) {

    for(int o = 0; o < 2; o++) {
      FILE *out = outs[o];
      fprintf(out, "Phasing families... %d / %d",
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
    theFam->initPhase();
    if (CmdLineOpts::verbose) {
      fprintf(log, "Phasing family with parents %s and %s\n",
	      parentIds[0], parentIds[1]);
    }

    FILE *resultsFiles[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
    for(int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
      bool shouldPhase = openFilesToWrite(filename, resultsFiles, chrIdx,
			  parentIds, famIdLen,
			  /*totalFileLen=*/ prefixLen + famSpecificFileLen,
			  allocFilenameLen, log);
      shouldPhase = shouldPhase || jsonParFile || jsonFile;

      if (!shouldPhase) {
	for(int o = 0; o < 2; o++) {
	  FILE *out = outs[o];
	  fprintf(out, "Unable to write any output for chr %s: will not phase it\n",
		  Marker::getChromName(chrIdx));
	}
	continue;
      }

      if (CmdLineOpts::verbose) {
	fprintf(log, "  Chromosome %s:\n", Marker::getChromName(chrIdx));
      }

      Phaser::run(theFam, chrIdx, log);

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
    if (CmdLineOpts::vcfOutput && resultsFiles[5]) {
      theFam->printPhasedVCF(resultsFiles[5], "HAPI v" VERSION_NUMBER);
      fclose(resultsFiles[5]);
    }
    if (CmdLineOpts::pedOutput && resultsFiles[3] && resultsFiles[4]) {
      theFam->printPhasedPed(resultsFiles[3]);
      Marker::printMapFile(resultsFiles[4]);
      fclose(resultsFiles[3]);
      fclose(resultsFiles[4]);
    }
    if (CmdLineOpts::jsonParOutput && jsonParFile) {
      if (f > 0)
	fprintf(jsonParFile, ","); // need comma between each JSON entry
      theFam->printHapJson(jsonParFile, /*withChildren=*/ false);
    }
    if (CmdLineOpts::jsonOutput && jsonFile) {
      if (f > 0)
	fprintf(jsonFile, ","); // need comma between each JSON entry
      theFam->printHapJson(jsonFile, /*withChildren=*/ true);
    }
    theFam->deletePhase();

    numFinished++;
  }

  // Finished phasing all families!
  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "Phasing families... done.");
  }
  printf("               \n");
  fprintf(log, "\n");

  fclose(log);

  if (jsonParFile) {
    fprintf(jsonParFile, "}\n");
    fclose(jsonParFile);
  }
  if (jsonFile) {
    fprintf(jsonFile, "}\n");
    fclose(jsonFile);
  }

//  Marker::cleanUp();
//  PersonIO<PersonBulk>::cleanUp();
//  NuclearFamily::cleanUp();
//  delete [] filename;
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

// Open and print the meta data for JSON haplotype output
FILE *setupJsonOutput(const char *jsonOutName, char *filename, FILE **outs,
		      bool noOtherJson) {
  FILE *out = NULL;

  sprintf(filename, "%s/%s", CmdLineOpts::outPrefix, jsonOutName);
  bool okToWrite = checkIfFileExists(filename, CmdLineOpts::forceWrite);
  if (okToWrite) {
    out = fopen(filename, "w");
    if (!out) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      perror("open");
      fprintf(outs[1], "ERROR: couldn't open %s for writing!\n", filename);
      fprintf(outs[1], "open: %s\n", strerror(errno));
    }
  }

  if (!out) {
    if (noOtherJson) {
      bool nonJsonOutput = CmdLineOpts::txtOutput || CmdLineOpts::pedOutput ||
			   CmdLineOpts::vcfOutput || CmdLineOpts::ivOutput ||
			   CmdLineOpts::detectCO;
      if (!nonJsonOutput) {
	// No output possible: exit
	for(int o = 0; o < 2; o++) {
	  FILE *out = outs[o];
	  fprintf(out, "Unable to write any output: will not phase\n");
	}
	exit(5);
      }
    }
    else
      return NULL;
  }
  else {
    // Marker ids:
    fprintf(out, "{\"marker\":[");
    fprintf(out, "\"%s\"", Marker::getMarker(0)->getName());
    for(int m = 1; m < Marker::getNumMarkers(); m++)
      fprintf(out, ",\"%s\"", Marker::getMarker(m)->getName());
    fprintf(out, "],");

    // Physical positions:
    fprintf(out, "\"physpos\":[");
    fprintf(out, "%d", Marker::getMarker(0)->getPhysPos());
    for(int m = 1; m < Marker::getNumMarkers(); m++)
      fprintf(out, ",%d", Marker::getMarker(m)->getPhysPos());
    fprintf(out, "],");

    // Chromosome names:
    fprintf(out, "\"chr\":[");
    fprintf(out, "\"%s\"", Marker::getMarker(0)->getChromName());
    for(int m = 1; m < Marker::getNumMarkers(); m++)
      fprintf(out, ",\"%s\"", Marker::getMarker(m)->getChromName());
    fprintf(out, "],");

    // Chromosome starts:
    fprintf(out, "\"chrstr\":{");
    fprintf(out, "\"%s\":0", Marker::getMarker(0)->getChromName());
    for(int chrIdx = 1; chrIdx < Marker::getNumChroms(); chrIdx++)
      fprintf(out, ",\"%s\":%d", Marker::getChromName(chrIdx),
	      Marker::getFirstMarkerNum(chrIdx));
    fprintf(out, "},");
  }

  return out;
}

// Get list of families to be phased based on minimum number of children and
// parents with data
void getFamiliesToBePhased(FILE * log, dynarray<NuclearFamily *> &toBePhased) {
  std::vector<int> counts[3]; // for 0, 1, 2 parents

  for(int p = 0; p < 3; p++)
    // note: current limit is 32 children (in phaser.cc)
    for(int c = 0; c <= 32; c++)
      counts[p].push_back(0);

  for(NuclearFamily::fam_ht_iter iter = NuclearFamily::familyIter();
			       iter != NuclearFamily::familyIterEnd(); iter++) {
    NuclearFamily *theFam = iter->second;
    if (theFam->numChildren() > 1) {
      bool shouldPhase = true; // initially assume

      int childrenWithData = 0;
      for(int c = 0; c < theFam->numChildren(); c++)
	if (theFam->_children[c]->hasData())
	  childrenWithData++;

      if (childrenWithData < CmdLineOpts::minNumChildrenData)
	shouldPhase = false;

      if (childrenWithData > 32) {
	fflush(stdout);
	fprintf(stderr, "\n");
	fprintf(stderr, "ERROR: cannot phase more than 32 children in a family.\n");
	fprintf(stderr, "       changing to >64 bits for inheritance vectors would fix this\n");
	exit(9);
      }

      int parentsWithData = 0;
      if (theFam->_parents->first->hasData())
	parentsWithData++;
      if (theFam->_parents->second->hasData())
	parentsWithData++;

      if (CmdLineOpts::minNumParentsData > 0 &&
			      parentsWithData < CmdLineOpts::minNumParentsData)
	shouldPhase = false;

      counts[parentsWithData][childrenWithData]++;

//      if (parentsWithData == 2 && childrenWithData >= 3)
//	fprintf(stderr, "%s %s %d\n", theFam->_parents->first->getId(),
//		theFam->_parents->second->getId(), childrenWithData);

      if (shouldPhase)
	toBePhased.append(theFam);
    }
  }

  bool headerPrinted = false;
  for(int p = 0; p < 3; p++) {
    for(int c = 0; c <= 32; c++) {
      if (counts[p][c] > 0) {
	if (!headerPrinted) {
	  fprintf(log, "\nN_parents\tN_child\tCount\n");
	  fprintf(log, "------------------------------------------\n");
	  headerPrinted = true;
	}
	fprintf(log, "%d\t%d\t%d\n", p, c, counts[p][c]);
      }
    }
  }
  fprintf(log, "\n");
}


// Attempts to open the files to be printed to, and alerts the user if any
// exist. Returns true if at least one output file exists
bool openFilesToWrite(char *&filename, FILE *resultsFiles[6], int chrIdx,
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
  const char *prefix[6] = { "hap-", "iv-", "co-", "", "", "" };
  const char *ext[6] = { ".txt", ".csv", "", ".ped", ".map", ".vcf" };
  const int  wantType[6] = { CmdLineOpts::txtOutput, CmdLineOpts::ivOutput,
			     CmdLineOpts::detectCO, CmdLineOpts::pedOutput,
			     CmdLineOpts::pedOutput, CmdLineOpts::vcfOutput };
  // Do we want separate files for each chromsome? Not for all file types
  const bool separateChrs[6] = { true, true, true, false, false, false };
  for(int o = 0; o < 6; o++) {
    if (!separateChrs[o] && chrIdx != 0)
      // Will only open the single file corresponding to this file type once,
      // so skip now that we're past the first chromosome
      continue;

    if (wantType[o]) {
      if (famIdLen == 0)
	sprintf(filename, "%s/%s%s-%s%s%s%s", CmdLineOpts::outPrefix,
		prefix[o], parentIds[0], parentIds[1],
		(separateChrs[o]) ? "." : "",
		(separateChrs[o]) ? chrName : "", ext[o]);
      else
	sprintf(filename, "%s/%s%.*s-%s-%s%s%s%s", CmdLineOpts::outPrefix,
		prefix[o], famIdLen, parentIds[0],
		&parentIds[0][famIdLen+1], &parentIds[1][famIdLen+1],
		(separateChrs[o]) ? "." : "",
		(separateChrs[o]) ? chrName : "", ext[o]);

      bool okToWrite = checkIfFileExists(filename, CmdLineOpts::forceWrite);
      if (okToWrite) {
	resultsFiles[o] = fopen(filename, "w");
	if(!resultsFiles[o]) {
	  fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
	  perror("open");
	  fprintf(log, "ERROR: couldn't open %s for writing!\n", filename);
	  fprintf(log, "open: %s\n", strerror(errno));
	}
      }
      else
	resultsFiles[o] = NULL;
    }
  }
  // Any of these files being writable gives us an output:
  haveAnOutput = resultsFiles[0] || resultsFiles[1] || resultsFiles[2] ||
								resultsFiles[5];
  // need both the ped and map files writable to have an output
  haveAnOutput = haveAnOutput || (resultsFiles[3] && resultsFiles[4]);

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
