// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <sys/stat.h>
#include <string.h>
#include <zlib.h>
#include <genetio/marker.h>
#include <genetio/personbulk.h>
#include <genetio/personio.h>
#include <genetio/util.h>
#include "cmdlineopts.h"

void checkIfFileExists(char *filename, bool printWarning);
void printPhaseType(FILE *out);

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  char filename[FILENAME_LEN];

  // Ensure that we'll be able to print to the output file at the end:
  if (strlen(CmdLineOpts::outFile) >FILENAME_LEN - 8) {//8 chars for .phgeno\0
    fprintf(stderr, "ERROR: output filename too long!");
    exit(1);
  }
  if (CmdLineOpts::vcfOutput) {
    sprintf(filename, "%s.vcf.gz", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
  }
  else if (CmdLineOpts::useImpute2Format) {
    // Check whether the .haps output file exists:
    sprintf(filename, "%s.haps.gz", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
    // Check whether the .sample output file exists:
    sprintf(filename, "%s.sample", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
  }
  else {
    // Check whether the phgeno output file exists:
    sprintf(filename, "%s.phgeno.gz", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
    // Check whether the phind output file exists:
    sprintf(filename, "%s.phind", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
    // Check whether the phsnp output file exists:
    sprintf(filename, "%s.phsnp", CmdLineOpts::outFile);
    checkIfFileExists(filename, CmdLineOpts::forceWrite);
  }

  // open the .log file for writing
  sprintf(filename, "%s.log", CmdLineOpts::outFile);
  FILE *log = fopen(filename, "w");
  if (!log) {
    fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
    exit(1);
  }

  //////////////////////////////////////////////////////////////////////////
  // Output files don't already exist, can go forward!

  // Print status to stdout and the log
  FILE *outs[2] = {stdout, log };
  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "\n");
    fprintf(out, "HAPI v%s!    (Released %s)\n\n", VERSION_NUMBER,
	    RELEASE_DATE);

    if (CmdLineOpts::vcfInput) {
      fprintf(out, "VCF input file:\t  %s\n", CmdLineOpts::genoFile);
    }
    else {
      fprintf(out, "Genotype file:\t  %s\n", CmdLineOpts::genoFile);
      fprintf(out, "SNP file:\t  %s\n", CmdLineOpts::markerFile);
      fprintf(out, "Individual file:  %s\n", CmdLineOpts::indFile);
    }
    fprintf(out, "Output prefix:\t  %s\n", CmdLineOpts::outFile);

    fprintf(out, "\n");

    if (CmdLineOpts::onlyChr != NULL) {
      printf("Chromosome:\t  %s\n", CmdLineOpts::onlyChr);
    }
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

  // open the .phgeno file for writing before phasing, even though we won't
  // write it until phasing completes.  This ensures we have write
  // permissions for this most important file (others can be generated by the
  // user from the input if they fail to write later).
  gzFile gzout;
  if (CmdLineOpts::useImpute2Format) {
    sprintf(filename, "%s.haps.gz", CmdLineOpts::outFile);
    gzout = gzopen(filename, "w");
    if (!gzout) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      exit(1);
    }
  }
  else {
    sprintf(filename, "%s.phgeno.gz", CmdLineOpts::outFile);
    gzout = gzopen(filename, "w");
    if (!gzout) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      exit(1);
    }
  }

  printf("\nPhasing families with two or more children out of %lu families... ",
	 PersonBulk::numFamilies());

  // Phase!
  PersonBulk::fam_ht_iter iter = PersonBulk::familyIter();
  for( ; iter != PersonBulk::familyIterEnd(); iter++) {
    dynarray<PersonBulk*> *children = iter->second;
    if (children->length() > 1) {
      // phase the current family:
      // require at least two children as trios are better to phase in a
      // population context
      // TODO:
//      Phaser::run(iter->first, children);
    }
  }

  // Finished phasing all families!
  printf("done.\n");

  // TODO: after phasing results are stored, remove this exit call
  exit(0);

  mult_printf(outs, "\nPrinting... ");
  fflush(stdout);

  if (CmdLineOpts::useImpute2Format) {
    PersonIO<PersonBulk>::printGzImpute2Haps(gzout);
    gzclose(gzout);

    // Print sample file:
    sprintf(filename, "%s.sample", CmdLineOpts::outFile);
    FILE *out = fopen(filename, "w");
    if (!out) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      perror("open");
    }
    else {
      PersonIO<PersonBulk>::printImpute2SampleFile(out);
      fclose(out);
    }
  }
  else {
    // Print final haplotypes to phgeno file (opened above):
    PersonIO<PersonBulk>::printGzEigenstratPhased(gzout);
    gzclose(gzout);

    // Print phind file:
    sprintf(filename, "%s.phind", CmdLineOpts::outFile);
    FILE *out = fopen(filename, "w");
    if (!out) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      perror("open");
    }
    else {
      PersonIO<PersonBulk>::printPhasedIndFile(out);
      fclose(out);
    }
    // Print phsnp file:
    sprintf(filename, "%s.phsnp", CmdLineOpts::outFile);
    out = fopen(filename, "w");
    if (!out) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      perror("open");
    }
    else {
      Marker::printSNPFile(out);
      fclose(out);
    }
  }

  for(int o = 0; o < 2; o++) {
    FILE *out = outs[o];
    fprintf(out, "done.\n");
  }

  fclose(log);
}

void checkIfFileExists(char *filename, bool printWarning) {
  struct stat buffer;
  if (stat( filename, &buffer) == 0) {
    if (printWarning) {
      fprintf(stderr, "WARNING: output filename %s exists; --force set, so overwriting\n",
	      filename);
    }
    else {
      fprintf(stderr, "ERROR: output filename %s exists; won't overwrite so dying...\n",
	      filename);
      exit(1);
    }
  }
}
