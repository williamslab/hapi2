// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <genetio/nuclearfamily.h>
#include <genetio/dynarray.h>
#include <list>

#ifndef ANALYSIS_H
#define ANALYSIS_H

struct Recomb {
  Recomb(uint8_t c, int p, int r) {
    child = c; numObs = 1;
    prevInfIdx = p; recombMarker = r;
  }
  uint8_t child;   // which child index received recombination?
  uint8_t numObs;  // how many observed recombined?
  int prevInfIdx;  // index into _inform
  int recombMarker;
};

class Analysis {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void findCOs(NuclearFamily *theFam, FILE *out, int chrIdx);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void printCO(std::list<Recomb>::iterator record, FILE *out, int p,
			NuclearFamily *theFam, const char *chrName,
			int firstMarker);

    // List storing marker index of informative positions for the two parents
    // (indexed) that are upstream of the current position in an analysis
    static dynarray<int> _informMarkers[2];

    // Parallel array (technically arrays for the first index) with
    // <_informMarkers> that stores the <phase.ambigMiss> value corresponding
    // to the given informative marker
    static dynarray<uint64_t> _ambigMiss[2];

    // Number of consistent informative markers observed for each parent/child
    // combination. This is only relevant at the beginning of the chromosome
    // where we are trying to establish what the background haplotype is. When
    // all children have greater than CmdLineOpt::edgeCO we
    // use a simple boolean to indicate things are trustworthy
    static dynarray<int> _numInformForChild[2];

    // Records of recombinations encountered that may or may not have enough
    // evidence to be called a crossover
    static std::list<Recomb> _recombs[2];
};

#endif // ANALYSIS_H
