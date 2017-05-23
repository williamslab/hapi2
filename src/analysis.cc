// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <genetio/marker.h>
#include "analysis.h"
#include "cmdlineopts.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<int>      Analysis::_informMarkers[2];
dynarray<uint64_t> Analysis::_ambigMiss[2];
dynarray<int>      Analysis::_numInformForChild[2];
std::list<Recomb>  Analysis::_recombs[2];


// Detect COs in <theFam>
void Analysis::findCOs(NuclearFamily *theFam, FILE *out, int chrIdx) {
  int numChildren = theFam->_children.length();

  // Make space for each child
  for(int p = 0; p < 2; p++) _numInformForChild[p].resize(numChildren);

  const char *chrName = Marker::getChromName(chrIdx);

  // init <_numInformForChild>
  for(int p = 0; p < 2; p++) {
    for(int c = 0; c < numChildren; c++) {
      _numInformForChild[p][c] = 0;
    }
  }
  std::list<Recomb>::iterator itr;
  bool informSeen = false;
  uint64_t prevIV = 0;
  // confidently know starting haplotypes for each parent?
  bool allChildrenSolid[2] = { false, false };

  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);
  for(int marker = firstMarker; marker <= lastMarker; marker++) {
    const PhaseVals &phase = theFam->getPhase(marker);

    if (phase.status == PHASE_OK) {
      // In various loops below, must examine whichever parent(s)
      // this marker is informative for. To do a for loop, need to figure out
      // which parents to to include in the range; either same as the
      // starting parent, or, if the marker is partly informative (PI) --
      // i.e., if hetParent == 2, we analyze both:
      int isPI = phase.hetParent >> 1;
      int startP = phase.hetParent & 1;
      int endP = isPI * 1 + (1-isPI) * startP;

      if (informSeen) {
	uint64_t recombs = prevIV ^ phase.iv;

	for(int p = startP; p <= endP; p++) {
	  for(itr = _recombs[p].begin(); itr != _recombs[p].end(); ) {
	    // check if the recombination (back to the original) is present
	    // here. If so, erase the record
	    uint64_t thisRecombVal = 1ul << (2 * itr->child + p);
	    if (recombs & thisRecombVal) {
	      // Recombined back! Not considered a crossover -- erase
	      itr = _recombs[p].erase(itr);
	      // Also no need to inspect this recombination further below
	      recombs -= thisRecombVal;
	      continue;
	    }
	    // Does the child have data here? Is missing if the first bit of
	    // the two bits allotted to it in <phase.ambigMiss> is 1
	    int childBit = 2 * itr->child;
	    if ( ((phase.ambigMiss >> childBit) & 1) == 0 ) { // non-missing?
	      if (itr->numObs + 1 == CmdLineOpts::detectCO) {
		// Valid crossover! print the info
		printCO(itr, out, p, theFam, chrName, firstMarker);

		// done with this recomb record
		itr = _recombs[p].erase(itr);
		continue;
	      }
	      else
		itr->numObs++;
	    }
	    ++itr;
	  }
	}

	while (recombs > 0) {
	  // At least one child recombined
	  // Figure out which one and from which parent
	  // Note: ambiguous recombinations will in general be non-crossovers
	  // (if real). TODO: ideally report them as ambiguous

	  // Determine the bit index of the lowest order recombination event:
	  int recombBit = ffsll(recombs) - 1;
//	  assert(recombBit >= 0);

	  int recombChild = recombBit / 2;  // integer division gives child
	  int recombParent = recombBit % 2; // even is parent 0, odd 1
	  // Should only observe recombinations at markers that are
	  // heterozygous for the parent in which the recombination occurred
	  assert(phase.hetParent == recombParent || phase.hetParent == 2);

	  // Remove this recombination from further consideration:
	  recombs -= 1 << recombBit;

	  if (!allChildrenSolid[recombParent] &&
	      _numInformForChild[recombParent][recombChild] <
							  CmdLineOpts::edgeCO) {
	    // Need at least CmdLineOpts::edgeCO consecutive markers in
	    // agreement about which haplotype a child has received before we
	    // can call a crossover. We skip this crossover and reset the
	    // count for this child
	    // TODO: document this
	    _numInformForChild[recombParent][recombChild] = 0;
	    continue;
	  }

	  // Observed a recombination. Now must ensure that there are at least
	  // CmdLineOpts::detectCO downstream informative markers observed
	  // with this status before printing the event. Store away record:
	  _recombs[recombParent].emplace_front(recombChild,
		  /*prevInform=*/_informMarkers[recombParent].length() - 1,
		  /*recombMarker=*/ marker);

	}
      }


      // Update various state for next marker
      for(int p = startP; p <= endP; p++) {
	_informMarkers[p].append(marker);
	_ambigMiss[p].append(phase.ambigMiss);
	if (!allChildrenSolid[p]) {
	  bool haveNonSolid = false; // for updating allChildrenSolid
	  // increment counts of children's informative 
	  for(int c = 0; c < numChildren; c++) {
	    // Note: we reset the count to 0 above whenever a child recombines
	    // early on on the chromosome (i.e., before allChildrenSolid[p])
	    if ( ((phase.ambigMiss >> (2 * c)) & 1) == 0 ) { // non-missing?
	      _numInformForChild[p][c]++;
	      if (_numInformForChild[p][c] < CmdLineOpts::edgeCO)
		haveNonSolid = true;
	    }
	    else
	      haveNonSolid = true;
	  }
	  allChildrenSolid[p] = !haveNonSolid;
	}
      }
      prevIV = phase.iv;
      informSeen = true;
    }
  }

  // Finally report any recombinations at the end of the chromosome. Any with
  // more than CmdLineOpts::edgeCO informative sites observed suffice
  for(int p = 0; p < 2; p++) {
    for(itr = _recombs[p].begin(); itr != _recombs[p].end(); itr++) {
      if (itr->numObs >= CmdLineOpts::edgeCO) {
	// Valid crossover! print the info
	printCO(itr, out, p, theFam, chrName, firstMarker);
      }
    }
  }

  // Reset
  for(int p = 0; p < 2; p++) {
    _recombs[p].clear();
    _informMarkers[p].clear();
    _ambigMiss[p].clear();
  }
}

// Helper function for printing detected CO events
void Analysis::printCO(std::list<Recomb>::iterator record, FILE *out, int p,
		       NuclearFamily *theFam, const char *chrName,
		       int firstMarker) {
  PersonBulk *parents[2] = { theFam->_parents->first, theFam->_parents->second};
  int numChildren = theFam->_children.length();

  // First narrow down where the previous informative marker was -- child must
  // have non-missing data
  int prevInfIdx = record->prevInfIdx;
  int childBit = 2 * record->child;
  while( (_ambigMiss[p][prevInfIdx] >> childBit) & 1 ) {//missing?
    prevInfIdx--;
    // To have a recomb record child must have several
    // informative markers in succession so this must not ever
    // go off the beginning of the chromosome
//    assert(prevInfIdx >= 0);
  }
  int prevInform = _informMarkers[p][prevInfIdx];
  fprintf(out, "CO: %s %d %d %s %s %d %s %d %s\n",
	  parents[p]->getId(), p, numChildren,
	  theFam->_children[ record->child ]->getId(), chrName,
	  prevInform - firstMarker, Marker::getMarker(prevInform)->getName(),
	  record->recombMarker - firstMarker,
	  Marker::getMarker(record->recombMarker)->getName());
}
