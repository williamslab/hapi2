// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <genetio/marker.h>
#include "analysis.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<int>     Analysis::_informMarkers[2];
dynarray<int>     Analysis::_numInformForChild[2];
std::list<Recomb> Analysis::_recombs[2];


// Detect COs in <theFam>
void Analysis::findCOs(NuclearFamily *theFam, FILE *out) {
  PersonBulk *parents[2] = { theFam->_parents->first, theFam->_parents->second};
  int numChildren = theFam->_children.length();

  // Make space for each child
  for(int p = 0; p < 2; p++) _numInformForChild[p].resize(numChildren);

  int numChrs = Marker::getNumChroms();
  for(int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
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
	if (informSeen) {
	  uint64_t recombs = prevIV ^ phase.iv;

	  // Must examine/update recombination records for whichever parent(s)
	  // this marker is informative for. To do a for loop, need to figure
	  // out which parents to to include in the range; either same as the
	  // starting parent, or, if the marker is partly informative (PI) --
	  // i.e., if hetParent == 2, we analyze both:
	  int isPI = phase.hetParent >> 1;
	  int startP = phase.hetParent & 1;
	  int endP = isPI * 1 + (1-isPI) * startP;
	  for(int p = startP; p <= endP; p++) {
	    for(itr = _recombs[p].begin(); itr != _recombs[p].end(); ) {
	      // check if the recombination (back to the original) is present
	      // here. If so, erase the record
	      uint64_t thisRecombVal = 1ul << (2 * itr->child + p);
	      if (recombs & thisRecombVal) {
		// Recombined back! Not considered a crossover -- erase
		itr = _recombs[p].erase(itr);
		// Also no need to inspect this recombination here
		recombs -= thisRecombVal;
		continue;
	      }
	      // Does the child have data here? If so, increment the obs count
	      uint8_t geno = theFam->_children[itr->child]->getBitGeno(marker);
	      if (geno != G_MISS) {
		// TODO! make this 10 value a user-specified option
		if (itr->numObs + 1 == 10) {
		  // Valid crossover! print the info and erase the record

		  // First narrow down where the previous informative marker was
		  // -- child must have non-missing data
		  int prevInfIdx = itr->prevInfIdx;
		  int prevInform = _informMarkers[p][prevInfIdx];
		  while(theFam->_children[itr->child]->getBitGeno(prevInform) ==
								       G_MISS) {
		    prevInfIdx--;
		    // To have a recomb record child must have several
		    // informative markers in succession so this must not ever
		    // go off the beginning of the chromosome
		    assert(prevInfIdx >= 0);
		    prevInform = _informMarkers[p][prevInfIdx];
		  }
		  // TODO: don't subtract firstMarker
		  fprintf(out, "CO: %s %d %d %s %s %d %s %d %s\n",
			  parents[p]->getId(), p, numChildren,
			  theFam->_children[ itr->child ]->getId(), chrName,
			  prevInform - firstMarker,
			  Marker::getMarker(prevInform)->getName(),
			  itr->recombMarker - firstMarker,
			  Marker::getMarker(itr->recombMarker)->getName());

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
	    // TODO! deal with ambiguous bits

	    // Determine the bit index of the lowest order recombination event:
	    int recombBit = ffsll(recombs) - 1;
	    assert(recombBit >= 0);

	    int recombChild = recombBit / 2;  // integer division gives child
	    int recombParent = recombBit % 2; // even is parent 0, odd 1
	    // Should only observe recombinations at markers that are
	    // heterozygous for the parent in which the recombination occurred
	    assert(phase.hetParent == recombParent || phase.hetParent == 2);

	    // Remove this recombination from further consideration:
	    recombs -= 1 << recombBit;

	    // TODO: use CmdLineOpts here (search for 10 in comments and in
	    // the .h file)
	    if (!allChildrenSolid[recombParent] &&
			  _numInformForChild[recombParent][recombChild] < 10) {
	      // Need at least 10 consecutive markers in agreement about which
	      // haplotype a child has received before we can call a crossover.
	      // We skip this crossover and reset the count for this chlid
	      _numInformForChild[recombParent][recombChild] = 0;
	      continue;
	    }

	    // Observed a recombination. Ensure that there are at least 10
	    // informative markers observed with this status before printing.
	    // Store away record:
	    _recombs[recombParent].emplace_front(recombChild,
		  /*prevInform=*/_informMarkers[recombParent].length() - 1,
		  /*recombMarker=*/ marker);

	  }
	}


	// Update various state for next marker
	// Similar code above to determine bounds on parent loop (see comment
	// there)
	int isPI = phase.hetParent >> 1;
	int startP = phase.hetParent & 1;
	int endP = isPI * 1 + (1-isPI) * startP;
	for(int p = startP; p <= endP; p++) {
	  _informMarkers[p].append(marker);
	  if (!allChildrenSolid[p]) {
	    bool haveNonSolid = false; // for updating allChildrenSolid
	    // increment counts of children's informative 
	    for(int c = 0; c < numChildren; c++) {
	      // Note: we reset the count to 0 whenever a child recombines early
	      // on on the chromosome (i.e., before allChildrenSolid[p])
	      if (theFam->_children[c]->getBitGeno(marker) != G_MISS) {
		_numInformForChild[p][c]++;
		if (_numInformForChild[p][c] < 10)
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

    // Reset
    for(int p = 0; p < 2; p++) {
      _recombs[p].clear();
      _informMarkers[p].clear();
    }
  }
}
