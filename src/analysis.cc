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
dynarray<int>      Analysis::_numInformForChild[2];
std::list<Recomb>  Analysis::_recombs[2];


// Detect COs in <theFam>, on chromosome with index <chrIdx>, and print
// detected COs to <out>
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
  uint64_t prevIVUnassigned = UINT64_MAX;
  // confidently know starting haplotypes for each parent?
  bool allChildrenSolid[2] = { false, false };

  // Bits that correspond to either or both parents in IV values:
  uint64_t allBitsSet = ~0ul; // initially: fewer depending on numChildren
  if (numChildren < 32)
    allBitsSet &= (1ul << (2*numChildren)) - 1;
  uint64_t parBits[3];
  parBits[0] = 0x5555555555555555 & allBitsSet;
  parBits[1] = 0xAAAAAAAAAAAAAAAA & allBitsSet;
  parBits[2] = allBitsSet;

  // print header
  fprintf(out, "#event parent_id parent_sex num_child recipient_id chrom marker_num_before marker_id_before marker_num_after marker_id_after\n");

  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);
  for(int marker = firstMarker; marker <= lastMarker; marker++) {
    const PhaseVals &phase = theFam->getPhase(marker);

    uint8_t ambigParHet = phase.ambigParHet | (1 << phase.hetParent);

    // Have an informative marker so long as the status is PHASE_OK.
    // Also need for the marker to be conssitently informative for one or the
    // other parents (or both): if ambigParHet has the two lowest order bits
    // set, it could be inforamtive for either parent. As such, in that case,
    // it's ambiguous and we treat it as informative for neither (skip it)
    if (phase.status == PHASE_OK && (ambigParHet & 3) != 3) {
      // Are both parents potentially heterozygous? Yes if either the bit for
      // both parents heterozygous is set to 1 or if both the other two bits
      // are set to 1.
      uint8_t bothParHet = (ambigParHet >> 2) |
			    (((ambigParHet & 2) >> 1) & (ambigParHet & 1));
      // For determining which set of prior recombinations to examine, which
      // parent index should we start from? 0 if <bothParHet>, but otherwise
      // it suffices to shift by 1 which will give 0 if parent 0 only is
      // heterozygous and 1 if only parent 1 is.
      int startP = (1 - bothParHet) * (ambigParHet >> 1);
      // What is the upper bound of the loop over informative parents? Should
      // either be 1 if <bothParHet>. Otherwise, for only one parent, we only
      // go through the loop once, so <startP> == <endP>.
      int endP = bothParHet + (1 - bothParHet) * startP;

      // TODO: remove
      if (ambigParHet & 4) {
	assert(startP == 0 && endP == 1);
      }

      // Which IV bits are uninformative at this marker? Any that are in
      // ivFlippable and any corresponding to children that have mssing data
      uint64_t curUninformIV = phase.ivFlippable |
					  ((phase.ambigMiss & parBits[0]) * 3);
      // Also not informative for IV values corresponding to any
      // non-heterozygous parent
      if (!bothParHet) { // TODO: optimize
	uint8_t informPar = startP;
	uint8_t otherPar = informPar ^ 1;
	curUninformIV |= parBits[otherPar];
      }
      if (informSeen) {
	// which haplotypes recombined? Will only inspect those where the
	// IV value is fixed in the state -- that is, bits where ivFlippable ==0
	// -- and where the IV has an assigned value -- that is, we've
	// encountered a state where the corresponding bit in ivFlippable == 0
	uint64_t recombs = (prevIV ^ phase.iv) &
					    ~(curUninformIV | prevIVUnassigned);

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
	    // non-missing and non-flippable?
	    if (((phase.ambigMiss >> childBit) & 1) == 0 && // non-missing?
	        ((phase.ivFlippable >>(childBit+p)) & 1) == 0) {//non-flippable?
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
//	  assert(phase.hetParent == recombParent || phase.hetParent == 2);
	  assert((ambigParHet & (1 << recombParent)) |
						      (ambigParHet & (1 << 2)));

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
	if (!allChildrenSolid[p]) {
	  bool haveNonSolid = false; // for updating allChildrenSolid
	  // increment counts of children's informative 
	  for(int c = 0; c < numChildren; c++) {
	    // Note: we reset the count to 0 above whenever a child recombines
	    // early on on the chromosome (i.e., before allChildrenSolid[p])
	    if (((phase.ambigMiss >> (2 * c)) & 1) == 0 && // non-missing?
		((phase.ivFlippable >> (2 * c + p)) & 1) == 0) {//non-flippable?
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
      // We'll propagate forward the IV value for bits for which the current
      // marker is uninformative. Also must update unassigned bits: any for
      // which the current marker is uninformative remain unassigned.
      prevIV &= curUninformIV;
      prevIV |= phase.iv & ~curUninformIV;
      prevIVUnassigned &= curUninformIV;

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
  }
}

// Helper function for printing detected CO events
void Analysis::printCO(std::list<Recomb>::iterator record, FILE *out, int p,
		       NuclearFamily *theFam, const char *chrName,
		       int firstMarker) {
  PersonBulk *parents[2] = { theFam->_parents->first, theFam->_parents->second};
  int numChildren = theFam->_children.length();

  // First narrow down where the previous informative marker was -- child must
  // have non-missing data and their inheritance vector must be non-flippable
  int prevInfIdx = record->prevInfIdx;
  int childBit = 2 * record->child;
  int parMask = 1 << p; // TODO: comment
  while(prevInfIdx >= 0) {
    int marker = _informMarkers[p][prevInfIdx];
    const PhaseVals &phase = theFam->getPhase(marker);
    uint64_t ambigMiss = phase.ambigMiss;
    uint64_t ivFlippable = phase.ivFlippable;
    if ( !( ((ambigMiss >> childBit) & 1) | // missing?
	    ((ivFlippable >> childBit) & parMask) ) ) { //flippable?
      uint8_t ambigParHet = phase.ambigParHet | (1 << phase.hetParent);
      assert((ambigParHet & (1 << p)) | (ambigParHet & (1 << 2)));
      break;
    }
    prevInfIdx--;
  }

  if (prevInfIdx < 0) {
    fprintf(out, "CO: %s %d %d %s %s %d %s %d %s\n",
	    parents[p]->getId(), p, numChildren,
	    theFam->_children[ record->child ]->getId(), chrName,
	    -1, "NA", record->recombMarker - firstMarker,
	    Marker::getMarker(record->recombMarker)->getName());
  }
  else {
    int prevInform = _informMarkers[p][prevInfIdx];
    fprintf(out, "CO: %s %d %d %s %s %d %s %d %s\n",
	    parents[p]->getId(), p, numChildren,
	    theFam->_children[ record->child ]->getId(), chrName,
	    prevInform - firstMarker, Marker::getMarker(prevInform)->getName(),
	    record->recombMarker - firstMarker,
	    Marker::getMarker(record->recombMarker)->getName());
  }
}
