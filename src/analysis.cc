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
      // Determine which parents this marker is necessarily informative for. If
      // the two lowest order bits are set then it's heterozygous for only one
      // parent (though there maybe an ambiguous way in which it can be
      // heterozygous for both).
      uint8_t oneParHetBits = ambigParHet & 3;
      // The following is 1 if either of the two lowest order bits are set and 0
      // otherwise. Thus it is 0 if the only possible heterozygous parent status
      // is both parents het.
      uint8_t haveOneParHet = ((oneParHetBits & 2) >> 1) ^ (oneParHetBits & 1);
      // start from 0 for bit 0 set; 1 for bit 1 set; and 0 if both parents are
      // heterozygous (when niether is set)
      int startP = oneParHetBits >> 1;
      // upper bound -- strictly greater than -- for the parent values that are
      // heterozygous is 1 for bit 0 set; 2 for bit 1 set; and 2 if both parents
      // are heterozgyous
      int boundP = haveOneParHet * oneParHetBits + (1 - haveOneParHet) * 2;

      if (informSeen) {
	// which haplotypes recombined? Will only inspect those where the
	// IV value is fixed in the state -- that is, bits where ivFlippable ==0
	// -- and where the IV has an assigned value -- that is, we've
	// encountered a state where the corresponding bit in ivFlippable == 0
	uint64_t recombs = (prevIV ^ phase.iv) &
					~(phase.ivFlippable | prevIVUnassigned);

	for(int p = startP; p < boundP; p++) {
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
      for(int p = startP; p < boundP; p++) {
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
      // We'll propagate forward the IV value for bits that are flippable at the
      // current marker. Also update unassigned bits: they are those for which
      // all previous markers have had ivFlippable = 1
      prevIV &= phase.ivFlippable;
      prevIV |= phase.iv & ~phase.ivFlippable;
      prevIVUnassigned &= phase.ivFlippable;

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
