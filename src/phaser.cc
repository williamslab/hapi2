// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <bitset>
#include "phaser.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray< dynarray<State*> > Phaser::_hmm;
dynarray<int> Phaser::_hmmMarker;
Phaser::state_ht Phaser::_stateHash;
uint64_t Phaser::_parBits[3];
uint64_t Phaser::_flips[4];
uint64_t Phaser::_ambigFlips[4];


// TODO: remove? At least move to a class if we keep this
void printGeno(FILE *out, int marker, uint8_t bitGeno);


void Phaser::run(PersonBulk::par_pair parents, dynarray<PersonBulk*> &children,
		 int chrIdx) {
  // Get first and last marker numbers for the chromosome to be phased:
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);

  if (children.length() > 32) {
    fflush(stdout);
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: cannot phase more than 32 children in a family.\n");
    fprintf(stderr, "       changing to >64 bits for inheritance vectors would fix this\n");
    exit(9);
  }

  init(children.length());

  dynarray<State> partialStates;

  // build states/phase each marker
  for(int m = firstMarker; m <= lastMarker; m++) {
    uint8_t parentData, parentGenoTypes, childGenoTypes;
    // Each index corresponds to a genotype (see the Geno enumerated type).
    // The two bits corresponding to each child are set to 1 for the index
    // of its genotype.
    // Index 4 (the last value) stores the raw genotype data in PLINK format
    uint64_t childrenData[5];

    ///////////////////////////////////////////////////////////////////////////
    // Step 0: get the data for this marker
    getFamilyData(parents, children, m, parentData, parentGenoTypes,
		  childrenData, childGenoTypes);


    ///////////////////////////////////////////////////////////////////////////
    // Step 1: Determine marker type and check for Mendelian errors
    int mt = getMarkerType(parentGenoTypes, childGenoTypes);
    assert(mt > 0 && (mt & ( ~((1 << MT_N_TYPES) -1) )) == 0);

    if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
      assert(mt == (1 << MT_ERROR) || mt == (1 << MT_AMBIG));
      // TODO: indicate that the marker is erroneous / ambiguous
      continue;
    }

    if (mt == (1 << MT_UN)) {
      // TODO: if the children are heterozygous, phase when parent data
      // available; if not, choose an arbitrary phase and set the status somehow
      continue;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step 2: For each consistent marker type, deduce the inheritance vector
    // bits where possible and generate partial states that only have bits
    // defined for the informative parent(s)
    makePartialStates(partialStates, mt, parentData, childrenData);

    ///////////////////////////////////////////////////////////////////////////
    // Step 3: Using states for the previous informative marker, fill in
    // unknown inheritance vector bits by propagation and generate full states.
    // Simultaneously calculate the minimum recombination counts for each state
    // and identify which previous state(s) produce that minimum count (i.e.,
    // perform count-based Viterbi calculation)
    makeFullStates(partialStates, /*marker=*/ m,
		   /*childrenMiss=*/ childrenData);


    // Clean up for next iteration
    partialStates.clear();
  }
}

// Do initial setup of values used throughout phasing the current chromosome
void Phaser::init(int numChildren) {
  _hmm.clear();
  _hmmMarker.clear();
  _stateHash.set_empty_key(NULL);

  // alternating bits set starting with bit 0 then bit 2, ...
  uint64_t allBitsSet = ~0ul; // initially: fewer depending on numChildren
  if (numChildren < 32)
    allBitsSet &= (1 << (2*numChildren)) - 1;
  _parBits[0] = 0x5555555555555555 & allBitsSet;
  _parBits[1] = 0xAAAAAAAAAAAAAAAA & allBitsSet;
  _parBits[2] = allBitsSet;

  // For flipping bits in all children. See comment at declaration.
  // Storing these in an indexed array avoids the need to branch based on
  // the index value. That is, the necessary values are already calculated and
  // stored, but lookup in an array is more efficient than branching.
  // TODO: is there a way to just use _parBits?
  _flips[0] = 0;
  _flips[1] = _parBits[0];
  _flips[2] = _parBits[1];
  _flips[3] = _parBits[2];
  // For flipping bits in all ambiguous children. See lookupState()
  _ambigFlips[0] = _ambigFlips[3] = 0;
  _ambigFlips[1] = _ambigFlips[2] = _parBits[1];
}

// Looks up and stores the genotype values for the parents and children
// For speedy marker type detection, uses bit representation in <*GenoTypes>
// variables to indicate which genotypes the parents and children have.
void Phaser::getFamilyData(PersonBulk::par_pair parents,
			   dynarray<PersonBulk*> &children, int marker,
			   uint8_t &parentData, uint8_t &parentGenoTypes,
			   uint64_t childrenData[5], uint8_t &childGenoTypes) {
  // Get data for dad (first 2 bits) and mom (bits 3-4)
  uint8_t dadData = parents->first->getBitGeno(marker);
  uint8_t momData = parents->second->getBitGeno(marker);
  parentData = dadData + (momData << 2);

  // Which of the genotypes are present in the parents/children? There are
  // four possible genotypes, and the first four bits will have a value of 1
  // iff the corresponding genotype value is present in at least one child
  parentGenoTypes = (1 << dadData) | (1 << momData);
  childGenoTypes = 0;

  for(int g = 0; g < 4; g++)
    childrenData[g] = 0; // note: limited to 32 children (checked above)

  int numChildren = children.length();
  for(int c = 0; c < numChildren; c++) {
    uint8_t curChildData = children[c]->getBitGeno(marker);
    childrenData[ curChildData ] += 3 << (c*2);
    childrenData[4] += curChildData << (c*2);
    childGenoTypes |= 1 << curChildData; // observed genotype <curChildData>
  }
}

// Determines what type of marker this is using data for the parents if present
// or based on the observed genotype values for the the children when one or
// both parent's data are missing
int Phaser::getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes) {
  // Only valid values for parentGenoTypes are between 1 and 12 (excluding 7
  // and 11 which are caught below)
  assert(parentGenoTypes >= 1 && parentGenoTypes <= 12);
  assert(childGenoTypes >= 1 && childGenoTypes <= 15);

  // stores the genotype of the non-missing parent if present; if both are
  // missing, stores 1, the value for missing data
  int missingType = -1;

  switch (parentGenoTypes) {
    //////////////////////////////////////////////////////////////////////////
    // Parents are both homozygous
    case (1 << G_HOM0):  // both homozygous for 0
      if (childGenoTypes & ((1 << G_HOM1) | (1 << G_HET)) )
	// Mendelian error -- invalid bits set: only hom for 0 and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;
      break;
    case (1 << G_HOM1):  // both homozygous for 1
      if (childGenoTypes & ((1 << G_HOM0) | (1 << G_HET)) )
	// Mendelian error -- invalid bits set: only hom for 1 and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;
      break;
    case ((1 << G_HOM0) | (1 << G_HOM1)):  // one homozygous for 0, other 1
      if (childGenoTypes & ((1 << G_HOM0) | (1 << G_HOM1)))
	// Mendelian error -- invalid bits set: only het and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;
      break;

    //////////////////////////////////////////////////////////////////////////
    // One parent homozygous the other heterozygous
    case ((1 << G_HOM0) | (1 << G_HET)): // one homozygous for 0, other het
      if (childGenoTypes & (1 << G_HOM1))
	// Mendelian error -- child is homozygous for allele not present in
	// the homozygous parent
	return 1 << MT_ERROR;
      else
	// one parent heterozygous, other homozygous: informative for one parent
	return 1 << MT_FI_1;
      break;
    case ((1 << G_HOM1) | (1 << G_HET)): // one homozygous for 1, other het
      if (childGenoTypes & (1 << G_HOM0))
	// Mendelian error -- child is homozygous for allele not present in
	// the homozygous parent
	return 1 << MT_ERROR;
      else
	// one parent heterozygous, other homozygous: informative for one parent
	return 1 << MT_FI_1;
      break;

    //////////////////////////////////////////////////////////////////////////
    // Both parents heterozygous
    case (1 << G_HET):  // both heterozygous
      // both heterozygous for same alleles (by virtue of biallelic only data):
      // partly informative
      //
      // No Mendelian errors possible in this case
      return 1 << MT_PI;

    //////////////////////////////////////////////////////////////////////////
    // Both parents missing
    case (1 << G_MISS): // both missing
      // Further analysis done below using observed children genotypes
      missingType = G_MISS;
      break;

    //////////////////////////////////////////////////////////////////////////
    // One parent missing, other non-missing
    // Further analysis done below using observed children genotypes
    // (including checks for Mendelian errors)
    case ((1 << G_HOM0) | (1 << G_MISS)): // one homozygous for 0, other missing
      missingType = G_HOM0;
      break;
    case ((1 << G_HOM1) | (1 << G_MISS)): // one homozygous for 1, other missing
      missingType = G_HOM1;
      break;
    case ((1 << G_HET) | (1 << G_MISS)):  // one heterozygous, other missing
      missingType = G_HET;
      break;

    default:
      fprintf(stderr, "ERROR: got impossible parent genotype value %d\n",
	      parentGenoTypes);
      exit(5);
      break;
  }

  assert(missingType >= G_HOM0 && missingType <= G_HOM1);

  switch (childGenoTypes) {
    case ((1 << G_HOM0) | (1 << G_HOM1)):  // both homozygous types observed
    case ((1 << G_HOM0) | (1 << G_HOM1) | (1 << G_MISS)): // above and missing
    // both homozgyous types and heterozygous type observed
    case ((1 << G_HOM0) | (1 << G_HOM1) | (1 << G_HET)):
    case 15: // all types observed (= to just above as well as missing)
      // if children have both homozygous types, parents must both be
      // heterozygous: partly informative
      //
      // Ensure this accords with what we know about the parent genotypes:
      switch(missingType) {
	case G_MISS: // both missing data -- no information to add -- consistent
	case G_HET:  // known parent is heterozgyous -- consistent
	  return 1 << MT_PI;
	case G_HOM0: // one parent homozygous for 0
	case G_HOM1: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for
	  // allele that is not present in the parent. Since both homozygous
	  // types are present, this is violated.
	  return 1 << MT_ERROR;
      }
      break;
    case (1 << G_HOM0):                    // only homozygous for 0 observed
    case ((1 << G_HOM0) | (1 << G_MISS)):  // above and missing type observed
      // could in principle be anything; see if there is information on one
      // of the parents that constrains the type:
      switch (missingType) {
	case G_MISS: // both missing data -- could be any type:
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  return (1 << MT_UN) | (1 << MT_FI_1);
	case G_HOM1: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
      }
      break;
    case (1 << G_HOM1):                   // only homozygous for 1 observed
    case ((1 << G_HOM1) | (1 << G_MISS)): // above and missing type observed
      // could in principle be anything; see if there is information on one
      // of the parents that constrains the type:
      switch (missingType) {
	case G_MISS: // both missing data -- could be any type:
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case G_HOM1: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  return (1 << MT_UN) | (1 << MT_FI_1);
      }
      break;
    case ((1 << G_HOM0) | (1 << G_HET)): // homozygous 0 and het types observed
    case ((1 << G_HOM0) | (1 << G_HET) | (1 << G_MISS)): // above and missing
      // could be informative for one parent or partly informative; see if there
      // is information on one of the parents that constrains the type:
      switch (missingType) {
	case G_MISS: // both missing data -- no information to add
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  return 1 << MT_FI_1;
	case G_HOM1: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
      }
      break;
    case ((1 << G_HOM1) | (1 << G_HET)): // homozygous 1 and het types observed
    case ((1 << G_HOM1) | (1 << G_HET) | (1 << G_MISS)): // above and missing
      // could be informative for one parent or partly informative; see if there
      // is information on one of the parents that constrains the type:
      switch (missingType) {
	case G_MISS: // both missing data -- no information to add
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case G_HOM1: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  return 1 << MT_FI_1;
      }
      break;
    case (1 << G_HET):                   // only heterozygous types observed
    case ((1 << G_HET) | (1 << G_MISS)): // above and missing type observed
      // This is ambiguous unless we have a parent that is homozygous
      switch (missingType) {
	case G_MISS: // both missing data -- no information to add
	  return 1 << MT_AMBIG;
	case G_HET: // one heterozygous parent
	  // no hints about phasing with all observed samples heterozygous
	  return 1 << MT_AMBIG;
	case G_HOM0: // one parent homozygous for 0
	case G_HOM1: // one parent homozygous for 1
	  // other parent must be heterozygous, so fully informative for that
	  // parent
	  return 1 << MT_FI_1;
      }
      break;
    case (1 << G_MISS): // all children missing data
      // no real information
      return 1 << MT_AMBIG;
    default:
      fprintf(stderr, "ERROR: got impossible parent genotype value %d\n",
	      parentGenoTypes);
      exit(5);
      break;
  }

  return -1; // shouldn't happen
}

// TODO: remove this
//    for(int t = MT_FI_1; t < MT_N_TYPES;  t++) {
//      if (mt & (1 << t)) {
//	switch(t) {
//	  case MT_FI_1:
//	    fprintf(stderr, " FI");
//	    break;
//	  case MT_PI:
//	    fprintf(stderr, " PI");
//	    break;
//	  case MT_UN:
//	    fprintf(stderr, " UN");
//	    break;
//	  case MT_AMBIG:
//	    fprintf(stderr, " AM");
//	    break;
//	  case MT_ERROR:
//	    fprintf(stderr, " ER");
//	    break;
//	}
//      }
//    }
//    fprintf(stderr, "\n");

// Given the marker types <markerTypes> that are consistent with the family
// data, generates partial states indicating, where known, which allele
// each parent transmitted to the children.
void Phaser::makePartialStates(dynarray<State> &partialStates,
			       int markerTypes, uint8_t parentData,
			       const uint64_t childrenData[5]) {
  // Fully informative for one parent:
  if (markerTypes & (1 << MT_FI_1)) {
    makePartialFI1States(partialStates, parentData, childrenData);
  }

  // Partly informative
  if (markerTypes & (1 << MT_PI)) {
    makePartialPIStates(partialStates, parentData, childrenData);
  }

  // Uninformative
  if (markerTypes & (1 << MT_UN)) {
    // TODO: in fact, a state with this genotype for the parents will always
    // have 0 recombinations observed. Perhaps we need a way of evaluating the
    // evidence for or against this case? In some sense should assume the
    // marker is not this and use surrounding markers to revert to this case
    // only when necessary
  }
}

// Helper for makePartialStates(): applicable to fully informative for one
// parent markers (or the states that correspond to this possibility when the
// parent's genotypes aren't fully known)
void Phaser::makePartialFI1States(dynarray<State> &partialStates,
				  uint8_t parentData,
				  const uint64_t childrenData[5]) {
  uint8_t startPar, endPar;

  // If both parents are missing, we'll make states corresponding to each
  // being heterozygous (other homozygous), consistent with the ambiguity
  if (parentData == (G_MISS << 2) + (G_MISS)) {
    startPar = 0;
    endPar = 1;
  }
  // parent 0 het: either we have data for parent 0 as a het, or we have data
  // for parent 1 as homozygous (or both)
  else if ((parentData & 3) == G_HET || ((parentData >> 2) & 3) == G_HOM0 ||
	   ((parentData >> 2) & 3) == G_HOM1) {
    startPar = endPar = 0;
  }
  // parent 1 het:
  else {
    assert(((parentData >> 2) & 3) == G_HET || (parentData & 3) == G_HOM0 ||
	   (parentData & 3) == G_HOM1);
    startPar = endPar = 1;
  }

  for(uint8_t hetPar = startPar; hetPar <= endPar; hetPar++) {
    uint8_t homPar = hetPar ^ 1;

    // Add space for a new (partial) state. The use of addEmpty() allows us to
    // get a reference to the state stored here and enables direct writing to
    // the stored state (no need to copy values in)
    int newIndex = partialStates.length();
    partialStates.addEmpty();
    State &newState = partialStates[newIndex];

    // Initially assumes the phase is assigned such that allele 0 is on
    // haplotype 0 and allele 1 is on haplotype 1. More specifically, the
    // heterozygous children are assumed to receive haplotype 1 and the
    // homozygous children haplotype 2
    newState.iv = _parBits[hetPar] & childrenData[ G_HET ];
    // No ambiguous bits (though a missing data child could have ambiguous bits
    // propagated from a previous marker, but that happens during full state
    // construction
    newState.ambig = 0;
    // No information about transmission from homozygous parent and for
    // children that are missing data:
    newState.unassigned = _parBits[homPar] | childrenData[ G_MISS ];
    newState.hetParent = hetPar;
    // Won't assign <parentPhase> field as all partial states have a value of
    // 0 here and it's never used
    //newState.parentPhase = 0;

    assert((newState.iv & newState.unassigned) == 0ul);
  }
}

// Helper for makePartialStates(): applicable to partly informative markers (or
// the states that correspond to this possibility when the parent's genotypes
// aren't fully known)
void Phaser::makePartialPIStates(dynarray<State> &partialStates,
				 uint8_t parentData,
				 const uint64_t childrenData[5]) {
  // Add space for a new (partial) state. The use of addEmpty() allows us to
  // get a reference to the state stored here and enables direct writing to
  // the stored state (no need to copy values in)
  int newIndex = partialStates.length();
  partialStates.addEmpty();
  State &newState = partialStates[newIndex];

  // Initially assumes the phase is assigned such that allele 0 is on
  // haplotype 0 and allele 1 is on haplotype 1. In particular, the genotypes
  // of the children (excluding the missing data children) already encode the
  // inheritance vector assuming this phase assignment in the parents.
  // This also assumes that all heterozygous children received haplotype 0
  // from parent 0 and haplotype 1 from parent 1, but this will be changed
  // when making full states.
  newState.iv = childrenData[ 4 ] & ( ~childrenData[ G_MISS ] );
  // Potentially ambiguous bits are those where the child is heterozygous:
  newState.ambig = childrenData[ G_HET ];
  // Have full transmission information except for children with missing data:
  newState.unassigned = childrenData[ G_MISS ];
  newState.hetParent = 2; // both parents heterozygous
  // Won't assign <parentPhase> field as all partial states have a value of
  // 0 here and it's never used
  //newState.parentPhase = 0;

  assert((newState.iv & newState.unassigned) == 0ul);
}

// Using the states at the previous marker where present and <partialStates>,
// generates all full states. Stores these as the last value in <_hmm> and
// stores the marker index <marker> that these states apply to in <_hmmMarker>.
void Phaser::makeFullStates(const dynarray<State> &partialStates, int marker,
			    const uint64_t childrenData[5]) {
  // The new entry in <_hmm> corresponds to <marker>:
  _hmmMarker.append(marker);

  if (_hmm.length() == 0) {
    // Partial states are equivalent to full states when there are no previous
    // states
    _hmm.addEmpty();
    int len = partialStates.length();
    for(int i = 0; i < len; i++) {
      State *newState = new State(partialStates[i]);
      _hmm[0].append(newState);
    }
    return;
  }

  int prevHMMIndex = _hmm.length() - 1;
  _hmm.addEmpty();
  dynarray<State*> &prevStates = _hmm[prevHMMIndex];

  // TODO: memory leak: the keys are heap allocated
  _stateHash.clear_no_resize();

  uint16_t numPrev = prevStates.length();
  for(uint16_t prevIndex = 0; prevIndex < numPrev; prevIndex++) {
    const State *prevState = prevStates[prevIndex];

    int numPartial = partialStates.length();
    for(int curInd = 0; curInd < numPartial; curInd++) {
      const State &curPartial = partialStates[curInd];

      // (1) For each (current) partial state, determine the full state that
      // <prevState> maps to. The <iv> and <unassigned> values are determined
      // based the <iv> and <unassigned> fields in <prevState> and <curPartial>.
      // <curPartial.unassigned> has bits set to 1 for the parent that is
      // homozygous/uninformative or both bits set to 1 for a child that is
      // missing data. There is no information in <curPartial> about the
      // haplotype transmissions for these cases and so we propagate the
      // inheritance vector values from the previous inheritance vector.
      uint64_t fullIV = curPartial.iv | (prevState->iv & curPartial.unassigned);
      uint64_t fullUnassigned = prevState->unassigned & curPartial.unassigned;
      uint64_t fullAmbig; // assigned below
      uint64_t hetChildBitFlip = 0;  // only applicable to MT_PI states
      uint64_t newAmbigChildren = 0; // only applicable to MT_PI states

      // Which haplotypes recombined? Note that the current <iv> value assumes
      // a certain phase for the parents and may have more recombinations than
      // another possibility. We consider all the possibilities below.
      // The following gives us enough information to determine how to handle
      // the heterozygous children at MT_PI markers.
      // Note: because (fullIV & curPartial.unassigned) ==
      //                                (prevState->iv & curPartial.unassigned),
      // we know that no recombinations occur at those bits.
      // We must avoid counting recombinations for bits that were unassigned
      // in <prevState>, however.
      uint64_t recombs = (prevState->iv ^ fullIV) & (~prevState->unassigned);

      // (2) if curPartial is partly informative, for heterozygous children:
      //     (a) When two recombinations occur relative to the previous marker
      //         and the child was unambiguous previously, flip both
      //         corresponding bits in the inheritance vector. HAPI avoids a
      //         state space explosion for heterozygous children at MT_PI
      //         markers by a combination of only modeling states with 0
      //         recombinations instead of 2 when possible and (b),(c) below.
      //     (b) When one recombination occurs relative to the previous marker,
      //         leave the iv value as assigned and set the child as ambiguous.
      //     Note: ambiguous bits are only set for:
      //      (i)  children that are newly ambiguous: i.e, those that are
      //           heterozygous and exhibit one recombination relative to the
      //           previous marker
      //      (ii) or children that are heterozygous and were ambiguous at the
      //           previous marker (this status propagates forward until an
      //           unambiguous marker, including the child being homozygous at
      //           MT_PI)
      if (curPartial.hetParent == 2) { // MT_PI state
	// determines which previously unambiguous children that are
	// heterozygous at this position also have two recombinations relative
	// to <prevState>
	uint64_t recombsPrevUnambig = recombs & (~prevState->ambig);
	calcHetChildPIVals(recombsPrevUnambig, childrenData, hetChildBitFlip,
			   newAmbigChildren);
	fullIV ^= hetChildBitFlip;
	fullAmbig = newAmbigChildren |
	      ((childrenData[G_MISS] | childrenData[G_HET]) & prevState->ambig);
	// Remove the recombinations from both parents from <recombs>
	recombs ^= hetChildBitFlip;
      }
      else {
	// Fully informative marker: only ambiguous bits are those where a
	// child was ambiguous in the previous state and missing data at
	// this marker
	fullAmbig = childrenData[G_MISS] & prevState->ambig;
      }

      // (3) Look up or create a full state with equivalent <iv> and <ambig>
      // values to <fullIV> and <fullAmbig>, and determine if <prevState>
      // yields fewer recombinations for these states. If so, update the
      // necessary values in the state.
      // Also examines an alternate phase type which may or may not map to
      // an equivalent state. If not, looks up that value, if so, compares the
      // recombinations for the two possibilities separately.
      // TODO: optimize
      uint8_t parPhaseFlip = (curPartial.hetParent < 2) ? 1 : 3;
      updateStates(fullIV, fullAmbig, fullUnassigned, recombs,
		   curPartial.hetParent, /*initParPhase=default phase=*/ 0,
		   parPhaseFlip, prevIndex, prevState->minRecomb,
		   childrenData);

      // TODO! call updateStates() a second time when hetParent == 2
      // That is, modify <fullIV> and <fullAmbig> for 1 parent away
    }
  }
  // TODO: explore the optimization in which we find the state(s) with minimum
  // recombinations and compute the difference between all other states and
  // this state/s. If those states have more recombinations than the numbers of
  // differences between them, remove those states.
}

// For MT_PI states, determines:
//   <hetChildBitFlip>  : which <iv> values need to be flipped in order to
//                        avoid having a recombination from both parents
//   <newAmbigChildren> : which <iv> values are newly ambiguous in this state
//                        due to being heterozgyous and having one recombination
//                        relative to <prevState>
// Further details are in the comments of makeFullStates() before this
// method gets called.
void Phaser::calcHetChildPIVals(const uint64_t recombsPrevUnambig,
				const uint64_t childrenData[5],
				uint64_t &hetChildBitFlip,
				uint64_t &newAmbigChildren) {
  uint64_t recombsUnambigHetChild = recombsPrevUnambig & childrenData[G_HET];
  uint64_t parRecomb[2];
  for(int p = 0; p < 2; p++)
    parRecomb[p] = recombsUnambigHetChild & _parBits[p];

  /////////////////////////////////////////////////////////////////////////////
  // (1) flip as needed for children that have two recombinations

  // For every two bits, this sum will only be 2 decimal = 10 binary iff the
  // parent 0 bit and the parent 1 bit are 1:
  uint64_t sumRecomb = parRecomb[0] + (parRecomb[1] >> 1);
  // Extract the higher order bit of the two bits for each child. The following
  // will give a 01 binary value in the given child's two bits if both
  // recombinations occurred
  uint64_t bothRecomb = (sumRecomb & parRecomb[1]) >> 1;

  // Can simply multiply by 3 decimal = 11 binary and the corresponding two
  // bits will be set for any child with both haplotypes recombined.
  // We only flip the iv values for children that were unambiguous at the
  // previous marker:
  hetChildBitFlip = bothRecomb * 3;

  /////////////////////////////////////////////////////////////////////////////
  // (2) identify children that are newly ambiguous

  // If <sumRecomb> has a value of 01 for a given child, there's only 1
  // recombination, which is ambiguous for heterozygous children (the only
  // recombinations we're examining now), so:
  uint64_t oneRecomb = sumRecomb & parRecomb[0];

  newAmbigChildren = oneRecomb * 3;
}

// Look up or create a full state with equivalent <iv> and <ambig> values to
// <fullIV> and <fullAmbig>, and determine if <prevState> yields fewer
// recombinations for these states. If so, update the necessary values in the
// state.
// Also examines an alternate phase type which may or may not map to an
// equivalent state. If not, looks up that value, if so, compares the
// recombinations for the two possibilities separately.
void Phaser::updateStates(uint64_t fullIV, uint64_t fullAmbig,
			  uint64_t fullUnassigned, uint64_t recombs,
			  uint8_t hetParent, uint8_t initParPhase,
			  uint8_t parPhaseFlip, uint16_t prevIndex,
			  uint16_t prevMinRecomb,
			  const uint64_t childrenData[5]) {
  // How many iterations of the loop? See various comments below.
  int numIter = 2;

  // How many recombinations for the initial phase assignment?
  size_t numRecombs = popcount(recombs);
  uint8_t curParPhase = initParPhase;

  // First, decide how many iterations of the loop below. If the opposite
  // phase assignment yields an equivalent state, we only need to loop once.
  // Only way the opposite phase assignment is equivalent is if there is
  // no missing data (we propagate <iv> values from the previous marker for
  // these individuals so changing the parent's phase at the current marker
  // won't change anything.
  // Also, for PI type states, must have all heterozygous markers be ambiguous.
  // Otherwise, since flipping will necessitate fixing the unambiguous children
  // to match the previous state (to avoid recombinations from both parents),
  // the resulting state will not be equivalent.
  if (childrenData[G_MISS] == 0 && (hetParent < 2 ||
					    childrenData[G_HET] == fullAmbig)) {
    numIter = 1;

    // Determine which phase assignment has minimal recombinations before the
    // loop below.

    // Don't flip ambiguous bits: these should be the same as for the original
    // option.
    recombs ^= _parBits[ hetParent ] & (~fullAmbig);
    size_t curCount = popcount(recombs);
    if (curCount < numRecombs) {
      numRecombs = curCount;
      // TODO: move this comment?
      // Note: when curPartial.hetParent == 2, we actually flip both alleles
      // This is first because _parBits[2] has all bits set, but more
      // importantly because flipping only one bit will in generalhave
      // different ambiguous status whereas flipping both bits yields a IV
      // value with the same ambiguous status.
      curParPhase ^= parPhaseFlip;
    }
    else if (curCount == numRecombs) {
      // TODO: ambiguous locally. Indicate in State value.
    }
  }


  // Examine the 1 or 2 needed phase assignments
  for(int i = 0; i < numIter; i++) {
    // First look up or create a full state with equivalent <iv> and <ambig>
    // values to <fullIV> and <fullAmbig>
    // TODO: what about using <unassigned> to distinguish values? Only
    // relevant near the beginning of a chromosome, but complex given that
    // some markers might be either FI or PI: that value effectively
    // distinguishes between them.
    State *theState = lookupState(fullIV, fullAmbig);
    // TODO: add a check for prevState->minRecomb > theState.minRecomb?

    int totalRecombs = prevMinRecomb + numRecombs;
    if (totalRecombs < theState->minRecomb) {
      theState->iv = fullIV;
      theState->ambig = fullAmbig;
      theState->unassigned = fullUnassigned;
      theState->prevState = prevIndex;
      theState->minRecomb = totalRecombs;
      theState->hetParent = hetParent;
      theState->parentPhase = curParPhase;
    }
    else if (totalRecombs == theState->minRecomb) {
      // TODO: ambiguous for previous state; need way to represent
    }

    if (numIter == 2) {
      // Update the various values as needed for the inverted phase in the next
      // iteration.

      // Children that have missing data should not be flipped
      // Also, for both parent het markers, heterozygous children's bits do
      // not get flipped: they're either unambigous and constrained by
      // <prevState->iv> or their ambiguous and should remain the same to stick
      // with the convention used in <_stateHash>.
      // Note: (hetParent >> 1) == 1 iff hetParent == 2. It is 0 for the other
      //       values.
      uint64_t noFlipBits = childrenData[G_MISS] |
				      (childrenData[G_HET] * (hetParent >> 1));
      uint64_t flipVal = _parBits[ hetParent ] & (~noFlipBits);
      fullIV ^= flipVal;
      recombs ^= flipVal;
      numRecombs = popcount(recombs);
      curParPhase ^= parPhaseFlip;
      // <fullAmbig> and <fullUnassigned> don't change
    }
  }
}

// Given <iv> and <ambig>, first convert <iv> to the equivalent canonical value
// (see comment just inside method) and then do a hash table lookup for the
// State. If it does not exist, create it. Either way, return it to the caller.
State * Phaser::lookupState(const uint64_t iv, const uint64_t ambig) {
  // As inheritance vectors have four equivalent values, we have a fixed key
  // defined via the following convention:
  // (1) The lowest order two bits for which <iv> is ambiguous, the value
  //     should be 00.
  // (2) The <iv> values for ambiguous bits are either the default of 10 binary
  //     or 00 if satisfying (1) involves inverting the inheritance value of
  //     only one parent. (Note that for ambiguous bits, <iv> values of 10 == 01
  //     and 00 == 11, so it suffices to only consider the case of one parent's
  //     bits flipped).

  uint64_t unambig = ~ambig;

  // Get the fixed key value
  // First determine the <iv> value for the lowest order unambiguous child:
  int lowOrderChildBit = ffsll(unambig) - 1;
  // Markers that have all heterozygous children are ambiguous and not
  // considered at this stage, so there must be at least one unambiguous child:
  assert(lowOrderChildBit >= 0);

  // The genotype of the child tells us what bits we need to flip: the exact
  // bits that are assigned 1 need to be flipped in all unambiguous children.
  uint8_t flipType = (iv >> lowOrderChildBit) & 3;

  // Conveniently, we've got _flips indexed by the 4 possible flip types with
  // the values to flip assigned in each child:
  uint64_t lookupIV = iv ^ (_flips[flipType] & unambig);
  // And we've done something analogous for ambiguous bits; must flip them too:
  lookupIV ^= _ambigFlips[flipType] & ambig;

  // Do the lookup
  iv_ambig_real theKey(lookupIV, ambig);
  state_ht_iter it = _stateHash.find( &theKey );
  if (it == _stateHash.end()) {
    // need to create state
    State *newState = new State;
    newState->minRecomb = UINT16_MAX;
    iv_ambig newStateKey = new iv_ambig_real(lookupIV, ambig);
    _stateHash[ newStateKey ] = newState;
    int curHMMIndex = _hmm.length() - 1;
    _hmm[curHMMIndex].append(newState);
    // caller will assign necessary values
    return newState;
  }
  else {
    return it->second;
  }
}

void printGeno(FILE *out, int marker, uint8_t bitGeno) {
  const char *alleles = Marker::getMarker(marker)->getAlleleStr();
  switch(bitGeno) {
    case G_HOM0:
      fprintf(out, "%c/%c", alleles[0], alleles[0]);
      break;
    case G_MISS:
      fprintf(out, "0/0");
      break;
    case G_HET:
      fprintf(out, "%c/%c", alleles[0], alleles[2]);
      break;
    case G_HOM1:
      fprintf(out, "%c/%c", alleles[2], alleles[2]);
      break;
  }
}
