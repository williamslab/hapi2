// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <bitset>
#include <genetio/util.h>
#include "phaser.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray< dynarray<State*> > Phaser::_hmm;
dynarray<int>      Phaser::_hmmMarker;
dynarray<uint16_t> Phaser::_minStates;
Phaser::state_ht Phaser::_stateHash;
uint64_t Phaser::_parBits[3];
uint64_t Phaser::_flips[4];
uint64_t Phaser::_ambigFlips[4];


void Phaser::run(NuclearFamily *theFam, int chrIdx) {
  // Get first and last marker numbers for the chromosome to be phased:
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);

  int numChildren = theFam->numChildren();
  if (numChildren > 32) {
    fflush(stdout);
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: cannot phase more than 32 children in a family.\n");
    fprintf(stderr, "       changing to >64 bits for inheritance vectors would fix this\n");
    exit(9);
  }

  // TODO: I think the clear() calls are not necessary
  _hmm.clear();
  _hmmMarker.clear();

  parBitsInit(numChildren);

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
    getFamilyData(theFam, m, parentData, parentGenoTypes, childrenData,
		  childGenoTypes);

    ///////////////////////////////////////////////////////////////////////////
    // Step 1: Determine marker type and check for Mendelian errors
    int mt = getMarkerType(parentGenoTypes, childGenoTypes);
    assert(mt > 0 && (mt & ( ~((1 << MT_N_TYPES) -1) )) == 0);

    if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
      // should only be one of the above:
      assert(mt == (1 << MT_ERROR) || mt == (1 << MT_AMBIG));
      switch(mt) {
	case 1 << MT_ERROR:
	  theFam->setStatus(/*marker=*/ m, PHASE_ERROR, parentData,
			    childrenData[4]);
	  break;
	case 1 << MT_AMBIG:
	  theFam->setStatus(/*marker=*/ m, PHASE_AMBIG, parentData,
			    childrenData[4]);
	  break;
      }
      continue;
    }

    if (mt == (1 << MT_UN)) {
      // TODO: this should only happen when we have data for both parents. In
      // that case, when the children are heterozygous, we can phase them and
      // should do so.
      theFam->setStatus(/*marker=*/ m, PHASE_UNINFORM, parentData,
			childrenData[4]);
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
    makeFullStates(partialStates, /*marker=*/ m, childrenData);

    // Clean up for next marker
    partialStates.clear();
    for (state_ht_iter it = _stateHash.begin(); it != _stateHash.end(); it++) {
      delete it->first;
    }
    _stateHash.clear_no_resize();
  }

  /////////////////////////////////////////////////////////////////////////////
  // Step 4: HMM analysis finished. Back trace and assign phase!
  backtrace(theFam);
}

// Do initial setup of values used throughout phasing the current chromosome
void Phaser::parBitsInit(int numChildren) {
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
void Phaser::getFamilyData(NuclearFamily *theFam, int marker,
			   uint8_t &parentData, uint8_t &parentGenoTypes,
			   uint64_t childrenData[5], uint8_t &childGenoTypes) {
  NuclearFamily::par_pair parents = theFam->_parents;
  dynarray<PersonBulk*> &children = theFam->_children;

  // Get data for dad (first 2 bits) and mom (bits 3-4)
  uint8_t dadData = parents->first->getBitGeno(marker);
  uint8_t momData = parents->second->getBitGeno(marker);
  parentData = dadData + (momData << 2);

  // Which of the genotypes are present in the parents/children? There are
  // four possible genotypes, and the first four bits will have a value of 1
  // iff the corresponding genotype value is present in at least one child
  parentGenoTypes = (1 << dadData) | (1 << momData);
  childGenoTypes = 0;

  for(int g = 0; g < 5; g++)
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
      // No Mendelian errors possible in this case, but if all children are
      // het/missing (or a combination thereof), the marker is ambiguous
      switch(childGenoTypes) {
	case 1 << G_HET:
	case 1 << G_MISS:
	case (1 << G_HET) | (1 << G_MISS):
	  return 1 << MT_AMBIG;
	default:
	  return 1 << MT_PI;
      }

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
  // TODO: getMarkerType() can infer what this genotype must be for missing
  // data cases. Get it to do so
  uint8_t homParGeno; // homozygous parent genotype

  // If both parents are missing, we'll make states corresponding to each
  // being heterozygous (other homozygous), consistent with the ambiguity
  if (parentData == (G_MISS << 2) + (G_MISS)) {
    startPar = 0;
    endPar = 1;
    homParGeno = G_MISS;
  }
  // parent 0 het: either we have data for parent 0 as a het, or we have data
  // for parent 1 as homozygous (or both)
  else if ((parentData & 3) == G_HET || (parentData >> 2) == G_HOM0 ||
	   (parentData >> 2) == G_HOM1) {
    startPar = endPar = 0;
    homParGeno = parentData >> 2;
  }
  // parent 1 het:
  else {
    assert((parentData >> 2) == G_HET || (parentData & 3) == G_HOM0 ||
	   (parentData & 3) == G_HOM1);
    startPar = endPar = 1;
    homParGeno = parentData & 3;
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
    newState.homParentGeno = homParGeno;
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
  // Won't assign <homParentGeno> as it's meaningless when both parents are het
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
    assert(len <= 3); // can only have 2 forms of FI partial states and one PI
    for(int i = 0; i < len; i++) {
      if (i == 1 && partialStates[1].hetParent == 1 &&
						partialStates[0].hetParent == 0)
	// TODO: not sure we want this: perhaps if we know there's no parent
	// data anywhere on the length of the chromosome
	// For the very first marker, if we're not sure which parent is
	// heterozygous, arbitrarily pick parent 0 as such. Otherwise, since the
	// two states are equivalent but with opposite parent labels, there will
	// be two equal paths through the HMM.
	continue;
      State *newState = new State(partialStates[i]);
      newState->minRecomb = 0;
      newState->parentPhase = 0;
      _hmm[0].append(newState);
    }
    return;
  }

  int prevHMMIndex = _hmm.length() - 1;
  _hmm.addEmpty();
  dynarray<State*> &prevStates = _hmm[prevHMMIndex];

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
      // only applicable to MT_PI states:
      // Relative to the default phase for the parents, which bits recombined
      // in the heterozygous children that were not ambiguous at the previous
      // marker? Recombination bits may be 00, 01, 10, 11 binary:
      uint64_t unambigHetRecombs[4];
      uint64_t propagateAmbig;

      // Which haplotypes recombined? Note that the current <iv> value assumes
      // a certain phase for the parents and may have more recombinations than
      // another possibility. We consider other possibilities below.
      // The following gives us enough information to determine how to handle
      // the heterozygous children at MT_PI markers.
      // Note: because (fullIV & curPartial.unassigned) ==
      //                                (prevState->iv & curPartial.unassigned),
      // we know that no recombinations occur at <curPartial.unassigned> bits.
      // We must avoid counting recombinations for bits that were unassigned
      // in <prevState>, however: differences between <prevState->iv> and
      // <fullIV> are not meaningful since those bits weren't assigned in
      // <prevState>.
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
      // The following gives 1 for <hetParent> == 2, 0 for <hetParent> == 0,1
      uint8_t isPI = curPartial.hetParent >> 1;
      if (isPI) { // MT_PI state
	// determines which previously unambiguous children that are
	// heterozygous at this position also have two recombinations relative
	// to <prevState>
	uint64_t unambigHets = childrenData[G_HET] & (~prevState->ambig);
	calcHetChildPIVals(recombs, unambigHets, unambigHetRecombs);
	// Remove the recombinations from both parents from <recombs>
	fullIV ^= unambigHetRecombs[3];
	recombs ^= unambigHetRecombs[3];
	propagateAmbig =
		(childrenData[G_MISS] | childrenData[G_HET]) & prevState->ambig;
	fullAmbig = unambigHetRecombs[1] | unambigHetRecombs[2] |propagateAmbig;
      }
      else {
	// Fully informative marker: only ambiguous bits are those where a
	// child was ambiguous in the previous state and missing data at
	// this marker
	fullAmbig = childrenData[G_MISS] & prevState->ambig;
      }

      // (3) Look up or create a full state with equivalent <iv> and <ambig>
      // values to <fullIV> and <fullAmbig>, and determine if <prevState>
      // yields fewer recombinations for these states than the currently stored
      // previous state (if any). If so, update the necessary values in the
      // state.
      // Also examines an alternate phase type which may or may not map to
      // an equivalent state. If not, looks up that value, if so, compares the
      // recombinations for the two possibilities separately.
      //
      // What type of phase is the alternative? For MT_PI states, the
      // alternative, which must have the same ambiguous bits, is different for
      // both parents, so the type is 3 decimal == 11 binary (as opposed to the
      // default of 0).
      // See just below for what this code does:
      uint8_t parPhaseFlip = isPI * 3 + (1 - isPI);
      // Above equivalent to the line below but has no branching
      //uint8_t parPhaseFlip = (isPI) ? 3 : 1;
      updateStates(fullIV, fullAmbig, fullUnassigned, recombs,
		   curPartial.hetParent, curPartial.homParentGeno,
		   /*initParPhase=default phase=*/ 0, parPhaseFlip, prevIndex,
		   prevState->minRecomb, childrenData);

      // For MT_PI states, have 1 or 2 more states to examine:
      if (isPI) {
	// Above considered default and state with both parents' phase inverted.
	// Now consider the two states in which each parent alone has inverted
	// phase.

	///////////////////////////////////////////////////////////////////////
	// Generate the <fullIV> and <fullAmbig> values for parent 0 flipped:
	// (Note: <fullUnassigned> does not change)
	flipPIVals(fullIV, fullAmbig, childrenData, propagateAmbig,
		   unambigHetRecombs);
	recombs = (prevState->iv ^ fullIV) & (~prevState->unassigned);

	///////////////////////////////////////////////////////////////////////
	// Now ready to look up or create full states with the <fullIV> and
	// <fullAmbig> values, etc.
	updateStates(fullIV, fullAmbig, fullUnassigned, recombs,
		     /*curPartial.hetParent=*/2, /*homParentGeno=*/G_MISS,
		     /*initParPhase=parent 0 flip=*/ 1, /*parPhaseFlip=*/ 3,
		     prevIndex, prevState->minRecomb, childrenData);
      }
    }
  }
  // TODO: explore the optimization in which we find the state(s) with minimum
  // recombinations and compute the difference between all other states and
  // this state/s. If those states have more recombinations than the numbers of
  // differences between them, remove those states.
}

// For MT_PI states, determines:
//   <bothParRecomb>  : which heterozygous children have recombinations from
//			both parents. These bits will be flipped.
//   <oneParRecomb[2]> : which heterozygous children have a recombination only
//			 from parent 0 or parent 1 (the index of the array)
// Further details are in the comments of makeFullStates() before this
// method gets called.
void Phaser::calcHetChildPIVals(const uint64_t recombs,
				const uint64_t unambigHets,
				uint64_t unambigHetRecombs[4]) {
  uint64_t recombsUnambigHetChild = recombs & unambigHets;
  uint64_t parRecomb[2];
  for(int p = 0; p < 2; p++)
    parRecomb[p] = recombsUnambigHetChild & _parBits[p];

  // set both bits for children that inherit a recombination from the given
  // parents
  parRecomb[0] |= parRecomb[0] << 1;
  parRecomb[1] |= parRecomb[1] >> 1;

  // unambiguous het children that exhibit recombinations from both parents:
  unambigHetRecombs[3] = parRecomb[0] & parRecomb[1];

  // unambiguous het children with a recombination from only one parent:
  for(int r = 1; r <= 2; r++)
    unambigHetRecombs[r] = parRecomb[r-1] ^ unambigHetRecombs[3];

  // unambiguous het children with no recombinations:
  unambigHetRecombs[0] = unambigHets - (parRecomb[0] | parRecomb[1]);
}

// For MT_PI states, flipping a single parent's phase is complex for
// heterozygous children. Here we generate the <fullIV>, <fullAmbig>, and
// <recomb> values wherein parent 0's phase is flipped relative to the default.
// We base this on the values passed in as arguments, including the <fullIV>
// and <fullAmbig> values generated for the default phase.
void Phaser::flipPIVals(uint64_t &fullIV, uint64_t &fullAmbig,
			const uint64_t childrenData[5], uint64_t propagateAmbig,
			const uint64_t unambigHetRecombs[4]) {
  ///////////////////////////////////////////////////////////////////////
  // Get the flipped <fullIV> value:

  // <iv> value is the same as the default for all missing data children:
  uint64_t newFullIV = fullIV & childrenData[G_MISS];
  // flip the transmitted haplotype for parent 0 for all homozygous
  // children
  newFullIV |= (fullIV ^ _parBits[0]) &
				  (childrenData[G_HOM0] | childrenData[G_HOM1]);
  // For heterozygous children (only remaining children to be defined), the
  // value in <newFullIV> is 00. This is the value we want for all ambiguous
  // children -- it is the canonical value for looking up in <_stateHash> for
  // this phase type (and for parent 1 flipped phase).
  // Thus, we only need to modify heterozygous children's IV values if they are
  // non-ambiguous
  // Default <iv> value was 10 binary. Any children that recombined on parent 1
  // under that setup now have no recombinations and have the correct <iv>
  // value. Children that recombined on parent 0 now have 2 recombinations and
  // need to be flipped to 11 binary. Recombination on parent 0 is 01 binary:
  newFullIV ^= unambigHetRecombs[1];
  // Done!
  fullIV = newFullIV;

  // Children that either recombined on both parents or did not recombine at
  // all previously are ambiguous when the phase one parent is flipped, so:
  uint64_t newAmbigChildren = unambigHetRecombs[0] | unambigHetRecombs[3];

  fullAmbig = newAmbigChildren | propagateAmbig;
}

// Look up or create a full state with equivalent <iv> and <ambig> values to
// <fullIV> and <fullAmbig>, and determine if <prevState> yields fewer
// recombinations for these states than the currently stored previous state (if
// any). If so, update the necessary values in the state.
// Also examines an alternate phase type which may or may not map to an
// equivalent state. If they are equivalent, compares the two possibilities
// to determine which has fewer recombinations. If they are not equivalent,
// looks up or creates the other full state and does the same analysis relative
// to the optimality of <prevState>.
void Phaser::updateStates(uint64_t fullIV, uint64_t fullAmbig,
			  uint64_t fullUnassigned, uint64_t recombs,
			  uint8_t hetParent, uint8_t homParentGeno,
			  uint8_t initParPhase, uint8_t parPhaseFlip,
			  uint16_t prevIndex, uint16_t prevMinRecomb,
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
  // missing data children so changing the parent's phase at the current marker
  // will modify some children's <iv> values but not the missing data ones).
  // Also, for PI type states, must have all heterozygous markers be ambiguous.
  // Otherwise, since flipping will necessitate fixing the unambiguous children
  // to match the previous state (to avoid recombinations from both parents),
  // the resulting state will not be equivalent (akin to missing data children).
  if (childrenData[G_MISS] == 0 && (hetParent < 2 ||
					    childrenData[G_HET] == fullAmbig)) {
    numIter = 1;

    // Determine which phase assignment has minimal recombinations before the
    // loop below.

    // Don't flip ambiguous bits: these should be the same as for the original
    // option to stay with the convention used in the hash (see lookupState()).
    recombs ^= _parBits[ hetParent ] & (~fullAmbig);
    size_t curCount = popcount(recombs);
    if (curCount < numRecombs) {
      numRecombs = curCount;
      curParPhase ^= parPhaseFlip;
      fullIV ^= _parBits[ hetParent ] & (~fullAmbig);
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
      theState->homParentGeno = homParentGeno;
      theState->parentPhase = curParPhase;
    }
    else if (totalRecombs == theState->minRecomb) {
      // TODO: ambiguous for previous state; need way to represent
    }

    if (i == 0 && numIter == 2) {
      // Update the various values as needed for the inverted phase in the next
      // iteration.

      // Children that have missing data should not be flipped
      // Also, for both parent het markers, heterozygous children's bits do
      // not get flipped: they're either unambigous and constrained by
      // <prevState->iv> or they're ambiguous and should remain the same to
      // stick with the convention used in <_stateHash>.
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
  // (1) In the lowest order two bits for which <iv> is unambiguous, the value
  //     should be 00.
  // (2) The <iv> values for ambiguous bits are either the default of 10 binary
  //     or 00 if satisfying (1) involves inverting the inheritance value of
  //     only one parent. (Note that for ambiguous bits, <iv> values of 10 can
  //     be flipped to 01 without increasing the numbers of recombinations.
  //     Likewise 00 can be flipped to 11. So it suffices to only consider only
  //     two ambiguous values: 10 and 00.)

  uint64_t unambig = _parBits[2] - ambig;

  //////////////////////////////////////////////////////////////////////////
  // Get the canonical key value

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

  //////////////////////////////////////////////////////////////////////////
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

// Back traces and minimum recombinant phase using the states in <_hmm>.
void Phaser::backtrace(NuclearFamily *theFam) {
  int lastIndex = _hmm.length() - 1;

  // Find the state with minimum recombinations:
  findMinStates(lastIndex);

  // TODO: this ignores ambiguities at last index
  State *curState = _hmm[lastIndex][ _minStates[0] ];
  State *prevState = NULL;
  for(int hmmIndex = lastIndex; hmmIndex >= 0; hmmIndex--) {
    theFam->setPhase(_hmmMarker[hmmIndex], curState->iv, curState->ambig,
		     curState->hetParent, curState->homParentGeno,
		     curState->parentPhase);
    // TODO: memory leak: bunch of State*s being thrown away here
    _hmm[hmmIndex].clear();


    // In the previous state, resolve ambiguous <iv> values and propagate
    // backward any <iv> values that were unassigned in that state
    if (hmmIndex - 1 >= 0) {
      prevState = _hmm[hmmIndex-1][ curState->prevState ];

      // First propagate backward any <iv> values that were unassigned
      // previously
      prevState->iv |= curState->iv & prevState->unassigned;

      // Note: ambiguities that remain in <curState> will not give information
      // about resolving such in <prevState>
      uint64_t ambigToResolve = prevState->ambig & (~curState->ambig);
      uint64_t ambigRecombs = (curState->iv ^ prevState->iv) & ambigToResolve;
      uint64_t parRecomb[2];
      for(int p = 0; p < 2; p++)
	parRecomb[p] = ambigRecombs & _parBits[p];
      // set both bits for children that inherit a recombination from the given
      // parents
      parRecomb[0] |= parRecomb[0] << 1;
      parRecomb[1] |= parRecomb[1] >> 1;

      // Invert the phase for ambiguous assignments that recombine on both
      // homologs relative to <curState>
      uint64_t bothRecomb = parRecomb[0] & parRecomb[1];
      prevState->iv ^= bothRecomb;

      uint64_t noRecomb = ambigToResolve - (parRecomb[0] | parRecomb[1]);

      // children with 0 recombinations (after flipping <bothRecomb>) are no
      // longer ambiguous.
      // children that are ambiguous in <prevState> _and_ recombine on one
      // homolog relative to <curState> truly have ambiguous phase.
      prevState->ambig -= bothRecomb | noRecomb;
    }

    // ready for next iteration
    curState = prevState;
  }
}

// Given an index <hmmIndex> into <_hmm>, does a linear search to find the
// states with minimum numbers of recombinations. Stores these in <_minStates>.
void Phaser::findMinStates(int hmmIndex) {
  uint16_t minLastMarker = UINT16_MAX;

  int numStates = _hmm[hmmIndex].length();
  for(int i = 0; i < numStates; i++) {
    uint16_t curRecomb = _hmm[hmmIndex][i]->minRecomb;
    if (curRecomb < minLastMarker) {
      minLastMarker = curRecomb;
      _minStates.clear();
      _minStates.append(i);
    }
    else if (curRecomb == minLastMarker) {
      _minStates.append(i);
    }
  }
}
