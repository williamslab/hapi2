// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <bitset>
#include <genetio/util.h>
#include "phaser.h"
#include "cmdlineopts.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray< dynarray<State*> >    Phaser::_hmm;
dynarray<int>                   Phaser::_hmmMarker;
dynarray< dynarray<uint32_t> >  Phaser::_ambigPrevLists;
dynarray< std::pair<uint8_t,uint64_t> > Phaser::_genos;
Phaser::state_ht                Phaser::_stateHash;
Phaser::state_idx_set           *Phaser::_curIdxSet;
Phaser::state_idx_set           *Phaser::_prevIdxSet;
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

  // Ready storage containers to analyze this chromosome
  _hmm.clear();
  _hmmMarker.clear();
  _genos.clear();
  for(int i = 0; i < _ambigPrevLists.length(); i++) {
    _ambigPrevLists[i].clear();
  }
  _ambigPrevLists.clear();
  parBitsInit(numChildren);

  bool bothParMissing = !(theFam->_parents->first->hasData() ||
					   theFam->_parents->second->hasData());

  dynarray<State> partialStates;

  // build states/phase each marker
  for(int m = firstMarker; m <= lastMarker; m++) {
    uint8_t parentData, parentGenoTypes, homParGeno, childGenoTypes;
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
    int mt = getMarkerType(parentGenoTypes, childGenoTypes, homParGeno);
    assert(mt > 0 && (mt & ~((1 << MT_N_TYPES) -1) ) == 0);

    if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
      // should only be one of the above:
      assert(mt == (1 << MT_ERROR) || mt == (1 << MT_AMBIG));
      switch(mt) {
	case 1 << MT_ERROR:
	  theFam->setStatus(/*marker=*/ m, PHASE_ERROR, parentData,
			    childrenData[4], childrenData[G_MISS] &_parBits[0]);
	  break;
	case 1 << MT_AMBIG:
	  theFam->setStatus(/*marker=*/ m, PHASE_AMBIG, parentData,
			    childrenData[4], childrenData[G_MISS] &_parBits[0]);
	  break;
      }
      continue;
    }

    if (mt & (1 << MT_UN)) {
      // TODO: When we the children are heterozygous, phase with parent's
      // data. Note that if we are missing data for both parents and all
      // the children are heterozygous, the marker is ambiguous and handled
      // just above.
      theFam->setStatus(/*marker=*/ m, PHASE_UNINFORM, parentData,
			childrenData[4], childrenData[G_MISS] &_parBits[0]);
      continue;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step 2: For each consistent marker type, deduce the inheritance vector
    // bits where possible and generate partial states that only have bits
    // defined for the informative parent(s)
    makePartialStates(partialStates, mt, parentData, homParGeno, childrenData);

    ///////////////////////////////////////////////////////////////////////////
    // Step 3: Using states for the previous informative marker, fill in
    // unknown inheritance vector bits by propagation and generate full states.
    // Simultaneously calculate the minimum recombination counts for each state
    // and identify which previous state(s) produce that minimum count (i.e.,
    // perform count-based Viterbi calculation)
    makeFullStates(partialStates, /*marker=*/ m, childrenData, bothParMissing,
		   numChildren);

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
    allBitsSet &= (1ul << (2*numChildren)) - 1;
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
    childrenData[ curChildData ] += 3ul << (c*2);
    childrenData[4] += ((uint64_t) curChildData) << (c*2);
    childGenoTypes |= 1ul << curChildData; // observed genotype <curChildData>
  }
}

// Determines what type of marker this is using data for the parents if present
// or based on the observed genotype values for the the children when one or
// both parent's data are missing
int Phaser::getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes,
			  uint8_t &homParGeno) {
  // Only valid values for parentGenoTypes are between 1 and 12 (excluding 7
  // and 11 which are caught below)
  assert(parentGenoTypes >= 1 && parentGenoTypes <= 12);
  assert(childGenoTypes >= 1 && childGenoTypes <= 15);

  if (childGenoTypes == (1 << G_MISS)) {
    // check no data for all children case first: is ambiguous
    return 1 << MT_AMBIG;
  }

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
      else {
	// one parent heterozygous, other homozygous: informative for one parent
	homParGeno = G_HOM0;
	return 1 << MT_FI_1;
      }
      break;
    case ((1 << G_HOM1) | (1 << G_HET)): // one homozygous for 1, other het
      if (childGenoTypes & (1 << G_HOM0))
	// Mendelian error -- child is homozygous for allele not present in
	// the homozygous parent
	return 1 << MT_ERROR;
      else {
	// one parent heterozygous, other homozygous: informative for one parent
	homParGeno = G_HOM1;
	return 1 << MT_FI_1;
      }
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
	  homParGeno = G_HOM0;
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  homParGeno = G_HOM0;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  homParGeno = G_MISS; // not really sure of other parent's genotype
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
	  homParGeno = G_HOM1;
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  homParGeno = G_HOM1;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case G_HOM1: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  homParGeno = G_MISS; // not really sure of other parent's genotype
	  return (1 << MT_UN) | (1 << MT_FI_1);
      }
      break;
    case ((1 << G_HOM0) | (1 << G_HET)): // homozygous 0 and het types observed
    case ((1 << G_HOM0) | (1 << G_HET) | (1 << G_MISS)): // above and missing
      // could be informative for one parent or partly informative; see if there
      // is information on one of the parents that constrains the type:
      switch (missingType) {
	case G_MISS: // both missing data -- no information to add
	  homParGeno = G_HOM0;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  homParGeno = G_HOM0;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  homParGeno = G_HOM0;
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
	  homParGeno = G_HOM1;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HET: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  homParGeno = G_HOM1;
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case G_HOM0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case G_HOM1: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  homParGeno = G_HOM1;
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
	  // other parent either heterozygous or homozygous for other allele
	  homParGeno = missingType;
	  return (1 << MT_UN) | (1 << MT_FI_1);
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
			       uint8_t homParGeno,
			       const uint64_t childrenData[5]) {
  int lastIndex = _genos.length();
  _genos.addEmpty();
  _genos[lastIndex].first = parentData;
  _genos[lastIndex].second = childrenData[4];

  // Fully informative for one parent:
  if (markerTypes & (1 << MT_FI_1)) {
    makePartialFI1States(partialStates, parentData, homParGeno, childrenData);
  }

  // Partly informative
  if (markerTypes & (1 << MT_PI)) {
    makePartialPIStates(partialStates, parentData, childrenData);
  }
}

// Helper for makePartialStates(): applicable to fully informative for one
// parent markers (or the states that correspond to this possibility when the
// parent's genotypes aren't fully known)
void Phaser::makePartialFI1States(dynarray<State> &partialStates,
				  uint8_t parentData, uint8_t homParGeno,
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
  else if ((parentData & 3) == G_HET || (parentData >> 2) == G_HOM0 ||
	   (parentData >> 2) == G_HOM1) {
    startPar = endPar = 0;
  }
  // parent 1 het:
  else {
    assert((parentData >> 2) == G_HET || (parentData & 3) == G_HOM0 ||
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
    // haplotype 0 and allele 1 is on haplotype 1. When the homozygous parent
    // is homozygous for allele 0, the heterozygous children will have received
    // haplotype 1 from the heterozygous parent. However, if the homozygous
    // parent is homozygous for allele 1, the heterozygous children will have
    // received haplotype 0 from the heterozygous parent, and the homozygous
    // children received allele 1 from this parent. So:
    // Assuming the following below:
    static_assert(G_HOM1 == 3 && G_HOM0 == 0, "genotype encoding non-standard");
    // Following always holds; commented out to improve efficiency:
//    assert(homParGeno == G_HOM0 || homParGeno == G_HOM1);
    uint8_t homParAll1 = homParGeno & 1; // binary: homozy parent have allele 1?
    // Which genotype in the children received allele 1 from het parent?
    uint8_t childAll1Geno = homParAll1 * G_HOM1 + (1 - homParAll1) * G_HET;
    newState.iv = _parBits[hetPar] & childrenData[ childAll1Geno ];
    // No ambiguous bits (though a missing data child could have ambiguous bits
    // propagated from a previous marker, but that happens during full state
    // construction)
    newState.ambig = 0;
    // No information about transmission from homozygous parent and for
    // children that are missing data:
    newState.unassigned = _parBits[homPar] | childrenData[ G_MISS ];
    newState.hetParent = hetPar;
    newState.homParentGeno = homParGeno;
    // Won't assign <parentPhase> field as all partial states have a value of
    // 0 here and it's never used
    //newState.parentPhase = 0;
    // Won't assign <arbitraryPar> field: will set this in makeFullStates()

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
  newState.iv = childrenData[ 4 ] & ~childrenData[ G_MISS ];
  // Potentially ambiguous bits are those where the child is heterozygous.
  // If this is the first state on this chromosome, these are ambiguous and
  // their phase is resolved by later markers.
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
			    const uint64_t childrenData[5],
			    bool bothParMissing, int numChildren) {
  // The new entry in <_hmm> corresponds to <marker>:
  _hmmMarker.append(marker);

  if (_hmm.length() == 0) {
    // Partial states are equivalent to full states when there are no previous
    // states
    _hmm.addEmpty();
    int len = partialStates.length();
    assert(len <= 3); // can only have 2 forms of FI partial states and one PI
    for(int i = 0; i < len; i++) {
      if (i == 1 && partialStates[1].hetParent == 1
		 && partialStates[0].hetParent == 0
		 && bothParMissing)
	// For the very first marker, if we're not sure which parent is
	// heterozygous, and both are missing data, arbitrarily pick parent
	// 0 as such. Otherwise, since the two states are equivalent but
	// with opposite parent labels, there will be two equal paths
	// through the HMM.
	// TODO: document this behavior
	continue;
      State *newState = new State(partialStates[i]);
      newState->minRecomb = 0;
      newState->parentPhase = 0;
      newState->ambigParPhase = 0;
      newState->ambigPrev = 0;
      newState->ambigParHet = 0;
      newState->error = 0;
      newState->arbitraryPar = 1;
      _hmm[0].append(newState);
    }
    return;
  }

  int prevHMMIndex = _hmm.length() - 1;
  _hmm.addEmpty();
  dynarray<State*> &prevStates = _hmm[prevHMMIndex];

  // Minimum and maximum recombination counts: for applying the optimization in
  // rmBadStatesCheckErrorFlag().
  // Note that the maximum value may be stale and may not be present in any
  // state, but we track it as an optimization and avoid going through the
  // states when this value is too low to bother.
  uint16_t minMaxRec[2] = { UINT16_MAX, 0 };

  uint32_t numPrev = prevStates.length();
  int numPartial = partialStates.length();
  for(uint32_t prevIdx = 0; prevIdx < numPrev; prevIdx++) {
    const State *prevState = prevStates[prevIdx];

    // See comment above the isIVambigPar() method. We deal with
    // <IVambigPar> == 1 below and in updateStates()
    uint8_t IVambigPar = isIVambigPar(prevState, bothParMissing);

    for(int curIdx = 0; curIdx < numPartial; curIdx++) {
      const State &curPartial = partialStates[curIdx];
      if (IVambigPar && curPartial.hetParent == 1)
	// When both parents are missing data, initial PI states do not actually
	// distinguish the two parents and so later FI markers will be ambiguous
	// for which parent is heterozygous. Here we arbitrarily pick parent 0
	// as the first heterozygous parent in these cases and omit states with
	// parent 1 het.
	// See the code at the beginning of this function for a similar choice
	// at the first first informative marker.
	continue;
      mapPrevToFull(prevState, prevIdx, curPartial, minMaxRec, childrenData,
		    IVambigPar);
    }
  }

  // Introduce error states in prevStates if doing so saves at least
  // <CmdLineOpts::max1MarkerRecomb> recombinations
  // TODO: may be able to optimize. If there a states at both prev and cur that
  // introduce 0 recombinations relative to the most recent one, can we be
  // certain that the marker won't need error states?
  if (CmdLineOpts::max1MarkerRecomb > 0 && prevHMMIndex - 1 >= 0) {
    dynarray<State*> &back2States = _hmm[prevHMMIndex - 1];

    int64_t numBack2 = back2States.length();
    for(int64_t back2Idx = 0; back2Idx < numBack2; back2Idx++) {
      State *back2State = back2States[back2Idx];

      // See comment above the isIVambigPar() method. We deal with
      // <IVambigPar> == 1 below and in updateStates()
      uint8_t IVambigPar = isIVambigPar(back2State, bothParMissing);

      for(int curIdx = 0; curIdx < numPartial; curIdx++) {
	const State &curPartial = partialStates[curIdx];
	if (IVambigPar && curPartial.hetParent == 1)
	  continue; // See comment above in equivalent non-error code
	// Set back2Idx negative and shift by 1 so that the index 0 is also
	// negative
	mapPrevToFull(back2State, /*prevIdx=negative=>error*/ -back2Idx-1,
		      curPartial, minMaxRec, childrenData, IVambigPar);
      }
    }
  }

  rmBadStatesCheckErrorFlag(_hmm[ prevHMMIndex+1 ], minMaxRec, numChildren);
}

// If the IV values transmitted to each child by the two parents are either all
// identical or all opposite, then when <bothParMissing>, which parent is
// heterozygous (for FI markers) will be ambiguous. Additionally, a
// recombination at a PI state that follows these ambiguous markers will have
// two parent phase assignments (opposite each other) that will produce
// equivalent numbers of recombinations.
uint8_t Phaser::isIVambigPar(const State *state, bool bothParMissing) {
  uint64_t IVparDiff = (state->iv & _parBits[0]) ^
					      ((state->iv & _parBits[1]) >> 1);
  // Omit differences at standard ambiguous positions (not only bit 0 in each
  // child is set if there's a difference)
  uint64_t stdAmbig = (state->ambig & _parBits[1]) >> 1;
  IVparDiff &= ~stdAmbig;
  if (bothParMissing && (IVparDiff == 0 ||
					IVparDiff == (_parBits[0] & ~stdAmbig)))
    return 1;

  return 0;
}

// Do the work of mapping a previous state to all states at the current as
// stored in <curPartial> and generate states as needed at the current marker
// (these are stored in _hmm).
// If <prevIdx> is negative, it is the index for a state that is two markers
// previous, and is introduced in order to detect erroneous markers evidenced
// by large numbers of recombinations at a single marker (in this case, the
// immediately previous marker; it will be skipped over and marked as an error
// if by mapping a state from two markers back to the current state with a
// penalty term there are overall fewer recombinations).
void Phaser::mapPrevToFull(const State *prevState, int64_t prevIdx,
			   const State &curPartial, uint16_t minMaxRec[2],
			   const uint64_t childrenData[5],
			   uint8_t IVambigPar) {
  // (1) For each (current) partial state, map to an initial set of full state
  // values from <prevState>. The full <iv> and <unassigned> values are
  // determined based the corresponding fields in <prevState> and <curPartial>.
  // <curPartial.unassigned> has bits set to 1 for the parent that is
  // homozygous/uninformative or both bits set to 1 for a child that is missing
  // data. There is no information in <curPartial> about the haplotype
  // transmissions for these cases and so we propagate the inheritance vector
  // values from the previous inheritance vector.
  uint64_t fullIV = curPartial.iv | (prevState->iv & curPartial.unassigned);
  uint64_t fullUnassigned = prevState->unassigned & curPartial.unassigned;
  uint64_t fullAmbig; // assigned below
  /////////////////////////////////////////////////////////////////////////////
  // only applicable to MT_PI states and defined below:
  // Relative to the default phase for the parents, which bits recombined
  // in the heterozygous children that were not ambiguous at the previous
  // marker? Recombination bits may be 00, 01, 10, 11 binary:
  uint64_t unambigHetRecombs[4];
  // Indicates the children have <prevState->iv> values that were
  // unassigned for parent 0/1 (but not both)
  uint64_t childPrevUnassigned[2];
  uint64_t propagateAmbig;
  uint64_t defaultPhaseHasRecomb;

  // Which haplotypes recombined? Note that the current <iv> value assumes a
  // certain phase for the parents and may have more recombinations than
  // another possibility. We consider other possibilities below in
  // updateStates() and flipPIVals().
  // The following gives us enough information to determine how to handle the
  // heterozygous children at MT_PI markers.
  // Note: because (fullIV & curPartial.unassigned) ==
  //                                	(prevState->iv & curPartial.unassigned),
  // we know that no recombinations occur at <curPartial.unassigned> bits.
  // We must avoid counting recombinations for bits that were unassigned
  // in <prevState>, however: such differences between <prevState->iv> and
  // <fullIV> are not meaningful since those bits weren't assigned in
  // <prevState>.
  uint64_t recombs = (prevState->iv ^ fullIV) & ~prevState->unassigned;
  // set both bits for children that inherit a recombination from the given
  // parent
  uint64_t parRecombs[2] = { (recombs & _parBits[0]) * 3,
			     ((recombs & _parBits[1]) >> 1) * 3 };
  // These variables are used to determine whether we need to explore all the
  // different possible phase types for parents. Typically the answer is yes,
  // but at the beginning of the chromosome, if an informative marker for one
  // of the parents hasn't yet been seen, then the two possible phase types for
  // the first marker that is informative for that parent will not in fact
  // differ in their numbers of recombinations. To avoid indicating these as
  // ambiguous, we set these values and check them in updateStates() and below
  // for PI states.
  bool hetParentUndefined, oneParentUndefined;

  // (2) if curPartial is partly informative, for heterozygous children:
  //     (a) When two recombinations occur relative to the previous marker and
  //         the child was unambiguous previously, flip both corresponding bits
  //         in the inheritance vector. HAPI avoids a state space explosion for
  //         heterozygous children at MT_PI markers by a combination of only
  //         modeling states with 0 recombinations instead of 2 when possible
  //         and (b),(c) below.
  //     (b) When one recombination occurs relative to the previous marker,
  //         leave the iv value as assigned and set the child as ambiguous.
  //     Note: ambiguous bits are only set for:
  //      (i)  children that are newly ambiguous: i.e, those that are
  //           heterozygous and exhibit one recombination relative to the
  //           previous marker
  //      (ii) or children that are heterozygous or missing and were ambiguous
  //           at the previous marker (this status propagates forward until an
  //           unambiguous marker, including the child being homozygous at
  //           MT_PI)
  // ... Also deals with complexities around children that were unassigned in
  //     <prevState>. These produce ambig1 type ambiguities described in
  //     handlePI().
  // The following gives 1 for <hetParent> == 2, 0 for <hetParent> == 0,1
  uint8_t isPI = curPartial.hetParent >> 1;
  if (isPI) { // MT_PI state
    handlePI(prevState, fullIV, fullAmbig, recombs, parRecombs, propagateAmbig,
	     defaultPhaseHasRecomb, childPrevUnassigned, unambigHetRecombs,
	     childrenData);

    bool parUndefined[2];
    for (int p = 0; p < 2; p++) {
      uint64_t parBits = _parBits[p];
      parUndefined[p] = (prevState->unassigned & parBits) == parBits;
    }
    oneParentUndefined = parUndefined[0] || parUndefined[1];
    hetParentUndefined = false; // not applicable to PI states
  }
  else {
    // Fully informative marker: only ambiguous bits are those where a child
    // was ambiguous in the previous state and missing data at this marker
    fullAmbig = childrenData[G_MISS] & prevState->ambig;

    uint64_t hetParBits = _parBits[curPartial.hetParent];
    hetParentUndefined = (prevState->unassigned & hetParBits) == hetParBits;
    oneParentUndefined = false; // not applicable to FI states
  }

  // (3) As needed, remove apparent recombinations from <iv> values that were
  // ambiguous in the previous state
  // Which children were ambiguous in the previous state but not here?
  // Also various values relating to ambig1 type ambiguous values
  uint64_t stdAmbigOnlyPrev, ambig1PrevInfo, ambig1Unassigned;
  fixRecombFromAmbig(fullIV, recombs, parRecombs, isPI,
		     /*ambigOnlyPrev=*/ prevState->ambig & ~fullAmbig,
		     curPartial.hetParent, stdAmbigOnlyPrev, ambig1PrevInfo,
		     ambig1Unassigned);

  // (4) Look up or create a full state with equivalent <iv> and <ambig> values
  // to <fullIV> and <fullAmbig>, and determine if <prevState> yields fewer
  // recombinations for these states than the currently stored previous state
  // (if any). If so, update the necessary values in the state.
  // Also examines an alternate phase type which may or may not map to an
  // equivalent state. If not, looks up that value, if so, compares the
  // recombinations for the two possibilities separately.
  //
  // What type of phase is the alternative? For MT_PI states, the alternative,
  // which must have the same ambiguous bits, is 3 decimal == 11 binary (vs.
  // default of 0).
  // See just below for what this code does:
  uint8_t altPhaseType = isPI * 3 + (1 - isPI);
  // Above equivalent to the line below but has no branching
  //uint8_t altPhaseType = (isPI) ? 3 : 1;
  updateStates(fullIV, fullAmbig, fullUnassigned, ambig1Unassigned, recombs,
	       prevState->unassigned, stdAmbigOnlyPrev, ambig1PrevInfo,
	       curPartial.hetParent, curPartial.homParentGeno,
	       /*initParPhase=default phase=*/ 0, altPhaseType, prevIdx,
	       prevState->minRecomb, prevState->error, IVambigPar,
	       minMaxRec, hetParentUndefined, childrenData);

  // For MT_PI states, have 1 or 2 more states to examine:
  // Exception is if one of the parents doesn't have any transmitted haplotypes
  // defined; in that case, it suffices to consider only the two possible
  // states defined above, since the transmissions for the unassigned parent
  // are trivially ambiguous
  if (isPI && !oneParentUndefined) {
    // Above considered default and state with both parents' phase inverted.
    // Now consider the two states in which each parent alone has inverted
    // phase.

    ///////////////////////////////////////////////////////////////////////////
    // Generate the <fullIV> and <fullAmbig> values for parent 0 flipped:
    // (Note: <fullUnassigned> does not change)
    flipPIVals(fullIV, fullAmbig, childrenData, propagateAmbig,
	       unambigHetRecombs, childPrevUnassigned, defaultPhaseHasRecomb);
    recombs = (prevState->iv ^ fullIV) & ~prevState->unassigned;
    parRecombs[0] = (recombs & _parBits[0]) * 3;
    parRecombs[1] = ((recombs & _parBits[1]) >> 1) * 3;

    // TODO: potentially can optimize by storing information obtained in the
    // first call to fixRecombFromAmbig()
    fixRecombFromAmbig(fullIV, recombs, parRecombs, /*isPI=*/ 1,
		       /*ambigOnlyPrev=*/ prevState->ambig & ~fullAmbig,
		       curPartial.hetParent, stdAmbigOnlyPrev,
		       ambig1PrevInfo, ambig1Unassigned);

    // Following is not failing, so comment out for efficiency:
//    assert((recombs & childrenData[G_HET] & (childPrevUnassigned[0] |
//						 childPrevUnassigned[1])) == 0);

    ///////////////////////////////////////////////////////////////////////////
    // Now ready to look up or create full states with the <fullIV> and
    // <fullAmbig> values, etc.
    updateStates(fullIV, fullAmbig, fullUnassigned, ambig1Unassigned, recombs,
		 prevState->unassigned, stdAmbigOnlyPrev, ambig1PrevInfo,
		 /*curPartial.hetParent=*/2, /*homParentGeno=*/G_MISS,
		 /*initParPhase=parent 0 flip=*/1,
		 /*altPhaseType=parent 1 flip=*/ 2, prevIdx,
		 prevState->minRecomb, prevState->error, IVambigPar,
		 minMaxRec, /*hetParentUndefined=*/ false, childrenData);
  }
}

// See long comment in makeFullStates() -- handles heterozygous children and
// previously unassigned <iv> values at MT_PI states:
void Phaser::handlePI(const State *prevState, uint64_t &fullIV,
		      uint64_t &fullAmbig, uint64_t &recombs,
		      uint64_t parRecombs[2], uint64_t &propagateAmbig,
		      uint64_t &defaultPhaseHasRecomb,
		      uint64_t childPrevUnassigned[2],
		      uint64_t unambigHetRecombs[4],
		      const uint64_t childrenData[5]) {
  // Set both bits for any children that have unassigned values for each
  // parent in the previous state:
  childPrevUnassigned[0] = (prevState->unassigned & _parBits[0]) * 3;
  childPrevUnassigned[1] = ((prevState->unassigned & _parBits[1]) >> 1)*3;
  uint64_t bothPrevUnassigned = childPrevUnassigned[0] & childPrevUnassigned[1];
  uint64_t eitherPrevUnassigned = childPrevUnassigned[0]|childPrevUnassigned[1];
  // Want <childPrevUnassigned> to be specific to each parent, not both.
  // This comes up later in flipPIVals()
  for(int p = 0; p < 2; p++)
    childPrevUnassigned[p] -= bothPrevUnassigned;

  // In certain circumstances, we assign only one parent (parent 0) as
  // ambiguous. This is to deal with complexities around when only one parent
  // has an assigned <iv> value and a PI marker occurs in which a child is
  // heterozygous. Such a marker and heterozygous genotype gives information
  // about allelic transmissions from both parents, but switching to the
  // opposite <iv> will only incur one recombination relative to the one parent
  // whose <iv> value has been assigned, not two. So, as distinct from
  // standard ambiguous <iv> values, there are two <iv> values, one of which
  // gives 0 recombinations relative to the previous markers, the other of
  // which, despite being inverted for both parents, gives only 1 recombination.
  // We call these ambig1 and deal with them separately. Note that they are
  // similar to unambigous <iv> values in that once established (i.e., after
  // the single ambig bit is set) we seek to match the <iv> values to the prior
  // marker while propagating foward the ambig1 status. Thus, in particular,
  // they should be considered by the calcHetChildPIRecombs() method
  uint64_t prevStdAmbig = ((prevState->ambig & _parBits[1]) >> 1) * 3;

  // Determine the recombination values for heterozygous children that are
  // unambiguous and fully assigned in <prevState>. To avoid recombination on
  // both homologs on these children, we must sometimes flip the <fullIV>
  // values. Also, a single recombination is ambiguous and must be accounted
  // for. See below.
  uint64_t unambigAssignedHets = childrenData[G_HET] &
					~(prevStdAmbig | eitherPrevUnassigned);
  calcHetChildPIRecombs(parRecombs, unambigAssignedHets, unambigHetRecombs);

  // Remove the recombinations from both parents for het children (whose
  // <iv> values can be swapped)
  fullIV ^= unambigHetRecombs[3];
  recombs ^= unambigHetRecombs[3];
  for(int p = 0; p < 2; p++)
    parRecombs[p] -= unambigHetRecombs[3];

  // Determine which children are ambiguous:
  // (1) Some from the previous marker (when het/missing here):
  propagateAmbig =
		(childrenData[G_MISS] | childrenData[G_HET]) & prevState->ambig;
  // (2) Near the beginning of the chromosome, if a child has not had either
  //     <iv> value assigned previously and is het here, it is ambiguous. We
  //     can choose either phase assignment without affecting the number of
  //     recombinations.
  //   (We put this together with propagateAmbig so that it gets reused for the
  //   alternate phase types addressed in flipPIVals().)
  propagateAmbig |= childrenData[G_HET] & bothPrevUnassigned;
  // (3) Some newly ambiguous due to recombinations from one parent:
  uint64_t newAmbigChildren = unambigHetRecombs[1] | unambigHetRecombs[2];
  // (4) New 1ambig children:
  //     Children that are het with only one parent assigned previously are
  //     ambig1 type ambiguous (see above); mark them as such.
  //     Want only one bit set in <fullAmbig> for these children: arbitrarily
  //     use parent 0 bits.
  //   (We put this together with propagateAmbig so that it gets reused for the
  //   alternate phase types addressed in flipPIVals().)
  uint64_t onePrevUnassigned = eitherPrevUnassigned - bothPrevUnassigned;
  uint64_t hetChildPrevUnassign = childrenData[G_HET] & onePrevUnassigned;
  propagateAmbig |= hetChildPrevUnassign & _parBits[0]; // new ambig1 children

  fullAmbig = newAmbigChildren | propagateAmbig;

  // One final issue with unassigned <iv> values:
  // For het children with only one bit unassigned previously, we can flip both
  // <iv> values to avoid recombinations with respect to the parent that is
  // assigned in the previous marker and not incur a recombination from the
  // alternate (currently unassigned) parent.
  // Need to have both bits set for any children with a recombination from
  // either parent:
  defaultPhaseHasRecomb = parRecombs[0] | parRecombs[1];
  // Flip the <iv> value and remove any recombinations for these children:
  uint64_t toRemove = hetChildPrevUnassign & defaultPhaseHasRecomb;
  fullIV ^= toRemove;
  recombs &= ~toRemove;

  // Following is not failing, so comment out for efficiency:
//  assert(recombs == ((prevState->iv ^ fullIV) & ~prevState->unassigned &
//			~(childrenData[G_HET] &
//			   (childPrevUnassigned[0] | childPrevUnassigned[1]))));
}

// For MT_PI states, determines <unambigHetRecombs>, an array with both bits
// set when a heterozygous but unambiguous child has the indexed type of
// recombinations. For example, if the children's bits in <unambigHetRecombs[1]>
// are set, that child recombines on parent 0. Index 0 corresponds to no
// recombinations, 2 corresponds to parent 1 recombinations, 3 corresponds to
// both parents recombining.
void Phaser::calcHetChildPIRecombs(const uint64_t parRecombs[2],
				   const uint64_t unambigAssignedHets,
				   uint64_t unambigHetRecombs[4]) {
  uint64_t parRecombUnambigHet[2];
  for(int p = 0; p < 2; p++)
    parRecombUnambigHet[p] = parRecombs[p] & unambigAssignedHets;

  // unambiguous het children that exhibit recombinations from both parents:
  unambigHetRecombs[3] = parRecombUnambigHet[0] & parRecombUnambigHet[1];

  // unambiguous het children with a recombination from only one parent:
  for(int r = 1; r <= 2; r++)
    unambigHetRecombs[r] = parRecombUnambigHet[r-1] ^ unambigHetRecombs[3];

  // unambiguous het children with no recombinations:
  unambigHetRecombs[0] = unambigAssignedHets -
			      (parRecombUnambigHet[0] | parRecombUnambigHet[1]);
}

// For MT_PI states, flipping a single parent's phase is complex for
// heterozygous children. Here we generate the <fullIV>, <fullAmbig>, and
// <recomb> values wherein parent 0's phase is flipped relative to the default.
// We base this on the values passed in as arguments, including the <fullIV>
// and <fullAmbig> values generated for the default phase.
void Phaser::flipPIVals(uint64_t &fullIV, uint64_t &fullAmbig,
			const uint64_t childrenData[5], uint64_t propagateAmbig,
			const uint64_t unambigHetRecombs[4],
			const uint64_t childPrevUnassigned[2],
			uint64_t defaultPhaseHasRecomb) {
  ///////////////////////////////////////////////////////////////////////
  // Get the flipped <fullIV> value:

  // <iv> value is the same as the default for all missing data children:
  uint64_t newFullIV = fullIV & childrenData[G_MISS];
  // flip the transmitted haplotype for parent 0 for all homozygous
  // children
  newFullIV |= (fullIV ^ _parBits[0]) &
				  (childrenData[G_HOM0] | childrenData[G_HOM1]);
  // For heterozygous children (only remaining children to be defined), the
  // value in <newFullIV> is 00. This is the value we want for all standard
  // ambiguous children -- it is the canonical value for looking up in
  // <_stateHash> for this phase type (and for parent 1 flipped phase).
  // Thus, we only need to modify heterozygous children's IV values if they are
  // non-ambiguous.
  // Default <iv> value was 10 binary. Any children that recombined on parent 1
  // under that setup now have no recombinations and have the correct <iv>
  // value. Children that recombined on parent 0 now have 2 recombinations and
  // need to be flipped to 11 binary. Recombination on parent 0 is index 01:
  newFullIV ^= unambigHetRecombs[1];
  // For het children that had unassigned <prevState->iv> values, can flip the
  // <iv> so that no recombinations result. Since the default assignment has
  // <iv> value of 10 binary, and the value now is 00, if the default phase
  // has a recombination and parent 0 is unassigned, we'd like to have value 11
  // so that there's no recombination from parent 1. Alternatively, if the
  // default phase has no recombination and parent 1 is unassigned, we'd also
  // like a value of 00 so that there's no recombination from parent 0:
  uint64_t flipPartlyUnassigned = childrenData[G_HET] &
			    ((childPrevUnassigned[0] & ~defaultPhaseHasRecomb) |
			     (childPrevUnassigned[1] & defaultPhaseHasRecomb));
  newFullIV ^= flipPartlyUnassigned;
  // Done!
  fullIV = newFullIV;

  // Children that either recombined on both parents or did not recombine at
  // all previously are ambiguous when the phase one parent is flipped, so:
  uint64_t newAmbigChildren = unambigHetRecombs[0] | unambigHetRecombs[3];

  fullAmbig = newAmbigChildren | propagateAmbig;
}

// For children with an ambiguous <iv> value in the previous state,
// apparent recombinations to the current state are not recombinations since
// we can flip the <iv> value in the previous state. Here we flip any such
// recombinations to avoid counting them. The <iv> value in the previous state
// gets changed during back tracing.
void Phaser::fixRecombFromAmbig(uint64_t &fullIV, uint64_t &recombs,
				const uint64_t parRecombs[2], uint8_t isPI,
				uint64_t ambigOnlyPrev, uint8_t hetParent,
				uint64_t &stdAmbigOnlyPrev,
				uint64_t &ambig1PrevInfo,
				uint64_t &ambig1Unassigned) {
  // As detailed elsewhere, have ambig1 cases where only parent bit 0 is set
  // ambiguous that are distinct from standard ambiguous cases
  uint64_t parAmbigOnlyPrev[2] = { (ambigOnlyPrev & _parBits[0]) * 3,
				   ((ambigOnlyPrev & _parBits[1]) >> 1) * 3 };
  stdAmbigOnlyPrev = parAmbigOnlyPrev[1];
  uint64_t ambig1OnlyPrev = parAmbigOnlyPrev[0] - parAmbigOnlyPrev[1];
  uint64_t anyAmbigOnlyPrev = parAmbigOnlyPrev[0];

  // TODO: test optimization: don't call this function if ambigOnlyPrev == 0?

  // Get their recombinations specific to each parent:
  uint64_t parRecombsFromAmbig[2];
  for(int p = 0; p < 2; p++)
    parRecombsFromAmbig[p] = anyAmbigOnlyPrev & parRecombs[p];

  // Two general cases:
  // (1) if the current state is PI (i.e., <isPI> == 1), we will only flip when
  // there's a recombination on both parents
  // The <toRemove> assignment below is only non-zero if <isPI> == 1:
  uint64_t bothParRecomb = parRecombsFromAmbig[0] & parRecombsFromAmbig[1];
  // When only one parent is marked ambiguous, this is the ambig1 case and
  // flipping the previous <iv> yields 1 recombination (not 0 or 2), so we
  // only remove one recombination (arbitrarily on parent 1) in that case:
  uint64_t toRemove = isPI * (bothParRecomb & ~(ambig1OnlyPrev & _parBits[0]));
  // (2) If the current state is FI (i.e., <isPI> == 0), we will only flip when
  // there's a recombination for the parent that is heterozygous here (note that
  // the homozygous parent won't have observed recombinations anyway). This only
  // applies for standard ambiguous values; ambig1 is below.
  // The <toRemove> assignment is only non-zero if <isPI> == 0:
  uint8_t hetParIdx = (1 - isPI) * hetParent;
  uint64_t hetParRecombBothAmbig = parRecombsFromAmbig[ hetParIdx ] &
							      stdAmbigOnlyPrev;
  toRemove += (1 - isPI) * (hetParRecombBothAmbig & _parBits[hetParIdx]);
  // Can/will flip <prevState->iv> during back tracing, so don't count
  // recombinations now
  recombs ^= toRemove;

  // The ambig1 type ambiguous values have either 0 or 1 recombinations:
  // recombinations from both parents become 1 recombination at an earlier
  // marker (see above). The current <recombs> values applies for the current
  // assignment (in the caller) of <fullIV>, but in updateStates(), we examine
  // two phase possibilities. It suffices when <isPI> == 0, to simply flip
  // 0 to 1 and 1 to 0 for the heterozygous parent, and 2 recombinations are
  // impossible.
  // When <isPI> == 1, 2 recombinations can occur and, for values with no
  // recombination on either parent or 2 recombinations, we need to produce
  // 0 or 1 recombinations and then to invert these values for the alternate
  // phase type. We do this using <ambig1PrevInfo>: it indicates which
  // <recomb> bits should stay the same -- so we flip bits that are set 0
  // in this value. In general, if there is a recombination only on parent 0
  // or only on parent 1, it suffices _not_ to flip <recomb> since both phase
  // possibilities produce 1 recombination. Thus those bits stay as is and are
  // set to 11 binary below. When both recombine or neither do, we must invert
  // the recombination assignment made (or not made) above in
  // <toRemove = isPI * ...>. We chose there to have a recombination on parent 0
  // so the assignment below sets (via a ~ that makes this confusing)
  // the parent 0 bit in <ambigPrevInfo> to 0 for children that show
  // recombinations here from both parents or from neither parent. Thus that
  // bit will be flipped later in updateStates().
  uint64_t noRecomb = ~(parRecombsFromAmbig[0] | parRecombsFromAmbig[1]);
  ambig1PrevInfo = isPI * (ambig1OnlyPrev & ~((bothParRecomb | noRecomb) &
								  _parBits[0]));

  // When we change the number of recombinations, we have decided that the
  // previous state should have inverted <iv>. In this case, when <isPI> == 0,
  // the <iv> values corresponding to the homozygous parent are no longer
  // reflective of the value that will occur once this flipping occurs.
  // Therefore, we must flip the corresponding bits in <fullIV>:
  uint8_t homozyParent = (1 - isPI) * (1 - hetParent);
  fullIV ^= (1 - isPI) * (hetParRecombBothAmbig & _parBits[ homozyParent ]);

  // For ambig1 values transitioning to FI markers with a recombination from one
  // parent, in fact, either <iv> value assignment at the previous marker --
  // which affects the <iv> value for the homozygous parent -- incurs 0
  // additional recombinations (beyond the 1 for transitioning to this current
  // state). This is analogous to the <iv> value being unassigned, and we make
  // use of this fact to prevent the need of tacking this completely ambiguous
  // (but only for one parent) case.
  // Note that parent's transmission is considered assigned if there are 0
  // recombinations (since a later difference would recombine relative to the
  // <iv> value assigned at the previous marker. But when we consider the
  // alternate phase type in updateStates(), the children that recombined/didn't
  // will be inverted. We therefore use <ambig1PrevInfo> to store which
  // ambig1Unassigned bits to flip later in updateStates():
  ambig1PrevInfo += (1 - isPI) * ambig1OnlyPrev & _parBits[homozyParent];
  ambig1Unassigned = (1 - isPI) * (ambig1PrevInfo &
					      parRecombsFromAmbig[ hetParIdx ]);
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
			  uint64_t fullUnassigned, uint64_t ambig1Unassigned,
			  uint64_t recombs, uint64_t prevUnassigned,
			  uint64_t stdAmbigOnlyPrev, uint64_t ambig1PrevInfo,
			  uint8_t hetParent, uint8_t homParentGeno,
			  uint8_t initParPhase, uint8_t altPhaseType,
			  int64_t prevIndex, uint16_t prevMinRecomb,
			  uint8_t prevError, uint8_t IVambigPar,
			  uint16_t minMaxRec[2], bool hetParentUndefined,
			  const uint64_t childrenData[5]) {
  // How many iterations of the loop? See various comments below.
  int numIter = 2;

  // How many recombinations for the initial phase assignment?
  size_t numRecombs = popcount(recombs);
  uint8_t curParPhase = initParPhase;

  // Note: (hetParent >> 1) == 1 iff hetParent == 2. It is 0 for the other
  //       values.
  uint8_t isPI = hetParent >> 1;

  // Is it possible to swap the phase of this state and obtain the same
  // number of recombinations from the previous state? Decided below
  uint8_t ambigLocal = 0;

  // First, decide how many iterations of the loop below. If the opposite
  // phase assignment yields an equivalent state, we only need to loop once.
  // Only way the opposite phase assignment is equivalent is if:
  // (1) There is no missing data (we propagate <iv> values from the previous
  // marker for missing data children so changing the parent's phase at the
  // current marker will modify some children's <iv> values but not the missing
  // data ones).
  // An exception to (1) occurs at PI states -- if the child both has missing
  // data and is ambiguous, the IV value will be equivalently ambiguous for
  // both states. (That's not the case for FI states since the ambiguous values
  // the are equivalent to one another are flipped for both parents. At FI
  // states we consider only flipping one parent and so the [fixed] ambig value
  // will be different relative to the other children that are flipped.)
  // (2) Similarly, there have to be no children that were only ambiguous
  // in the previous state and not here. To avoid recombinations relative to
  // the previous state the previous <iv> value can be flipped. When this
  // happens, the propagated <iv> value differs between the two possible states.
  // An exception to (2) occurs at PI states -- when a child is homozygous
  // there's nothing to propagate so the two states will be equivalent for such
  // children.
  // (3) For PI type states, must have all heterozygous markers be (standard)
  // ambiguous. Otherwise, since flipping will necessitate fixing the
  // unambiguous children to match the previous state (to avoid recombinations
  // from both parents), the resulting state will not be equivalent (akin to
  // missing data children).
  if (hetParentUndefined) {
    // also when the heterozygous parent wasn't defined in the previous state
    // we need only consider one possibility. Otherwise the two possibilities
    // will be trivially ambiguous
    numIter = 1;
  }
  else if ( ((1-isPI) && childrenData[G_MISS] == 0 && stdAmbigOnlyPrev == 0 &&
						      ambig1Unassigned == 0) ||
       (isPI && (stdAmbigOnlyPrev & childrenData[G_HET]) == 0
	     && (childrenData[G_MISS] | childrenData[G_HET]) == fullAmbig)) {
    numIter = 1;

    // Determine which phase assignment has minimal recombinations before the
    // loop below.

    // Don't flip ambiguous bits: these should be the same as for the original
    // option to stay with the convention used in the hash (see lookupState()).
    // Also don't flip recombs when the previous state was ambig1. Earlier
    // code removed recombs that could be addressed by flipping that ambiguous
    // <iv> value and that change applies equally to both phase possibilities
    // here.
    // Note we exclude any cases where the previous state had standard ambig
    // values in the condition above (see comment).
    // See comment below in the analogous operation regarding <ambig1PrevInfo>.
    // Finally, don't inadvertently introduce recombinations by flipping values
    // that aren't assigned in the previous state.
    recombs ^= _parBits[ hetParent ] & ~(fullAmbig | (isPI * ambig1PrevInfo) |
								prevUnassigned);
    size_t curCount = popcount(recombs);
    if (curCount < numRecombs) {
      numRecombs = curCount;
      curParPhase = altPhaseType;
      fullIV ^= _parBits[ hetParent ] & ~fullAmbig;
      // No need to excute this as we ensure that (1-isPI) * ambig1PrevInfo == 0
//      ambig1Unassigned ^= (1 - isPI) * ambig1PrevInfo;
    }
    else if (curCount == numRecombs) {
      // Either of the two phase types give the same number of recombinations.
      // Will indicate this in the state below.
      ambigLocal = 1;
    }
  }

  State *lastState = NULL;

  // Examine the 1 or 2 needed phase assignments
  for(int i = 0; i < numIter; i++) {
    // First look up or create a full state with equivalent <iv> and <ambig>
    // values to <fullIV> and <fullAmbig>
    State *theState = lookupState(fullIV, fullAmbig,
				  fullUnassigned | ambig1Unassigned);
    // TODO: add a check for prevState->minRecomb > theState.minRecomb?

    // Should never map to the same state as the previous iteration (this should
    // be caught in the code above that sets numIter = 1)
    assert(lastState == NULL || theState != lastState);
    lastState = theState;

    int totalRecombs = prevMinRecomb + numRecombs;
    if (prevIndex < 0)
      totalRecombs += CmdLineOpts::max1MarkerRecomb; // penalty for error states

    // Does the current previous state lead to minimum recombinations for
    // <theState>? (Could be either a new minimum or an ambiguous one)
    bool theStateUpdated = false;

    // TODO: document the error case (prefer paths without errors)
    // Is the previous state <prevIndex> the new minimally recombinant path to
    // <theState>? Certainly if it produces fewer recombinations.
    // Also if it has equal numbers of recombinations and <theState> is either
    // an error state or has an error state in its previous state path and the
    // proposed new state has no error in its path nor is it an error state.
    if (totalRecombs < theState->minRecomb ||
	(totalRecombs == theState->minRecomb &&
	 theState->error > 0 && prevError == 0 && prevIndex >= 0)) {
      theStateUpdated = true;

      theState->iv = fullIV;
      // next two already assigned in lookupState():
//      theState->ambig = fullAmbig;
//      theState->unassigned = fullUnassigned | ambig1Unassigned;
      theState->minRecomb = totalRecombs;
      theState->hetParent = hetParent;
      theState->homParentGeno = homParentGeno;
      theState->parentPhase = curParPhase;
      theState->ambigParPhase = ambigLocal * (1 << altPhaseType);
      theState->ambigPrev = 0;
      theState->ambigParHet = 0;
      theState->arbitraryPar = (1-isPI) * IVambigPar;
      if (prevIndex >= 0) {
	theState->prevState = prevIndex;
	// Propagate error state information: if the path of states leading to
	// this one has an error in it, indicate this by having
	// <theState->error> == 2. So we'll copy forward the value of previous
	// value of 2 or convert a previous value of 1 into 2:
	// The following produces 0 if prevError == 0 and 2 otherwise:
	theState->error = (prevError & 2) | ((prevError & 1) << 1);
      }
      else { // erroneous previous state
	theState->prevState = -(prevIndex + 1);
	theState->error = 1;
      }

      if (totalRecombs < minMaxRec[0])
	minMaxRec[0] = totalRecombs;
      else if (totalRecombs > minMaxRec[1])
	minMaxRec[1] = totalRecombs;
    }
    // Is the proposed state ambiguous with the current minimum?
    // Requires that they produce equal numbers of recombinations and:
    // Either the new state is not an error or both it and <theState> are
    // marked as errors. Ensure that they are the same types of errors as well:
    // the state is either a brand new error for all previous paths leading to
    // it or has an error somewhere in its previous path. Ideally would like
    // to do something slightly more elaborate than this as this will
    // arbitrarily pick the location of the error state with no indication
    // about the ambiguity.  The effect of this is to prefer paths without
    // errors in the case of ambiguities.
    else if (totalRecombs == theState->minRecomb &&
	     (prevIndex >= 0 || (theState->error & 1)) &&
	     (prevError == 0 || (theState->error & 2))) {
      theStateUpdated = true;

      // binary indicator for whether the two <hetParent> values differ
      uint8_t hetParIsDiff;
      theState->ambigParHet |= calcAmbigParHetBits(theState->hetParent,
						   hetParent,
						   hetParIsDiff);

      // ambiguous: could get here from >1 previous state:
      // Note: the order of adding error states in makeFullStates() is such
      //       that you can't have <theState->error> set and prevIndex not
      //       also be an error.
      theState->ambigParPhase |= ambigLocal * (1 - hetParIsDiff) *
				    ((1 << curParPhase) | (1 << altPhaseType));

      const uint32_t MAX_IDX_IN_STATE = 1 << 6;

      // Have an ambiguous previous state -- add to list if one exists or
      // create one
      uint32_t prevState = (prevIndex >= 0) ? prevIndex : -(prevIndex + 1);
      ////////////////////////////////////////////////////////////////////////
      // See comment at declaration of State::ambigPrev
      // TODO: refactor: put in its own function
      switch (theState->ambigPrev) {
	case 2: // already storing list, so simply append
	  {
	    // Can happen that <prevState> is already in the list; check whether
	    // that's the case (will necessarily be the last element since we
	    // consider prevStates sequentially in makeFullStates())
	    int lastElement = _ambigPrevLists[theState->prevState].length() - 1;
	    if (_ambigPrevLists[theState->prevState][lastElement] != prevState)
	      _ambigPrevLists[ theState->prevState ].append(prevState);
	    else
	      // only ambiguous parent phase if the heterozygous parent matches
	      // the one listed
	      theState->ambigParPhase |= (1 - hetParIsDiff) *(1 << curParPhase);
	  }
	  break;

	case 1: // storing ambiguous values in series of 6 bits in <prevState>
	  {
	    bool inserted = false;
	    const uint32_t mask = (1 << 6) - 1;

	    if (prevState < MAX_IDX_IN_STATE) {
	      // <prevState> fits in 6 bits; simultaneously check; if there's
	      // room for this new value; find where to put it; and ensure the
	      // same value doesn't occur twice

	      // first 6 bits shouldn't be the same as <prevState>
	      uint32_t curPrevs = theState->prevState;
	      assert((curPrevs & mask) != prevState);
	      curPrevs >>= 6;
	      for(int shift = 6; shift < 30; shift += 6) {
		if (curPrevs == 0) {
		  // found slot to insert: do so
		  theState->prevState |= prevState << shift;
		  inserted = true;
		  break;
		}
		if ((curPrevs & mask) == prevState) {
		  theState->ambigParPhase |= (1 - hetParIsDiff) *
							    (1 << curParPhase);
		  inserted = true;
		  break;
		}
		curPrevs >>= 6;
	      }
	    }

	    if (!inserted) {
	      // Must store all ambiguous values in a list
	      theState->ambigPrev = 2;
	      uint32_t curPrevs = theState->prevState;
	      uint32_t ambigIndex = _ambigPrevLists.length();
	      theState->prevState = ambigIndex;
	      _ambigPrevLists.addEmpty();
	      for(int shift = 0; shift < 30; shift += 6) {
		_ambigPrevLists[ambigIndex].append( curPrevs & mask );
		curPrevs >>= 6;
		if (curPrevs == 0)
		  break; // got to last previous index
	      }
	      _ambigPrevLists[ambigIndex].append(prevState);
	    }
	  }
	  break;

	case 0: // previously unambiguous
	  {
	    uint32_t curPrev = theState->prevState;
	    if (curPrev == prevState) {
	      theState->ambigParPhase |= (1 - hetParIsDiff) *(1 << curParPhase);
	    }
	    else if (curPrev < MAX_IDX_IN_STATE &&
						prevState < MAX_IDX_IN_STATE) {
	      // Both values fit in 6 bits, set ambiguity and add the new
	      // previous state
	      theState->ambigPrev = 1;
	      theState->prevState |= prevState << 6;
	    }
	    else {
	      theState->ambigPrev = 2;
	      uint32_t ambigIndex = _ambigPrevLists.length();
	      theState->prevState = ambigIndex;
	      _ambigPrevLists.addEmpty();
	      _ambigPrevLists[ambigIndex].append(curPrev);
	      _ambigPrevLists[ambigIndex].append(prevState);
	    }
	  }
	  break;

	default:
	  fprintf(stderr, "ERROR: got impossible ambigPrev value %d\n",
		  theState->ambigPrev);
	  break;
      }

      // indicate that the state has an arbitrary parent assignment if
      // any previous states that lead to it are of this class
      theState->arbitraryPar |= (1-isPI) * IVambigPar;

      // Sanity checks:
      assert(theState->unassigned == (fullUnassigned | ambig1Unassigned));
      // TODO: probably change this: have the assertion hold & update comment
      // This doesn't hold -- if two previous states lead to different parent
      // phase values, to avoid extreme amounts of ambiguity, we don't store all
      // alternatives here, but do store all possible previous states. Then
      // during back tracing, we follow all previous paths and ultimately find
      // the source of the ambiguity.
//      assert(theState->parentPhase == curParPhase ||
//	      (hetParIsDiff | (theState->ambigParPhase & (1 << curParPhase))));
    }

    if (i == 0 && numIter == 2) {
      // Update the various values as needed for the inverted phase in the next
      // iteration.

      // Children that have missing data should not be flipped
      // Also, for both parent het markers, heterozygous children's bits do
      // not get flipped: they're either unambigous and constrained by
      // <prevState->iv> or they're ambiguous and should remain the same to
      // stick with the convention used in <_stateHash>.
      uint64_t noFlipBits = childrenData[G_MISS] | (childrenData[G_HET] * isPI);
      uint64_t flipVal = _parBits[ hetParent ] & ~noFlipBits;
      // When the child was ambiguous previously and <isPI> == 0, we need to
      // flip the propagated <iv> values from the homozygous parent, so:
      uint8_t homozyParent = (1 - isPI) * (1 - hetParent);
      flipVal |= (1 - isPI) * (_parBits[ homozyParent ] & stdAmbigOnlyPrev);
      fullIV ^= flipVal;
      // Don't flip recombs when the previous state was ambiguous. Earlier
      // code removed recombs that could be addressed by flipping that ambiguous
      // <iv> value and that change applies equally to both phase possibilities
      // here.
      // Also see long comment in fixRecombFromAmbig() regarding ambig1 prev
      // states when <isPI>. The gist of it is that we sometimes keep the
      // recombination count for these at 1, otherwise we flip between 0 and1
      // via <ambig1PrevInfo>.
      // Finally, don't inadvertently introduce recombinations by flipping
      // values that aren't assigned in the previous state.
      recombs ^= flipVal & ~(stdAmbigOnlyPrev | (isPI * ambig1PrevInfo) |
								prevUnassigned);
      size_t oldNumRecombs = numRecombs;
      numRecombs = popcount(recombs);
      if ((isPI*IVambigPar) && theStateUpdated && oldNumRecombs == numRecombs) {
	theState->arbitraryPar = 1;
	break;
      }
      curParPhase = altPhaseType;

      // ambig1 values affect fullUnassigned. The value inverts depending on the
      // phase type at MT_FI markers. Note that the unassigned values are for
      // the homozygous parent and only for children that were ambig1 type
      // ambiguous at the previous marker.
      ambig1Unassigned ^= (1 - isPI) * ambig1PrevInfo;
      // <fullAmbig> doesn't change
    }
  }
}

// Given <iv> and <ambig>, first convert <iv> to the equivalent canonical value
// (see comment just inside method) and then do a hash table lookup for the
// State. If it does not exist, create it. Either way, return it to the caller.
State * Phaser::lookupState(const uint64_t iv, const uint64_t allAmbig,
			    const uint64_t unassigned) {
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

  // Get bits for children that are have an standard ambig value. Those with
  // ambig1 type values have similarities with unambiguous cases: their phase
  // is fixed relative to a previous value, and only after there's a
  // recombination does anything to do with the ambig1 value come in. Perhaps
  // more importantly, whereas with standard ambiguous values there are two
  // underlying IV values (two pairs of equivalent IVs 00 == 11 and 01 == 10 for
  // ambig IVs), the ambig1 type values have the standard four. A final way
  // to think about this is that missing data children have values propagated
  // from a previous marker, but states with equivalent IV values as defined
  // here need not think about the fact that the child is missing data: however
  // a state gets produced, if it's equivalent to some other state in terms of
  // IVs, it will produce the same number of recombinations downstream. That's
  // the case for these ambig1 type values, with the priviso of course that the
  // ambig1 type ambig bits must be equivalent for the equivalence to hold.
  uint64_t ambigStd = ((allAmbig & _parBits[1]) >> 1) * 3;

  uint64_t unambig = _parBits[2] - ambigStd;

  //////////////////////////////////////////////////////////////////////////
  // Get the canonical key value

  // First determine the <iv> value for the lowest order unambiguous child:
  int lowOrderChildBit = ffsll(unambig) - 1;
  // Markers that have all heterozygous children are totally ambiguous and not
  // considered at this stage, so there must be at least one unambiguous child:
  assert(lowOrderChildBit >= 0);

  // The genotype of the child tells us what bits we need to flip: the exact
  // bits that are assigned 1 need to be flipped in all unambiguous children.
  uint8_t flipType = (iv >> lowOrderChildBit) & 3;

  // Conveniently, we've got _flips indexed by the 4 possible flip types with
  // the values to flip assigned in each child:
  uint64_t lookupIV = iv ^ (_flips[flipType] & unambig);
  // And we've done something analogous for ambiguous bits; must flip them too:
  lookupIV ^= _ambigFlips[flipType] & ambigStd;

  //////////////////////////////////////////////////////////////////////////
  // Do the lookup

  iv_ambig_real theKey(lookupIV, allAmbig, unassigned);
  state_ht_iter it = _stateHash.find( &theKey );
  if (it == _stateHash.end()) {
    // need to create state
    State *newState = new State;
    newState->ambig = allAmbig;
    newState->unassigned = unassigned;
    newState->minRecomb = UINT16_MAX;
    iv_ambig newStateKey = new iv_ambig_real(lookupIV, allAmbig, unassigned);
    _stateHash[ newStateKey ] = newState;
    int curHMMIndex = _hmm.length() - 1;
    _hmm[curHMMIndex].append(newState);
    // caller will assign other necessary values
    return newState;
  }
  else {
    return it->second;
  }
}

// This method does two things:
// 1. It identifies states that are proveably suboptimal -- having more
// recombinations than are necessary to transition to them from at least one
// other state.
// 2. If all the states have <State::error> == 2, then all paths have an error
// in them and this indicator isn't necessary to track. In fact, doing so will
// interfere with decisions about states that are equivalent, so we clear the
// flag when all are 2.
void Phaser::rmBadStatesCheckErrorFlag(dynarray<State*> &curStates,
				       uint16_t minMaxRec[2], int numChildren) {
  bool allError2 = true; // assume all states have <error> == 2 initially
  int shiftToState = -1;

  // TODO: potential optimizations:
  // 1. If we know there are no added errors this time and that we cleared the
  // errors in the last iteration, could avoid checking for them
  // 2. If we know something about the <IV> of the minimum state, may be able
  // to throw out more states -- 2 * numChildren is potentially conservative

  if (minMaxRec[1] >= minMaxRec[0] + 2 * numChildren) {
    for(int i = 0; i < curStates.length(); i++) {
      if (curStates[i]->minRecomb >= minMaxRec[0] + 2 * numChildren) {
	// safe to discard state <i>
	if (shiftToState < 0) {
	  // Will bump states up, overwritting this index
	  shiftToState = i;
	}
      }
      else if (shiftToState >= 0) {
	curStates[shiftToState] = curStates[i];
	shiftToState++;
      }

      if (curStates[i]->error != 2)
	// one state that is either an error state (newly erroneous) or not an
	// error at all: can't clear error flag
	allError2 = false;
    }

    // Lastly must update for the fact that the list of states is now smaller
    if (shiftToState >= 0)
      curStates.resize(shiftToState);
  }
  else {
    // no states to remove; still must to check on clearing the error flag, but
    // if there's nothing to clear then we can optimize by returning
    for(int i = 0; i < curStates.length(); i++) {
      if (curStates[i]->error != 2)
	// one state that is either an error state (newly erroneous) or not an
	// error at all: can't clear error flag
	return;
    }
  }

  if (allError2) {
    // All states have an error status of 2, can set them all to 0: there's no
    // difference between them
    for(int i = 0; i < curStates.length(); i++) {
      curStates[i]->error = 0;
    }
  }
}

// TODO: go through and comment this and possibly refactor.
// Back traces and minimum recombinant phase using the states in <_hmm>.
void Phaser::backtrace(NuclearFamily *theFam) {
  int lastIndex = _hmm.length() - 1;

  uint32_t curStateIdx = findMinStates(_hmm[lastIndex]);
  uint32_t prevStateIdx = UINT32_MAX;
  // Number of recombinations in <curState> relative to <prevState> below
  uint8_t numRecombs; // TODO: optimization: do we want this?
  for(int hmmIndex = lastIndex; hmmIndex >= 0; hmmIndex--) {
    int curHmmIndex = hmmIndex; // not redundant: errors modify hmmIndex
    State *curState = _hmm[curHmmIndex][curStateIdx];

    // Analogous to <State::ambigParPhase> -- bits 0-3 indicate which parental
    // phase values can produce minimum recombinations.
    // We assign this based on the states indicated in <curIdxSet>
    uint8_t curAmbigParPhase = curState->ambigParPhase;
    uint8_t curAmbigParHet = curState->ambigParHet;
    uint8_t curArbitraryPar = curState->arbitraryPar;
    // Indicates which inheritance vector values differ among any ambiguous
    // states that are in <_curIdxSet>. The true IV value is ambiguous in that
    // case.
    uint64_t ivFlippable = 0;

    // Find all the possible parent phase types that are have equal and minimal
    // numbers of recombinations at this marker and append their previous states
    // to <_prevIdxSet>.
    for(state_set_iter it = _curIdxSet->begin(); it != _curIdxSet->end(); it++){
      uint32_t stateIdx = *it;
      if (stateIdx == curStateIdx)
	// don't treat the same state as <curState> as ambiguous
	continue;

      State *ambigState = _hmm[curHmmIndex][stateIdx];

      ivFlippable |= curState->iv ^ ambigState->iv;

      // binary indicator for whether the two states have different
      // <hetParent> values
      uint8_t hetParIsDiff;
      curAmbigParHet |= calcAmbigParHetBits(curState->hetParent,
					    ambigState->hetParent,
					    hetParIsDiff);
      curAmbigParHet |= ambigState->ambigParHet;

      // only track differences in parent phase assignments if the two
      // <hetParent> values are the same.  Otherwise the values of the
      // parent phase assignments are not comparable.
      curAmbigParPhase |= (1 - hetParIsDiff) *
		   ((1 << ambigState->parentPhase) | ambigState->ambigParPhase);

      curArbitraryPar |= ambigState->arbitraryPar;

      // Add the previous state index(es) to <_prevIdxSet> so long as the
      // ambiguous state has the same error status as the current state.
      // They are allowed in fact to have different error statuses (i.e., 0 and
      // 2) so long as one of them isn't 1: if one has a value of 1 and the
      // other doesn't, the previous state indexes they each reference are to
      // different markers.
      if (ambigState->error != curState->error &&
	  ((curState->error & 1) || (ambigState->error & 1)))
	// error values reference two markers previous not one
	// TODO: want to find a way to allow this? I think we just need one more
	//       set storing the ambiguous state indexes two back
	continue;

      int prevHmmIndex = hmmIndex - 1 - (curState->error & 1);
      if (!ambigState->ambigPrev) {
	uint32_t thePrevIdx = ambigState->prevState;
	_prevIdxSet->insert(thePrevIdx); // only one
	// Note: we call this here and in collectAmbigPrevIdxs() to facilitate
	// properly calculating <ivFlippable>. Without this, states that have
	// ambiguous IV values can be such that both <IV> values are uncertain
	// when only one is (or even that the opposite parent is uncertain
	// relative to the true one [at least in principle this could happen])
	propagateBackIV(ambigState,
			/*prevState=*/ _hmm[prevHmmIndex][thePrevIdx]);
      }
      else {
	uint32_t dontcare;
	collectAmbigPrevIdxs(ambigState, _hmm[prevHmmIndex], dontcare);
      }
    }

    // Remove any ambig par phase types that have equivalent phase to <curState>
    // Sometimes states differ in the IV values the propagate forward from the
    // previous markers, and the true ambiguity lies in those earlier markers
    curAmbigParPhase &= ~(1 << curState->parentPhase);

    // In the previous state, resolve ambiguous <iv> values and propagate
    // backward any <iv> values that were unassigned in that state
    if (hmmIndex - 1 >= 0) {
      // Note: curState->error == 2 means the path has an error in it. We only
      // need deal with curState->error == 1:
      if (curState->error & 1) {
	// Set error for immediately previous marker.
	uint64_t childrenData = _genos[curHmmIndex].second;
	uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
	theFam->setStatus(/*marker=*/ _hmmMarker[hmmIndex-1], PHASE_ERR_RECOMB,
			  _genos[hmmIndex-1].first, childrenData, missing);
	// <curState->prevState> references a state two indexes back
	deleteStates(_hmm[hmmIndex-1]);
	hmmIndex--;
      }

      if (!curState->ambigPrev) {
	// only one previous state
	prevStateIdx = curState->prevState;

	numRecombs = propagateBackIV(curState,
				  /*prevState=*/_hmm[hmmIndex-1][prevStateIdx]);
      }
      else {
	// multiple previous states: add all but the first to <_prevIdxSet>
	// the first will be the new previous state
	numRecombs = collectAmbigPrevIdxs(curState, _hmm[hmmIndex-1],
					  prevStateIdx);
      }
    }
    else
      numRecombs = 0;

    // Missing genotype value is 01, not any other, so:
    uint8_t parentData = _genos[curHmmIndex].first;
    uint8_t parMissing = (parentData & 5) & ~((parentData & 10) >> 1);
    uint64_t childrenData = _genos[curHmmIndex].second;
    uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
    theFam->setPhase(_hmmMarker[curHmmIndex], curState->iv,
		     curState->ambig & _parBits[1], missing, ivFlippable,
		     parMissing, curState->hetParent, curState->homParentGeno,
		     curState->parentPhase, numRecombs, curAmbigParHet,
		     curAmbigParPhase, curArbitraryPar);
    deleteStates(_hmm[curHmmIndex]);

    // ready for next iteration -- make prev into cur for states and state sets
    curStateIdx = prevStateIdx;
    state_idx_set *tmp = _curIdxSet;
    _curIdxSet = _prevIdxSet;
    _prevIdxSet = tmp;
    _prevIdxSet->clear();
  }
}

// Given a list of states, does a linear search to find the states with minimum
// numbers of recombinations. Returns the index of the first state encountered
// that is minimal and stores any other states that are tied for minimum in
// <_curIdxSet>.
uint32_t Phaser::findMinStates(dynarray<State*> &theStates) {
  uint16_t minLastMarker = UINT16_MAX;
  // The state we will designate as the minimum even as there may ambiguity
  // stored in <_stateIdxSet>
  uint32_t minStateIdx = UINT32_MAX;

  int numStates = theStates.length();
  for(int i = 0; i < numStates; i++) {
    uint16_t curRecomb = theStates[i]->minRecomb;
    if (curRecomb < minLastMarker) {
      // new minimum: clear any that were equivalent to the previous minimum
      _curIdxSet->clear();
      minLastMarker = curRecomb;
      minStateIdx = i;
    }
    else if (curRecomb == minLastMarker) {
      _curIdxSet->insert(i);
    }
  }

  return minStateIdx;
}

// Free/delete the states stored in <theStates> and clear it
void Phaser::deleteStates(dynarray<State*> &theStates) {
  int numStates = theStates.length();
  for(int i = 0; i < numStates; i++) {
    delete theStates[i];
  }
  theStates.clear();
}

// Given <curState> that is reachable with minimum recombinations via
// <prevState>, propagates back information about the inheritance vector values
// contained in <curState>. This assigns previously-<unassigned> bits in
// <prevState> and makes many previously ambiguous <iv> values unambiguous. The
// knowledge of the <iv> values at the later state allows the choice between the
// two possible <iv> values which are opposite one another and therefore differ
// in the transmitted homolog from both parents. Ambiguities only remain if
// there are recombinations on both sides of the ambiguous <iv> value (such as
// for a putative non-crossover/gene conversion) in which case the recombining
// parent is truly ambiguous.
uint8_t Phaser::propagateBackIV(State *curState, State *prevState) {
  // TODO: from some contexts we don't need to do quite this much work. Make
  // a separate function? May not be worth it.
  // First propagate backward any <iv> values that were unassigned
  // previously
  prevState->iv = (prevState->iv & ~prevState->unassigned) |
					(curState->iv & prevState->unassigned);

  // Note: ambiguities that remain in <curState> will not give information
  // about resolving such in <prevState>
  uint64_t ambigToResolve = prevState->ambig & ~curState->ambig;
  // Get both bits set for children that have the two different types of
  // ambiguities (standard and ambig1):
  uint64_t prevAnyAmbig = (ambigToResolve & _parBits[0]) * 3;
  uint64_t prevStdAmbig = ((ambigToResolve & _parBits[1]) >> 1) * 3;
  assert((prevAnyAmbig & prevStdAmbig) == prevStdAmbig);
  uint64_t prevAmbig1 = prevAnyAmbig - prevStdAmbig;

  uint64_t ambigRecombs = (curState->iv ^ prevState->iv) & prevAnyAmbig;
  // set both bits for children that inherit a recombination from the given
  // parent
  uint64_t parRecombs[2] = { (ambigRecombs & _parBits[0]) * 3,
			    ((ambigRecombs & _parBits[1]) >> 1) * 3 };

  // Invert the phase for ambiguous assignments that recombine on both homologs
  // relative to <curState>
  uint64_t bothRecomb = parRecombs[0] & parRecombs[1];
  prevState->iv ^= bothRecomb;

  // children that are ambiguous in <prevState> _and_ recombine on one homolog
  // relative to <curState> truly have ambiguous phase.
  // Ambig1 ambiguities are always resolved; the issue with them has to do with
  // counting of recombinations when there were undefined <iv> values before
  // encountering a PI state in which the relevant sample is heterozygous.
  // TODO: optimization -- test this:
//  uint64_t oneRecomb = parRecombs[0] ^ parRecombs[1];
//  prevState->ambig = prevStdAmbig & oneRecomb;
  // old less efficient version of updating <prevState->ambig>
  uint64_t noRecomb = ambigToResolve - (parRecombs[0] | parRecombs[1]);
  prevState->ambig &= ~(bothRecomb | noRecomb | prevAmbig1);

  // Return number of recombinations
  // Note: this number may be off by 1 for any ambig1 values. If we end up
  // flipping both bits for ambig1 values (via <bothRecomb> above), the <iv>
  // values will look like there are 0 recombinations relative to <prevState>,
  // but <curState->minRecomb> will encode 1 recombination.  This 1
  // recombination occurs earlier at the establishment of the ambig1 <iv>
  // value.
  return curState->minRecomb - prevState->minRecomb -
			  (curState->error & 1) * CmdLineOpts::max1MarkerRecomb;
}

// Assumes that curState->ambigPrev > 0 (i.e., there is some ambiguity).
// Inserts all indexes referenced in <curState->prevState> into
// <_prevIdxSet>. Also calls propagateBackIV() for all such previous states
// which is needed to properly calculate <ivFlippable> in backtrace().
// <prevStates> is the list of states at the previous marker referenced by
// <curState->prevState>.
//
// There are effectively two return values from this method. One is the
// parameter <thePrevStateIdx>, which is the previous state that will be
// considered the "main" state for assigning phase. The second is the actual
// return value which is the number of recombinations between <curState> and
// the previous state with index <thePrevStateIdx>.
uint8_t Phaser::collectAmbigPrevIdxs(State *curState,
				     dynarray<State*> &prevStates,
				     uint32_t &thePrevStateIdx) {
  const uint32_t mask = (1 << 6) - 1;
  uint8_t maxNumRecombs;

  switch (curState->ambigPrev) {
    // We assume the caller deals with this case:
//    case 0: // only one
//      _prevIdxSet->insert( curState->prevState );
//      break;

    case 1: // multiple previous states in <prevState>: add all
      {
	uint32_t prevs = curState->prevState;
	// Set the first index as the (default) previous state
	thePrevStateIdx = prevs & mask;
	_prevIdxSet->insert( thePrevStateIdx );
	prevs >>= 6;
	maxNumRecombs = propagateBackIV(curState, prevStates[ thePrevStateIdx]);
	for(int shift = 6; shift < 30; shift += 6) {
	  uint32_t curPrevIdx = prevs & mask;
	  _prevIdxSet->insert( curPrevIdx );
	  uint8_t curNumRecombs = propagateBackIV(curState,
						  prevStates[ curPrevIdx ]);
	  if (curNumRecombs > maxNumRecombs) {
	    maxNumRecombs = curNumRecombs;
	    thePrevStateIdx = curPrevIdx;
	  }
	  prevs >>= 6;
	  if (prevs == 0)
	    break; // got to last previous index
	}
      }
      break;

    case 2: // previous states in <_ambigPrevLists>: add all
      {
	dynarray<uint32_t> &prevs = _ambigPrevLists[curState->prevState];
	thePrevStateIdx = prevs[0];
	_prevIdxSet->insert( thePrevStateIdx );
	maxNumRecombs = propagateBackIV(curState, prevStates[ thePrevStateIdx]);
	for(int i = 1; i < prevs.length(); i++) {
	  uint32_t curPrevIdx = prevs[i];
	  _prevIdxSet->insert( curPrevIdx );
	  uint8_t curNumRecombs = propagateBackIV(curState,
						  prevStates[ curPrevIdx ]);
	  if (curNumRecombs > maxNumRecombs) {
	    maxNumRecombs = curNumRecombs;
	    thePrevStateIdx = curPrevIdx;
	  }
	}
      }
      break;

    default:
      fprintf(stderr, "ERROR: got impossible ambigPrev value %d\n",
	      curState->ambigPrev);
      exit(5);
      break;
  }

  return maxNumRecombs;
}
