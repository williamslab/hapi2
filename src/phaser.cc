// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include "phaser.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray< dynarray<State> > Phaser::_hmm;
dynarray<int> Phaser::_hmmMarker;
uint64_t Phaser::_allBitsSet;
uint64_t Phaser::_parBits[2];


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

    partialStates.clear();
  }
}

// Do initial setup of values used throughout phasing the current chromosome
void Phaser::init(int numChildren) {
  _hmm.clear();
  _hmmMarker.clear();

  // alternating bits set starting with bit 0 then bit 2, ...
  _allBitsSet = ~0ul; // initially: fewer depending on numChildren
  if (numChildren < 32)
    _allBitsSet &= (1 << (2*numChildren)) - 1;
  _parBits[0] = 0x5555555555555555 & _allBitsSet;
  _parBits[1] = 0xAAAAAAAAAAAAAAAA & _allBitsSet;
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

  // TODO: remove
  assert( 1 == (1 << G_HOM0));
  assert(12 == ((1 << G_HOM1) | (1 << G_HET)));
  assert(14 == ((1 << G_HOM1) | (1 << G_HET) | (1 << G_MISS)));
  assert( 5 == ((1 << G_HOM0) | (1 << G_HET)));
  assert( 7 == ((1 << G_HOM0) | (1 << G_HET) | (1 << G_MISS)));
  assert( 9 == ((1 << G_HOM0) | (1 << G_HOM1)));
  assert( 8 == (1 << G_HOM1));
  assert( 4 == (1 << G_HET));
  assert( 2 == (1 << G_MISS));
  assert( 3 == ((1 << G_HOM0) | (1 << G_MISS)));
  assert(10 == ((1 << G_HOM1) | (1 << G_MISS)));
  assert( 6 == ((1 << G_HET) | (1 << G_MISS)));
  assert(11 == ((1 << G_HOM0) | (1 << G_HOM1) | (1 << G_MISS)));
  assert(13 == ((1 << G_HOM0) | (1 << G_HOM1) | (1 << G_HET)));
  assert(1 == G_MISS);
  assert(2 == G_HET);
  assert(0 == G_HOM0);
  assert(3 == G_HOM1);

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
			       uint64_t childrenData[5]) {
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
				  uint8_t parentData,uint64_t childrenData[5]) {
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
    // No information about transmission from homozygous parent and for
    // children that are missing data:
    newState.unassigned = _parBits[homPar] | childrenData[ G_MISS ];
  }
}

// Helper for makePartialStates(): applicable to partly informative markers (or
// the states that correspond to this possibility when the parent's genotypes
// aren't fully known)
void Phaser::makePartialPIStates(dynarray<State> &partialStates,
				 uint8_t parentData, uint64_t childrenData[5]) {
  // Add space for a new (partial) state. The use of addEmpty() allows us to
  // get a reference to the state stored here and enables direct writing to
  // the stored state (no need to copy values in)
  int newIndex = partialStates.length();
  partialStates.addEmpty();
  State &newState = partialStates[newIndex];

  // Initially assumes the phase is assigned such that allele 0 is on
  // haplotype 0 and allele 1 is on haplotype 1. In particular, the genotypes
  // of the children (excluding the missing data children) already encode the
  // inheritance vector assuming this phase assignment in the parents
  // TODO: want to store a value that is the inverse of missing data?
  newState.iv = childrenData[ 4 ] & ( ~childrenData[ G_MISS ] );
  // Have full transmission information except for children with missing data:
  newState.unassigned = childrenData[ G_MISS ];
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
