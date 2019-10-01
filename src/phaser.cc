// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <bitset>
#include <genetio/util.h>
#include <float.h>
#include "phaser.h"
#include "cmdlineopts.h"

////////////////////////////////////////////////////////////////////////////////
// declare/initialize static members
dynarray< dynarray<State*> >    Phaser::_hmm;
dynarray<int>                   Phaser::_hmmMarker;
dynarray< dynarray<uint32_t> >  Phaser::_ambigPrevLists;
dynarray< std::pair<uint8_t,uint64_t> > Phaser::_genos;
Phaser::state_ht                Phaser::_stateHash;
Phaser::BT_state_set            *Phaser::_curBTAmbigSet;
Phaser::BT_state_set            *Phaser::_prevBTAmbigSet;
PhaseMethod                     Phaser::_phaseMethod;
uint64_t Phaser::_parBits[3];
uint64_t Phaser::_flips[4];
uint64_t Phaser::_ambigFlips[4];
int      Phaser::_lastInformMarker;
int      Phaser::_curMarker;

extern uint8_t swap01Phase[4][16];
extern uint8_t swap2Phase[4][16];

void Phaser::run(NuclearFamily *theFam, int chrIdx, FILE *log) {
  // TODO! remove this and make it a command line option
  _phaseMethod = PHASE_MINREC;

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

  // bit vector with bit 0: dad, bit 1: mom
  uint8_t missingPar = 0;
  if (!theFam->_parents->first->hasData())
    missingPar += 1;
  if (!theFam->_parents->second->hasData())
    missingPar += 2;

  dynarray<State> partialStates;
  _lastInformMarker = -1;

  // build states/phase each marker
  for(_curMarker = firstMarker; _curMarker <= lastMarker; _curMarker++) {
    uint8_t parentData, parentGenoTypes, homParGeno, childGenoTypes;
    // Each index corresponds to a genotype (see the Geno enumerated type).
    // The two bits corresponding to each child are set to 1 for the index
    // of its genotype.
    // Index 4 (the last value) stores the raw genotype data in PLINK format
    uint64_t childrenData[5];
    int numMissChildren = 0; // how many children are missing data?

    ///////////////////////////////////////////////////////////////////////////
    // Step 0: get the data for this marker
    getFamilyData(theFam, _curMarker, parentData, parentGenoTypes, childrenData,
		  childGenoTypes, numMissChildren);

    ///////////////////////////////////////////////////////////////////////////
    // Step 1: Determine marker type and check for Mendelian errors
    int mt = getMarkerType(parentGenoTypes, childGenoTypes, homParGeno);
    assert(mt > 0 && (mt & ~((1 << MT_N_TYPES) -1) ) == 0);

    if (CmdLineOpts::verbose)
      printMarkerType(mt, log);

    if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
      // should only be one of the above:
      assert(mt == (1 << MT_ERROR) || mt == (1 << MT_AMBIG));
      switch(mt) {
	case 1 << MT_ERROR:
	  theFam->setStatus(_curMarker, PHASE_ERROR, parentData,
			    childrenData[4], childrenData[G_MISS] &_parBits[0]);
	  break;
	case 1 << MT_AMBIG:
	  theFam->setStatus(_curMarker, PHASE_AMBIG, parentData,
			    childrenData[4], childrenData[G_MISS] &_parBits[0]);
	  break;
      }
      continue;
    }

    if (mt & (1 << MT_UN)) {
      uint8_t swapHetChildren = 0;
      if ((parentData & 3) == G_HOM1 || ((parentData >> 2) & 3) == G_HOM0)
	// The default way to print phase of heterozygous children at
	// uninformative markers is with allele 0 transmitted by dad and 1
	// transmitted by mom. When Dad is homozygous for allele 1 or (when
	// he is missing) Mom is homozygous for allele 0, this order needs to
	// be swapped:
	swapHetChildren = 1;
      theFam->setUninform(_curMarker, parentData, childrenData[4],
			  childrenData[G_MISS] & _parBits[0], homParGeno,
			  swapHetChildren);
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
    makeFullStates(partialStates, _curMarker, firstMarker, childrenData,
		   missingPar, numChildren, numMissChildren);

    // Clean up for next marker
    partialStates.clear();
    for (state_ht_iter it = _stateHash.begin(); it != _stateHash.end(); it++) {
      delete it->first;
    }
    _stateHash.clear_no_resize();
    _lastInformMarker = _curMarker;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Step 4: HMM analysis finished. Back trace and assign phase!
  if (CmdLineOpts::verbose) {
    fprintf(log, "  Back tracing... ");
    fflush(log);
  }
  backtrace(theFam, missingPar, firstMarker, lastMarker);
  if (CmdLineOpts::verbose) {
    fprintf(log, "done\n");
    fflush(log);
  }
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
			   uint64_t childrenData[5], uint8_t &childGenoTypes,
			   int &numMissChildren) {
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
    // following is 1 when the genotype is 1 (missing), and 0 otherwise
    numMissChildren += (curChildData & 1) & (~curChildData >> 1);
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
	// both parents homozygous: uninformative marker
	return 1 << MT_UN;
      break;
    case (1 << G_HOM1):  // both homozygous for 1
      if (childGenoTypes & ((1 << G_HOM0) | (1 << G_HET)) )
	// Mendelian error -- invalid bits set: only hom for 1 and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homozygous: uninformative marker
	return 1 << MT_UN;
      break;
    case ((1 << G_HOM0) | (1 << G_HOM1)):  // one homozygous for 0, other 1
      if (childGenoTypes & ((1 << G_HOM0) | (1 << G_HOM1)))
	// Mendelian error -- invalid bits set: only het and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homozygous: uninformative marker
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
//	  homParGeno = G_MISS; // not really sure of other parent's genotype
	  // though the parent's genotype is potentially heterozygous, we have
	  // no evidence that that is the case here. As such, in the caller, if
	  // the marker type includes MT_UN, we set the status to uninformative.
	  // We also wish to impute, as much as possible, any missing genotypes
	  // in the parent. Here, we know the missing parent must have
	  // transmitted allele 0 to the children. As such, we impute it as
	  // carrying that allele, and later (in backtrace()), we determine
	  // whether both haplotypes were transmitted or only one and we print
	  // the known information accordingly. IMPUTING
	  homParGeno = G_HOM0;
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
//	  homParGeno = G_MISS; // not really sure of other parent's genotype
	  // see comment above marked IMPUTING
	  homParGeno = G_HOM1;
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
	  // gives inverse homozygous genotype relative to <missingType>:
	  // see comment above marked IMPUTING
	  homParGeno = G_HOM1 - missingType;
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

// For use with verbose mode, prints the marker type to the log
void Phaser::printMarkerType(int mt, FILE *log) {
  fprintf(log, "    Marker %d:", _curMarker);
  if (mt == 1 << MT_ERROR)
    fprintf(log, " Mendelian error\n");
  else if (mt == 1 << MT_AMBIG)
    fprintf(log, " Ambiguous (all individuals heterozyogus or missing)\n");
  else if (mt & (1 << MT_UN))
    fprintf(log, " Uninformative (neither parent heterozygous)\n");
  else {
    int FI1 = 1 << MT_FI_1;
    int PI = 1 << MT_PI;
    int both = FI1 | PI;
    if (mt & both)
      fprintf(log, " One OR both parents heterozygous\n");
    else if (mt & FI1)
      fprintf(log, " One parent heterozygous\n");
    else if (mt & PI)
      fprintf(log, " Both parents heteterozygous\n");
  }
  fflush(log);
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
  newState.homParentGeno = G_MISS;
  // Won't assign <parentPhase> field as all partial states have a value of
  // 0 here and it's never used
  //newState.parentPhase = 0;

  assert((newState.iv & newState.unassigned) == 0ul);
}

// Using the states at the previous marker where present and <partialStates>,
// generates all full states. Stores these as the last value in <_hmm> and
// stores the marker index <marker> that these states apply to in <_hmmMarker>.
void Phaser::makeFullStates(const dynarray<State> &partialStates, int marker,
			    int firstMarker, const uint64_t childrenData[5],
			    uint8_t missingPar, int numChildren,
			    int numMissChildren) {
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
		 && missingPar == 3)
	// For the very first marker, if we're not sure which parent is
	// heterozygous, and both are missing data, arbitrarily pick parent
	// 0 as such. Otherwise, since the two states are equivalent but
	// with opposite parent labels, there will be two equal paths
	// through the HMM.
	// TODO: document this behavior
	continue;
      State *newState = new State(partialStates[i]);
      newState->ambigPrev = 0;
      newState->minRecomb = 0;
      uint8_t isPI = newState->hetParent >> 1;
      // PI states are informative for both parents (both are heteterozygous),
      // so there are 0 informative markers since seeing the non-het parent if
      // <newState> is PI: (1 - isPI) below.
      // If the state is not PI, AND if the homozygous parent
      // (1 - newState->hetParent) is missing, we count the number of markers
      // since the beginning of the chromosome:
      bool homParMissing =
		  (1 - isPI) && ((1 << (1 - newState->hetParent)) & missingPar);
      newState->numMarkersSinceNonHetPar =
			(1 - isPI) * homParMissing * (marker - firstMarker + 1);
      // PI states are heterozygous for both parents, so count is 1 when <isPI>
      // and zero otherwise (is a one het par marker in that case)
      newState->numMarkersSinceOneHetPar = isPI;
      // initially probability of log(1) == 0; note that in principle this
      // should be log(1/N), where N is the number of initial states, but really
      // all that encodes is that they have equal probability, so setting this
      // to 0 suffices:
      newState->maxLikelihood = 0;
      newState->maxPrevRecomb = 0;
      newState->ambigParHet = 1 << newState->hetParent;
      newState->parentPhase = 0;
      // 1 << 0, i.e., 1 << newState->parentPhase
      newState->ambigParPhase = (1 << 0) << (2 * newState->hetParent);
      // arbitrary if both parents are missing:
      newState->arbitraryPar = missingPar == 3;
      newState->error = 0;
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
    uint8_t IVambigPar = isIVambigPar(prevState, missingPar);

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
		    IVambigPar, missingPar,
		    /*numDataChildren=*/numChildren - numMissChildren,
		    marker - _hmmMarker[ prevHMMIndex ]);
    }
  }

  // Introduce error states in prevStates if doing so saves at least
  // <CmdLineOpts::max1MarkerRecomb> recombinations
  // TODO: may be able to optimize. If there a states at both prev and cur that
  // introduce 0 recombinations relative to the most recent one, can we be
  // certain that the marker won't need error states?
  // TODO! (3) likelihood-based error detection:
  if (CmdLineOpts::max1MarkerRecomb > 0 && prevHMMIndex - 1 >= 0) {
    dynarray<State*> &back2States = _hmm[prevHMMIndex - 1];

    int64_t numBack2 = back2States.length();
    for(int64_t back2Idx = 0; back2Idx < numBack2; back2Idx++) {
      State *back2State = back2States[back2Idx];

      // See comment above the isIVambigPar() method. We deal with
      // <IVambigPar> == 1 below and in updateStates()
      uint8_t IVambigPar = isIVambigPar(back2State, missingPar);

      for(int curIdx = 0; curIdx < numPartial; curIdx++) {
	const State &curPartial = partialStates[curIdx];
	if (IVambigPar && curPartial.hetParent == 1)
	  continue; // See comment above in equivalent non-error code
	// Set back2Idx negative and shift by 1 so that the index 0 is also
	// negative
	mapPrevToFull(back2State, /*prevIdx=negative=>error*/ -back2Idx-1,
		      curPartial, minMaxRec, childrenData, IVambigPar,
		      missingPar,
		      /*numDataChildren=*/numChildren - numMissChildren,
		      marker - _hmmMarker[ prevHMMIndex - 1 ]);
      }
    }
  }

  rmBadStatesCheckErrorFlag(_hmm[ prevHMMIndex+1 ], minMaxRec, numChildren);
}

// If the IV values transmitted to each child by the two parents are either all
// identical or all opposite, then when <missingPar> == 3 (both missing), which
// parent is heterozygous (for FI markers) will be ambiguous. Additionally, a
// recombination at a PI state that follows these ambiguous markers will have
// two parent phase assignments (opposite each other) that will produce
// equivalent numbers of recombinations.
uint8_t Phaser::isIVambigPar(const State *state, uint8_t missingPar) {
  uint64_t IVparDiff = (state->iv & _parBits[0]) ^
					      ((state->iv & _parBits[1]) >> 1);
  // Omit differences at standard ambiguous positions (not only bit 0 in each
  // child is set if there's a difference)
  uint64_t stdAmbig = (state->ambig & _parBits[1]) >> 1;
  IVparDiff &= ~stdAmbig;
  if (missingPar == 3 && (IVparDiff == 0 ||
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
			   const uint64_t childrenData[5], uint8_t IVambigPar,
			   uint8_t missingPar, int numDataChildren,
			   int numMarkersSincePrev) {
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

  // The following gives 1 for <hetParent> == 2, 0 for <hetParent> == 0,1
  uint8_t isPI = curPartial.hetParent >> 1;

  // (2) Check for penalties and do associated book keeping. See long comment
  //     above checkPenalty()
  // stores the IV for the relevant parent having only one haplotype transmitted
  // enables checks for whether the penalty should be applied in updateStates()
  uint64_t penalty = 0;
  // For determining whether to trigger the two types of penalties:
  // Number of markers since the last one heterozygous for whichever parent is
  // homozygous here (for non-PI states)
  int16_t numMarkersSinceNonHetPar = 0;
  // Number of markers since the last one that is heterozygous for only one
  // parent
  int16_t numMarkersSinceOneHetPar = 0;
  if (missingPar > 0) {
    checkPenalty(prevState, curPartial, isPI, missingPar, numMarkersSincePrev,
		 penalty, numMarkersSinceNonHetPar, numMarkersSinceOneHetPar);
  }

  // (3) if curPartial is partly informative, for heterozygous children:
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

  // These variables are used to determine whether we need to explore all the
  // different possible phase types for parents. Typically the answer is yes,
  // but at the beginning of the chromosome, if an informative marker for one
  // of the parents hasn't yet been seen, then the two possible phase types for
  // the first marker that is informative for that parent will not in fact
  // differ in their numbers of recombinations. To avoid indicating these as
  // ambiguous, we set these values and check them in updateStates() and below
  // for PI states.
  bool hetParentUndefined, oneParentUndefined;

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

  // (4) As needed, remove apparent recombinations from <iv> values that were
  // ambiguous in the previous state
  // Which children were ambiguous in the previous state but not here?
  // Also various values relating to ambig1 type ambiguous values
  uint64_t stdAmbigOnlyPrev, ambig1PrevInfo, ambig1Unassigned;
  fixRecombFromAmbig(fullIV, recombs, parRecombs, isPI,
		     /*ambigOnlyPrev=*/ prevState->ambig & ~fullAmbig,
		     curPartial.hetParent, stdAmbigOnlyPrev, ambig1PrevInfo,
		     ambig1Unassigned);

  // (5) Look up or create a full state with equivalent <iv> and <ambig> values
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
	       penalty, prevState, stdAmbigOnlyPrev, ambig1PrevInfo,
	       curPartial.hetParent, curPartial.homParentGeno,
	       /*initParPhase=default phase=*/ 0, altPhaseType, prevIdx,
	       IVambigPar, minMaxRec, hetParentUndefined, childrenData,
	       numDataChildren, numMarkersSinceNonHetPar,
	       numMarkersSinceOneHetPar);

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
		 penalty, prevState, stdAmbigOnlyPrev,
		 ambig1PrevInfo, /*curPartial.hetParent=*/2,
		 /*homParentGeno=*/G_MISS, /*initParPhase=parent 0 flip=*/1,
		 /*altPhaseType=parent 1 flip=*/ 2, prevIdx,
		 IVambigPar, minMaxRec, /*hetParentUndefined=*/ false,
		 childrenData, numDataChildren, numMarkersSinceNonHetPar,
		 numMarkersSinceOneHetPar);
  }
}

// For dealing with regions in which a parent for which we do not have data
// transmitted only one of their two haplotypes. These will be silent and
// manifest as a long stretch of markers wherein only one parent is
// heterozygous. To ensure that we find the right IV when we eventually do
// encounter an informative marker for this parent, we detect this scenario
// (including requiring the IV to be one recombination away from having only
// one haplotype at the most recent informative marker for this parent) and
// impose a penalty for IVs that keep the same transmission status for the
// child that had the outlier haplotype at the most recent informative marker.
//
// The penalty occurs when transitioning to a state that is heterozygous for
// the indicated parent. We will be performing such a transition whenever the
// current state's heterozygous parent is different from the previous state
// (which was, for the <numMarkersSinceNonHetPar> value to be large,
// necessarily heterozygous for the other parent only).
//
// A second and related case occurs when we have data for one parent and the
// inheritance vector gets switched such that the transmitted haplotypes for
// the parent with data are associated with the other parent. In this case,
// sites that are heterozygous for the data parent will be heterozygous for
// both parents (since the children's alleles will be misphased). We penalize
// such both parent het states after a sufficient number of them occurs in
// sequence.
void Phaser::checkPenalty(const State *prevState, const State &curPartial,
			  uint8_t isPI, uint8_t missingPar,
			  int numMarkersSincePrev, uint64_t &penalty,
			  int16_t &numMarkersSinceNonHetPar,
			  int16_t &numMarkersSinceOneHetPar) {
  bool haveAllButOne = false; // all but one haplotype?
  // want to analyze parent 0 if bit 0 in missingPar is 1; otherwise start from
  // parent 1 (result of this is either 0 or 1):
  int firstP = 1 - (missingPar & 1);
  // want to analyze parent 2 if bit 1 in missingPar is 1; otherwise shouldn't
  // iterate to 1 (result of this is either 1 or 2):
  int limitP = (missingPar >> 1) + 1;
  for(int p = firstP; !isPI && !haveAllButOne && p < limitP; p++) {
    if (curPartial.hetParent == p)
      // parent is heterozygous in this state: no need to check
      continue;

    uint64_t allButOneIV = 0;

    // The inheritance vector with all but one haplotype being the same:
    uint64_t parHap = prevState->iv & _parBits[p];
    uint64_t oppParHap = parHap ^ _parBits[p];

    if (prevState->unassigned & _parBits[p]) {   // Some unassigned bits for <p>
      uint64_t assignedParHap = parHap & ~prevState->unassigned;
      // If all assigned bits are 0 or 1:
      if (assignedParHap == 0 ||
	  assignedParHap == (_parBits[p] & ~prevState->unassigned)) {
	haveAllButOne = true;
	allButOneIV = prevState->unassigned & _parBits[p];
      }
    }
    else {   // All bits assigned for <p>
      // Only a single bit set in parHap or oppParHap?
      // This corresponds to either all but one of the children or all children
      // except one receiving the same haplotype.
      if (parHap && !(parHap & (parHap - 1))) {
	assert(!haveAllButOne); // if violated, we need to track two penalties
	haveAllButOne = true;
	allButOneIV = parHap;
      }
      else if (oppParHap && !(oppParHap & (oppParHap - 1))) {
	assert(!haveAllButOne); // if violated, we need to track two penalties
	haveAllButOne = true;
	allButOneIV = oppParHap;
      }
    }

    if (!haveAllButOne)
      // want numMarkersSinceNonhetPar == 0 (which it is in the caller) if this
      // is true for both parents
      continue;

    if (prevState->numMarkersSinceNonHetPar < 0) {
      // have already applied the penalty previously; check whether the
      // heterozygous parent matches the previous state
      if (curPartial.hetParent == prevState->hetParent ||
	  ((prevState->ambigParHet >> curPartial.hetParent) & 1) == 1) {
	// same parent is heterozygous, possibly apply a second penalty
	numMarkersSinceNonHetPar = prevState->numMarkersSinceNonHetPar -
							    numMarkersSincePrev;
	if (numMarkersSinceNonHetPar < -CmdLineOpts::oneHapTransThreshold &&
	    prevState->numMarkersSinceNonHetPar >=
					 -CmdLineOpts::oneHapTransThreshold) {
	  // Has been 2x threshold distance: apply another penalty (which should
	  // force a recombination)
	  penalty = allButOneIV;
	}
      }
      else {
	// penalty was applied for a different heterozygous parent
	if (curPartial.hetParent < 2 && (curPartial.hetParent & missingPar))
	  // one parent is heterozygous, but different than the previous state;
	  // reset the count:
	  numMarkersSinceNonHetPar = numMarkersSincePrev;
      }
    }
    else {
      if (prevState->numMarkersSinceNonHetPar + numMarkersSincePrev >
					    CmdLineOpts::oneHapTransThreshold) {
	// above threshold => apply penalty:
	penalty = allButOneIV;
      }

      if (curPartial.hetParent < 2) {
	if (curPartial.hetParent == prevState->hetParent ||
	    ((prevState->ambigParHet >> curPartial.hetParent) & 1) == 1) {
	  if (penalty)
	    // penalty applied: make count negative
	    numMarkersSinceNonHetPar = -(prevState->numMarkersSinceNonHetPar +
					 numMarkersSincePrev -
					 CmdLineOpts::oneHapTransThreshold);
	  else
	    // accumulate markers that we have encountered since the last
	    // heterozygote:
	    numMarkersSinceNonHetPar = prevState->numMarkersSinceNonHetPar +
							    numMarkersSincePrev;
	}
	else if (curPartial.hetParent & missingPar) {
	  // one parent is heterozygous, but different than the previous state;
	  // reset the count:
	  numMarkersSinceNonHetPar = numMarkersSincePrev;
	}
      }
      // else: numMarkersSinceNonHetPar == 0 (assigned in caller)
    }
  }

  if (isPI && missingPar < 3 && prevState->numMarkersSinceOneHetPar >= 0) {
    // Case where we have data for one parent; attempt to detect switch in
    // which parent is which

    // Limit the impact of long stretches of uninformative markers on this
    // penalty: count a stretch as no more than 10 markers:
    numMarkersSinceOneHetPar = prevState->numMarkersSinceOneHetPar +
						  min(numMarkersSincePrev, 10);
    // NOTE: it may not be necessary, and we currently don't check, but
    //       this penalty should really only apply if the commented out
    //       conditional below holds. It's complicated by ambiguous IVs
    //       and the fact that the four phase options for this PI state
    //       behave differently with respect to this check
//    uint8_t dataPar = 2 - missingPar; // the parent with data
//    uint8_t dataParIV = fullIV & _parBits[dataPar];
//    if (!( dataParIV == _parBits[ dataPar ] || dataParIV == 0 )) {
//      // If the data parent transmitted the same haplotype to all children,
//      // then all sites where the children are heterozygous will be PI, so
//      // there shouldn't be a penalty. Here, we've inverted those cases and
//      // there should be a penalty:
//    }
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
  // (4) New ambig1 children:
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
			  uint64_t recombs, uint64_t penalty,
			  const State *prevState, uint64_t stdAmbigOnlyPrev,
			  uint64_t ambig1PrevInfo, uint8_t hetParent,
			  uint8_t homParentGeno, uint8_t initParPhase,
			  uint8_t altPhaseType, int64_t prevIndex,
			  uint8_t IVambigPar, uint16_t minMaxRec[2],
			  bool hetParentUndefined,
			  const uint64_t childrenData[5], int numDataChildren,
			  int16_t numMarkersSinceNonHetPar,
			  int16_t numMarkersSinceOneHetPar) {
  // How many iterations of the loop? See various comments below.
  int numIter = 2;

  // How many recombinations for the initial phase assignment?
  // Element 0 stores total;
  // When needed (below), element 1 stores the recombination count for parent 1
  size_t numRecombs[2];
  numRecombs[0] = popcount(recombs);
  uint8_t curParPhase = initParPhase;

  // Is the penalty due to the fact that the IV values corresponding to a
  // parent are unassigned?
  uint64_t unassignedPenalty = 0;
  if (fullUnassigned && (penalty & fullUnassigned) == fullUnassigned) {
    assert(penalty == fullUnassigned);
    unassignedPenalty = penalty;
    penalty = 0;
  }

  // Note: (hetParent >> 1) == 1 iff hetParent == 2. It is 0 for the other
  //       values.
  uint8_t isPI = hetParent >> 1;

  // Is it possible to swap the phase of this state and obtain the same
  // number of recombinations from the previous state? Decided below
  uint8_t ambigLocal = 0;

  // For cases where a missing data parent transmitted only a single haplotype:
  // (which is only detectable by examining the number of states where that
  // parent is not inferred to be heterozygous)
  //
  // how many penalty "recombinations" should we apply?
  int8_t penaltyCount = 0;
  const uint8_t THE_PENALTY = 1; // if there is a penalty
  if (penalty) { // are the IV values in the penalty scenario?
    for(int p = 0; p < 2; p++) {
      if ((~recombs) & penalty & _parBits[p])
	// no recombination in the one child that differs for a missing data
	// parent that hasn't had an informative marker for a long stretch:
	// apply penalty
	penaltyCount += THE_PENALTY;
    }
  }

  // For maximum likelihood phasing, get the relevant (log) probabilities and
  // the likelihood of transitioning to this state given the current parent
  // phase configuration:
  float localLikehood = -FLT_MAX, mapDist;
  float noRecombProb = -FLT_MAX, recombProb = -FLT_MAX;
  if (_phaseMethod == PHASE_MAXLIKE) {
    // TODO! this (currently non-functional) code doesn't account for the
    //       missing parent no het penalty and/or the error penalty
    mapDist = Marker::getMarker(_curMarker)->getMapPos() -
			 Marker::getMarker(_lastInformMarker)->getMapPos();
    noRecombProb = -mapDist;
    // TODO: optimization -- store table with log probabilities?
    recombProb = log(mapDist) + noRecombProb;

    if (isPI) {
      numRecombs[1] = popcount(recombs & _parBits[1]);
      // TODO! (2) do parent-specific calculations
      localLikehood = numRecombs[0] * recombProb +
			      (numDataChildren - numRecombs[0]) * noRecombProb;
    }
    else
      localLikehood = numRecombs[0] * recombProb +
			      (numDataChildren - numRecombs[0]) * noRecombProb;
  }

  // First, decide how many iterations of the loop below. If the opposite
  // phase assignment yields an equivalent state, we only need to loop once.
  // Only way the opposite phase assignment is equivalent is if:
  // (1) There is no missing data for children that have an assigned IV (we
  // propagate <iv> values from the previous marker for missing data children
  // so changing the parent's phase at the current marker will modify some
  // children's <iv> values but not the missing data ones that have an assigned
  // IV value).
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
  else if ( ((1-isPI) &&
	     (childrenData[G_MISS] & _parBits[hetParent] & ~fullUnassigned) == 0
	     && stdAmbigOnlyPrev == 0 && ambig1Unassigned == 0) ||
       (isPI && (stdAmbigOnlyPrev & childrenData[G_HET]) == 0
	     && ((childrenData[G_MISS] & ~fullUnassigned) |
					  childrenData[G_HET]) == fullAmbig)) {
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
							prevState->unassigned);
    size_t curCount[2]; // current count of number of recombinations
    curCount[0] = popcount(recombs);

    int8_t curPenalty = 0;
    if (penalty) { // are the IV values in the penalty scenario?
      for(int p = 0; p < 2; p++) {
	if ((~recombs) & penalty & _parBits[p])
	  // no recombination in the one child that differs for a missing data
	  // parent that hasn't had an informative marker for a long stretch:
	  // apply penalty
	  curPenalty += THE_PENALTY;
      }
    }

    if (_phaseMethod == PHASE_MINREC) { // minimum recombinant
      if (curCount[0] + curPenalty < numRecombs[0] + penaltyCount) {
	numRecombs[0] = curCount[0];
	curParPhase = altPhaseType;
	penaltyCount = curPenalty;
	fullIV ^= _parBits[ hetParent ] & ~fullAmbig;
	// No need to execute this as we ensured that
	// (1-isPI) * ambig1PrevInfo == 0
//	ambig1Unassigned ^= (1 - isPI) * ambig1PrevInfo;
      }
      else if (curCount[0] + curPenalty == numRecombs[0] + penaltyCount) {
	// Either of the two phase types give the same number of recombinations.
	// Will indicate this in the state below.
	ambigLocal = 1;
      }
    }
    else { // maximum likelihood
      // TODO! this (currently non-functional) code doesn't account for the
      //       missing parent no het penalty and/or the error penalty
      float curLocLikehood;
      if (isPI) {
	curCount[1] = popcount(recombs & _parBits[1]);
	// TODO! (2) do parent-specific calculations
	curLocLikehood = curCount[0] * recombProb +
				(numDataChildren - curCount[0]) * noRecombProb;
      }
      else
	curLocLikehood = curCount[0] * recombProb +
				(numDataChildren - curCount[0]) * noRecombProb;

      if (curLocLikehood > localLikehood) {
	localLikehood = curLocLikehood;
	curParPhase = altPhaseType;
	fullIV ^= _parBits[ hetParent ] & ~fullAmbig;
      }
      else if (curLocLikehood == localLikehood) {
	// Either of the two phase types give the same likelihood.
	// Will indicate this in the state below.
	ambigLocal = 1;
      }
    }
  }

  // When there is a penalty due to many markers being unassigned for a given
  // parent, we will generate two states: one "standard" one that preserves
  // the unassigned IV values -- and includes a penalty.
  // The second is a state that has all the IV values for the unassigned
  // parent the same -- thus being a scenario with only one haplotype
  // transmitted by the parent.
  int nUpenaltyIters = 1;
  if (unassignedPenalty)
    nUpenaltyIters = 2;

  State *lastState = NULL;

  // Examine the 1 or 2 needed phase assignments
  for(int i = 0; i < numIter; i++) {

    // State(s) for this iteration (over variable <i>). There's usually only
    // one of these (index 0), but if uPenaltyIters == 2, there will be two
    State *theState[2];

    // Does the current previous state lead to minimum recombinations for
    // <theState>? (Could be either a new minimum or an ambiguous one)
    // Assigned below
    bool theStateUpdated[2];

    for(int upenaltyIter = 0; upenaltyIter < nUpenaltyIters; upenaltyIter++) {

      uint64_t iterUnassigned = fullUnassigned;
      uint64_t iterIV = fullIV;
      if (upenaltyIter == 1) {
	iterUnassigned = 0;
	// Which parent transmitted only one haplotype? Since fullUnassigned
	// must be only set for one of the parents, if the following boolean is
	// true (1), the one hap transmission parent is parent 1. Otherwise,
	// it's parent 0.
	bool oneHapParent = fullUnassigned & _parBits[1];
	// for the assigned bits, are any non-zero? If so, we'll make the
	// IV with all children assigned the same value all 1s for those
	// children; otherwise all 0s:
	bool anyNonZero = iterIV & _parBits[ oneHapParent ];
	iterIV |= anyNonZero * fullUnassigned;
      }

      // Which bit corresponds to the lowest order unambiguous child?
      int lowOrderChildBit;
      // First look up or create a full state with equivalent <iv> and <ambig>
      // values to <fullIV> (here <iterIV>) and <fullAmbig>
      theState[upenaltyIter] = lookupState(iterIV, fullAmbig,
					   iterUnassigned | ambig1Unassigned,
					   lowOrderChildBit);
      // TODO: test optimization: add a check for prevState->minRecomb >
      // theState.minRecomb (or analogous comparison for likelihood)

      if (upenaltyIter == 0) {
	// only need check this for standard state (i.e., <upenaltyIter> == 0);
	// the unassigned penalty state could be the same between the two.

	// Should never map to the same state as the previous iteration over <i>
	// (this should be caught in the code above that sets numIter = 1)
	assert(lastState == NULL || theState[upenaltyIter] != lastState);
	lastState = theState[upenaltyIter];
      }

      int totalRecombs = 0;
      float totalLikehood = -FLT_MAX;
      if (_phaseMethod == PHASE_MINREC) {
	totalRecombs = prevState->minRecomb + numRecombs[0];

	if (prevIndex < 0) // apply penalty for error states
	  totalRecombs += CmdLineOpts::max1MarkerRecomb;

	// apply any penalty for only one haplotype transmission
	totalRecombs += penaltyCount;

	if (nUpenaltyIters == 2 && upenaltyIter == 0)
	  // special case one haplotype transmission penalty at the beginning of
	  // the chromosome (where the bits for the corresponding parent will be
	  // unassigned
	  totalRecombs += THE_PENALTY;

	if (numMarkersSinceOneHetPar > CmdLineOpts::bothParHetThreshold) {
	  // apply penalty for children's IVs switching from one parent to other
	  totalRecombs += 1;
	  numMarkersSinceOneHetPar = -1; // penalty applied
	}
      }
      else
	// TODO! this (currently non-functional) code doesn't account for the
	//       missing parent no het penalty and/or the error penalty
	totalLikehood = prevState->maxLikelihood + localLikehood;

      theStateUpdated[upenaltyIter] = checkMinRecomb(iterIV, iterUnassigned,
				       ambig1Unassigned, theState[upenaltyIter],
				       prevState, hetParent,
				       homParentGeno, curParPhase,
				       altPhaseType, ambigLocal,
				       prevIndex, IVambigPar,
				       minMaxRec,
				       numMarkersSinceNonHetPar,
				       numMarkersSinceOneHetPar,
				       totalRecombs, totalLikehood,
				       numRecombs, lowOrderChildBit);
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
							prevState->unassigned);
      size_t oldNumRecombs = numRecombs[0];
      numRecombs[0] = popcount(recombs);

      // For next state: reset and update penalty "recombinations":
      penaltyCount = 0;
      if (penalty) { // are the IV values in the penalty scenario?
	for(int p = 0; p < 2; p++) {
	  if ((~recombs) & penalty & _parBits[p])
	    // no recombination in the one child that differs for a missing data
	    // parent that hasn't had an informative marker for a long stretch:
	    // apply penalty
	    penaltyCount += THE_PENALTY;
	}
      }

      bool equalStates;
      if (_phaseMethod == PHASE_MINREC)
	equalStates = oldNumRecombs == numRecombs[0];
      else {
	float oldLocalLikehood = localLikehood;
	if (isPI) {
	  numRecombs[1] = popcount(recombs & _parBits[1]);
	  // TODO! (2) do parent-specific calculations
	  localLikehood = numRecombs[0] * recombProb +
			      (numDataChildren - numRecombs[0]) * noRecombProb;
	}
	else
	  localLikehood = numRecombs[0] * recombProb +
			      (numDataChildren - numRecombs[0]) * noRecombProb;
	equalStates = oldLocalLikehood == localLikehood;
      }

      if (isPI && IVambigPar && equalStates) {
	for(int upi = 0; upi < nUpenaltyIters; upi++)
	  // TODO: below right? Or only want one theStateUpdated?
	  if (theStateUpdated[upi])
	    theState[upi]->arbitraryPar = 1;
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
			    const uint64_t unassigned,
			    int &lowOrderChildBit) {
  // As inheritance vectors have four equivalent values, we have a fixed key
  // defined via the following convention:
  // (1) In the lowest order two bits for which <iv> is unambiguous and
  //     assigned, the value should be 00.
  // (2) The <iv> values for ambiguous bits are either the default of 10 binary
  //     or 00, with the latter value occurring if satisfying (1) involves
  //     inverting the inheritance value of only one parent. (Note that for
  //     ambiguous bits, <iv> values of 10 can be flipped to 01 without
  //     increasing the numbers of recombinations. Likewise 00 can be flipped
  //     to 11. So it suffices to only consider only two ambiguous values: 10
  //     and 00.)
  // (3) Any unassigned IV values are 0.

  uint64_t unambig, ambigStd;
  lowOrderChildBit = lowOrderUnambigUnassignedBit(allAmbig, unassigned,
						  unambig, ambigStd);

  // The genotype of the child tells us what bits we need to flip: the exact
  // bits that are assigned 1 need to be flipped in all unambiguous children.
  uint8_t flipType = (iv >> lowOrderChildBit) & 3;

  //////////////////////////////////////////////////////////////////////////
  // Get the canonical key value

  // Conveniently, we've got _flips indexed by the 4 possible flip types with
  // the values to flip assigned in each child:
  uint64_t lookupIV = iv ^ (_flips[flipType] & unambig);
  // And we've done something analogous for ambiguous bits; must flip them too:
  lookupIV ^= _ambigFlips[flipType] & ambigStd;

  // Lastly, ensure that any unassigned IV values are 0:
  lookupIV &= ~unassigned;

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
    newState->maxLikelihood = -FLT_MAX;
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

int Phaser::lowOrderUnambigUnassignedBit(const uint64_t allAmbig,
					 const uint64_t unassigned,
					 uint64_t &unambig,
					 uint64_t &ambigStd) {
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
  // the case for these ambig1 type values, with the proviso of course that the
  // ambig1 type ambig bits must be equivalent for the equivalence to hold.
  ambigStd = ((allAmbig & _parBits[1]) >> 1) * 3;

  unambig = _parBits[2] - ambigStd;

  // which bits are unassigned for each parent? See the loop: must ignore
  // scenarios where all children have unassigned IV from one parent.
  uint64_t parUnassigned[2];
  for(int p = 0; p < 2; p++) {
    parUnassigned[p] = unassigned & _parBits[p];
    // is this parent fully unassigned?
    uint8_t parFullUnassign = parUnassigned[p] == _parBits[p];
    // if so, zero out bits: can't omit all children from being used to map to
    // the canonical value
    parUnassigned[p] *= (1 - parFullUnassign);
  }
  // bit 0 set for all children with one (bit 0/1) or both IV values unassigned
  uint64_t unassignedBit0 = (parUnassigned[0] | (parUnassigned[1] >> 1)) &
								    _parBits[0];
  uint64_t assignedChildren = (_parBits[0] - unassignedBit0) * 3;
  uint64_t toFlipBits = unambig & assignedChildren;

  // First determine the <iv> value for the lowest order unambig/assigned child:
  int lowOrderChildBit = ffsll(toFlipBits) - 1;

  // Markers that have all heterozygous children are totally ambiguous and not
  // considered at this stage, so there must be at least one unambiguous child:
  assert(lowOrderChildBit >= 0 && lowOrderChildBit % 2 == 0);

  return lowOrderChildBit;
}

// Check whether the proposed state values lead to minimum or equal numbers
// of recombinations compared to the current values for the state. In either
// case, update the state values accordingly.
bool Phaser::checkMinRecomb(uint64_t fullIV, uint64_t fullUnassigned,
			    uint64_t ambig1Unassigned, State *theState,
			    const State *prevState, uint8_t hetParent,
			    uint8_t homParentGeno, uint8_t curParPhase,
			    uint8_t altPhaseType, uint8_t ambigLocal,
			    int64_t prevIndex, uint8_t IVambigPar,
			    uint16_t minMaxRec[2],
			    int16_t numMarkersSinceNonHetPar,
			    int16_t numMarkersSinceOneHetPar,
			    int totalRecombs, float totalLikehood,
			    size_t numRecombs[2], int lowOrderChildBit) {
  // Does the current previous state lead to minimum recombinations for
  // <theState>? (Could be either a new minimum or an ambiguous one)
  bool theStateUpdated = false;

  // Is the path via <prevIndex> an error? It is if prevIndex is an error,
  // which would make the current state an error if it is used as the new
  // previous state
  bool prevPathError = prevState->error != 0 || prevIndex < 0;
  // Is the current state an error state and/or does it include an error
  // state in its previous state?
  bool curIsError = theState->error > 0;

  uint8_t isPI = hetParent >> 1;

  // TODO! (3) Error case for max likelihood? Relevant to this if statement
  // and following else statement

  // TODO: document the error case (prefer paths without errors)
  // Is the previous state <prevIndex> the new minimally recombinant path to
  // <theState>? Certainly if it produces fewer recombinations.
  // Also if it has equal numbers of recombinations and <theState> is either
  // an error state or has an error state in its previous state path and the
  // proposed new state has no error in its path nor is it an error state.
  if ((_phaseMethod == PHASE_MINREC &&
	(totalRecombs < theState->minRecomb ||
	 (totalRecombs == theState->minRecomb &&
					     curIsError && !prevPathError)) )
	|| (_phaseMethod == PHASE_MAXLIKE &&
				    totalLikehood > theState->maxLikelihood) ) {
    theStateUpdated = true;

    theState->iv = fullIV;
    // next two already assigned in lookupState():
//    theState->ambig = fullAmbig;
//    theState->unassigned = fullUnassigned | ambig1Unassigned;
    if (prevIndex >= 0) {
      theState->prevState = prevIndex;
      // Propagate error state information: if the path of states leading to
      // this one has an error in it, indicate this by having
      // <theState->error> == 2. So we'll copy forward the value of previous
      // value of 2 or convert a previous value of 1 into 2:
      // The following produces 0 if <prevState->error> == 0 and 2 otherwise:
      theState->error = (prevState->error & 2) | ((prevState->error & 1) << 1);
    }
    else { // erroneous previous state
      theState->prevState = -(prevIndex + 1);
      theState->error = 1;
    }
    theState->ambigPrev = 0;
    theState->minRecomb = totalRecombs;
    theState->numMarkersSinceNonHetPar = numMarkersSinceNonHetPar;
    theState->numMarkersSinceOneHetPar = numMarkersSinceOneHetPar;
    theState->maxLikelihood = totalLikehood;
    theState->maxPrevRecomb = numRecombs[0];
    theState->hetParent = hetParent;
    theState->ambigParHet = 1 << hetParent;
    theState->homParentGeno = homParentGeno;
    theState->parentPhase = curParPhase;
    uint8_t parPhaseBits = (1 << curParPhase) |
					    (ambigLocal * (1 << altPhaseType));
    theState->ambigParPhase = parPhaseBits << (2 * hetParent);
    theState->arbitraryPar = (1-isPI) * IVambigPar;

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
  // TODO: can we optimize the conditions here?
  // TODO! need a tolerance for equality comparison of likelihoods, need to
  //       test this before if statement above (epsilon greater isn't really
  //       greater)
  else if ((_phaseMethod == PHASE_MINREC &&
	    totalRecombs == theState->minRecomb &&
	    (prevIndex >= 0 || (theState->error & 1)) &&
	    (prevState->error == 0 || (theState->error & 2))) ||
	   (_phaseMethod == PHASE_MAXLIKE &&
	    totalLikehood == theState->maxLikelihood) ) {
    theStateUpdated = true;

    bool newBestPrev = false;
    if (numRecombs[0] > theState->maxPrevRecomb) {
      // We prefer to back trace to previous states that have the most
      // recombinations relative to the current state. We track the largest
      // number of these local recombinations and, when the current max is
      // exceeded, we update the values below to enable back tracing to use
      // those relevant to the chosen previous state
      newBestPrev = true;
      // Because the canonical swap/phase type of <theState> may differ from
      // the state that it's being changed to (i.e., the old values are an
      // ambiguous alternative), we must swap <theState>'s values so that
      // they match the new ones.  That way, bits set in e.g.,
      // <theState->ambigParPhase> will have a consistent meaning.
      uint8_t flipType = ((theState->iv ^ fullIV) >> lowOrderChildBit) & 3;
      theState->iv = fullIV;
      theState->maxPrevRecomb = numRecombs[0];
      theState->hetParent = hetParent;
      theState->parentPhase = curParPhase;
      
      theState->ambigParPhase =
	swap01Phase[flipType][ theState->ambigParPhase & 15 ] |
	(swap2Phase[flipType][ theState->ambigParPhase >> 4 ] << 4);
    }
    else {
      // See the note in the prevoius branch above <flipType>; here we flip the
      // values that will be newly added (as ambiguous alternatives) to be
      // consistent with those in <theState>
      // those in <theState>.
      uint8_t flipType = ((theState->iv ^ fullIV) >> lowOrderChildBit) & 3;
      // if hetParent < 2, parentPhase flips differently: only have two
      // possibilities: flipped (1) or not (0):
      flipType = isPI * flipType + (1 - isPI) * ((flipType >> hetParent) & 1);
      curParPhase ^= flipType;
      altPhaseType ^= flipType;
    }

    if (theState->homParentGeno == G_MISS)
      theState->homParentGeno = homParentGeno;
    assert(hetParent == 2 || theState->homParentGeno == homParentGeno);

    theState->numMarkersSinceNonHetPar =
	      min(theState->numMarkersSinceNonHetPar, numMarkersSinceNonHetPar);
    theState->numMarkersSinceOneHetPar =
	      min(theState->numMarkersSinceOneHetPar, numMarkersSinceOneHetPar);

    theState->ambigParHet |= 1 << hetParent;
    uint8_t parPhaseBits = (1 << curParPhase) |
					    (ambigLocal * (1 << altPhaseType));
    theState->ambigParPhase |= parPhaseBits << (2 * hetParent);

    // indicate that the state has an arbitrary parent assignment if
    // any previous states that lead to it are of this class
    theState->arbitraryPar |= (1-isPI) * IVambigPar;

    updateAmbigPrev(theState, prevIndex, newBestPrev);

    // Sanity check:
    assert(theState->unassigned == (fullUnassigned | ambig1Unassigned));
  }

  return theStateUpdated;
}

// When there are multiple previous states that can optimally reach <theState>,
// add the corresponding previous state to the list of these previous states
// if such a list exists; otherwise create one
void Phaser::updateAmbigPrev(State *theState, int64_t prevIndex,
			     bool newBestPrev) {
  const uint32_t MAX_IDX_IN_STATE = 1 << 6;

  uint32_t prevState = (prevIndex >= 0) ? prevIndex : -(prevIndex + 1);
  ////////////////////////////////////////////////////////////////////////
  // See comment at declaration of State::ambigPrev
  switch (theState->ambigPrev) {
    case 2: // already storing list, so simply append
      {
	// Can happen that <prevState> is already in the list; check whether
	// that's the case (will necessarily be either the last element or the
	// first since we consider prevStates sequentially in makeFullStates())
	// TODO: optimization: use set not list?
	dynarray<uint32_t> &thePrevs = _ambigPrevLists[theState->prevState];
	int lastElement = thePrevs.length()-1;
	if (newBestPrev) {
	  uint32_t oldFirst = thePrevs[0];
	  if (oldFirst != prevState) {
	    thePrevs[0] = prevState;
	    if (thePrevs[lastElement] == prevState)
	      thePrevs[lastElement] = oldFirst;
	    else
	      thePrevs.append(oldFirst);
	  }
	}
	else {
	  if (thePrevs[lastElement] != prevState && thePrevs[0] !=prevState)
	    thePrevs.append(prevState);
	}
      }
      break;

    case 1: // storing ambiguous values in series of 6 bits in <prevState>
      {
	bool inserted = false;
	const uint32_t mask = MAX_IDX_IN_STATE - 1;

	// Is the value small enough to fit in 6 bits and is there room in
	// <theState->prevState>? The latter is the case so long as the end bit
	// (see below in case 0 for info on this) is not located at bit 30
	if (prevState < MAX_IDX_IN_STATE && theState->prevState < (1 << 30)) {
	  // <prevState> fits in 6 bits; simultaneously: check if there's room
	  // for this new value; find where to put it; and ensure the same
	  // value doesn't occur twice

	  uint32_t curPrevs = theState->prevState;
	  for(int shift = 0; shift < 30; shift += 6) {
	    if (curPrevs-1 == 0) { // subtract end val 1 (see case 0 below)
	      if (newBestPrev) // <prevState> best: beginning of list
		theState->prevState = (theState->prevState << 6) |prevState;
	      else { // append to end
		// Take original <theState->prevState> value, subtract off end
		// value (1 << shift) and append the new <prevState> to be
		// added along with the end value (1<<6), shifted in <shift>
		// bits:
		theState->prevState = (theState->prevState - (1 << shift)) |
					      ((prevState | (1 << 6)) << shift);
	      }
	      inserted = true;
	      break;
	    }
	    if ((curPrevs & mask) == prevState) {
	      // Already present. For the purposes of the code below, this
	      // means it's <inserted>
	      inserted = true;
	      if (newBestPrev && shift > 0) {
		// <prevState> best: move to list beginning
		//
		// Pull out the lower order bits / elements before the current
		// location of <prevState>
		uint32_t lowOrderBits = theState->prevState & ((1 << shift) -1);
		// and the higher order bits / elements after <prevState>
		uint32_t highOrderBits = theState->prevState &
							~((1 << (shift+6)) - 1);
		theState->prevState = highOrderBits | (lowOrderBits << 6) |
								      prevState;
	      }
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
	  dynarray<uint32_t> &thePrevs = _ambigPrevLists[ambigIndex];

	  if (newBestPrev) // <prevState> best: beginning of list
	    thePrevs.append(prevState);
	  // <curPrevs>-1 == 0 when at the end of the list; see case 0 below
	  for(int shift = 0; shift < 30 && curPrevs-1 > 0; shift += 6) {
	    thePrevs.append( curPrevs & mask );
	    curPrevs >>= 6;
	  }
	  if (!newBestPrev)
	    thePrevs.append(prevState);
	}
      }
      break;

    case 0: // previously unambiguous
      {
	uint32_t curPrev = theState->prevState;
	if (curPrev == prevState) {
	  // same state: do nothing
	}
	else if (curPrev < MAX_IDX_IN_STATE && prevState < MAX_IDX_IN_STATE) {
	  // Both values fit in 6 bits, set ambiguity and add the new previous
	  // state
	  theState->ambigPrev = 1;
	  if (newBestPrev)
	    // best previous state: make it the first (lowest order) entry so
	    // that it's used to backtrace
	    theState->prevState = (theState->prevState << 6) | prevState;
	  else
	    theState->prevState |= prevState << 6;
	  // We append a single bit to the end of the list of previous this
	  // allows us to easily determine whether, for example, a 0 value is
	  // an index since a series of 6 bits that equal 0 but that are not
	  // the uppermost bits in <prevState> must be an index. This is
	  // therefore an indicator of the end of the list:
	  theState->prevState |= 1 << 12;
	}
	else {
	  theState->ambigPrev = 2;
	  uint32_t ambigIndex = _ambigPrevLists.length();
	  theState->prevState = ambigIndex;
	  _ambigPrevLists.addEmpty();
	  dynarray<uint32_t> &thePrevs = _ambigPrevLists[ambigIndex];
	  if (newBestPrev) {
	    // <prevState> best: put at beginning of list
	    thePrevs.append(prevState);
	    thePrevs.append(curPrev);
	  }
	  else {
	    thePrevs.append(curPrev);
	    thePrevs.append(prevState);
	  }
	}
      }
      break;

    default:
      fprintf(stderr, "ERROR: got impossible ambigPrev value %d\n",
	      theState->ambigPrev);
      break;
  }
}

// This method does two things:
// TODO! (3) likelihood-based approach to this:
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
	delete curStates[i];
	if (shiftToState < 0) {
	  // Will bump states up, overwritting this index
	  shiftToState = i;
	}
      }
      else {
	if (shiftToState >= 0) {
	  curStates[shiftToState] = curStates[i];
	  shiftToState++;
	}

	if (curStates[i]->error != 2)
	  // one state that is either an error state (newly erroneous) or not an
	  // error at all: can't clear error flag
	  allError2 = false;
      }
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
void Phaser::backtrace(NuclearFamily *theFam, uint8_t missingPar,
		       int chrFirstMarker, int chrLastMarker) {
  // TODO! maximum likelihood-based back tracing code

  int lastIndex = _hmm.length() - 1;

  if (lastIndex < 0)
    // no data: done (can happen if the children are missing data almost
    // everywhere and the marker is ambiguous)
    return;

  uint32_t curStateIdx = findMinStates(_hmm[lastIndex]);
  uint32_t prevStateIdx = UINT32_MAX;

  // For when both parents are missing data:
  // Are we in a run of arbitrary/ambiguous parent assignments that begins at
  // the end of the chromosome? Due to ambiguous recombinations, it is possible
  // for the states to be such that the parent assignments are completely
  // ambiguous. When this happens, we can select one of the parent assignments
  // at the end of the chromosome and trace back until the point where the
  // assignments are no longer ambiguous. That state is the point at which
  // the arbitrary parent assignment starts.
  bool parArbitraryRun = false;
  if (_curBTAmbigSet->size() == 1 && missingPar == 3) {
    state_set_iter it = _curBTAmbigSet->begin();

    // Parent heterozygosity is necessarily ambiguous when:
    // - For children whose IV values differ, the transmitted haplotypes from
    //   both parents differ.
    // - If two different states' IV values match, the parent haplotype
    //   transmissions are opposite each other for the individual child. That is
    //   each such child received AB or BA. See <ivValsExpOpp>.
    //   Note this doesn't apply to ambiguous IV values stored in the ambig
    //   field of the state. These necessarily aren't informative about the
    //   parents as they can take on two different values.
    // - Ambig IV values match between the states
    // It's non-trivial to argue for ambiguity here, but: if all children are
    // AB and BA, the marker is certainly ambiguous as to which parent is
    // heterozygous. When a recombination happens from a marker where all
    // children are AB or BA, it can be attributed to either parent. This yields
    // either AA or BB and these two IV values satisfy the first condition
    // above. Ambiguous IV values are similar to the first case: the child has
    // one of two different IV values that are consistent with two possible
    // phase assignments in the parents.

    uint64_t curIV = _hmm[lastIndex][curStateIdx]->iv;
    uint64_t ambigIV = it->iv;
    uint64_t ivDiffBits = curIV ^ ambigIV;
    // Which bits do we expect to have opposite values for?
    uint64_t ivExpectOppositeBits = ~(ivDiffBits | it->ambig);
    uint64_t ivValsExpOpp = ambigIV & ivExpectOppositeBits;
    parArbitraryRun = _hmm[lastIndex][curStateIdx]->ambig == it->ambig &&
	      ((ivDiffBits >> 1) & _parBits[0]) == (ivDiffBits & _parBits[0]) &&
	      (((ivValsExpOpp >> 1) ^ ivValsExpOpp) & _parBits[0]) ==
					  (_parBits[0] & ivExpectOppositeBits);

  }

  // Should we assign the phase with status reflecting the potential for one
  // haplotype transmission only? This is used to liberally make such
  // assignments at the beginning and end of a chromosome where the IV is
  // potentially consistent with one hap trans
  bool assignOneHapTrans = false;
  uint8_t oneHapTransHetParent;
  uint8_t oneHapTransHap;
  uint64_t oneHapTransOutlierIV;

  // Number of recombinations in <curState> relative to <prevState> below
  uint8_t numRecombs; // TODO: optimization: do we want this?
  int lastAssignedMarker = chrLastMarker + 1; // technically not assigned yet
  uint64_t lastAssignedIV = 0, lastIVFlip = 0;
  uint64_t lastIVparDiff = 0;
  bool lastIVSet = false;
  for(int hmmIndex = lastIndex; hmmIndex >= 0; hmmIndex--) {
    int curHmmIndex = hmmIndex; // not redundant: errors modify hmmIndex
    State *curState = _hmm[curHmmIndex][curStateIdx];

    uint8_t curAmbigParPhase = curState->ambigParPhase;
    uint8_t curAmbigParHet = curState->ambigParHet;
    uint8_t curArbitraryPar = curState->arbitraryPar;
    // Indicates which inheritance vector values differ among any ambiguous
    // states that are in <_curIdxSet>. The true IV value is ambiguous in that
    // case.
    uint64_t ivFlippable = 0;

    if (!assignOneHapTrans &&
	_hmm[curHmmIndex][curStateIdx]->numMarkersSinceNonHetPar > 0 &&
	(curHmmIndex == lastIndex || // always assign one hap trans at end ...
	 // ... and beginning of a chromosome
	 _hmm[curHmmIndex][curStateIdx]->numMarkersSinceNonHetPar ==
				_hmmMarker[curHmmIndex] - chrFirstMarker + 1)) {
      uint8_t hetParent = _hmm[curHmmIndex][curStateIdx]->hetParent;

      // detect which haplotype was transmitted:
      uint8_t homParent = 1 - hetParent;
      uint64_t parHap =_hmm[curHmmIndex][curStateIdx]->iv & _parBits[homParent];
      uint64_t oppParHap = parHap ^ _parBits[homParent];

      if ((parHap & (parHap - 1)) == 0) { // test for 1 bit set
	assignOneHapTrans = true;
	oneHapTransHetParent = hetParent;
	oneHapTransOutlierIV = parHap;
	// transmitted haplotype is 0
	oneHapTransHap = 0;
      }
      else if ((oppParHap & (oppParHap - 1)) == 0) { // test for all but 1 set
	assignOneHapTrans = true;
	oneHapTransHetParent = hetParent;
	oneHapTransOutlierIV = oppParHap;
	// transmitted haplotype is 1
	oneHapTransHap = 1;
      }
      else {
	// Note: because the IVs are unassigned initially, it's not possible
	// to check for the one hap trans pattern at the initial markers.
	// Therefore, the <State::numMarkersSinceNonHetPar> value will be
	// non-zero for missing data parents that are homozygous and haven't yet
	// had their IVs assigned
	assert(_hmmMarker[curHmmIndex] - chrFirstMarker <=
					    CmdLineOpts::oneHapTransThreshold);
      }
    }

    if (assignOneHapTrans) {
      uint8_t hetParent = _hmm[curHmmIndex][curStateIdx]->hetParent;

      // has it ended?
      if (_hmm[curHmmIndex][curStateIdx]->numMarkersSinceNonHetPar <= 0 ||
	  oneHapTransHetParent != hetParent)
	assignOneHapTrans = false;
      else {
	// can do this by (1) setting <ivFlippable> for the child that has the
	// outlier haplotype ...
	ivFlippable |= oneHapTransOutlierIV;

	// ... (2) indicating the site may be heterozygous for both parents ...
	curAmbigParHet |= 1 << 2;

	// ... and (3) assigning the phase type for the both parent het case:
	// the phase type for the both parent het determines which haplotype is
	// set to be printed (other is set missing); we want this to be
	// <transHap>.
	// The homozygous genotype is either 0 or 3, and we'll divide by 3 to
	// get a boolean:
	uint8_t homGenoKind = _hmm[curHmmIndex][curStateIdx]->homParentGeno / 3;
	// If the phase type is 0, when the genotype is heterozygous,
	// haplotype 0 is allele 0, and haplotype 1 is allele 1. Whichever
	// haplotype is the same as homGenoKind will be printed (non-missing).
	// So, if <homGenoKind> == 0 and <transHap> == 0, we want the phase type
	// to be the same as <transHap>. Same if <homGenoKind> == 0 and
	// <transHap> == 1.
	// The reverse is true when <homGenoKind> == 1. So it suffices to
	// xor:
	// Note: must shift the phase types for <homParent> and <hetParent> to
	// their respective positions in the 2-bit <phaseType>.
	uint8_t homParent = 1 - hetParent;
	uint8_t phaseType = ((oneHapTransHap ^ homGenoKind) << homParent) |
		    (_hmm[curHmmIndex][curStateIdx]->parentPhase << hetParent);
	// (Note: shifting by 4 bits to get to the initial bit that stores the
	// both parent het phase types)
	curAmbigParPhase |= 1 << (4 + phaseType);

      }
    }

    uint8_t curHomParGeno = curState->homParentGeno;
    // above only valid for one parent het states or states with ambiguous
    // parent heterozygosity
    bool homParGenoAssigned = curState->hetParent < 2 ||
						      (curAmbigParHet & 3) != 0;
    uint8_t ambigHomParGeno = 0;

    // Find all the possible parent phase types that are have equal and minimal
    // numbers of recombinations at this marker and append their previous states
    // to <_prevIdxSet>.
    for(state_set_iter it = _curBTAmbigSet->begin(); it !=_curBTAmbigSet->end();
									it++) {
      State *ambigState = _hmm[curHmmIndex][ it->stateIdx ];
      uint64_t ivDiffBits = curState->iv ^ it->iv;

      if (parArbitraryRun && (ivDiffBits > 0 || curState->ambig != it->ambig)) {
	// Only remains arbitrary if the following conditions hold (see longer
	// comment above when parArbtitrary run is first set true)
	// (1) All differing bits different for both parents?
	// (2) All non-diff (same) bits the opposite of each other in the IV?
	// (3) The ambig IV values (curState->ambig and it->ambig) match
	uint64_t ivExpectOppositeBits = ~(ivDiffBits | it->ambig);
	uint64_t ivValsExpOpp = it->iv & ivExpectOppositeBits;
	parArbitraryRun = curState->ambig == it->ambig &&
	      ((ivDiffBits >> 1) & _parBits[0]) == (ivDiffBits & _parBits[0]) &&
	      (((ivValsExpOpp >> 1) ^ ivValsExpOpp) & _parBits[0]) ==
					  (_parBits[0] & ivExpectOppositeBits);
	if (!parArbitraryRun)
	  curArbitraryPar = 1;
      }

      if (!parArbitraryRun) {
	ivFlippable |= ivDiffBits;

	curAmbigParHet |= ambigState->ambigParHet;
	curAmbigParPhase |= it->ambigParPhase;
	curArbitraryPar |= ambigState->arbitraryPar;
      }

      if (ambigState->hetParent < 2 || (ambigState->ambigParHet & 3) != 0) {
	assert(ambigState->homParentGeno != G_MISS);
	if (!homParGenoAssigned) {
	  curHomParGeno = ambigState->homParentGeno;
	  homParGenoAssigned = true;
	}
	else
	  ambigHomParGeno |= curHomParGeno ^ ambigState->homParentGeno;
      }

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
      if (prevHmmIndex >= 0) {
	BT_ambig_info dontcare1;
	uint8_t dontcare2;
	collectAmbigPrevIdxs(it->iv, it->ambig, ambigState->ambigPrev,
			     ambigState->prevState, _hmm[prevHmmIndex],
			     dontcare1, dontcare2);
      }
    }

    // Remove any ambig par phase types that have equivalent phase to <curState>
    curAmbigParPhase -= (1 << curState->parentPhase) << 2 * curState->hetParent;
    // Remove any ambig par het types that are equivalent to <curState>
    curAmbigParHet -= 1 << curState->hetParent;

    // In the previous state, resolve ambiguous <iv> values and propagate
    // backward any <iv> values that were unassigned in that state
    if (hmmIndex - 1 >= 0) {
      // Note: curState->error == 2 means the path has an error in it. We only
      // need deal with curState->error == 1:
      if (curState->error & 1) {
	// Set error for immediately previous marker.
	uint64_t childrenData = _genos[hmmIndex-1].second;
	uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
	theFam->setStatus(/*marker=*/ _hmmMarker[hmmIndex-1], PHASE_ERR_RECOMB,
			  _genos[hmmIndex-1].first, childrenData, missing);
	// <curState->prevState> references a state two indexes back
	deleteStates(_hmm[hmmIndex-1]);
	hmmIndex--;
      }

      BT_ambig_info prevStateInfo;
      uint8_t numAmbig1Recombs;
      collectAmbigPrevIdxs(curState->iv, curState->ambig, curState->ambigPrev,
			   curState->prevState, _hmm[hmmIndex-1],
			   prevStateInfo, numAmbig1Recombs);

      prevStateIdx = prevStateInfo.stateIdx;
      State *prevState = _hmm[hmmIndex-1][prevStateIdx];
      prevState->iv = prevStateInfo.iv;
      prevState->ambig = prevStateInfo.ambig;
      // TODO: optimization: remove entry for prevStateInfo from _prevBTAmbig

      // Both the following are equivalent, but now we do something different
      // as noted next
//      numRecombs = curState->minRecomb - prevState->minRecomb -
//			  (curState->error & 1) * CmdLineOpts::max1MarkerRecomb;
//      numRecombs = curState->maxPrevRecomb;
      // Changes in IV when a child is ambig1 induce a recombination that is
      // counted at a later marker than it originally occurred. The
      // <numAmbig1Recombs> counts the number of such recombinations that must
      // be attributed to an earlier marker to be consistent with the changes
      // in IV values:
      prevState->maxPrevRecomb += numAmbig1Recombs;
      numRecombs = curState->maxPrevRecomb - numAmbig1Recombs;
    }
    else
      numRecombs = 0;

    // Missing genotype value is 01, not any other, so:
    uint8_t parentData = _genos[curHmmIndex].first;
    uint8_t parMissing = (parentData & 5) & ~((parentData & 10) >> 1);
    uint64_t childrenData = _genos[curHmmIndex].second;
    uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
    // if <curHomParGeno> is non-missing iff <homParGenoAssigned>
    assert((curHomParGeno != G_MISS) == homParGenoAssigned);
    // either we have a homozygous parent genotype OR the site is PI with no
    // ambiguity in parent heterozygosity
    assert(curHomParGeno != G_MISS ||
	   (curState->hetParent == 2 && (curAmbigParHet & 3) == 0));
    assert(!ambigHomParGeno);
    theFam->setPhase(_hmmMarker[curHmmIndex], curState->iv,
		     curState->ambig & _parBits[1], missing, ivFlippable,
		     parMissing, curState->hetParent, curHomParGeno,
		     curState->parentPhase, numRecombs, curAmbigParHet,
		     curAmbigParPhase, curArbitraryPar);

    //////////////////////////////////////////////////////////////////////////
    // For all markers between the next informative one and the current one,
    // assign (1) which parent haplotypes were _un_transmitted, and
    // (2) for any ambiguous sites, indicate whether they are necessarily
    // homozygous:

    // Determine which parent haplotypes were _un_transmitted. The first two
    // bits are parent 0's haplotype 0 and 1, and the second two bits are
    // parent 1's haplotype 0 and 1. If the corresponding bit is set to 1, the
    // haplotype was _not_ transmitted
    uint64_t untrans;
    // TODO: need the flanking informative markers for each parent do this right
    // We use the IV values at the flanking informative markers to determine
    // this. <ivFlippable> values are uncertain and so those values will not
    // figure into which haplotypes the parents transmitted.
    // start with <ivFlippable> at the current and subsequent marker:
    uint64_t uncertainIV = lastIVFlip | ivFlippable;
    if (lastIVSet)
      // if we have a meaningful IV for the subsequent marker, indicate that the
      // IV values are uncertain for any IV values that are flipped between them
      uncertainIV |= curState->iv ^ lastAssignedIV;

    uint64_t curIVparDiff = (curState->iv & _parBits[0]) ^
					    ((curState->iv & _parBits[1]) >> 1);

    // make the transmitted haplotype assignments for all markers between the
    // current marker and the subsequent informative marker
    int curMarker = _hmmMarker[curHmmIndex];
    int startIndex = curMarker + 1;
    if (curHmmIndex == 0)
      // TODO: bug: what if there's an error at the first position and the true
      // last marker is at curHmmIndex == 1?
      startIndex = chrFirstMarker;
    for(int m = startIndex; m < lastAssignedMarker; m++) {
      if (m == curMarker)
	continue; // Don't set untrans for the current marker (is PHASE_OK)

      untrans = 0;

      const PhaseVals &phase = theFam->getPhase(m);
      assert(phase.status != PHASE_OK);
      // If the child is missing data, it should not figure into the transmitted
      // haplotypes regardless of the imputed IV value at this uninformative
      // marker. In cases, for example, where all but one child is missing data
      // each parent will only have transmitted one haplotype and we will have
      // limited information about their overall genotype status. Note that we
      // do this work of determining the <untrans> value to enable imputation
      // of the parent's genotype, but we can't impute based on children whose
      // data is missing.
      // Missing status in the 0th bit per child
      uint64_t thisMarkMiss = phase.ambigMiss & _parBits[0];

      for(int p = 0; p < 2; p++) {
	// align the missing bit with the relevant <curState->iv> and
	// <uncertainIV> values
	uint64_t shiftedMiss = thisMarkMiss << p;
	// (curState->iv & _parBits[p]) == 0, the first haplotype was
	// transmitted. To detect when it was _un_transmitted, we use
	// ~curState->iv and check when that value is 0
	if ((~curState->iv & _parBits[p] & ~(uncertainIV | shiftedMiss)) == 0)
	  // first haplotype for this parent untransmitted:
	  untrans |= 1 << (2 * p); // use 1 for first haplotype
	// (curState->iv & _parBits[p]) == 0 implies the second haplotype was
	// untransmitted, so:
	if ((curState->iv & _parBits[p] & ~(uncertainIV | shiftedMiss)) == 0)
	  untrans |= 2 << (2 * p); // use 2 for second haplotype
      }

      // Below, wish only to consider the IV values where the children are
      // non-missing and where their IV value is certain.
      // Because the <curIVparDiff> value is only defined for the _parBits[0]
      // values, here we OR in both uncertainIV and (uncertainIV >> 1),
      // capturing both parent IV values. Since this is a mask, the _parBits[1]
      // values don't change the result below.
      uint64_t ivToInclude = ~(uncertainIV | thisMarkMiss | (uncertainIV >> 1));

      // At PHASE_AMBIG sites, for many IV values, we can infer that the parents
      // are in fact homozygous (for opposite alleles). The IV values where this
      // is not the case are those where both parent transmitted the same
      // pattern (or the exact opposite pattern).
      // This is related to isIVambigPar().
      // Are we uncertain about homozygous status?
      bool curUncertainHomozyAtAmbig = (curIVparDiff & ivToInclude) == 0 ||
			    ( (curIVparDiff ^ _parBits[0]) & ivToInclude ) == 0;
      bool lastUncertainHomozyAtAmbig = (lastIVparDiff & ivToInclude) == 0 ||
			    ( (lastIVparDiff ^ _parBits[0]) & ivToInclude ) ==0;
      // Must both be homozygous if we are _not_ uncertain for cur and last:
      bool bothParHomozyAtAmbig = 1 -
		      (curUncertainHomozyAtAmbig | lastUncertainHomozyAtAmbig);

      theFam->setUntransPar(m, untrans, bothParHomozyAtAmbig);
    }

    // To detect untrans, we exclude samples that recombine between informative
    // positions. Because the IV values are set even for uninformative values
    // at prior markers, we'd like to propagate back the IV values here. Because
    // of the way we construct states, the recombined value is always guaranteed
    // to be later on the chromosome. Thus, we propagate back IV values for
    // the uninformative parent (if one is homozygous) and for any children that
    // are missing data. This will enable the <uncertainIV> value above to
    // identify IV values that do differ between (informative) markers relative
    // to the data for each child.
    // TODO: what if the heterozygous status of the parent is ambiguous?
    // TODO: what about ambiguous states? Should we be propagating back IV
    // values in the loop toward the beginning of this method in order to
    // calculate IVflippable, etc.? In some sense that is what we're calculating
    // here: which IV values are potentially flippable and therefore shouldn't
    // factor into which parent's haplotypes were transmitted.
    uint8_t hetParent = curState->hetParent;
    uint8_t oneParHomozy = 1 - (hetParent >> 1);
    uint64_t propagateMask = (missing * 3) |
					  (oneParHomozy * _parBits[hetParent]);
    if (lastIVSet)
      lastAssignedIV = (curState->iv & ~propagateMask) |
					       (lastAssignedIV & propagateMask);
    else
      lastAssignedIV = curState->iv;
    lastIVSet = true;
    lastIVFlip = ivFlippable;
    lastIVparDiff = curIVparDiff;


    deleteStates(_hmm[curHmmIndex]);

    // ready for next iteration -- make prev into cur for states and state sets
    curStateIdx = prevStateIdx;
    BT_state_set *tmp = _curBTAmbigSet;
    _curBTAmbigSet = _prevBTAmbigSet;
    _prevBTAmbigSet = tmp;
    _prevBTAmbigSet->clear();

    lastAssignedMarker = _hmmMarker[curHmmIndex];
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
      _curBTAmbigSet->clear();
      minLastMarker = curRecomb;
      minStateIdx = i;
    }
    else if (curRecomb == minLastMarker) {
      BT_ambig_info values(i, theStates[i]->iv, theStates[i]->ambig,
			   theStates[i]->parentPhase,
			   theStates[i]->ambigParPhase);
      _curBTAmbigSet->insert(values);
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

// Given <curIV,curAmbig> -- state values that are reachable with minimum
// recombinations via <prevState>, determines what the previous <iv> and <ambig>
// values would be if the <cur*> values applied to the current state. This
// works by propagating back information about the inheritance vector values
// contained in <cur*>. This assigns previously-<unassigned> bits in
// <prevState> and makes many previously ambiguous <iv> values unambiguous. The
// knowledge of the <iv> values at the later state allows the choice between
// the two possible <iv> values which are opposite one another and therefore
// differ in the transmitted homolog from both parents. Ambiguities only remain
// if there are recombinations on both sides of the ambiguous <iv> value (such
// as for a putative non-crossover/gene conversion) in which case which parent
// recombined is truly ambiguous.
// Returns (via parameters) the updated <iv> and <ambig> values: <newPrevIV> and
// <newPrevAmbig>.
void Phaser::propagateBackIV(uint64_t curIV,uint64_t curAmbig, State *prevState,
			     uint64_t &newPrevIV, uint64_t &newPrevAmbig,
			     uint8_t &newPrevParPhase,
			     uint8_t &newPrevAmbigParPhase,
			     uint8_t &numAmbig1Recombs) {
  // TODO: potential optimization: don't execute if <prevState->unassigned> == 0
  // and <prevState->ambig> == 0

  uint64_t unambig, ambigStd;
  int lowOrderChildBit = lowOrderUnambigUnassignedBit(prevState->ambig,
						      prevState->unassigned,
						      unambig, ambigStd);

  uint8_t flipType = ((curIV ^ prevState->iv) >> lowOrderChildBit) & 3;

  // Will swapping lead only to the same IV in prev as in cur? Only swap if so.
  // Note that the unassigned bits should not factor into this
  uint8_t shouldSwap =
		((prevState->iv ^ curIV) & unambig & ~prevState->unassigned) ==
			  (_flips[flipType] & unambig & ~prevState->unassigned);
  flipType *= shouldSwap; // if 0, flipType has no effect
  newPrevIV = prevState->iv ^ _flips[flipType];
  
  newPrevAmbigParPhase =
	swap01Phase[flipType][ prevState->ambigParPhase & 15 ] |
	(swap2Phase[flipType][ prevState->ambigParPhase >> 4 ] << 4);
  // if hetParent < 2, parentPhase flips differently
  uint8_t isPI = prevState->hetParent >> 1;
  flipType = isPI * flipType +
			  (1 - isPI) * ((flipType >> prevState->hetParent) & 1);
  newPrevParPhase = prevState->parentPhase ^ flipType;

  // Now propagate backward any <iv> values that were unassigned previously
  newPrevIV = (newPrevIV & ~prevState->unassigned) |
					(curIV & prevState->unassigned);

  // Note: ambiguities that remain in <curAmbig> will not give information
  // about resolving such in <prevState>
  uint64_t ambigToResolve = prevState->ambig & ~curAmbig;
  // Get both bits set for children that have the two different types of
  // ambiguities (standard and ambig1):
  uint64_t prevAnyAmbig = (ambigToResolve & _parBits[0]) * 3;
  uint64_t prevStdAmbig = ((ambigToResolve & _parBits[1]) >> 1) * 3;
  assert((prevAnyAmbig & prevStdAmbig) == prevStdAmbig);
  uint64_t prevAmbig1 = prevAnyAmbig - prevStdAmbig;

  uint64_t ambigRecombs = (curIV ^ newPrevIV) & prevAnyAmbig;
  // set both bits for children that inherit a recombination from the given
  // parent
  uint64_t parRecombs[2] = { (ambigRecombs & _parBits[0]) * 3,
			    ((ambigRecombs & _parBits[1]) >> 1) * 3 };

  // Invert the phase for ambiguous assignments that recombine on both homologs
  // relative to <cur*>
  uint64_t bothRecomb = parRecombs[0] & parRecombs[1];
  newPrevIV ^= bothRecomb;

  // children that are ambiguous in <prevState> _and_ recombine on one homolog
  // relative to <cur*> truly have ambiguous phase.
  // Ambig1 ambiguities are always resolved; the issue with them has to do with
  // counting of recombinations when there were undefined <iv> values before
  // encountering a PI state in which the relevant sample is heterozygous.
  // TODO: optimization -- test this:
//  uint64_t oneRecomb = parRecombs[0] ^ parRecombs[1];
//  prevState->ambig = prevStdAmbig & oneRecomb;
  // old less efficient version of updating <prevState->ambig>
  uint64_t noRecomb = ambigToResolve & ~(parRecombs[0] | parRecombs[1]);
  newPrevAmbig = prevState->ambig & ~(bothRecomb | noRecomb | prevAmbig1);

  // count the number of ambig1 recombinations that were resolved; we will
  // ensure that the real location of the crossover gets this count, not
  // some later marker
  numAmbig1Recombs = popcount(bothRecomb & prevAmbig1 & _parBits[0]);
}

// Inserts BT_ambig_info objects for previous states in <curState->prevState>
// into <_prevBTAmbigSet>. Also calls propagateBackIV() for all such previous
// states which is needed to properly calculate <ivFlippable> in backtrace().
// Without this call, states that have ambiguous IV values can be such that both
// <IV> values are indicated as uncertain when only one is (or even that the
// opposite parent is uncertain relative to the true one [at least in principle
// this could happen])
//
// <prevStates> is the list of states at the previous marker that are
// referenced by <curState->prevState>.
//
// Returns via reference <thePrevInfo>, which provides the index of the
// preferred previous state for back tracing and the updated (via
// propagateBackIV() <iv> and <ambig> values for that state). Note that
// these values are only useful when we're considering the state that will be
// assigned to the current marker. It is ignored in other contexts (other
// ambiguous states at the current marker)
void Phaser::collectAmbigPrevIdxs(uint64_t curIV, uint64_t curAmbig,
				  uint8_t ambigPrev, uint32_t prevState,
				  const dynarray<State*> &prevStates,
				  BT_ambig_info &thePrevInfo,
				  uint8_t &numAmbig1Recombs) {
  switch (ambigPrev) {
    case 0: // only one
      {
	uint32_t thePrevStateIdx = prevState;
	uint64_t prevIV, prevAmbig;
	uint8_t prevParPhase, prevAmbigParPhase;
	propagateBackIV(curIV, curAmbig, prevStates[thePrevStateIdx],
			prevIV, prevAmbig, prevParPhase, prevAmbigParPhase,
			numAmbig1Recombs);
	thePrevInfo.stateIdx = thePrevStateIdx;
	thePrevInfo.iv = prevIV;
	thePrevInfo.ambig = prevAmbig;
	thePrevInfo.parPhase = prevParPhase;
	thePrevInfo.ambigParPhase = prevAmbigParPhase;
	_prevBTAmbigSet->insert( thePrevInfo );
      }
      break;

    case 1: // multiple previous states in <prevState>: add all
      {
	const uint32_t mask = (1 << 6) - 1;
	uint32_t prevs = prevState;

	// Set the first index as the previous state according to our convention
	uint32_t thePrevStateIdx = prevs & mask;
	uint64_t prevIV, prevAmbig;
	uint8_t prevParPhase, prevAmbigParPhase;
	propagateBackIV(curIV, curAmbig, prevStates[ thePrevStateIdx ],
			prevIV, prevAmbig, prevParPhase, prevAmbigParPhase,
			numAmbig1Recombs);
	thePrevInfo.stateIdx = thePrevStateIdx;
	thePrevInfo.iv = prevIV;
	thePrevInfo.ambig = prevAmbig;
	thePrevInfo.parPhase = prevParPhase;
	thePrevInfo.ambigParPhase = prevAmbigParPhase;
	_prevBTAmbigSet->insert( thePrevInfo );
	prevs >>= 6;

	// <prevs>-1 == 0 when at the end of the list. Described elsewhere when
	// ambigous prevous states are dealt with.
	for(int shift = 6; shift < 30 && prevs-1 > 0; shift += 6) {
	  uint32_t curPrevIdx = prevs & mask;
	  uint8_t dontcare; // numAmbig1Recombs will default to first prev state
	  propagateBackIV(curIV, curAmbig, prevStates[ curPrevIdx ],
			  prevIV, prevAmbig, prevParPhase, prevAmbigParPhase,
			  dontcare);
	  BT_ambig_info values(curPrevIdx, prevIV, prevAmbig, prevParPhase,
			       prevAmbigParPhase);
	  _prevBTAmbigSet->insert( values );
	  prevs >>= 6;
	}
      }
      break;

    case 2: // previous states in <_ambigPrevLists>: add all
      {
	dynarray<uint32_t> &prevs = _ambigPrevLists[prevState];
	uint32_t thePrevStateIdx = prevs[0];
	uint64_t prevIV, prevAmbig;
	uint8_t prevParPhase, prevAmbigParPhase;
	propagateBackIV(curIV, curAmbig, prevStates[ thePrevStateIdx ],
			prevIV, prevAmbig, prevParPhase, prevAmbigParPhase,
			numAmbig1Recombs);
	thePrevInfo.stateIdx = thePrevStateIdx;
	thePrevInfo.iv = prevIV;
	thePrevInfo.ambig = prevAmbig;
	thePrevInfo.parPhase = prevParPhase;
	thePrevInfo.ambigParPhase = prevAmbigParPhase;
	_prevBTAmbigSet->insert( thePrevInfo );

	for(int i = 1; i < prevs.length(); i++) {
	  uint32_t curPrevIdx = prevs[i];
	  uint8_t dontcare; // numAmbig1Recombs will default to first prev state
	  propagateBackIV(curIV, curAmbig, prevStates[ curPrevIdx ],
			  prevIV, prevAmbig, prevParPhase, prevAmbigParPhase,
			  dontcare);
	  BT_ambig_info values(curPrevIdx, prevIV, prevAmbig, prevParPhase,
			       prevAmbigParPhase);
	  _prevBTAmbigSet->insert( values );
	}
      }
      break;

    default:
      fprintf(stderr, "ERROR: got impossible ambigPrev value %d\n",
	      ambigPrev);
      exit(5);
      break;
  }
}

// For debugging:
void Phaser::tracePrior(int hmmIndex, int stateIdx) {
  if (hmmIndex > _hmm.length())
    return;

  do {
    State *state = _hmm[ hmmIndex ][ stateIdx ];
    printf("iv = %ld, ambig = %ld, unassigned = %ld, minRecomb = %d, numSinceNonHet = %d\n",
	   state->iv, state->ambig, state->unassigned, state->minRecomb,
	   state->numMarkersSinceNonHetPar);
    printf("  hmmIndex = %d, state = %d  (marker = %d)\n",
	   hmmIndex, stateIdx, _hmmMarker[hmmIndex]);

    hmmIndex--;
    hmmIndex -= state->error & 1; // skip immediately previous marker
    stateIdx = state->prevState;
    if (state->ambigPrev != 0) {
      printf("  multiple previous states: (marker = %d)\n",
	     _hmmMarker[hmmIndex]);
      if (state->ambigPrev == 1) {
	const uint32_t mask = (1 << 6) - 1;
	uint32_t prevs = state->prevState;
	printf("  ");
	for(int shift = 0; shift < 30 && prevs-1 > 0; shift += 6) {
	  uint32_t curPrevIdx = prevs & mask;
	  printf("  %d,%d", hmmIndex, curPrevIdx);
	  prevs >>= 6;
	}
	printf("\n");
      }
      else {
	printf("    ambigPrev = %d, prevState = %d\n",
	       state->ambigPrev, state->prevState);
      }
      break;
    }
  } while(hmmIndex >= 0);
  fflush(stdout);
}
