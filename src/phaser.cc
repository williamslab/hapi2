// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include <cmath>
#include <bit>
#include <float.h>
#include "phaser.h"
#include "cmdlineopts.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// declare/initialize static members
dynarray< dynarray<State*> >    Phaser::_hmm;
dynarray<int>                   Phaser::_hmmMarker;
dynarray<int>                   Phaser::_parInformMarkers[2];
dynarray< dynarray<uint32_t> >  Phaser::_ambigPrevLists;
std::deque< std::pair<int, uint32_t> >  Phaser::_prevErrorStates;
dynarray< std::pair<uint8_t,uint64_t> > Phaser::_genos;
Phaser::state_ht                Phaser::_stateHash;
Phaser::BT_state_set            *Phaser::_curBTAmbigSet;
Phaser::BT_state_set            *Phaser::_prevBTAmbigSet;
Phaser::uint64_set              Phaser::_prevCanonIVsBT;
PhaseMethod                     Phaser::_phaseMethod;
uint8_t  Phaser::_parentData;
uint64_t Phaser::_childrenData[5];
uint64_t Phaser::_parBits[3];
uint64_t Phaser::_childSexes[2];
uint64_t Phaser::_flips[4];
uint64_t Phaser::_ambigFlips[4];
int      Phaser::_curMarker;
int      Phaser::_lastInformMarker;
int      Phaser::_lastRealInformIndex;
int      Phaser::_lastForceInformMarker;
int      Phaser::_lastForceInformIndex;
int16_t  Phaser::_maxNumMarkersSinceHetPar[2];
int      Phaser::_lastMarkerInMaxSince;
bool     Phaser::_onChrX;
uint8_t  Phaser::_missingPar;
uint8_t  Phaser::_firstMissP;
uint8_t  Phaser::_limitMissP;

extern uint8_t swap01Phase[4][16];
extern uint8_t swap2Phase[4][16];

void Phaser::run(NuclearFamily *theFam, int chrIdx, FILE *log) {
  _onChrX = Marker::isChromX(chrIdx);
  _phaseMethod = PHASE_MINREC; // TODO! remove this and make command line option

  // Ready storage containers/various state to analyze this chromosome
  initPhaseState(theFam);

  dynarray<State> partialStates;
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);

  // PHASE! Build states/calculate corresponding scores for each marker
  for(_curMarker = firstMarker; _curMarker <= lastMarker; _curMarker++) {
    uint8_t parentGenoTypes, homParGeno, childGenoTypes;
    int numMissChildren = 0; // how many children are missing data?

    ///////////////////////////////////////////////////////////////////////////
    // Step 0: get the data for this marker
    getFamilyData(theFam, parentGenoTypes, childGenoTypes, numMissChildren);

    ///////////////////////////////////////////////////////////////////////////
    // Step 1: Determine marker type; handle trivial, ambiguous, or erroneous
    // markers; and check for the need to force an informative marker
    //
    // Should we add a forced informative marker? Set to which parent if >= 0.
    int8_t forceInform = -1;
    int mt = getMarkerType_prelimAnalyses(theFam, parentGenoTypes,
					  childGenoTypes, homParGeno,
					  forceInform, firstMarker, log);
    if (mt < 0)
      continue; // marker processing complete; no further phasing needed

    ///////////////////////////////////////////////////////////////////////////
    // Step 2: For each consistent marker type, deduce the inheritance vector
    // bits where possible and generate partial states that only have bits
    // defined for the informative parent(s)
    makePartialStates(partialStates, mt, homParGeno);

    ///////////////////////////////////////////////////////////////////////////
    // Step 3: Using states for the previous informative marker, fill in
    // unknown inheritance vector bits by propagation and generate full states.
    // Simultaneously calculate the minimum recombination counts for each state
    // and identify which previous state(s) produce that minimum count (i.e.,
    // perform count-based Viterbi calculation)
    makeFullStates(partialStates, firstMarker, theFam->numChildren(),
		   numMissChildren, forceInform);

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
    fprintf(log, "  Back tracing...\n");
    fflush(log);
  }
  backtrace(theFam, firstMarker, lastMarker, log);
  if (CmdLineOpts::verbose) {
    fprintf(log, "done\n");
    fflush(log);
  }
}

// Do initial setup of values used throughout phasing the current chromosome
// if on the X chromosome, assign <_childSexes>
void Phaser::initPhaseState(NuclearFamily *theFam) {
  int numChildren = theFam->numChildren();
  if (numChildren > 32) {
    fflush(stdout);
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: cannot phase more than 32 children in a family.\n");
    fprintf(stderr, "       changing to >64 bits for inheritance vectors would fix this\n");
    exit(9);
  }

  _hmm.clear();
  _hmmMarker.clear();
  _parInformMarkers[0].clear();
  _parInformMarkers[1].clear();
  _prevErrorStates.clear();
  _genos.clear();
  for(int i = 0; i < _ambigPrevLists.length(); i++) {
    _ambigPrevLists[i].clear();
  }
  _ambigPrevLists.clear();

  _lastInformMarker = -1;
  _lastRealInformIndex = -1;
  _lastForceInformMarker = -1;
  _lastForceInformIndex = -1;
  _maxNumMarkersSinceHetPar[0] = -1;
  _maxNumMarkersSinceHetPar[1] = -1;
  _lastMarkerInMaxSince = -1;

  // assign bits to indicate which parents are missing; bit 0: dad, bit 1: mom
  _missingPar = 0;
  if (!theFam->_parents->first->hasData())
    _missingPar |= 1;
  if (!theFam->_parents->second->hasData())
    _missingPar |= 2;

  // force parents to be missing according to command-line options:
  _missingPar |= CmdLineOpts::forceMissingParBits;

  // for looping over missing parents for one hap trans assignment
  // analyze parent 0 if bit 0 in _missingPar is 1; otherwise start looping
  // from parent 1 (result of this is either 0 or 1):
  _firstMissP = 1 - (_missingPar & 1);
  // don't analyze parent 0 (the dad) on the X chromosome: only has one hap
  // the following will flip 0 to 1 iff _firstMissP == 0 and _onChrX == 1
  _firstMissP ^= (~_firstMissP) & _onChrX;
  // analyze parent 1 if bit 1 in _missingPar is 1; otherwise shouldn't iterate
  // to 1 (result of this is either 1 or 2):
  _limitMissP = (_missingPar >> 1) + 1;

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

  // when phasing chrX, get inheritance vector bits corresponding to the male
  // and female children. Need these to deal with the differences between their
  // inheritance patterns.
  _childSexes[0] = _childSexes[1] = 0;
  if (_onChrX) {
    int numChildren = theFam->_children.length();
    for(int c = 0; c < numChildren; c++) {
      if (theFam->_children[c]->getSex() == 'U') {
	fprintf(stderr, "ERROR: attempt to phase X chromosome, but child %s has unknown sex\n",
	    theFam->_children[c]->getId());
	exit(18);
      }
      uint8_t index = (theFam->_children[c]->getSex() == 'M') ? 0 : 1;
      _childSexes[index] += 3ul << (c*2);
    }
  }
}

// Looks up and stores the genotype values for the parents and children.
// For speedy marker type detection, uses bit representation in <*GenoTypes>
// variables to indicate which genotypes the parents and children have.
void Phaser::getFamilyData(NuclearFamily *theFam, uint8_t &parentGenoTypes,
			   uint8_t &childGenoTypes, int &numMissChildren) {
  NuclearFamily::par_pair parents = theFam->_parents;
  dynarray<PersonBulk*> &children = theFam->_children;

  // Get data for dad (bits 0 and 1) and mom (bits 2 and 3)
  uint8_t dadData = (_missingPar & 1) ? G_MISS :
				       parents->first->getBitGeno(_curMarker);
  uint8_t momData = (_missingPar & 2) ? G_MISS :
				       parents->second->getBitGeno(_curMarker);
  if (_onChrX && dadData == G_HET) // set males to missing instead of het on X
    dadData = G_MISS;
  _parentData = dadData + (momData << 2);

  // Which of the genotypes are present in the parents/children? There are
  // four possible genotypes, and the first four bits will have a value of 1
  // iff the corresponding genotype value is present in at least one child
  parentGenoTypes = (1 << dadData) | (1 << momData);
  childGenoTypes = 0;

  for(int g = 0; g < 5; g++)
    _childrenData[g] = 0; // note: limited to 32 children (checked above)

  int numChildren = children.length();
  for(int c = 0; c < numChildren; c++) {
    uint8_t curChildData = children[c]->getBitGeno(_curMarker);
    if (_onChrX && children[c]->getSex() == 'M' && curChildData == G_HET)
      curChildData = G_MISS; // set males to missing instead of het on X
    _childrenData[ curChildData ] += 3ul << (c*2);
    _childrenData[4] += ((uint64_t) curChildData) << (c*2);
    childGenoTypes |= 1ul << curChildData; // observed genotype <curChildData>
    // following is 1 when the genotype is 1 (missing), and 0 otherwise
    numMissChildren += (curChildData & 1) & (~curChildData >> 1);
  }
}

// Determines the type of the current marker (informative for one parent,
// partly informative, uninformative, error, ambiguous) and processes the
// trivial marker types (error, ambiguous, and uninformative).
// Also checks for the need to force in an informative marker, which is needed
// when (on chrX) the Mom is missing data and she has transmitted the same
// chromosome to all the children.
int Phaser::getMarkerType_prelimAnalyses(NuclearFamily *theFam,
					 uint8_t parentGenoTypes,
					 uint8_t childGenoTypes,
					 uint8_t &homParGeno,
					 int8_t &forceInform,
					 int firstMarker,
					 FILE *log) {
  // First get the marker type:
  int mt;
  bool specialXMT = false; // see getMarkerTypeX() code for this special type
  uint8_t imputeUninfHomPar = G_MISS;
  homParGeno = G_MISS;
  if (!_onChrX) // autosomal?
    mt = getMarkerTypeAuto(parentGenoTypes, childGenoTypes, homParGeno,
			   imputeUninfHomPar);
  else // X
    mt = getMarkerTypeX(childGenoTypes, homParGeno, imputeUninfHomPar,
			specialXMT);
  assert(mt > 0 && (mt & ~((1 << MT_N_TYPES) -1) ) == 0);

  if (CmdLineOpts::verbose)
    printMarkerType(mt, log);

  // Can immediately handle erroneous or ambiguous markers:
  if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
    // should only be one of the above:
    assert(mt == (1 << MT_ERROR) || mt == (1 << MT_AMBIG));
    switch(mt) {
      case 1 << MT_ERROR:
	theFam->setStatus(_curMarker, PHASE_ERROR, _parentData,
			  _childrenData[4], _childrenData[G_MISS] &_parBits[0]);
	break;
      case 1 << MT_AMBIG:
	theFam->setStatus(_curMarker, PHASE_AMBIG, _parentData,
			  _childrenData[4], _childrenData[G_MISS] &_parBits[0]);
	break;
    }
    return -1; // no further phasing needed
  }

  // Handle uninformative markers:
  if (mt & (1 << MT_UN)) {
    // On chrX, when Mom is missing data, long stretches without any
    // informative markers often mean that she has transmitted only one
    // haplotype to the children. We'll force an informative marker after a
    // sufficiently long stretch so that the inheritance vector reflects the
    // one haplotype transmission status.
    // TODO: what if this is <specialXMT>?
    if (_onChrX && (_missingPar & 2) && (mt & (1 << MT_FI_1)) &&
	_curMarker - firstMarker >= CmdLineOpts::forceInformInit) {
      if (checkForceInform()) {
	// will add a new marker later, and its index is this:
	_lastForceInformIndex = _hmmMarker.length();
	_lastForceInformMarker = _curMarker;
	if (specialXMT || ((_parentData & 3) != G_MISS))
	  // going to use this as a fully informative marker; need to store
	  // dad's genotype in <homParGeno>. Occassionally (specifically in
	  // the cases that this if statement specifies) <homParGeno> will
	  // store what Mom's homozygous genotype would be if she is homozygous.
	  // Get dad's homozygous genotype (or may be missing)
	  homParGeno = _parentData & 3;
	return mt;
      }
    }
    // On the autosomes, long stretches without any informative markers for a
    // missing parent often mean that that parent has transmitted only one
    // haplotype to the children. Force an informative marker after
    // sufficiently long stretch.
    // Check for the need to force an informative marker if either
    // _lastRealInformIndex < 0 -- i.e., there has not been a forced
    // informative marker since the most recent "real" one -- or if there have
    // been at least CmdLineOpts::forceInformSeparation markers since the last
    // forced one.
    else if (!_onChrX && _missingPar > 0 &&
	     (_lastRealInformIndex < 0 ||
	      (_curMarker - _lastInformMarker) >=
					  CmdLineOpts::forceInformSeparation) &&
	     // need it to be possible for this marker to be heterozygous for
	     // the OHT parent
	     ( (mt & 1 << MT_PI) || (mt & (1 << MT_FI_1)) ) &&
	     _curMarker - firstMarker >= CmdLineOpts::forceInformInit) {
      int numMarkersSinceInform = _curMarker - _lastInformMarker;
      int lastHmmIdx = _hmm.length() - 1;
      if (lastHmmIdx >= 0 && _lastMarkerInMaxSince < lastHmmIdx) {
	// need to calculate _maxNumMarkersSinceHetPar for the most recent
	// marker
	_maxNumMarkersSinceHetPar[0] = _maxNumMarkersSinceHetPar[1] = -1;
	for(uint16_t i = 0; i < _hmm[lastHmmIdx].length(); i++) {
	  State *theState = _hmm[lastHmmIdx][i];
	  if ( _missingPar == 3 &&
				isIVambigPar(theState->iv, theState->ambig) ) {
	    if (
		((theState->iv & _parBits[0]) == 0 || (theState->iv & _parBits[0]) == _parBits[0]) &&
		((theState->iv & _parBits[1]) == 0 || (theState->iv & _parBits[1]) == _parBits[1])
	       )
	      // don't skip: this is a region in which both parents are OHT and
	      // we need to continue to add forced informative markers to
	      // prevent errors from interrupting the OHT region.
	      ;
	    else
	      // no need to force informative markers when the het parent is
	      // ambiguous
	      continue;
	  }
	  for (int p = _firstMissP; p < _limitMissP; p++) {
	    // _May_ want to force an informative marker; check the IV:
	    uint64_t parHap = theState->iv & _parBits[p];
	    uint64_t oppParHap = parHap ^ _parBits[p];
	    // 0, 1 or 2 bits set in parHap or oppParHap?
	    // * Only 1 corresponds to either all but one of the children or all
	    //   children except one receiving the same haplotype.
	    // * Also check for 2 since, in some cases, two children could
	    //   recombine near each other and produce a oneHapTrans scenario
	    uint8_t numTrans[2] = { (uint8_t) popcount(parHap),
				    (uint8_t) popcount(oppParHap) };
	    // which pattern (parHap/oppParHap) has min number of 1s (popcount)?
	    int minPat = (numTrans[0] <= numTrans[1]) ? 0 : 1;
	    if (numTrans[minPat] > 2)
	      continue; // too many children would need to recombine: no OHT

	    if (theState->numMarkersSinceNonHetPar[p] >
						_maxNumMarkersSinceHetPar[p])
	      _maxNumMarkersSinceHetPar[p] =
					  theState->numMarkersSinceNonHetPar[p];
	  }
	}
	_lastMarkerInMaxSince = lastHmmIdx;
      }

      // Decide whether to force an informative marker
      for (int p = _firstMissP; p < _limitMissP; p++) {
	if (_maxNumMarkersSinceHetPar[p] >= 0 &&
	    _maxNumMarkersSinceHetPar[p] + numMarkersSinceInform >=
						CmdLineOpts::forceInformInit) {
	  // Old version only assigned forceInform for one parent if only one
	  // parent had numMarkersSinceHetPar sufficiently large. This led to
	  // some bizarre state paths since there were ways to avoid getting
	  // any penalties. Always adding a forced informative marker of all
	  // the parent het types ensures that all paths get penalized.
	  forceInform = 2;
	  break;
	}
      }

      if (lastHmmIdx < 0)
	// no previous state => no informative marker though we've passed
	// CmdLineOpts::forceInformInit markers: should force one
	forceInform = 2;
      else if (forceInform >= 0 && _lastRealInformIndex < 0)
	_lastRealInformIndex = lastHmmIdx;
    }

    if (!specialXMT) { // standard uninformative
      uint8_t swapHetChildren = 0;
      if ((_parentData & 3) == G_HOM1 || ((_parentData >> 2) & 3) == G_HOM0)
	// The default way to print phase of heterozygous children at
	// uninformative markers is with allele 0 transmitted by dad and 1
	// transmitted by mom. When Dad is homozygous for allele 1 or (when
	// he is missing) Mom is homozygous for allele 0, this order needs to
	// be swapped:
	swapHetChildren = 1;
      if (imputeUninfHomPar != G_MISS)
	homParGeno = imputeUninfHomPar;
      theFam->setUninform(_curMarker, _parentData, _childrenData[4],
			  _childrenData[G_MISS] & _parBits[0], homParGeno,
			  swapHetChildren);
      if (forceInform >= 0) { // force an informative marker?
	if (_missingPar < 3) { // only one parent missing:
	  // Need <homParGeno> to be the genotype of the non-missing parent.
	  // Its current value may be <imputeUninfHomPar>, which applies only
	  // if the type is MT_UN. Forcing an informative marker is
	  // inconsistent with this, so we'll get the value directly from
	  // <_parentData>.
	  uint8_t missParent = _missingPar >> 1;
	  uint8_t nonMissParGeno;
	  if (missParent == 0)
	    nonMissParGeno = _parentData >> 2; // mom's genotype
	  else
	    nonMissParGeno = _parentData & 3;  // dad's genotype
	  if (nonMissParGeno != G_MISS)
	    homParGeno = nonMissParGeno;
	}
	return mt;
      }
    }
    else { // special X chromosome case
      uint8_t swapHetChildren = 0;
      if (homParGeno == G_HOM1) // homParGeno here is for dad
	// see case just above
	swapHetChildren = 1;
      theFam->setSpecialX(_curMarker, _parentData, _childrenData[4],
			  _childrenData[G_MISS] & _parBits[0], homParGeno,
			  swapHetChildren);
    }
    return -1; // no further phasing needed
  }

  if (!_onChrX && _missingPar > 0) {
    // track informative markers for each parent

    // marker should be informative to get here
    assert( (mt & ( (1 << MT_FI_1) | (1 << MT_PI) )) );

    uint8_t dadGeno = _parentData & 3;
    uint8_t momGeno = _parentData >> 2;
    if (dadGeno == G_HOM0 || dadGeno == G_HOM1)
      _parInformMarkers[1].append( _curMarker );
    else if (momGeno == G_HOM0 || momGeno == G_HOM1)
      _parInformMarkers[0].append( _curMarker );
    else if ( mt == (1 << MT_PI) ) {
      // Note: we only count PI markers as informative for both parents when
      // that's the only option. The reason is, if the marker _might_ be
      // FI for one parent, then the homozygous parent won't be informative
      _parInformMarkers[0].append( _curMarker );
      _parInformMarkers[1].append( _curMarker );
    }
  }

  return mt;
}

// Determines what type of marker this is using data for the parents if present
// or based on the observed genotype values for the the children when one or
// both parent's data are missing.
// For markers that are/may be MT_FI_1 type, also assigns <homParGeno> as the
// homozygous parent genotype
// In some circumstances, when a marker has the possibility of being both MT_UN
// and MT_FI_1, <imputeUninfHomPar> gives the homozgyous genotype of the
// uncertain parent (i.e., assuming it is MT_UN).
int Phaser::getMarkerTypeAuto(uint8_t parentGenoTypes, uint8_t childGenoTypes,
			      uint8_t &homParGeno, uint8_t &imputeUninfHomPar) {
  // Only valid values for parentGenoTypes are between 1 and 12 (excluding 7
  // and 11 which are caught below)
  assert(parentGenoTypes >= 1 && parentGenoTypes <= 12);
  assert(childGenoTypes >= 1 && childGenoTypes <= 15);

  // check no data for all children case first: is ambiguous
  if (childGenoTypes == (1 << G_MISS))
    return 1 << MT_AMBIG;

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
	  homParGeno = G_HOM0;
	  imputeUninfHomPar = G_HOM1;
	  return (1 << MT_UN) | (1 << MT_FI_1);
	case G_HOM1: // one parent homozygous for 1
	  homParGeno = G_HOM1;
	  imputeUninfHomPar = G_HOM0;
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

// For X chromosome sites: Determines what type of marker this is using data
// for the parents (if present) or based on the observed genotype values for
// the children when one or both parent's data are missing
//
// For markers that are/may be MT_FI_1 type, assigns <homParGeno> as dad's
// homozygous genotype
// When a marker has the possibility of being both MT_UN and MT_FI_1,
// <imputeUninfHomPar> gives the homozygous genotype of the mom (who is either
// heterozygous or homozygous). (Except for <specialXMT> described next.)
//
// Also determines whether this is a special marker type that could be fully
// informative or uninformative and where decisions about how to assign the
// parents' genotypes get deferred until backtracing.
int Phaser::getMarkerTypeX(uint8_t childGenoTypes, uint8_t &homParGeno,
			   uint8_t &imputeUninfHomPar, bool &specialXMT) {
  // ensure only first four bits set and that dad (lowest order 3 bits) is not
  // heterozygous (enforced in getParentData()):
  assert(_parentData <= 15 && (_parentData & 3) != G_HET);
  // ensure sons are not heterozygous (enforced in getParentData()):
  assert((_childrenData[G_HET] & _childSexes[0]) == 0);
  // ensure <childGenoTypes> is valid; this is bit vector representing which of
  // the four possible genotypes are assigned in the children
  assert(childGenoTypes >= 1 && childGenoTypes <= 15);

  // check no data for all children case first: is ambiguous
  if (childGenoTypes == (1 << G_MISS))
    return 1 << MT_AMBIG;

  uint8_t dadGeno = _parentData & 3;
  uint8_t oppDadHomType = (dadGeno == G_HOM0) ? G_HOM1 : G_HOM0;
  uint8_t momGeno = (_parentData >> 2);
  uint8_t oppMomHomType = (momGeno == G_HOM0) ? G_HOM1 : G_HOM0;
  switch (momGeno) {
    case G_HOM0:
    case G_HOM1:
      ////////////////////////////////////////////////////////////////////////
      // Mom homozygous
      if (dadGeno == momGeno) { // both homozyogus for same allele
	if (childGenoTypes & ((1 << oppMomHomType) | (1 << G_HET)) )
	  // Mendelian error -- invalid bits set: only hom for the same allele
	  // as the parents and missing are possible
	  return 1 << MT_ERROR;
	else
	  // both parents homozygous: uninformative marker
	  return 1 << MT_UN;
      }
      else { // dad either homozygous for a different allele than Mom or missing
	// at markers where the mom is homozygous, sons should not be
	// homozygous for the opposite allele:
	uint64_t maleChildHomOppMom =
				_childrenData[ oppMomHomType ] & _childSexes[0];
	if (maleChildHomOppMom > 0)
	  return 1 << MT_ERROR;

	if (dadGeno == oppMomHomType) { // one parent homozygous 0, other 1
	  // with both parents homozygous for different alleles, daughters
	  // should be heterozygous or missing:
	  uint64_t femaleChildHetGeno =
	      (_childrenData[ G_HET ] | _childrenData[G_MISS]) & _childSexes[1];
	  if (femaleChildHetGeno != _childSexes[1])
	    return 1 << MT_ERROR;
	}
	else if (dadGeno == G_MISS) {
	  // daughters must carry their dad's allele, and even though we don't
	  // know what it is, they should not be homozygous for opposite
	  // alleles:
	  uint64_t femaleChildHom0 = _childrenData[ G_HOM0 ] & _childSexes[1];
	  uint64_t femaleChildHom1 = _childrenData[ G_HOM1 ] & _childSexes[1];
	  if (femaleChildHom0 > 0 && femaleChildHom1 > 0)
	    return 1 << MT_ERROR;
	  else if (femaleChildHom0 > 0 || femaleChildHom1 > 0) {
	    if ((_childrenData[ G_HET ] & _childSexes[1]) > 0)
	      // Since Mom is homozygous, all female children should either be
	      // homozygous for the same genotype or all should be heterozygous
	      // (or could be missing). Having both types implies a Mendelian
	      // error
	      return 1 << MT_ERROR;
	    // impute dad:
	    homParGeno = (femaleChildHom0) ? G_HOM0 : G_HOM1;
	  }
	  else if ((_childrenData[ G_HET ] & _childSexes[1]) > 0)
	    // in this case, all the daughters are heterozygous; that means
	    // that dad must have the opposite homozygous genotype to Mom.
	    // impute him:
	    homParGeno = oppMomHomType;
	}

	// no error
	// both parents homozygous: uninformative marker
	return 1 << MT_UN;
      }

      break;

    case G_HET:
      //////////////////////////////////////////////////////////////////////////
      // Mom heterozygous

      // with Mom heterozygous, sons can be any homozygous (or missing)
      // genotype, and assertion above ensured they're not heterozygous, so
      // no further checks needed on them
      if (dadGeno != G_MISS) {
	// daughters must carry their dad's allele, so can't be oppDadHomType:
	uint64_t femaleChildHomOppDad =
			      _childrenData[ oppDadHomType ] & _childSexes[1];
	if (femaleChildHomOppDad != 0)
	  return 1 << MT_ERROR;

	// mom heterozgyous: informative for her
	// dad is always uninformative, and him being non-missing ensures this
	// marker is not phase ambiguous
	homParGeno = dadGeno;
	return 1 << MT_FI_1;
      }
      else {
	// dad is missing data

	// daughters must carry their dad's allele, and even though we don't
	// know what it is, they should not be homozygous for opposite alleles:
	uint64_t femaleChildHom0 = _childrenData[ G_HOM0 ] & _childSexes[1];
	uint64_t femaleChildHom1 = _childrenData[ G_HOM1 ] & _childSexes[1];
	if (femaleChildHom0 > 0 && femaleChildHom1 > 0)
	  return 1 << MT_ERROR;
	else if (femaleChildHom0 > 0 || femaleChildHom1 > 0)
	  // impute dad:
	  homParGeno = (femaleChildHom0) ? G_HOM0 : G_HOM1;
	else {
	  // dad's genotype is ambiguous; if there are no sons or they are all
	  // missing data, it's impossible to phase (Mom is het and all her
	  // daughters are, too, with no data for sons)
	  if ( (_childrenData[G_MISS] & _childSexes[0]) == _childSexes[0] )
	    return 1 << MT_AMBIG;
	}

	// mom heterozgyous: informative for her
	// dad is always uninformative, and we checked above that the marker is
	// not ambiguous
	return 1 << MT_FI_1;
      }

      break;

    case G_MISS:
      //////////////////////////////////////////////////////////////////////////
      // Mom is missing, will try to impute her

      // with Mom missing, she could be heterozygous, and in that case, sons
      // could be any homozygous (or missing) genotype. Assertion above ensured
      // the sons are not heterozygous, so no further checks needed on them as
      // far as Mendelian errors
      {
	// first, we'll take one step to try to infer Mom:
	// if sons or daughters differ in their genotype, Mom is heterozygous.
	// the following doesn't separate out sons and daughters, but we
	// don't need to: we'll check next for Mendelian errors, but if we
	// ignore that, the following condition is true if either (a) there are
	// sets of sons of opposite homozygous types or (b) there are daughters
	// of one type and sons of the other (or both [a] and [b]).
	bool momMustBeHet = _childrenData[ G_HOM0 ] > 0 &&
			    _childrenData[ G_HOM1 ] > 0;

	if (dadGeno != G_MISS) {
	  // daughters must carry their dad's allele, so can't be oppDadHomType:
	  uint64_t femaleChildHomOppDad =
			      _childrenData[ oppDadHomType ] & _childSexes[1];
	  if (femaleChildHomOppDad != 0)
	    return 1 << MT_ERROR;

	  // if there are heterozygous daughters and either daughters or sons
	  // that are homozygous for the same allele as Dad, Mom is
	  // heterozygous (the homozygous children inherited one allelic type
	  // and the het daughters necessarily got the other allele from the
	  // Mom)
	  momMustBeHet = momMustBeHet || ( _childrenData[ dadGeno ] > 0 &&
			 (_childrenData[ G_HET ] & _childSexes[1]) > 0 );
	  if (momMustBeHet) {
	    // mom heterozgyous: informative for her
	    // dad is always uninformative, and him being non-missing ensures
	    // this marker is not phase ambiguous
	    homParGeno = dadGeno;
	    return 1 << MT_FI_1;
	  }
	  else {
	    // if Mom is homozygous, she's homozygous for the same allele as
	    // the children are (note that if both homozygous types are present
	    // in the children, <momMustBeHet> is assigned above):
	    if (_childrenData[ G_HOM0 ] > 0)
	      imputeUninfHomPar = G_HOM0; // Mom's potential genotype
	    else if (_childrenData[ G_HOM1 ] > 0)
	      imputeUninfHomPar = G_HOM1; // Mom's potential genotype
	    else
	      // all children het; this implies there's > 0 daughters and they
	      // would have inherited their dad's allele. That means they
	      // inherited the opposite allele from their Mom, and if Mom is
	      // homozyogus, she's homozygous for the opposite allele to dad:
	      imputeUninfHomPar = oppDadHomType; // Mom's potential genotype
	    // Mom could be homozygous but may also be heterozygous with one of
	    // her alleles untransmitted:
	    homParGeno = dadGeno;
	    return (1 << MT_UN) | (1 << MT_FI_1);
	  }
	}
	else {
	  // both Mom and Dad missing

	  // daughters must carry their dad's allele, and even though we don't
	  // know what it is, they should not be homozygous for opposite
	  // alleles:
	  uint64_t femaleChildHom0 = _childrenData[ G_HOM0 ] & _childSexes[1];
	  uint64_t femaleChildHom1 = _childrenData[ G_HOM1 ] & _childSexes[1];
	  if (femaleChildHom0 > 0 && femaleChildHom1 > 0)
	    return 1 << MT_ERROR;
	  else if (femaleChildHom0 > 0 || femaleChildHom1 > 0)
	    // impute dad:
	    homParGeno = (femaleChildHom0) ? G_HOM0 : G_HOM1;
	  else {
	    // dad's genotype is ambiguous; if there are no sons or they are all
	    // missing data, it's impossible to phase (Mom may be het and all
	    // the daughters are het, too, with no data for sons)
	    if ( (_childrenData[G_MISS] & _childSexes[0]) == _childSexes[0] )
	      return 1 << MT_AMBIG;
	  }

	  // if there are both homozygous and heterozygous daughters, Mom is
	  // heterozygous
	  momMustBeHet = momMustBeHet ||
			 ( (femaleChildHom0 > 0 || femaleChildHom1 > 0) &&
			   (_childrenData[ G_HET ] & _childSexes[1]) > 0 );
	  if (momMustBeHet) {
	    // mom heterozgyous: informative for her
	    // dad is always uninformative, and we checked above that the
	    // marker is not ambiguous
	    return 1 << MT_FI_1;
	  }
	  else {
	    if (homParGeno == G_MISS && _childSexes[0] > 0) {
	      // Special X marker type here -- corner case which is:
	      // - all daughters are het (otherwise homParGeno != G_MISS)
	      // - at least one son is non-missing (otherwise MT_AMBIG)
	      // - all sons must have the same homozygous genotype or
	      //   <momMustBeHet> would be true
	      // => in this case, either:
	      //    a. the marker is FI
	      //    b. the marker is uninformative with Mom being the same
	      //       homozygous type as the son(s), and Dad being the
	      //       opposite homozygous type
	      // we'll initially assume that the site is uninformative and then
	      // during backtracing, check the surrounding IVs to see if it's
	      // possible for Mom to be FI without introducing more
	      // recombinations
	      // we'll store the homozygous type dad would have if the marker
	      // is uninformative (opposite the sons' homozygous type)
	      homParGeno = ( (_childrenData[ G_HOM0 ] & _childSexes[0]) > 0 ) ?
			   G_HOM1 : G_HOM0;
	      specialXMT = true;
	    }
	    // Mom could be homozygous but may also be heterozygous with one of
	    // her alleles untransmitted:
	    return (1 << MT_UN) | (1 << MT_FI_1);
	  }
	}
      }

      break;
  }

  return -1; // shouldn't happen
}

// For use with verbose mode, prints the marker type to the log
void Phaser::printMarkerType(int mt, FILE *log) {
  fprintf(log, "    Marker %d:", _curMarker);
  if (mt == 1 << MT_ERROR)
    fprintf(log, " Mendelian error");
  else if (mt == 1 << MT_AMBIG)
    fprintf(log, " Ambiguous");
  else if (mt & (1 << MT_UN))
    fprintf(log, " Uninformative");
  else {
    int FI1 = 1 << MT_FI_1;
    int PI = 1 << MT_PI;
    int both = FI1 | PI;
    if ((mt & both) == both)
      fprintf(log, " One OR both parents heterozygous");
    else if (mt & FI1)
      fprintf(log, " One parent heterozygous");
    else if (mt & PI)
      fprintf(log, " Both parents heteterozygous");
  }
  fprintf(log, " dad: ");
  printGeno((Geno) _parentData & 3, log);
  fprintf(log, " mom: ");
  printGeno((Geno) _parentData >> 2, log);
  fprintf(log, "\n");
  fflush(log);
}

// For use with verbose mode, prints a given genotype to the log
void Phaser::printGeno(uint8_t type, FILE *log) {
  switch(type) {
    case G_HOM0:
      fprintf(log, "0/0");
      break;
    case G_HET:
      fprintf(log, "0/1");
      break;
    case G_HOM1:
      fprintf(log, "1/1");
      break;
    case G_MISS:
      fprintf(log, "./.");
      break;
  }
}

// Check whether a forced informative marker should be placed at the current
// marker. Normally a marker that can be explained as uninformative will be
// treated as such, but on chrX, if Mom transmitted only one chromosome, there
// will not be any informative markers. To force HAPI to recognize that only
// one haplotype was transmitted, the follow code checks certain marker
// count-based thresholds, and will force in an informative marker (which will
// necessarily imply that Mom transmitted the same chromosome to all her
// children).
bool Phaser::checkForceInform() {
  // TODO: would like to check the IV, too
  if (_lastForceInformMarker < 0) { // no recent forced informative marker
    // which marker should we count from? We'll look back the maximum number of
    // informative markers that we tolerate in a force inform interval or to
    // the beginning of the chromosome if not enough markers have passed.
    int anchorInfMarker =
      (_hmmMarker.length() > CmdLineOpts::forceInformTolerance) ?
	_hmmMarker[_hmmMarker.length() - CmdLineOpts::forceInformTolerance] : 0;
    return /*forceInform=*/ _curMarker - anchorInfMarker >=
						  CmdLineOpts::forceInformInit;
  }
  else if (_curMarker - _lastForceInformMarker >=
						CmdLineOpts::forceInformInit) {
    // In a one hap trans region for Mom on chrX.
    // Because erroneous sites can interrupt such a region, we adopt a strategy
    // to ensure that these errors don't prevent a forced informative marker
    // from being added:
    // (1) if <= <forceInformTolerance> "informative" (likely erroneous)
    // markers have occurred since the last force informative marker, we'll use
    // the last force informative marker as the "anchor" for deciding when to
    // add another one: we'll add a new marker now.
    // (2) if more than <forceInformTolerance> but fewer than
    // <numInformToBreakForceInform>  "informative" (potentially erroneous)
    // markers have occurred since the last informative marker, we'll move the
    // "anchor" marker forward, ignoring <forceInformTolerance> markers and
    // potentially add a forceInformMarker later.
    // (3) If >= <numInformToBreakForceInform> "informative" markers have
    // occurred, we'll assume we've left the OHT region and stop tracking the
    // last force inform marker
    int numInformSinceLastForce = _hmmMarker.length() - 1 -
							  _lastForceInformIndex;
    if (numInformSinceLastForce <= CmdLineOpts::forceInformTolerance)
      return /*forceInform=*/ true;
    else if (numInformSinceLastForce <
				    CmdLineOpts::numInformToBreakForceInform) {
      _lastForceInformIndex++;
      _lastForceInformMarker = _hmmMarker[ _lastForceInformIndex ];
    }
    else
      _lastForceInformIndex = _lastForceInformMarker = -1;
  }

  // not forcing an informative marker (at least not at this marker)
  return false;
}

// Given the marker types <markerTypes> that are consistent with the family
// data at the current marker, generates partial states indicating, where
// known, which allele each parent transmitted to the children.
void Phaser::makePartialStates(dynarray<State> &partialStates,
			       int markerTypes, uint8_t homParGeno) {
  int lastIndex = _genos.length();
  _genos.addEmpty();
  _genos[lastIndex].first = _parentData;
  _genos[lastIndex].second = _childrenData[4];

  // Fully informative for one parent:
  if (markerTypes & (1 << MT_FI_1)) {
    if (homParGeno == G_MISS) {
      assert(_onChrX); // should only happen on X:
      // do we have any daughters? and are any non-missing? if so we can try
      // both possibilities for dad's genotype
      if (_childSexes[1] > 0 &&
		  (_childrenData[G_MISS] & _childSexes[1]) != _childSexes[1]) {
	// try both possibilities for dad's (homozygous) genotype
	makePartialFI1States(partialStates, /*homParGeno=*/ G_HOM0);
	makePartialFI1States(partialStates, /*homParGeno=*/ G_HOM1);
      }
      else {
	// no information to impute dad; use a placeholder genotype for him to
	// get the partial state made and then set his genotype to missing
	makePartialFI1States(partialStates, /*homParGeno=*/ G_HOM0);
	partialStates[ partialStates.length() - 1 ].homParentGeno = G_MISS;
      }
    }
    else
      makePartialFI1States(partialStates, homParGeno);
  }

  // Partly informative
  if (markerTypes & (1 << MT_PI)) {
    assert( !_onChrX );
    makePartialPIStates(partialStates);
  }
}

// Helper for makePartialStates(): applicable to fully informative for one
// parent markers (or the states that correspond to this possibility when the
// parent's genotypes aren't fully known)
void Phaser::makePartialFI1States(dynarray<State> &partialStates,
				  uint8_t homParGeno) {
  // which children are heterozygous or hom1 at this marker? store these.
  // this allows us to modify these values so that the autosomal code works
  // on chrX
  uint64_t childHetHom1[2] = { _childrenData[ G_HET ],_childrenData[ G_HOM1 ] };

  // First deal with chrX nuances:
  if (_onChrX) {
    // Use a hack to get the below autosomal code to work for the X chromosome
    // by making the sons diploid. This works by setting their genotypes such
    // that they carry their dad's allele -- i.e., making sons that are
    // homozygous for the opposite allele as dad heterozygous. (Sons that are
    // homozygous for the same allele as dad need not be changed: they'd
    // continue to have the same homozygous genotype if they were diploid.)
    // Start by getting those sons' bits:
    uint8_t oppDadHomType = (homParGeno == G_HOM0) ? G_HOM1 : G_HOM0;
    uint64_t maleChildHomOppDad =
				_childrenData[ oppDadHomType ] & _childSexes[0];
    // now we want to clear these sons' bits in their original homozygous
    // genotype; to avoid a conditional we make this a bit complex; our goal is
    //   _childrenData[ oppDadHomType ] ^= maleChildHomOppDad;
    // but we don't want to change _childrenData. First get an indicator of
    // whether oppDadHomType is G_HOM1 by getting its first bit:
    uint8_t oppDadIsHom1 = oppDadHomType & 1;
    // now we'll flip, and this will either flip the G_HOM1 value OR the G_HET
    // value. In the latter case, we're setting these sons as heterozygous,
    // since their bits in the G_HET vector will start as 0. Setting these sons
    // bits in the G_HET vector to 1 is something we want to do regardless, so
    // it won't hurt:
    childHetHom1[ oppDadIsHom1 ] ^= maleChildHomOppDad;
    // now, regardless of what happened above, we ensure that these sons bits
    // are set in the G_HET vector:
    childHetHom1[0] |= maleChildHomOppDad;
  }

  uint8_t startPar, endPar;

  // For chrX, only Mom is heterozygous
  if (_onChrX) {
    startPar = endPar = 1;
  }
  // If both parents are missing, we'll make states corresponding to each
  // being heterozygous (other homozygous), consistent with the ambiguity
  else if (_parentData == (G_MISS << 2) + (G_MISS)) {
    startPar = 0;
    endPar = 1;
  }
  // parent 0 het: either we have data for parent 0 as a het, or we have data
  // for parent 1 as homozygous (or both)
  else if ((_parentData & 3) == G_HET || (_parentData >> 2) == G_HOM0 ||
	   (_parentData >> 2) == G_HOM1) {
    startPar = endPar = 0;
  }
  // parent 1 het:
  else {
    assert((_parentData >> 2) == G_HET || (_parentData & 3) == G_HOM0 ||
	   (_parentData & 3) == G_HOM1);
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

    // Initially assumes the parent phase is assigned such that allele 0 is on
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
    // either gets the children that are G_HOM1 or those that are G_HET:
    newState.iv = _parBits[hetPar] & childHetHom1[ homParAll1 ];
    // No ambiguous bits (though a missing data child could have ambiguous bits
    // propagated from a previous marker, but that happens during full state
    // construction)
    newState.ambig = 0;
    // No information about transmission from homozygous parent and for
    // children that are missing data:
    // one exception on the homozygous parent being unassigned: if we're on the
    // X chromosome, we say that the father's bits are all assigned to
    // non-missing children
    newState.unassigned = ( (1-_onChrX) * _parBits[homPar] ) |
						      _childrenData[ G_MISS ];
    // on the X chromosome, females' paternal IV values is 0 and males are 1;
    // we use the latter when printing the (hemizygous) males
    newState.iv |= (_onChrX * ( _childSexes[0] & _parBits[0] )) &
						      ~_childrenData[ G_MISS ];
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
void Phaser::makePartialPIStates(dynarray<State> &partialStates) {
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
  newState.iv = _childrenData[ 4 ] & ~_childrenData[ G_MISS ];
  // Potentially ambiguous bits are those where the child is heterozygous.
  // If this is the first state on this chromosome, these are ambiguous and
  // their phase is resolved by later markers.
  newState.ambig = _childrenData[ G_HET ];
  // Have full transmission information except for children with missing data:
  newState.unassigned = _childrenData[ G_MISS ];
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
void Phaser::makeFullStates(const dynarray<State> &partialStates,
			    int firstMarker, int numChildren,
			    int numMissChildren, int8_t forceInform) {
  // make a new entry in <_hmm>
  // called "Abs" for "absolute" to distinguish from the <prevHMMIndex> field
  // in States
  int prevHMMIndexAbs = _hmm.length() - 1;
  _hmm.addEmpty();
  _hmmMarker.append(_curMarker); // new entry corresponds to <_curMarker>

  if (prevHMMIndexAbs < 0) { // no previous states
    addStatesNoPrev(partialStates, firstMarker);
    return;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Have previous states: add new states that transition from these

  // Because of forced informative markers, may need to transition from multiple
  // previous HMM indexes
  int stopHMMIndex = prevHMMIndexAbs;
  if (_lastRealInformIndex >= 0 && forceInform < 0)
    stopHMMIndex = _lastRealInformIndex;

  // Minimum and maximum recombination counts: for applying the optimization in
  // rmBadStatesCheckErrorFlag().
  float minMaxRec[2] = { FLT_MAX, 0 };

  // Count of the number of error states that get added by addStatesWithPrev()
  // below. We don't map these to the current marker (they are potential error
  // paths to the next and subsequent markers)
  int numPrevErrorStatesAdded = 0;

  // Typically, we transition from the states at the most recent informative
  // marker (whether forced or real) to the current marker. However, if the
  // previous marker was forced and the current marker is not, we also
  // transition from the most recent non-forced informative marker. This is
  // to allow for errors at the forced informative marker(s).
  for (int mapFromHMMIndex = prevHMMIndexAbs; true;
					      mapFromHMMIndex = stopHMMIndex) {
    // Only consider errors from "real" informative markers; stopHMMIndex is
    // "real" and if prevHMMIndexAbs != stopHMMIndex, prevHMMIndexAbs is not
    // real.
    addStatesWithPrev(partialStates, firstMarker, numChildren, numMissChildren,
		      forceInform, mapFromHMMIndex,
		      /*lastPrev=*/ mapFromHMMIndex == stopHMMIndex,
		      minMaxRec, numPrevErrorStatesAdded);

    if (mapFromHMMIndex == stopHMMIndex)
      break;
  }

  // reset "real" informative marker index
  if (_lastRealInformIndex >= 0 && forceInform < 0) {
    _lastRealInformIndex = -1;
  }
}

// helper function for makeFullStates() -- does the real work
//
// <lastPrev> indicates wheter this is the last previous HMM index we'll
// transition from; note that we may transition from multiple previous indexes
// if there are forced informative markers. Note that, when <lastPrev> is true
// the previous markers are at a "real" informative marker, while the others
// are for forced informative markers.
void Phaser::addStatesWithPrev(const dynarray<State> &partialStates,
			       int firstMarker, int numChildren,
			       int numMissChildren, int8_t forceInform,
			       int prevHMMIndexAbs, bool lastPrev,
			       float minMaxRec[2],
			       int &numPrevErrorStatesAdded) {
  // Final index is _hmm.length() - 1, but the caller ran _hmm.addEmpty(), so
  // the most recent index with states is -2:
  int mostRecentHMMIndex = _hmm.length() - 2;

  dynarray<State*> &prevStates = _hmm[prevHMMIndexAbs];
  uint32_t numPrev = prevStates.length();
  int numPartial = partialStates.length();

  //////////////////////////////////////////////////////////////////////////////
  // Map previous states to current states (i.e., at this marker) with no errors
  for(uint32_t prevIdx = 0; prevIdx < numPrev; prevIdx++) {
    if (prevHMMIndexAbs != mostRecentHMMIndex)
      // only need to map from the most recent marker; error-based transitions
      // happen below
      break;

    State *prevState = prevStates[prevIdx];

    // See comment above the isIVambigPar() method. We deal with
    // <IVambigPar> == 1 below and in updateStates()
    uint8_t IVambigPar = isIVambigPar(prevState->iv, prevState->ambig);
    // shouldn't force informative marker at site that is ambiguous for an
    // informative parent
    // Actually, we need to do this when both parents are in OHT region using
    // the condition above (where we don't skip a state)
    //assert(forceInform < 0 || !IVambigPar);

    // Does <prevState>, transition to any states with zero recombinations?
    // See code after next loop for why we track this
    bool zeroRecombsThisPrev = false; // initially assume

    for(int curIdx = 0; curIdx < numPartial; curIdx++) {
      const State &curPartial = partialStates[curIdx];
      if (forceInform >= 0) {
	if (curPartial.hetParent < 2 &&
	    (_missingPar & (1 << curPartial.hetParent)) == 0)
	  // have data for curPartial.hetParent: don't force informative marker
	  continue;
	else if (_missingPar == 3) {
	  if (curPartial.hetParent < 2) {
	    if (
		// only one haplotype transmitted for <curPartial.hetParent>
		((prevState->iv & _parBits[curPartial.hetParent]) == 0 ||
		 (prevState->iv & _parBits[curPartial.hetParent]) == _parBits[curPartial.hetParent]) &&
		// <curPartial.hetParent> has a more recent informative marker
		// than the other parent
		(prevState->numMarkersSinceNonHetPar[curPartial.hetParent] <
		 prevState->numMarkersSinceNonHetPar[1-curPartial.hetParent]) &&
		// <curPartial.hetParent>'s most recent informative marker is <
		// CmdLineOpts::forceInformInit markers ago, while the other
		// parent's is that far away
		(prevState->numMarkersSinceNonHetPar[curPartial.hetParent] <
		 CmdLineOpts::forceInformInit) &&
		(prevState->numMarkersSinceNonHetPar[1-curPartial.hetParent] >=
		 CmdLineOpts::forceInformInit)
	       ) {
	      // check other parent's potential for being forced
	      uint8_t otherHetParent = 1 - curPartial.hetParent;
	      uint64_t parHap = prevState->iv & _parBits[otherHetParent];
	      uint64_t oppParHap = parHap ^ _parBits[otherHetParent];
	      // 0, 1 or 2 bits set in parHap or oppParHap?
	      // * Only 1 corresponds to either all but one of the children or
	      //   all children except one receiving the same haplotype.
	      // * Also check for 2 since, in some cases, two children could
	      //   recombine near each other and produce a oneHapTrans scenario
	      uint8_t numTrans[2] = { (uint8_t) popcount(parHap),
				      (uint8_t) popcount(oppParHap) };
	      // which pattern (parHap/oppParHap) has min number of 1s
	      // (popcount)?
	      int minPat = (numTrans[0] <= numTrans[1]) ? 0 : 1;
	      if (numTrans[minPat] <= 2)
		// we want to force not curPartial.hetParent, but the other het
		// parent
		continue;
	    }
	  }
	}
      }

      if (IVambigPar && curPartial.hetParent == 1 && forceInform < 0)
	// When both parents are missing data, PI states do not actually
	// distinguish the two parents and so later FI markers will be
	// ambiguous for which parent is heterozygous. Here we arbitrarily pick
	// parent 0 as the heterozygous parent in these cases and omit states
	// with parent 1 het.
	// See the code at the beginning of this function for a similar choice
	// at the first informative marker.
	continue;
      mapPrevToFull(prevState,
		    /*prevHMMIndex=*/ mostRecentHMMIndex - prevHMMIndexAbs + 1,
		    prevIdx, /*error=*/ false, curPartial, minMaxRec,
		    IVambigPar,
		    /*numDataChildren=*/numChildren - numMissChildren,
		    _curMarker - _hmmMarker[ prevHMMIndexAbs ],
		    zeroRecombsThisPrev);
    }

    if (!zeroRecombsThisPrev || prevState->unassigned != 0) {
      // If the previous state transitioned with zero recombinations to a state
      // at the current marker, later paths should use that state always --
      // adding an error state will be suboptimal.
      // However, if the state transitioned with > 0 recombinations, there may
      // be a better path that skips the next state as an error and includes
      // <prevState>.
      // When <prevState> has unassigned haplotype transmissions, it will
      // not incur recombinations for those haplotypes, but may later recombine
      // more than it would if the next state were an error.
      // In either case, we'll consider error transitions from <prevState>.
      // This works by storing its indices in <_prevErrorStates>
      //
      // Note: we put these at the front so that states get transitioned from
      // in reverse order: states closer to the current one are preferred
      // (and we don't treat a state that's further back than some other as
      // ambiguous)
      _prevErrorStates.emplace_front(prevHMMIndexAbs, prevIdx);
      numPrevErrorStatesAdded++;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Introduce error states in prevStates if doing so saves at least
  // <CmdLineOpts::maxNoErrRecombs> recombinations
  // Only want to do this for "real" informative markers, not forced
  // informative markers, so we ensure <lastPrev> is true.
  if (CmdLineOpts::maxNoErrRecombs > 0 && lastPrev) {
    // Start <it> at begin() + <numPrevErrorStatesAdded> to skip that number of
    // states: these are for states added in the context of the current marker
    // and were just mapped from above (as non-error states)
    for(auto it = _prevErrorStates.begin() + numPrevErrorStatesAdded;
	     it != _prevErrorStates.end(); ) {
      int errorHMMIndex = it->first;

      // Only consider error paths that span fewer than a set number of markers
      // except for the immediately previous marker
      // Note that this interacts with the forced informative markers (see
      // checkForceInform()), and so we set the maximum distance to be a
      // meaningful number of sites more than this interval. (This is to ensure
      // that at least one forced informative marker has been added, presumably
      // at the current site, so that all previous markers are marked as
      // errors.)
      bool errorSpanTooLarge = _curMarker - _hmmMarker[ errorHMMIndex ] >
					  CmdLineOpts::forceInformInit + 50;
      // ... except we'll always allow for errors from the immediately previous
      // marker (meaning we'll allow transitions from one index before it):
      if (errorSpanTooLarge && errorHMMIndex != prevHMMIndexAbs -1) {
	it = _prevErrorStates.erase(it); // no need to consider this path later
	continue;
      }

      uint32_t errorPrevIdx = it->second;
      State *errorPathPrevState = _hmm[ errorHMMIndex ][ errorPrevIdx ];

      bool zeroRecombsThisPrev = false; // initially assume (see comments above)

      // See comment above the isIVambigPar() method. We deal with
      // <IVambigPar> == 1 below and in updateStates()
      uint8_t IVambigPar = isIVambigPar(errorPathPrevState->iv,
					errorPathPrevState->ambig);

      for(int curIdx = 0; curIdx < numPartial; curIdx++) {
	const State &curPartial = partialStates[curIdx];
	if (IVambigPar && curPartial.hetParent == 1)
	  continue; // See comment above in equivalent non-error code
	mapPrevToFull(errorPathPrevState,
		      /*prevHMMIndex=*/ mostRecentHMMIndex - errorHMMIndex + 1,
		      errorPrevIdx, /*error=*/ true, curPartial, minMaxRec,
		      IVambigPar,
		      /*numDataChildren=*/numChildren - numMissChildren,
		      _curMarker - _hmmMarker[ errorHMMIndex ],
		      zeroRecombsThisPrev);
      }

      // For next marker, need to erase the potential prev error states that
      // will be more than <CmdLineOpts::errorLength> informative markers
      // away, so here we erase those that are >= this distance
      // Alternatively, if we mapped this to the current marker with zero
      // recombinations, we should stop tracking the previous state:
      if (prevHMMIndexAbs - errorHMMIndex >= CmdLineOpts::errorLength ||
	  zeroRecombsThisPrev)
	it = _prevErrorStates.erase(it);
      else
	it++;
    }

    // Can have all the informative sites up to <_curMarker> be assigned as
    // erroneous. See comment above the <errorSpanTooLarge> variable for why
    // this bound is set at it is (interacts with forced informative markers)
    bool allowAllPrevSitesErrors =
		_hmmMarker.length() - 1 <= CmdLineOpts::errorLength &&
		(_curMarker - firstMarker) <= CmdLineOpts::forceInformInit + 50;
    if (allowAllPrevSitesErrors)
      addStatesNoPrev(partialStates, firstMarker, /*error=*/ true);
  }

  if (lastPrev)
    // only remove bad states, etc. when we've finished adding states to
    // the current HMM index:
    rmBadStatesCheckErrorFlag(_hmm[ mostRecentHMMIndex+1 ], minMaxRec,
			      numChildren);
}

// Add a state with no previous markers. This is for the very first informative
// marker on the chromosome and for adding a state that counts all the initial
// informative markers as errors
void Phaser::addStatesNoPrev(const dynarray<State> &partialStates,
			     int firstMarker, bool error) {
  int len = partialStates.length();
  assert(len <= 3); // can only have 2 forms of FI partial states and one PI
  for(int i = 0; i < len; i++) {
    if (i == 1 && partialStates[1].hetParent == 1
	       && partialStates[0].hetParent == 0
	       && _missingPar == 3)
      // For the very first marker, if we're not sure which parent is
      // heterozygous, and both are missing data, arbitrarily pick parent
      // 0 as such. Otherwise, since the two states are equivalent but
      // with opposite parent labels, there will be two equal paths
      // through the HMM.
      // TODO: document this behavior
      continue;

    // Partial states are equivalent to full states when there are no previous
    // states
    State *newState = new State(partialStates[i]);
    // how many indexes in the past to the previous state that maps here? This
    // index is *relative* to the current one, and when set to the current
    // number of indexes in the HMM, as below, implies that there are no
    // previous states to trace back to.
    newState->prevHMMIndex = _hmm.length();
    newState->ambigPrev = 0;

    // How many markers since each parent (0 or 1) passed by a heterozygous
    // marker?
    // First assume 0:
    newState->numMarkersSinceNonHetPar[0] =
                                    newState->numMarkersSinceNonHetPar[1] = 0;
    int numMarkersToCur = _curMarker - firstMarker + 1;
    // To match the behavior of other code, for any parent that is missing
    // data, we either
    // (a) assign <numMarkersToCur> if the parent is homozygous in the state
    // OR
    // (b) subtract 50 markers from <numMarkersToCur> if the parent is
    // heterozygous, ensuring it doesn't go negative:

    uint8_t isPI = newState->hetParent >> 1;
    for(int p = _firstMissP; p < _limitMissP; p++) {
      if (isPI || newState->hetParent == p)
        newState->numMarkersSinceNonHetPar[p] = max(0,
                                                    numMarkersToCur - 50);
      else
        newState->numMarkersSinceNonHetPar[p] = numMarkersToCur;
    }

    // Track the number of markers since the last one that is heterozygous for
    // only one parent (only need do this if only one parent is missing data)
    uint16_t numMarkersSinceOneHetPar = 0;
    if ( (_missingPar == 1 || _missingPar == 2) ) {
      // PI states are heterozygous for both parents, so count if <isPI> or
      // in a corner case (see comment near definition of <fakeOneHetPar> below)
      uint64_t fullIV = newState->iv;
      uint8_t hetParent = newState->hetParent;
      uint8_t nonHetPar = (1 - hetParent); // assuming !isPI
      bool fakeOneHetPar = ((_parentData >> (nonHetPar * 2)) & 3) == G_MISS &&
	    ( ( fullIV & ~_childrenData[G_MISS] & _parBits[nonHetPar]) == 0 ||
	      (~fullIV & ~_childrenData[G_MISS] & _parBits[nonHetPar]) == 0 );
      fakeOneHetPar = fakeOneHetPar && ( (_missingPar == (1 << hetParent)) ||
	  (( fullIV & _parBits[nonHetPar]) != 0 &&
	   (~fullIV & _parBits[nonHetPar]) != 0) );
      if (isPI || fakeOneHetPar)
	// Limit the impact of long stretches of uninformative markers on this
	// penalty: count a stretch as no more than 10 markers:
	numMarkersSinceOneHetPar = min(numMarkersToCur, 10);
    }
    newState->numMarkersSinceOneHetPar = numMarkersSinceOneHetPar;
    // initially probability of log(1) == 0; note that in principle this
    // should be log(1/N), where N is the number of initial states, but really
    // all that encodes is that they have equal probability, so setting this
    // to 0 is proportional:
    newState->maxLikelihood = 0;
    newState->maxPrevRecomb = 0;
    newState->ambigParHet = 1 << newState->hetParent;
    newState->parentPhase = 0;
    // 1 << 0, i.e., 1 << newState->parentPhase
    newState->ambigParPhase = (1 << 0) << (2 * newState->hetParent);
    // arbitrary if both parents are missing:
    newState->arbitraryPar = _missingPar == 3 && !_onChrX;
    if (error) {
      newState->error = 1;
      // see comment above State::minRecomb (phaser.h) for why this calculation:
      newState->minRecomb = CmdLineOpts::maxNoErrRecombs * 10 - 5;
      // Prefer paths that include a smaller number of markers indicated as
      // erroneous; do this by very slightly adjusting the error term depending
      // on the number of markers being spanned
      // Note: there are _hmm.length() - 1 markers before this one, and
      // _hmm.length() includes the current marker. We just applied a penalty
      // above, and we only want to further panelize if _hmm.length() > 2, so:
      newState->minRecomb += 1 * (_hmm.length() - 2);
      assert(_hmm.length() >= 2); // assumed above
    }
    else {
      newState->error = 0;
      newState->minRecomb = 0;
    }

    // TODO: optimization: in fact this may add two states with the same IV
    // ideally should lookup this state and only update if this one is better
    _hmm[ _hmm.length() - 1 ].append(newState);
  }
}

// If the IV values transmitted to each child by the two parents are either all
// identical or all opposite, then when <_missingPar> == 3 (both missing), which
// parent is heterozygous (for FI markers) will be ambiguous. Additionally, a
// recombination at a PI state that follows these ambiguous markers will have
// two parent phase assignments (opposite each other) that will produce
// equivalent numbers of recombinations.
uint8_t Phaser::isIVambigPar(uint64_t iv, uint64_t ambigIV,
			     uint64_t unassigned) {
  if (_onChrX) // only heterozygous parent on chrX is Mom, so no ambiguity
    return 0;
  uint64_t IVparDiff = (iv & _parBits[0]) ^ ((iv & _parBits[1]) >> 1);
  // Omit differences at standard ambiguous positions. This is relevant in two
  // contexts:
  // (1) at the beginning of a chromosome, PI states do not actually
  // distinguish the two parents and so later FI markers will be ambiguous for
  // which parent is heterozygous. The makeFullStates() function therefore uses
  // the return value from this function to decide whether to only introduce
  // FI for parent 0 states.
  // (2) similarly, just after a P region ends, a recombination at a PI state
  // will result in an ambiguous recombination that could be attributed to
  // either parent, and so later FI markers will be ambiguous for which parent
  // is heterozygous.
  // (note: only bit 0 in each child is set if there's a difference)
  uint64_t stdAmbig = (ambigIV & _parBits[1]) >> 1;
  // if <iv> hasn't been assigned for either parent (or both) for a given child,
  // we shouldn't factor in that part of the IV in determining whether the
  // IV is ambiguous in terms of parent assignments, so we collect all such
  // bits (putting them into position 0 as that's the only bit set).
  uint64_t unassignedBit0 = (unassigned & _parBits[0]) |
			    ((unassigned & _parBits[1]) >> 1);
  uint64_t ignoreMask = stdAmbig | unassignedBit0;
  IVparDiff &= ~ignoreMask;
  if (_missingPar == 3 && (IVparDiff == 0 ||
				      IVparDiff == (_parBits[0] & ~ignoreMask)))
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
// <prevHMMIndex> is the number of indexes previous to the current one that
// the previous state is from. When there are no errors, this is always 1.
void Phaser::mapPrevToFull(const State *prevState, uint8_t prevHMMIndex,
			   uint32_t prevIdx, bool error,
			   const State &curPartial, float minMaxRec[2],
			   uint8_t IVambigPar, int numDataChildren,
			   int numMarkersSincePrev, bool &zeroRecombsThisPrev) {
  uint8_t hetParent = curPartial.hetParent;

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
  uint8_t isPI = hetParent >> 1;

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
	     defaultPhaseHasRecomb, childPrevUnassigned, unambigHetRecombs);

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
    fullAmbig = _childrenData[G_MISS] & prevState->ambig;

    if (_missingPar == 3 && hetParent == 1 && prevState->unassigned > 0 &&
	isIVambigPar(fullIV, fullAmbig, fullUnassigned))
      // If one parent is in a one haplotype transmitted (OHT) state at the
      // _beginning_ of a chromosome, nothing constrains the code from assigning
      // the OHT parent to have the same IV as the other parent. In that case,
      // the IV of both parents will be equivalent, which is the ambigPar
      // condition (markers with 'P' status in output). Here we check if
      // this state will end up being a 'P' state, and, if it is assigned as
      // being heterozygous for parent 1, we omit it.
      // This parent 1 condition is arbitrary, but lines up with convention at
      // real 'P' markers/states: we only represent states that are heterozygous
      // for parent 0 at such states/markers.
      //
      // Note: if the truth is that both parents *did* transmit the same IV,
      // this will be detected at PI markers, and in fact, the marker *will* be
      // assigned as having P status.
      return;

    uint64_t hetParBits = _parBits[hetParent];
    hetParentUndefined = (prevState->unassigned & hetParBits) == hetParBits;
    oneParentUndefined = false; // not applicable to FI states
  }

  // What type of phase is the default? the alternative? For MT_PI states, the
  // alternative, which must have the same ambiguous bits, is 3 decimal == 11
  // binary (vs. default of 0). This corresponds to flipping both parents' phase
  // See just below for what this code does:
  uint8_t initPhaseType = 0;
  uint8_t altPhaseType = isPI * 3 + (1 - isPI);
  // Above equivalent to the line below but has no branching
  //uint8_t altPhaseType = (isPI) ? 3 : 1;

  // (3) As needed, remove apparent recombinations from <iv> values that were
  // ambiguous in the previous state
  // Which children were ambiguous in the previous state but not here?
  // Also various values relating to ambig1 type ambiguous values
  uint64_t stdAmbigOnlyPrev, ambig1PrevInfo, ambig1Unassigned;
  fixRecombFromAmbigIVambigPar(prevState, fullIV, fullAmbig, recombs, parRecombs,
		    isPI, IVambigPar,
		    /*ambigOnlyPrev=*/prevState->ambig & ~fullAmbig, hetParent,
		    stdAmbigOnlyPrev, ambig1PrevInfo, ambig1Unassigned,
		    initPhaseType, altPhaseType);

  // (4) Look up or create a full state with equivalent <iv> and <ambig> values
  // to <fullIV> and <fullAmbig>, and determine if <prevState> yields fewer
  // recombinations for these states than the currently stored previous state
  // (if any). If so, update the necessary values in the state.
  // Also examines an alternate phase type which may or may not map to an
  // equivalent state. If not, looks up that value, if so, compares the
  // recombinations for the two possibilities separately.
  updateStates(fullIV, fullAmbig, fullUnassigned, ambig1Unassigned, recombs,
	       prevState, stdAmbigOnlyPrev, ambig1PrevInfo, hetParent,
	       curPartial.homParentGeno, /*initParPhase=*/initPhaseType,
	       altPhaseType, prevHMMIndex, prevIdx, error, IVambigPar,
	       minMaxRec, hetParentUndefined, numDataChildren,
	       numMarkersSincePrev, zeroRecombsThisPrev);

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
    flipPIVals(fullIV, fullAmbig, propagateAmbig, unambigHetRecombs,
	       childPrevUnassigned, defaultPhaseHasRecomb);
    recombs = (prevState->iv ^ fullIV) & ~prevState->unassigned;
    parRecombs[0] = (recombs & _parBits[0]) * 3;
    parRecombs[1] = ((recombs & _parBits[1]) >> 1) * 3;

    initPhaseType = 1; // parent 0 flip
    altPhaseType = 2;  // parent 1 flip

    // TODO: potentially can optimize by storing information obtained in the
    // first call to fixRecombFromAmbig()
    fixRecombFromAmbigIVambigPar(prevState, fullIV, fullAmbig, recombs,
		      parRecombs, /*isPI=*/ 1, IVambigPar,
		      /*ambigOnlyPrev=*/ prevState->ambig & ~fullAmbig,
		      hetParent, stdAmbigOnlyPrev, ambig1PrevInfo,
		      ambig1Unassigned, initPhaseType, altPhaseType);

    // Following is not failing, so comment out for efficiency:
//    assert((recombs & _childrenData[G_HET] & (childPrevUnassigned[0] |
//						 childPrevUnassigned[1])) == 0);

    ///////////////////////////////////////////////////////////////////////////
    // Now ready to look up or create full states with the <fullIV> and
    // <fullAmbig> values, etc.
    updateStates(fullIV, fullAmbig, fullUnassigned, ambig1Unassigned, recombs,
		 prevState, stdAmbigOnlyPrev, ambig1PrevInfo, /*hetParent=*/2,
		 /*homParentGeno=*/G_MISS, /*initParPhase=*/ initPhaseType,
		 /*altPhaseType=*/ altPhaseType, prevHMMIndex, prevIdx,
		 error, IVambigPar, minMaxRec, /*hetParentUndefined=*/ false,
		 numDataChildren, numMarkersSincePrev, zeroRecombsThisPrev);
  }
}

// See long comment in makeFullStates() -- handles heterozygous children and
// previously unassigned <iv> values at MT_PI states:
void Phaser::handlePI(const State *prevState, uint64_t &fullIV,
		      uint64_t &fullAmbig, uint64_t &recombs,
		      uint64_t parRecombs[2], uint64_t &propagateAmbig,
		      uint64_t &defaultPhaseHasRecomb,
		      uint64_t childPrevUnassigned[2],
		      uint64_t unambigHetRecombs[4]) {
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

  // In certain circumstances, we assign only one parent (by convention parent
  // 0) as ambiguous. This is to deal with complexities around when only one
  // parent has an assigned <iv> value and a PI marker occurs in which a child
  // is heterozygous. Such a marker and heterozygous genotype gives information
  // about allelic transmissions from both parents, but switching to the
  // opposite <iv> will only incur one recombination relative to the one parent
  // whose <iv> value has been assigned, not two. So, as distinct from standard
  // ambiguous <iv> values, there are two <iv> values, one of which gives 0
  // recombinations relative to the previous markers, the other of which,
  // despite being inverted for both parents, gives only 1 recombination. We
  // call these ambig1 and deal with them separately. Note that they are
  // similar to unambigous <iv> values in that once established (i.e., after
  // the single ambig bit is set) we seek to match the <iv> values to the prior
  // marker while propagating forward the ambig1 status. Thus, in particular,
  // they should be considered by the calcHetChildPIRecombs() method
  uint64_t prevStdAmbig = ((prevState->ambig & _parBits[1]) >> 1) * 3;

  // Determine the recombination values for heterozygous children that are
  // unambiguous and fully assigned in <prevState>. To avoid recombination on
  // both homologs on these children, we must sometimes flip the <fullIV>
  // values. Also, a single recombination is ambiguous and must be accounted
  // for. See below.
  uint64_t unambigAssignedHets = _childrenData[G_HET] &
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
	      (_childrenData[G_MISS] | _childrenData[G_HET]) & prevState->ambig;
  // (2) Near the beginning of the chromosome, if a child has not had either
  //     <iv> value assigned previously and is het here, it is ambiguous. We
  //     can choose either phase assignment without affecting the number of
  //     recombinations.
  //   (We put this together with propagateAmbig so that it gets reused for the
  //   alternate phase types addressed in flipPIVals().)
  propagateAmbig |= _childrenData[G_HET] & bothPrevUnassigned;
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
  uint64_t hetChildPrevUnassign = _childrenData[G_HET] & onePrevUnassigned;
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
//			~(_childrenData[G_HET] &
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
			uint64_t propagateAmbig,
			const uint64_t unambigHetRecombs[4],
			const uint64_t childPrevUnassigned[2],
			uint64_t defaultPhaseHasRecomb) {
  ///////////////////////////////////////////////////////////////////////
  // Get the flipped <fullIV> value:

  // <iv> value is the same as the default for all missing data children:
  uint64_t newFullIV = fullIV & _childrenData[G_MISS];
  // flip the transmitted haplotype for parent 0 for all homozygous
  // children
  newFullIV |= (fullIV ^ _parBits[0]) &
				(_childrenData[G_HOM0] | _childrenData[G_HOM1]);
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
  uint64_t flipPartlyUnassigned = _childrenData[G_HET] &
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
void Phaser::fixRecombFromAmbigIVambigPar(const State *prevState,
			    uint64_t &fullIV, uint64_t fullAmbig,
			    uint64_t &recombs, const uint64_t parRecombs[2],
			    uint8_t isPI, uint8_t IVambigPar,
			    uint64_t ambigOnlyPrev, uint8_t hetParent,
			    uint64_t &stdAmbigOnlyPrev, uint64_t &ambig1PrevInfo,
			    uint64_t &ambig1Unassigned, uint8_t &initPhaseType,
			    uint8_t &altPhaseType) {
  // As detailed elsewhere, have ambig1 cases where only parent bit 0 is set
  // ambiguous that are distinct from standard ambiguous cases
  uint64_t parAmbigOnlyPrev[2] = { (ambigOnlyPrev & _parBits[0]) * 3,
				   ((ambigOnlyPrev & _parBits[1]) >> 1) * 3 };
  stdAmbigOnlyPrev = parAmbigOnlyPrev[1];
  uint64_t ambig1OnlyPrev = parAmbigOnlyPrev[0] - parAmbigOnlyPrev[1];
  uint64_t anyAmbigOnlyPrev = parAmbigOnlyPrev[0];

  // TODO: test optimization: don't call this function if ambigOnlyPrev == 0?
  //       actually also need !(isPI && IVambigPar) [see below]

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
  // <iv> value assigned at the previous marker). But when we consider the
  // alternate phase type in updateStates(), the children that recombined/didn't
  // will be inverted. We therefore use <ambig1PrevInfo> to store which
  // ambig1Unassigned bits to flip later in updateStates():
  ambig1PrevInfo += (1 - isPI) * ambig1OnlyPrev & _parBits[homozyParent];
  ambig1Unassigned = (1 - isPI) * (ambig1PrevInfo &
					      parRecombsFromAmbig[ hetParIdx ]);

  // When <IVambigPar> is true at a PI state, if two alternative phasings have
  // equal numbers of recombinations, need the initial phase to attribute more
  // recombinations to parent 0*. This is because we arbitrarily choose to only
  // track parent 0 het states when <IVambigPar> (i.e., when the heterozygous
  // parent is ambiguous).
  // * Actually only need consider the recombinations in homozygous children:
  //   the recombinations in children that are heterozygous can be attributed
  //   to either parent (are ambiguous).
  if (isPI && IVambigPar) {
    // Get recombinations under the alternate (flipped for both parents) phase
    // type:
    // code adapted from that in the latter part of updateStates() (where the IV
    // and other values are flipped to those for the alternate phase type)
    uint64_t noFlipBits = _childrenData[G_MISS] | _childrenData[G_HET];
    uint64_t flipVal = _parBits[2] & ~noFlipBits;
    uint64_t altRecombs = recombs ^
      (flipVal & ~(stdAmbigOnlyPrev | ambig1PrevInfo | prevState->unassigned));
    if (popcount(recombs) == popcount(altRecombs)) {
      // Ambiguous states! May want to flip the initial phase type:
      // Note: G_MISS won't recombine, so not G_HET suffices to get the G_HOM*
      // child recombinations
      uint64_t homChildRecombs = recombs & ~_childrenData[G_HET];
      uint8_t numHomChildRecombsPar0 = popcount(homChildRecombs & _parBits[0]);
      uint64_t altHomChildRecombs = altRecombs & ~_childrenData[G_HET];
      uint8_t numAltHomChildRecombsPar0 =
				    popcount(altHomChildRecombs & _parBits[0]);
      // the following is true when the alternate phase attributes more
      // recombinations to parent 0; since we only keep the parent 0 het state
      // when <IVambigPar>, using such a phase for the PI state means the two
      // can map to the same IV (this won't always happen, but when it can
      // happen, we prefer to have such a consistent IV), so we make the
      // alternate phase the default below
      bool useAltPhase = numAltHomChildRecombsPar0 > numHomChildRecombsPar0;

      if (!useAltPhase) {
	// One more scenario with the same motive (attributing recombinations
	// to parent 0) as just above. We prefer to match IV values at
	// previously ambiguous children (because code at the end of this
	// function ensures that such IVs get flipped when needed so that the
	// recombinations are attributed to parent 0*). But these recombinations
	// were erased above. We'll explicitly count them here and prefer the
	// phase that has fewer of them.
	// * Note: that code would have executed at a previous marker where the
	// child(ren) in question were heterozygous and ambiguously recombined.

	// The recombinations that were erased from (standard) ambiguous IVs
	// are those that recombined on both parents:
	uint8_t numHiddenRecombs = popcount(bothParRecomb);

	// The alternate phase will have the exact opposite pattern of
	// recombinations for <stdAmbigOnlyPrev> IV values. (This is because
	// such children are homozygous here and the alternate phase has both
	// parents flipped.) More specifically, if a child had both parents
	// recombined for the default phase type, neither will be recombined
	// for the alternate _and_ vice versa:
	uint64_t altHiddenRecombs = stdAmbigOnlyPrev &
			      ~(parRecombsFromAmbig[0] | parRecombsFromAmbig[1]);
	uint8_t numAltHiddenRecombs = popcount(altHiddenRecombs);
	useAltPhase = numHiddenRecombs > numAltHiddenRecombs;
      }

      if (useAltPhase) {
	uint64_t altIV = fullIV ^ flipVal;
	fullIV = altIV;
	// No changes to fullAmbig, fullUnassigned, or ambig1Unassigned needed
	recombs = altRecombs;
	// swap phase types:
	uint8_t tmp = initPhaseType;
	initPhaseType = altPhaseType;
	altPhaseType = tmp;
      }
    }

    // When transitioning from an IVambigPar state, we want to ensure that any
    // ambiguous recombinations get attributed to parent 0, so we flip the IV
    // for a child with any such recombinations.
    // In fact, the below also flips the IV if _both_ children recombine.
    // This ensures that the IV flip resulting from the ambiguous recombination
    // gets propagated to downstream PI states where the child remains
    // heterozygous.
    // We use _parBits[1] here because we want standard ambiguous children, not
    // ambig1 ambiguous children
    uint64_t stdAmbig = ((fullAmbig & _parBits[1]) >> 1) * 3;
    uint64_t parent1AmbigRecombs = (recombs & stdAmbig) & _parBits[1];
    uint64_t par1RecFlipVal = (parent1AmbigRecombs >> 1) * 3;
    fullIV ^= par1RecFlipVal;
    recombs ^= par1RecFlipVal;
  }
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
			  uint64_t recombs, const State *prevState,
			  uint64_t stdAmbigOnlyPrev, uint64_t ambig1PrevInfo,
			  uint8_t hetParent, uint8_t homParentGeno,
			  uint8_t initParPhase, uint8_t altPhaseType,
			  uint8_t prevHMMIndex, uint32_t prevIndex,
			  bool error, uint8_t IVambigPar, float minMaxRec[2],
			  bool hetParentUndefined, int numDataChildren,
			  int numMarkersSincePrev, bool &zeroRecombsThisPrev) {
  // How many iterations of the loop? See various comments below.
  int numIter = 2;

  // How many recombinations for the initial phase assignment?
  // Element 0 stores total;
  // When needed (below), element 1 stores the recombination count for parent 1
  size_t numRecombs[2];
  numRecombs[0] = popcount(recombs);
  uint8_t curParPhase = initParPhase;

  // Note: (hetParent >> 1) == 1 iff hetParent == 2. It is 0 for the other
  //       values.
  uint8_t isPI = hetParent >> 1;

  // Is it possible to swap the phase of this state and obtain the same
  // number of recombinations from the previous state? Decided below
  uint8_t ambigLocal = 0;

  // Track the number of markers since the last one that is heterozygous for
  // only one parent (only need do this if only one parent is missing data)
  uint16_t numMarkersSinceOneHetPar = 0;
  if (_missingPar == 1 || _missingPar == 2) {
    // Case where we have data for one parent; attempt to detect switch in
    // which parent is which: this manifests as long stretches where both
    // parents are heterozygous. Thus, we count markers below (in
    // <numMarkersSinceOneHetPar>) if <isPI> or in a corner case:

    // In addition to PI markers, at a FI1 state at a site that is (locally or
    // genome-wide) missing data for the het parent: if the IV has the non-het
    // parent transmitting all the same haplotype, a PI state and the FI1 state
    // will be indistinguishable. We don't want such a non-PI state to
    // interrupt a run of PI state counts.
    //
    // Thus we define <fakeOneHetPar> below according to the conditions in the
    // above paragraph: (a) the (<nonHetPar>) is missing, and (b) that
    // parent transmitted the same haplotype to all non-missing children.
    uint8_t nonHetPar = (1 - hetParent); // assuming !isPI
    bool fakeOneHetPar = ((_parentData >> (nonHetPar * 2)) & 3) == G_MISS &&
	    ( ( fullIV & ~_childrenData[G_MISS] & _parBits[nonHetPar]) == 0 ||
	      (~fullIV & ~_childrenData[G_MISS] & _parBits[nonHetPar]) == 0 );
    // ... also we don't want this case to trigger if we're in a OHT region for
    // the fully missing data parent. Thus we also require either that (a) the
    // fully missing data parent is heterozygous (implying we're not in a OHT
    // region, since in a OHT region the missing data parent has ambiguous
    // heterozygosity), or (b) that the non-het parent (necessarily fully
    // missing since this is an OR case) is _not_ in an OHT state once we
    // include _all_ children (even those who may be locally missing).
    fakeOneHetPar = fakeOneHetPar && ( (_missingPar == (1 << hetParent)) ||
	  (( fullIV & _parBits[nonHetPar]) != 0 &&
	   (~fullIV & _parBits[nonHetPar]) != 0) );
    if (isPI || fakeOneHetPar) {
      // Limit the impact of long stretches of uninformative markers on this
      // penalty: count a stretch as no more than 10 markers:
      numMarkersSinceOneHetPar = prevState->numMarkersSinceOneHetPar +
						  min(numMarkersSincePrev, 10);
      // NOTE: it may not be necessary, and we currently don't check, but
      //       this penalty should really only apply if the commented out
      //       conditional below holds. It's complicated by ambiguous IVs
      //       and the fact that the four phase options for this PI state
      //       behave differently with respect to this check
//      uint8_t dataPar = 2 - _missingPar; // the parent with data
//      uint8_t dataParIV = fullIV & _parBits[dataPar];
//      if (!( dataParIV == _parBits[ dataPar ] || dataParIV == 0 )) {
//        // If the data parent transmitted the same haplotype to all children,
//        // then all sites where the children are heterozygous will be PI, so
//        // there shouldn't be a penalty. Here, we've inverted those cases and
//        // there should be a penalty:
//      }
    }
    else {
      // don't just reset <numMarkersSinceOneHetPar> to 0 when we encounter one
      // non-PI state. That state could be the result of an error. Instead,
      // decrement by 10, capping at 0
      numMarkersSinceOneHetPar = max(0,
				     prevState->numMarkersSinceOneHetPar - 10);
    }
  }

  uint8_t ohpPenalty = 0;
  if (numMarkersSinceOneHetPar > CmdLineOpts::bothParHetThreshold) {
    ohpPenalty = 1;
    numMarkersSinceOneHetPar = 0; // penalty will be applied below
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
	     (_childrenData[G_MISS] & _parBits[hetParent] & ~fullUnassigned) ==0
	     && stdAmbigOnlyPrev == 0 && ambig1Unassigned == 0) ||
       (isPI && (stdAmbigOnlyPrev & _childrenData[G_HET]) == 0
	     && ((_childrenData[G_MISS] & ~fullUnassigned) |
					  _childrenData[G_HET]) == fullAmbig)) {
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

    if (_phaseMethod == PHASE_MINREC) { // minimum recombinant
      if (curCount[0] < numRecombs[0]) {
	numRecombs[0] = curCount[0];
	curParPhase = altPhaseType;
	fullIV ^= _parBits[ hetParent ] & ~fullAmbig;
	// No need to execute this as we ensured that
	// (1-isPI) * ambig1PrevInfo == 0
//	ambig1Unassigned ^= (1 - isPI) * ambig1PrevInfo;
      }
      else if (curCount[0] == numRecombs[0]) {
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

  State *lastState = NULL;

  // Examine the 1 or 2 needed phase assignments
  for(int i = 0; i < numIter; i++) {
    // Which bit corresponds to the lowest order unambiguous child?
    int lowOrderChildBit;

    // First look up or create a full state with equivalent <iv> and <ambig>
    // values to <fullIV> and <fullAmbig>
    State *theState = lookupState(fullIV, fullAmbig,
				  fullUnassigned | ambig1Unassigned,
				  lowOrderChildBit);
    // TODO: test optimization: add a check for prevState->minRecomb >
    // theState->minRecomb (or analogous comparison for likelihood)

    // Does the current previous state lead to minimum recombinations for
    // <theState>? (Could be either a new minimum or an ambiguous one)
    // Assigned below
    bool theStateUpdated;

    // Should never map to the same state as the previous iteration over <i>
    // (this should be caught in the code above that sets numIter = 1)
    assert(lastState == NULL || theState != lastState);
    lastState = theState;

    uint32_t totalRecombs = 0;
    float totalLikehood = -FLT_MAX;
    if (_phaseMethod == PHASE_MINREC) {
      // see comment above State::minRecomb (phaser.h) for why it's *10:
      totalRecombs = prevState->minRecomb + numRecombs[0] * 10;

      if (error) {
	// error state: apply penalty see comment above State::minRecomb
	// (phaser.h) for why this:
	totalRecombs += CmdLineOpts::maxNoErrRecombs * 10 - 5;
	// Prefer paths that include a smaller number of markers indicated as
	// erroneous; do this by very slightly adjusting the error term
	// depending on the number of markers being spanned. Note that
	// prevHMMIndex > 1 for all erroneous states (we're skipping at least
	// 1 marker; if we weren't, we wouldn't be in an error state and
	// prevHMMIndex == 1). We just applied a penalty for 1 marker, and we
	// only want to further penalize if prevHMMIndex > 2, so:
	totalRecombs += 1 * (prevHMMIndex - 2);
	assert(prevHMMIndex >= 2); // assumed above
      }

      // apply penalty for children's IVs switching from one parent to other
      totalRecombs += ohpPenalty * 10;
    }
    else
      // TODO! this (currently non-functional) code doesn't account for the
      //       missing parent no het penalty and/or the error penalty
      totalLikehood = prevState->maxLikelihood + localLikehood;

    theStateUpdated = checkMinUpdate(fullIV, fullUnassigned, ambig1Unassigned,
				     theState, prevState, hetParent,
				     homParentGeno, curParPhase, altPhaseType,
				     ambigLocal, prevHMMIndex, prevIndex, error,
				     IVambigPar, minMaxRec, numMarkersSincePrev,
				     numMarkersSinceOneHetPar,
				     totalRecombs, totalLikehood,
				     numRecombs[0], lowOrderChildBit);
    zeroRecombsThisPrev = zeroRecombsThisPrev || numRecombs[0] == 0;

    if (i == 0 && numIter == 2) {
      // Update the various values as needed for the inverted phase in the next
      // iteration.

      // Children that have missing data should not be flipped
      // Also, for both parent het markers, heterozygous children's bits do
      // not get flipped: they're either unambigous and constrained by
      // <prevState->iv> or they're ambiguous and should remain the same to
      // stick with the convention used in <_stateHash>.
      uint64_t noFlipBits = _childrenData[G_MISS] |
			    (_childrenData[G_HET] * isPI);
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
	if (theStateUpdated)
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
			    const uint64_t unassigned,
			    int &lowOrderChildBit) {
  uint64_t lookupIV = getCanonicalIV(iv, allAmbig, unassigned,
				     lowOrderChildBit);

  //////////////////////////////////////////////////////////////////////////
  // Do the lookup

  iv_ambig_real theKey(lookupIV, allAmbig, unassigned);
  state_ht_iter it = _stateHash.find( &theKey );
  if (it == _stateHash.end()) {
    // need to create state
    State *newState = new State;
    newState->ambig = allAmbig;
    newState->unassigned = unassigned;
    newState->minRecomb = UINT32_MAX;
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

// Returns the canonical inheritance vector value for <iv>. See the comment
// just inside the function for how this is defined.
uint64_t Phaser::getCanonicalIV(const uint64_t iv, const uint64_t allAmbig,
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
  //     XXX actually there is one exception: when transitioning out of an
  //     IVambigPar region to a PI state, we assign the ambiguous IVs so that
  //     the recombination gets attributed to parent 0. This is to stick with
  //     a convention adopted throughout that we start with parent 0 het states
  //     upon leaving IVambigPar (aka "P") regions.
  // (3) Any unassigned IV values are 0.

  uint64_t unambig, ambigStd;
  lowOrderChildBit = lowOrderUnambigAssignedBit(allAmbig, unassigned,
						unambig, ambigStd);

  // The genotype of the child tells us what bits we need to flip: the exact
  // bits that are assigned 1 need to be flipped in all unambiguous children.
  uint8_t flipType = (iv >> lowOrderChildBit) & 3;

  //////////////////////////////////////////////////////////////////////////
  // Get the canonical key value

  // Conveniently, we've got _flips indexed by the 4 possible flip types with
  // the values to flip assigned in each child:
  uint64_t canonIV = iv ^ (_flips[flipType] & unambig);
  // And we've done something analogous for ambiguous bits; must flip them too:
  canonIV ^= _ambigFlips[flipType] & ambigStd;

  // Lastly, ensure that any unassigned IV values are 0:
  canonIV &= ~unassigned;

  return canonIV;
}

int Phaser::lowOrderUnambigAssignedBit(const uint64_t allAmbig,
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
bool Phaser::checkMinUpdate(uint64_t fullIV, uint64_t fullUnassigned,
			    uint64_t ambig1Unassigned, State *theState,
			    const State *prevState, uint8_t hetParent,
			    uint8_t homParentGeno, uint8_t curParPhase,
			    uint8_t altPhaseType, uint8_t ambigLocal,
			    uint8_t prevHMMIndex, uint32_t prevIndex,
			    bool error, uint8_t IVambigPar, float minMaxRec[2],
			    int numMarkersSincePrev,
			    uint16_t numMarkersSinceOneHetPar,
			    uint32_t totalRecombs, float totalLikehood,
			    size_t numRecombs, int lowOrderChildBit) {
  // Does the current previous state lead to minimum recombinations for
  // <theState>? (Could be either a new minimum or an ambiguous one)
  bool theStateUpdated = false;

  uint8_t isPI = hetParent >> 1;

  // Returns -1 for the state defined by the arguments above being better
  // than <theState>, 0 for them being equal (leading to an ambiguity), and
  // 1 for <theState> being better
  int8_t whichOptimal = decideOptimalState(theState, prevState, prevHMMIndex,
					   error, totalRecombs, totalLikehood,
					   numRecombs);

  if (whichOptimal < 0) {
    theStateUpdated = true;

    theState->iv = fullIV;
    // next two already assigned in lookupState():
    //theState->ambig = fullAmbig;
    //theState->unassigned = fullUnassigned | ambig1Unassigned;
    theState->prevHMMIndex = prevHMMIndex;
    theState->prevState = prevIndex;
    if (error) // erroneous previous state
      theState->error = 1;
    else {
      // Propagate error state information: if the path of states leading to
      // this one has an error in it, indicate this by having
      // <theState->error> == 2. So we'll copy forward the value of previous
      // value of 2 or convert a previous value of 1 into 2:
      // The following produces 0 if <prevState->error> == 0 and 2 otherwise:
      theState->error = (prevState->error & 2) | ((prevState->error & 1) << 1);
    }
    theState->ambigPrev = 0;
    theState->minRecomb = totalRecombs;
    // assign <numMarkersSinceNonHetPar>:
    for (int p = _firstMissP; p < _limitMissP; p++) {
      // <p> is missing; update <theState->numMarkersSinceNonHetPar[p]>
      int totalNumSincePrev = numMarkersSincePrev;
      if (prevState->hetParent == 1 - p)
	// prevState is not heterozygous for <p>: sum prev count
	totalNumSincePrev += prevState->numMarkersSinceNonHetPar[p];
      theState->numMarkersSinceNonHetPar[p] = totalNumSincePrev;
    }
    theState->numMarkersSinceOneHetPar = numMarkersSinceOneHetPar;
    theState->maxLikelihood = totalLikehood;
    // see comment above State::minRecomb (phaser.h) for why it's *10:
    theState->maxPrevRecomb = numRecombs * 10;
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
  else if (whichOptimal == 0) {
    assert( theState->prevHMMIndex == prevHMMIndex );
    theStateUpdated = true;

    bool newBestPrev = false;
    if (numRecombs > theState->maxPrevRecomb) {
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
      // see comment above State::minRecomb (phaser.h) for why it's *10:
      theState->maxPrevRecomb = numRecombs * 10;
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
      uint8_t flipType = ((theState->iv ^ fullIV) >> lowOrderChildBit) & 3;
      // if hetParent < 2, parentPhase flips differently: only have two
      // possibilities: flipped (1) or not (0):
      flipType = isPI * flipType + (1 - isPI) * ((flipType >> hetParent) & 1);
      curParPhase ^= flipType;
      altPhaseType ^= flipType;
    }

    if (theState->homParentGeno == G_MISS && !_onChrX) {
      theState->homParentGeno = homParentGeno;
      // A reason that the homParentGeno might be missing even for an FI
      // state:
      // on the X chromosome, we try to impute Dad's genotype using the
      // daughters. Sometimes both possibilities remain and we consider states
      // of both <homParentGeno> values. Often one of these is better than
      // the other, but in some corner cases (see backtrace() for a detailed
      // example), there can be an ambiguity. We'll assign the genotype here
      // and, in later iterations, if the imputed genotype for Dad is
      // different from the current one, the else branch just below will set
      // him to be missing. Also see the comment above assertion below.
    }
    else if (homParentGeno != G_MISS &&
				      theState->homParentGeno != homParentGeno)
      theState->homParentGeno = G_MISS;
    // In general, we expect that, unless <isPI> (so <hetParent> == 2), the
    // <homParentGeno> should match across different ambiguous assignments of
    // states here. However, the X chromosome is a bit more involved, and as
    // needed, we'll set <homParentGeno> to missing (if there are conflicts)
    // and therefore the latter part of the assertion need not hold if <_onChrX>
    assert(_onChrX || hetParent == 2 ||
				      theState->homParentGeno == homParentGeno);

    // Update <numMarkersSinceNonHetPar>:
    // Have two values for <numMarkersSinceNonHetPar[p]>, the one in the state
    // and the one for the new path
    // Will be conservative and use the value closer to 0
    for (int p = _firstMissP; p < _limitMissP; p++) {
      if (prevState->hetParent == 1 - p) // prev not het for <p>: sum prev count
	theState->numMarkersSinceNonHetPar[p] =
		  min(prevState->numMarkersSinceNonHetPar[p] +
							  numMarkersSincePrev,
		      theState->numMarkersSinceNonHetPar[p]);
      else // prev is het for <p>: only accumulate <numMarkersSincePrev> for <p>
	theState->numMarkersSinceNonHetPar[p] =
				    min(numMarkersSincePrev,
					theState->numMarkersSinceNonHetPar[p]);
    }
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

  if (theStateUpdated && error) {
    // So that we don't get additional forced one hap trans states when:
    // (1) we're getting to this state via an error and
    // (2) a state was just added for this previous state (i.e., between the
    //     current marker and the <prevState>'s marker)
    // ... set the number of markers since non het parent to be
    // CmdLineOpts::forceInformSeparation less than CmdLineOpts::forceInformInit
    // That way, we can get another forced informative marker every
    // CmdLineOpts::forceInformSeparation markers.

    // the (first) error marker is always immediately after <prevHMMIndex>:
    int curHMMIndex = _hmmMarker.length() - 1;
    int errorHMMIndex = curHMMIndex - prevHMMIndex + 1;
    int errorInformMarker = _hmmMarker[errorHMMIndex];

    // We want to reduce numMarkersSinceNonHetPar to no more than the sum of:
    // (a) the distance to the error (don't accumulate since the error itself
    // isn't evidence for not having an informative marker) and (b)
    // (CmdLineOpts::forceInformInit - CmdLineOpts::forceInformSeparation) for
    // the reasons given at the top of this block.
    int numMarkersSinceError = _curMarker - errorInformMarker;
    int reducedNumSince = CmdLineOpts::forceInformInit -
			  CmdLineOpts::forceInformSeparation +
			  numMarkersSinceError;
    // We'll reduce both parents by the same amount and determine that amount
    // here:
    int reduceBy = 0;
    for (uint8_t p = _firstMissP; p < _limitMissP; p++) {
      if (reducedNumSince < theState->numMarkersSinceNonHetPar[p])
	reduceBy = max(reduceBy,
	            // we use this number so that when we do
	            // numMarkersSinceNonHetPar[p] - reduceBy below, we get
	            // numMarkersSinceNonHetPar[p] = reducedNumSince
		    theState->numMarkersSinceNonHetPar[p] - reducedNumSince);
    }
    for (uint8_t p = _firstMissP; p < _limitMissP; p++)
      theState->numMarkersSinceNonHetPar[p] -= reduceBy;
  }

  return theStateUpdated;
}

// Returns -1 if the proposed <prevState> is better than the current one stored
// in <theState>
// Returns 0 if they are ambiguous (equally optimal)
// Returns 1 if the proposed <prevState> is suboptimal
int8_t Phaser::decideOptimalState(State *theState, const State *prevState,
				  uint8_t prevHMMIndex, bool error,
				  uint32_t totalRecombs, float totalLikehood,
				  size_t newRecombs) {
  if (_phaseMethod == PHASE_MINREC) {
    if (totalRecombs < theState->minRecomb)
      return -1; // fewer recombinations: optimal
    else if (totalRecombs == theState->minRecomb) {
      // TODO: document that we prefer paths without errors
      // equal recombinations, so arguably ambiguous, BUT
      // we prefer paths without errors AND
      // when there are errors, we prefer paths that place recombinations
      // earlier (which can lead to fewer error markers)

      // Does the path via the previous state include an error? Yes if the
      // previous state indicates an error or if there's a new error
      bool prevPathError = prevState->error != 0 || error;
      // Is the current state an error state and/or does it include an error
      // state in its previous state?
      bool curIsError = theState->error > 0;

      if (!curIsError) {
	if (prevPathError)
	  return 1; // new has errors: is suboptimal
	else { // both non-error
	  // prefer paths that skip (forced informative) markers
	  if (prevHMMIndex > theState->prevHMMIndex)
	    return -1;
	  else if (prevHMMIndex == theState->prevHMMIndex)
	    return 0; // ambiguous
	  else
	    return 1;
	}
      }
      else { // curIsError
	if (!prevPathError)
	  return -1; // <theState> has errors, but new doesn't: is optimal
	else { // prevPathError
	  // brand new error in <theState>?
	  if (theState->error == 1) {
	    // need new error for this to be ambiguous
	    if (error) {
	      // prefer paths that skip (forced informative) markers
	      if (prevHMMIndex > theState->prevHMMIndex)
		return -1;
	      else if (prevHMMIndex == theState->prevHMMIndex)
		return 0; // ambiguous
	      else
		return 1;
	    }
	    else
	      return -1; // fewer error markers in new: is optimal
	  }
	  // else: older error in <theState>
	  if (error)
	    // brand new error in new prev; we'll arbitrarily make it optimal:
	    // we like recombinations happening earlier, and a new error path
	    // should typically have recombinations placed earlier than one
	    // with an error before
	    return -1;

	  // prefer paths that skip (forced informative) markers
	  if (prevHMMIndex > theState->prevHMMIndex)
	    return -1;
	  else if (prevHMMIndex == theState->prevHMMIndex)
	    return 0; // ambiguous
	  else
	    return 1;
	}
      }
    }
    else
      return 1;
  }
  else { // _phaseMethod == PHASE_MAXLIKE
    // TODO: (3) Error case for max likelihood? Relevant to this if statement
    // and following else statement
    // TODO: need a tolerance for equality comparison of likelihoods

    if (theState->maxLikelihood < totalLikehood)
      return -1;
    else if (theState->maxLikelihood == totalLikehood)
      return 0;
    else
      return 1;
  }
}

// When there are multiple previous states that can optimally reach <theState>,
// add the corresponding previous state to the list of these previous states
// if such a list exists; otherwise create one
void Phaser::updateAmbigPrev(State *theState, uint32_t prevState,
			     bool newBestPrev) {
  const uint32_t MAX_IDX_IN_STATE = 1 << 6;

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
	// is not located at bit 30 (see below in case 0 for info on end bit)
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
// Note that the maximum value may be stale and may not be present in any
// state, but we track it as an optimization and avoid going through the states
// when this value is too low to bother exploring a potential state path.
void Phaser::rmBadStatesCheckErrorFlag(dynarray<State*> &curStates,
				       float minMaxRec[2], int numChildren) {
  bool allError2 = true; // assume all states have <error> == 2 initially
  int shiftToState = -1;

  // TODO: potential optimizations:
  // 1. If we know there are no added errors this time and that we cleared the
  // errors in the last iteration, could avoid checking for them
  // 2. If we know something about the <IV> of the minimum state, may be able
  // to throw out more states -- 2 * numChildren is potentially conservative

  if (minMaxRec[1] >= minMaxRec[0] + 2 * numChildren) {
    for(int i = 0; i < curStates.length(); i++) {
      // see comment above State::minRecomb (phaser.h) for why *10:
      if (curStates[i]->minRecomb >= minMaxRec[0] + 2 * 10 * numChildren) {
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
void Phaser::backtrace(NuclearFamily *theFam, int chrFirstMarker,
		       int chrLastMarker, FILE *log) {
  int lastIndex = _hmm.length() - 1;

  if (lastIndex < 0)
    // no data: done (can happen if the children are missing data almost
    // everywhere and the markers are ambiguous)
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
  if (_curBTAmbigSet->size() == 1 && _missingPar == 3) {
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
  bool assignOneHapTrans[2] = { false, false };
  uint8_t oneHapTransHap[2] = { 0, 0 };
  uint64_t oneHapTransOutlierIV[2] = { 0, 0 };

  // To avoid including equivalent state paths where the parents are
  // flipped, we track the canonical IVs we've already applied at each state
  //
  // did we see an inverted IV (see below for when this is reset)
  bool sawAmbigInverse = false;
  // what is the most recent HMM index that has an inverted state?
  int mostRecentInvIdx = -1;

  // Number of recombinations in <curState> relative to <prevState> below
  uint16_t numRecombs; // TODO: optimization: do we want this?
  int lastAssignedMarker = chrLastMarker + 1; // technically not assigned yet
  uint64_t lastAssignedIV = 0, lastIVFlip = 0;
  uint64_t lastIVparDiff = 0;
  bool lastIVSet = false;
  // the most recent informative index per parent; assigned while accounting for
  // ambiguous parent heterozygosity
  int lastInformIndex[2];
  // which _hmm index has each parent fully assigned (actually should say not
  // fully missing: check is whether _all_ the IVs for a parent are unassigned
  // or not).
  // The name of this variable breaks with the convention above where last
  // means previous in the loop iteration. Here we mean first marker on the
  // chromosome
  int firstFullyAssignedIndex[2] = { -1, -1 };
  int nextHmmIndex;
  int prevHmmIndex = -1;
  for(int hmmIndex = lastIndex; hmmIndex >= 0; prevHmmIndex = hmmIndex,
					       hmmIndex = nextHmmIndex) {
    if (CmdLineOpts::verbose) {
      fprintf(log, "    Marker %d, hmmIndex %d:", _hmmMarker[ hmmIndex ],
	      hmmIndex);
      fflush(log);
    }
    nextHmmIndex = hmmIndex - 1; // modified when there are error states
    State *curState = _hmm[hmmIndex][curStateIdx];

    uint8_t curAmbigParPhase = curState->ambigParPhase;
    uint8_t curAmbigParHet = curState->ambigParHet;
    uint8_t curArbitraryPar = curState->arbitraryPar;
    // Indicates which inheritance vector values differ among any ambiguous
    // states that are in <_curIdxSet>. The true IV value is ambiguous in that
    // case.
    uint64_t ivFlippable = 0;

    uint8_t curHomParGeno = curState->homParentGeno;
    // above only valid for one parent het states or states with ambiguous
    // parent heterozygosity
    bool homParGenoAssigned = curState->hetParent < 2 ||
						      (curAmbigParHet & 3) != 0;
    // deal also with X chromosome corner case (see checkMinUpdate() and the
    // if (_onChrX) statement below):
    homParGenoAssigned = homParGenoAssigned &&
					  (!_onChrX || curHomParGeno != G_MISS);
    uint8_t ambigHomParGeno = 0;

    if (_missingPar == 3) {
      _prevCanonIVsBT.clear();
      int dontcare;
      uint64_t canonIV = getCanonicalIV(curState->iv, curState->ambig,
					curState->unassigned, dontcare);
      _prevCanonIVsBT.insert(canonIV);
    }

    // Find all the possible parent phase types that have equal and minimal
    // numbers of recombinations at this marker and append their previous states
    // to <_prevIdxSet>.
    int numAppliedStates = 0;
    for(state_set_iter it = _curBTAmbigSet->begin(); it !=_curBTAmbigSet->end();
									it++) {
      if (CmdLineOpts::verbose)
	fprintf(log, " %d", it->stateIdx);
      State *ambigState = _hmm[hmmIndex][ it->stateIdx ];
      // decide if it->iv has opposite parent assignments relative to any we
      // are applying
      bool isInverse = false;
      if (_missingPar == 3) {
	// invert which parent transmitted the IVs:
	uint64_t invertedIV = ((it->iv << 1) & _parBits[1]) |
			      ((it->iv >> 1) & _parBits[0]);
	int dontcare;
	uint64_t canonIV = getCanonicalIV(invertedIV, it->ambig,
					  /*unassigned=*/ 0, dontcare);
	auto it2 = _prevCanonIVsBT.find(canonIV);
	if (it2 != _prevCanonIVsBT.end())
	  isInverse = true;

	if (!isInverse) {
	  // we _will_ back traced to this state: store its canonical IV
	  canonIV = getCanonicalIV(it->iv, it->ambig, /*unassigned=*/ 0,
				   dontcare);
	  _prevCanonIVsBT.insert(canonIV);
	}
      }
      sawAmbigInverse = sawAmbigInverse || isInverse;
      if (isInverse)
	  mostRecentInvIdx = hmmIndex;
      if (!isInverse) {
	uint64_t ivDiffBits = curState->iv ^ it->iv;

	numAppliedStates++;

	for (uint8_t p = _firstMissP; p < _limitMissP; p++) {
	  if (firstFullyAssignedIndex[p] < 0 &&
					  ambigState->unassigned == _parBits[p])
	    firstFullyAssignedIndex[p] = prevHmmIndex;
	}

	bool thisStateParFlip = false;
	if (parArbitraryRun) {
	  // Only remains arbitrary if the following conditions hold (see longer
	  // comment above when parArbtitrary run is first set true)
	  // (1) All differing bits different for both parents?
	  // (2) All non-diff (same) bits the opposite of each other in the IV?
	  // (3) The ambig IV values (curState->ambig and it->ambig) match
	  uint64_t ivExpectOppositeBits = ~(ivDiffBits | it->ambig);
	  uint64_t ivValsExpOpp = it->iv & ivExpectOppositeBits;
	  thisStateParFlip = curState->ambig == it->ambig &&
	      ((ivDiffBits >> 1) & _parBits[0]) == (ivDiffBits & _parBits[0]) &&
	      (((ivValsExpOpp >> 1) ^ ivValsExpOpp) & _parBits[0]) ==
					  (_parBits[0] & ivExpectOppositeBits);
	  if (!thisStateParFlip &&
	      (ivDiffBits > 0 || curState->ambig != it->ambig)) {
	    parArbitraryRun = false;
	    curArbitraryPar = 1;
	  }
	}

	if (!thisStateParFlip) {
	  ivFlippable |= ivDiffBits;

	  curAmbigParHet |= ambigState->ambigParHet;
	  curAmbigParPhase |= it->ambigParPhase;
	  curArbitraryPar |= ambigState->arbitraryPar;
	}

	if (_onChrX) {
	  // On chrX, if dad is missing data, he is imputed using daughters; if,
	  // for example, there's only one daughter and she is heterozygous and
	  // _may_ have recombined relative to the previous marker, there's no
	  // information to tell us which allele is from Mom -- she could have
	  // transmitted either. Note that, this case -- which came up in
	  // simulated data -- arises when there is a recombination that becomes
	  // certain to have occurred once we encounter a marker where the
	  // recombined daughter is homozygous, but it may or may not have
	  // occurred at any markers between sites where she is homozygous.
	  if (homParGenoAssigned && ambigState->homParentGeno != curHomParGeno) {
	    curHomParGeno = G_MISS;
	    homParGenoAssigned = false;
	  }
	}
	else if (ambigState->hetParent < 2 || (ambigState->ambigParHet & 3) !=0) {
	  assert(ambigState->homParentGeno != G_MISS);
	  if (!homParGenoAssigned) {
	    curHomParGeno = ambigState->homParentGeno;
	    homParGenoAssigned = true;
	  }
	  else
	    ambigHomParGeno |= curHomParGeno ^ ambigState->homParentGeno;
	}
      }

      // Add the previous state index(es) to <_prevIdxSet> so long as the
      // ambiguous state has the same previous marker as the current one
      if (ambigState->prevHMMIndex != curState->prevHMMIndex)
	// TODO: want to find a way to allow this?
	continue;

      int prevHmmIndex = hmmIndex - curState->prevHMMIndex;
      if (prevHmmIndex >= 0) {
	BT_ambig_info dontcare1;
	uint8_t dontcare2;
	collectAmbigPrevIdxs(it->iv, it->ambig, ambigState->ambigPrev,
			     ambigState->prevState, _hmm[prevHmmIndex],
			     dontcare1, dontcare2);
      }
    }
    if (sawAmbigInverse) {
      if (numAppliedStates > 1)
	// set P flag: one of the ambiguous paths leads to inverted parents
	// at a later marker. Want to ensure that the arbitrarily chosen
	// path is indicated as arbitrary for the parent
	curArbitraryPar = 1;
      else if (mostRecentInvIdx > hmmIndex) {
	// stop run of P-ambig markers
	sawAmbigInverse = false;
	// ensure that the next marker has arbitrary parent choice (in some
	// cases it won't get marked above, but it definitely does have P
	// status)
	theFam->setArbitraryPar(_hmmMarker[mostRecentInvIdx]);
	mostRecentInvIdx = -1;
      }
      else if (mostRecentInvIdx < 0) {
	assert(false);
      }
    }


    // Decide whether to assign one hap trans for each missing data parent
    // TODO: can this assign one hap trans for dads on chrX?
    for (uint8_t p = _firstMissP; p < _limitMissP; p++) {
      if (!assignOneHapTrans[p]) {
	// search for potential OHT starting point

	if (hmmIndex == lastIndex || // always assign one hap trans at end ...
	    // ... and beginning of a chromosome (the beginning of the
	    // chromosome is completed after the backtracing loop)
	    firstFullyAssignedIndex[p] >= 0 ||
	    curState->prevHMMIndex == hmmIndex + 1) {
	  // detect which haplotype was transmitted:
	  uint64_t parHap = _hmm[hmmIndex][curStateIdx]->iv & _parBits[p];
	  uint64_t oppParHap = parHap ^ _parBits[p];

	  if ((parHap & (parHap - 1)) == 0) { // test for 1 bit set
	    assignOneHapTrans[p] = true;
	    oneHapTransOutlierIV[p] = parHap;
	    oneHapTransHap[p] = 0; // transmitted haplotype is 0
	  }
	  else if ((oppParHap & (oppParHap - 1)) == 0) { // test all but 1 set
	    assignOneHapTrans[p] = true;
	    oneHapTransOutlierIV[p] = oppParHap;
	    oneHapTransHap[p] = 1; // transmitted haplotype is 1
	  }
	}
      }

      // handle OHT for markers at the end of the chromosome (we handle the
      // beginning after the backtracing loop)
      if (assignOneHapTrans[p] && curState->prevHMMIndex != hmmIndex + 1 &&
	  firstFullyAssignedIndex[p] < 0) {
	uint8_t hetParent = _hmm[hmmIndex][curStateIdx]->hetParent;

	// has it ended?
	if (hetParent != 1 - p)
	  assignOneHapTrans[p] = false;
	else {
	  // can do this by (1) setting <ivFlippable> for the child that has
	  // the outlier haplotype ...
	  ivFlippable |= oneHapTransOutlierIV[p];

	  // ... (2) indicating the site may be heterozygous for both parents...
	  curAmbigParHet |= 1 << 2;

	  // ... and (3) assigning the phase type for the both parent het case:
	  // the phase type for the both parent het determines which haplotype
	  // is set to be printed (other is set missing); we want this to be
	  // <transHap>.
	  // The homozygous genotype is either 0 or 3, and we'll divide by 3 to
	  // get a boolean:
	  // Note: on the X chromosome, <homParentGeno> can be missing, which
	  // leads to phaseType == 0. This is OK: that parent (the dad) will
	  // have both alleles missing and his phase won't matter.
	  uint8_t homGenoKind = _hmm[hmmIndex][curStateIdx]->homParentGeno / 3;
	  // If the phase type is 0, when the genotype is heterozygous,
	  // haplotype 0 is allele 0, and haplotype 1 is allele 1. Whichever
	  // haplotype is the same as homGenoKind will be printed
	  // (non-missing).  So, if <homGenoKind> == 0 and <transHap> == 0, we
	  // want the phase type to be the same as <transHap>. Same if
	  // <homGenoKind> == 0 and <transHap> == 1.
	  // The reverse is true when <homGenoKind> == 1. So it suffices to
	  // xor:
	  // Note: must shift the phase types for <homParent> and <hetParent>
	  // to their respective positions in the 2-bit <phaseType>.
	  int8_t homParent = p;
	  int8_t phaseType = ((oneHapTransHap[p] ^ homGenoKind) << homParent) |
		      (_hmm[hmmIndex][curStateIdx]->parentPhase << hetParent);
	  // (Note: shifting by 4 bits to get to the initial bit that stores
	  // the both parent het phase types)
	  curAmbigParPhase |= 1 << (4 + phaseType);
	}
      }
    }

    // Remove any ambig par phase types that have equivalent phase to <curState>
    curAmbigParPhase -= (1 << curState->parentPhase) << 2 * curState->hetParent;
    // Remove any ambig par het types that are equivalent to <curState>
    curAmbigParHet -= 1 << curState->hetParent;

    // Note: curState->error == 2 means the path has an error in it. We only
    // need deal with curState->error == 1 (the previous marker is erroneous):
    if (curState->error & 1) {
      assert(curState->prevHMMIndex > 1);

      // Set error immediately for the skipped (error) markers
      for(uint8_t relIdx = 1; relIdx < curState->prevHMMIndex; relIdx++) {
	int theHMMIndex = hmmIndex - relIdx;
	if (theFam->getStatus(_hmmMarker[theHMMIndex]) != NUM_PHASE_STATUS) {
	  // this marker has already had its status set: it must be a forced
	  // informative marker and isn't truly an error. Leave the original
	  // status as is
	  deleteStates(_hmm[theHMMIndex]);
	  continue;
	}
	uint64_t childrenData = _genos[theHMMIndex].second;
	uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
	theFam->setStatus(/*marker=*/ _hmmMarker[theHMMIndex],
			  PHASE_ERR_RECOMB, _genos[theHMMIndex].first,
			  childrenData, missing);
	// <curState->prevState> references a state two indexes back
	deleteStates(_hmm[theHMMIndex]);
      }
      // skip the number the error markers above
      nextHmmIndex = hmmIndex - curState->prevHMMIndex;
    }
    else if (curState->prevHMMIndex > 1) {
      // skip some number of HMM indexes, freeing the memory associated with
      // them:
      // (these arose from forced informative markers; they're not necessarily
      // from this state path)
      for(uint8_t relIdx = 1; relIdx < curState->prevHMMIndex; relIdx++) {
	int theHMMIndex = hmmIndex - relIdx;
	deleteStates(_hmm[theHMMIndex]);
      }
      nextHmmIndex = hmmIndex - curState->prevHMMIndex;
    }

    // In the previous state, resolve ambiguous <iv> values and propagate
    // backward any <iv> values that were unassigned in that state
    if (hmmIndex - curState->prevHMMIndex >= 0) {
      BT_ambig_info prevStateInfo;
      uint8_t numAmbig1Recombs;
      int theHMMIndex = hmmIndex - curState->prevHMMIndex;
      collectAmbigPrevIdxs(curState->iv, curState->ambig, curState->ambigPrev,
			   curState->prevState, _hmm[theHMMIndex],
			   prevStateInfo, numAmbig1Recombs);

      prevStateIdx = prevStateInfo.stateIdx;
      State *prevState = _hmm[theHMMIndex][prevStateIdx];
      prevState->iv = prevStateInfo.iv;
      prevState->ambig = prevStateInfo.ambig;
      // TODO: optimization: remove entry for prevStateInfo from _prevBTAmbig

      // Both the following are equivalent, but now we do something different
      // as noted next
//      numRecombs = curState->minRecomb - prevState->minRecomb -
//			  (curState->error & 1) * CmdLineOpts::maxNoErrRecombs;
//      numRecombs = curState->maxPrevRecomb;
      // Changes in IV when a child is ambig1 induce a recombination that is
      // counted at a later marker than it originally occurred. The
      // <numAmbig1Recombs> counts the number of such recombinations that must
      // be attributed to an earlier marker to be consistent with the changes
      // in IV values:
      // see comment above State::minRecomb (phaser.h) for why it's *10:
      prevState->maxPrevRecomb += numAmbig1Recombs * 10;
      numRecombs = curState->maxPrevRecomb - numAmbig1Recombs * 10;
    }
    else
      numRecombs = 0;

    // Missing genotype value is 01 (not any other); 5 = 0101; 10 = 1010; so:
    uint8_t parentData = _genos[hmmIndex].first;
    uint8_t parMissing = (parentData & 5) & ~((parentData & 10) >> 1);
    uint64_t childrenData = _genos[hmmIndex].second;
    uint64_t missing = (childrenData & _parBits[0]) &
					  ~((childrenData & _parBits[1]) >> 1);
    // <curHomParGeno> is non-missing iff <homParGenoAssigned>
    assert((curHomParGeno != G_MISS) == homParGenoAssigned);
    // either we have a homozygous parent genotype OR the site is PI with no
    // ambiguity in parent heterozygosity OR this is chrX (with no data for
    // daughters)
    assert(curHomParGeno != G_MISS ||
	   (curState->hetParent == 2 && (curAmbigParHet & 3) == 0) || _onChrX);
    assert(!ambigHomParGeno);
    theFam->setPhase(_hmmMarker[hmmIndex], curState->iv,
		     curState->ambig & _parBits[1], missing, ivFlippable,
		     parMissing, curState->hetParent, curHomParGeno,
		     // comment above State::minRecomb (phaser.h) for why /10:
		     curState->parentPhase, numRecombs / 10, curAmbigParHet,
		     curAmbigParPhase, curArbitraryPar);

    if (curState->hetParent == 2 || curAmbigParHet)
      // informative for _both_ parents
      // _any_ ambiguity in the parent heterozygosity means that the marker may
      // be informative for each parent
      lastInformIndex[0] = lastInformIndex[1] = hmmIndex;
    else
      lastInformIndex[ curState->hetParent ] = hmmIndex;

    //////////////////////////////////////////////////////////////////////////
    // For all markers between the next informative one and the current one,
    // assign (1) which parent haplotypes were _un_transmitted, and
    // (2) for any ambiguous sites, indicate whether they are necessarily
    // homozygous:

    // make the transmitted haplotype assignments for all markers between the
    // current marker and the subsequent informative marker
    int curMarker = _hmmMarker[hmmIndex];
    uint64_t curIVparDiff =
      calcAndSetUntransPar(theFam, /*startMarker=*/ curMarker,
			   lastAssignedMarker, lastIVFlip, ivFlippable,
			   lastIVSet, curState, lastAssignedIV,
			   lastIVparDiff);

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
    uint8_t oneParHomozy = 1 - (curState->hetParent >> 1);
    uint8_t homParent = oneParHomozy * (1 - curState->hetParent);
    uint64_t propagateMask = (missing * 3) |
					  (oneParHomozy * _parBits[homParent]);
    if (lastIVSet)
      lastAssignedIV = (curState->iv & ~propagateMask) |
					       (lastAssignedIV & propagateMask);
    else
      lastAssignedIV = curState->iv;
    lastIVSet = true;
    lastIVFlip = ivFlippable;
    lastIVparDiff = curIVparDiff;
    lastAssignedMarker = _hmmMarker[hmmIndex];

    // About to finish back tracing; assign the untransmitted parent for the
    // first few markers (before the first informative one)
    if (hmmIndex - curState->prevHMMIndex < 0)
      calcAndSetUntransPar(theFam, /*startMarker=*/ chrFirstMarker,
			   lastAssignedMarker, lastIVFlip,
			   /*(first marker) ivFlippable=*/ 0, lastIVSet,
			   /*first marker state=first assigned state=*/curState,
			   lastAssignedIV, lastIVparDiff);

    deleteStates(_hmm[hmmIndex]);

    // ready for next iteration -- make prev into cur for states and state sets
    curStateIdx = prevStateIdx;
    BT_state_set *tmp = _curBTAmbigSet;
    _curBTAmbigSet = _prevBTAmbigSet;
    _prevBTAmbigSet = tmp;
    _prevBTAmbigSet->clear();

    if (CmdLineOpts::verbose)
      fprintf(log, "\n");
  }

  for (uint8_t p = _firstMissP; p < _limitMissP; p++) {
    if (assignOneHapTrans[p]) {
      int8_t homParent = p;
      int8_t hetParent = 1 - homParent;
      int8_t untransHap = 1 - oneHapTransHap[p];
      // have 4 untransParBits: first two for parent 0, second two for parent 1
      // we shift to the appropriate parent with << (2 * homParent);
      // either bit 0 or bit 1 of the given parent is untransmitted, and
      // haplotype 0 is bit 0, and hap 1 is bit 1, so (1 << untransHap) gives
      // this
      uint8_t untransParBits = (1 << untransHap) << (2 * homParent);
      int stopMarker = _hmmMarker[ lastInformIndex[ homParent ] ]; 
      if (firstFullyAssignedIndex[p] >= 0) {
	int alternativeStopMarker = _hmmMarker[ firstFullyAssignedIndex[p] ];
	assert(alternativeStopMarker >= stopMarker);
	stopMarker = alternativeStopMarker;
      }
      for(int m = chrFirstMarker; m < stopMarker; m++) {
	PhaseStatus curStat = theFam->getStatus(m);
	if (curStat == PHASE_UNINFORM)
	  theFam->setUntransPar(m, untransParBits);
	else if (curStat == PHASE_OK) {
	  const PhaseVals &phase = theFam->getPhase(m);
	  if (phase.hetParent == 2)
	    continue;
	  //assert(phase.hetParent == hetParent);
	  assert(phase.homParentGeno == G_HOM0 || phase.homParentGeno ==G_HOM1);
	  // See the earlier code that handles <assignOneHapTrans> during
	  // backtracing
	  uint8_t homGenoKind = phase.homParentGeno / 3;
	  uint8_t phaseType = ((oneHapTransHap[p] ^ homGenoKind) << homParent) |
				  (phase.parentPhase << hetParent);
	  theFam->updateOHTVals(m, /*ivFlippable |=*/ oneHapTransOutlierIV[p],
				/*ambigParHet |=*/ 1 << 2,
				/*ambigParPhase |=*/ 1 << (4 + phaseType));
	}
	// TODO: PHASE_X_SPECIAL?
	// TODO: for PHASE_AMBIG recalculate <bothParHomozyAtAmbig>?
      }
    }
  }
}

// Given a list of states, does a linear search to find the states with minimum
// numbers of recombinations. Returns the index of the first state encountered
// that is minimal and stores any other states that are tied for minimum in
// <_curIdxSet>.
uint32_t Phaser::findMinStates(dynarray<State*> &theStates) {
  uint32_t minLastMarker = UINT32_MAX;
  // The state we will designate as the minimum even as there may ambiguity
  // stored in <_stateIdxSet>
  uint32_t minStateIdx = UINT32_MAX;

  int numStates = theStates.length();
  for(int i = 0; i < numStates; i++) {
    uint32_t curRecomb = theStates[i]->minRecomb;
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
  int lowOrderChildBit = lowOrderUnambigAssignedBit(prevState->ambig,
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

// In regions where a parent has untransmitted chromosomes, this code determines
// which of the two is untransmitted and later code sets those alleles as
// missing. The determination is based on the IV of flanking informative markers
uint64_t Phaser::calcAndSetUntransPar(NuclearFamily *theFam, int startMarker,
				      int lastAssignedMarker,
				      uint64_t lastIVFlip, uint64_t ivFlippable,
				      bool lastIVSet, State *curState,
				      uint64_t lastAssignedIV,
				      uint64_t lastIVparDiff) {
  // Determine which parent haplotypes were _un_transmitted. The first two
  // bits are parent 0's haplotype 0 and 1, and the second two bits are
  // parent 1's haplotype 0 and 1. If the corresponding bit is set to 1, the
  // haplotype was _not_ transmitted
  uint8_t untrans;

  // TODO: need the flanking informative markers for each parent do this right
  // We use the IV values at the flanking informative markers to determine
  // this. <ivFlippable> values are uncertain and so those values will not
  // figure into which haplotypes the parents transmitted.
  // start with <ivFlippable> at the current and subsequent markers; we'll
  // modify this further below
  uint64_t uncertainIV = ivFlippable;

  uint64_t curIVparDiff = (curState->iv & _parBits[0]) ^
					    ((curState->iv & _parBits[1]) >> 1);

  if (curState->hetParent < 2)
    // Any IV values that differ between <curState->iv> and <lastAssignedIV>
    // are uncertain at the current marker (up through <lastAssignedMarker>) if
    // the parent transmitting those IV values is not heterozygous
    uncertainIV |= (curState->iv ^ lastAssignedIV) &
					      _parBits[1 - curState->hetParent];

  // for checking the special X marker/phase type case below
  uint64_t certainChildBits1[2];
  for(int sex = 0; sex < 2; sex++)
    certainChildBits1[sex] = _childSexes[sex] & _parBits[1] & ~uncertainIV;

  for(int m = startMarker; m < lastAssignedMarker; m++) {
    untrans = 0;

    if (m == startMarker + 1 && lastIVSet) {
      // although <startMarker> is only uncertain at <ivFlippable> values (it
      // is informative and its <ivFlippable> value applies exactly at that
      // site), subsequent markers are also uncertain at flippable values at
      // the physically subsequent marker -- since we're backtracing that's
      // the last visited site
      uncertainIV |= lastIVFlip;
      // if a child received a recombined haplotype, the location of the
      // recombination is uncertain, and which of the parent's haplotypes s/he
      // transmitted cannot be determined. As such, we won't include that
      // child's haplotype when determining whether a given haplotype was
      // transmitted.
      // We only apply this for startMarker + 1 and thereafter. The reason is
      // <startMarker> is informative and we *can* say which parent's haplotype
      // was transmitted there. That is, the recombination happened _after_
      // <startMarker>, not _at_ it. But just above the loop, we make the IVs
      // uncertain for the homozygous parent, if there is one.
      // (Corner case is when <startMarker> is <chrFirstMarker>, but there
      // curState->iv == lastAssignedIV and this next statement has no effect.)
      uncertainIV |= curState->iv ^ lastAssignedIV;
      // update which bits are certain to match the change in <uncertainIV>
      for(int sex = 0; sex < 2; sex++)
	certainChildBits1[sex] &= ~uncertainIV;
    }

    const PhaseVals &phase = theFam->getPhase(m);

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
    // is not the case are those where both parents transmitted the same
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
    bothParHomozyAtAmbig &= !_onChrX; // above case is only for autosomal SNPs

    if (_onChrX) {
      // on the X chromosome, for dad, there's only one haplotype, coded as 0
      // that's transmitted to the daughters; haplotype 1 should be set to be
      // the same as haplotype 0 -- i.e., dad is either fully missing or has
      // "both" alleles present (really only one)
      // we first clear bit 1 (13 == 1101 binary) and then OR with the value
      // in bit 0, shifted into bit 1
      untrans = (untrans & 13) | ( (untrans & 1) << 1);

      if (theFam->getStatus(m) == PHASE_X_SPECIAL) {
	// Check the condition for the special X chromosome marker/phase type:
	// these are probably just uninformative markers, but may in fact be
	// almost fully ambiguous

	uint64_t shiftedMiss = thisMarkMiss << 1; // shifted to Mom's IV bits
	bool sonsAllSameIV = false;
	uint8_t sonsMissing = 0;
	// when <sonsAllSameIV>, which haplotype (0/1) did the sons all inherit?
	uint8_t sonsIV = 0;
	// We examine only Mom's IV values for those children (sons just below)
	// where we are certain of the IV. First we check whether they're all
	// missing and then we determine whether they all have the same IV value
	if ( (certainChildBits1[0] & ~shiftedMiss) == 0 ) {
	  // all sons with definite IV values are missing, so the only
	  // question is whether the daughters all inherited the same allele
	  sonsMissing = 1;
	}
	else if ((curState->iv & certainChildBits1[0] & ~shiftedMiss) ==
					(certainChildBits1[0] & ~shiftedMiss)) {
	  sonsAllSameIV = true;
	  sonsIV = 1;
	}
	else if ((~curState->iv & certainChildBits1[0] & ~shiftedMiss) ==
					(certainChildBits1[0] & ~shiftedMiss)) {
	  sonsAllSameIV = true;
	  sonsIV = 0;
	}

	if (sonsAllSameIV || sonsMissing) { // first condition for ambiguity met
	  // if this site is ambiguous, we use the <untrans> variable to
	  // capture which parent we don't have information to correctly
	  // impute. In this form of ambiguity, dad is fully _un_imputable, so
	  // both his bits are set (3):
	  uint8_t setUntransIfAmbig = 3;
	  // ... and Mom can only be imputed for the haplotype transmitted to
	  // the certain sons:
	  // (Note: this is 3 if <sonsMissing> is true: in this case, Mom's
	  // haplotypes are unclear so she is set as fully missing)
	  uint8_t momSons = 3 - (1 - sonsMissing) * (1 << sonsIV);
	  setUntransIfAmbig |= momSons << 2;

	  // in order to get into the special case, the sons must all have the
	  // same genotype (or be missing) and therefore the same IV from Mom;
	  // now we need only check the daughters. First we check whether
	  // they're all missing and then we determine whether they (a) all have
	  // the same IV value and (b) have a different IV value than the sons.
	  // When that is the case, the marker is the special case ambiguous.
	  if ( (certainChildBits1[1] & ~shiftedMiss) == 0 )
	    // all daughters with definite IV values are missing, so we have
	    // an ambiguity
	    untrans |= setUntransIfAmbig;
	  else if ((curState->iv & certainChildBits1[1] & ~shiftedMiss) ==
					(certainChildBits1[1] & ~shiftedMiss)) {
	    uint8_t daughterIV = 1;
	    if (sonsMissing || sonsIV != daughterIV)
	      untrans |= setUntransIfAmbig;
	  }
	  else if ((~curState->iv & certainChildBits1[1] & ~shiftedMiss) ==
					(certainChildBits1[1] & ~shiftedMiss)) {
	    uint8_t daughterIV = 0;
	    if (sonsMissing || sonsIV != daughterIV)
	      untrans |= setUntransIfAmbig;
	  }
	}
      }
    }

    theFam->setUntransPar(m, untrans, bothParHomozyAtAmbig);
  }

  return curIVparDiff;
}

// For debugging:
void Phaser::tracePrior(int hmmIndex, int stateIdx) {
  if (hmmIndex > _hmm.length())
    return;

  do {
    State *state = _hmm[ hmmIndex ][ stateIdx ];
    printf("iv = %ld, ambig = %ld, unassigned = %ld, minRecomb = %.1f, numSinceOneHetPar = %d\n",
	   // comment above State::minRecomb (phaser.h) for why /10.0f:
	   state->iv, state->ambig, state->unassigned, state->minRecomb / 10.0f,
	   state->numMarkersSinceOneHetPar);
    printf("  hmmIndex = %d, state = %d  (marker = %d)\n",
	   hmmIndex, stateIdx, _hmmMarker[hmmIndex]);

    hmmIndex -= state->prevHMMIndex;
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

// For debugging:
void Phaser::printStates() {
  printStates(_hmm.length() - 1);
}

// For debugging:
void Phaser::printStates(int hmmIndex) {
  if (hmmIndex >= _hmm.length() || hmmIndex < 0)
    return;

  uint32_t minRecomb = UINT32_MAX;
  std::vector<int> minimal;

  int numStates = _hmm[hmmIndex].length();
  for(int idx = 0; idx < numStates; idx++) {
    State *state = _hmm[ hmmIndex ][ idx ];

    printf("iv = %ld, ambig = %ld, unassigned = %ld, minRecomb = %.1f, numSinceOneHetPar = %d\n",
	   // comment above State::minRecomb (phaser.h) for why /10.0f:
	   state->iv, state->ambig, state->unassigned, state->minRecomb / 10.0f,
	   state->numMarkersSinceOneHetPar);
    printf("  hmmIndex = %d, state = %d  (marker = %d)\n",
	   hmmIndex, idx, _hmmMarker[hmmIndex]);

    if (state->minRecomb < minRecomb) {
      minRecomb = state->minRecomb;
      minimal.clear();
      minimal.push_back(idx);
    }
    else if (state->minRecomb == minRecomb) {
      minimal.push_back(idx);
    }
  }

  printf("\n");
  // comment above State::minRecomb (phaser.h) for why /10.0f:
  printf("Overall minRecomb = %.1f, states:", minRecomb / 10.0f);
  for(auto it = minimal.begin(); it != minimal.end(); it++) {
    printf(" %d", *it);
  }
  printf("\n");
}
