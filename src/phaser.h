// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/nuclearfamily.h>
#include <genetio/personbulk.h>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <deque>

#ifndef PHASER_H
#define PHASER_H

// Phase methods: minimum recombinant or maximum likelihood
enum PhaseMethod {
  PHASE_MAXLIKE,
  PHASE_MINREC,
};

// Types of markers. We use bits set in an integer for markers that can
// have multiple types (occurs when one or both parents are missing data)
enum MarkerType {
  MT_FI_1,   // Fully informative for one parent
  MT_PI,     // Partly informative (both parents heterozygous for same alleles)
  MT_UN,     // Uninformative (both parents homozygous)
  MT_AMBIG,  // All children and observed parents are heterozygous or missing
  MT_ERROR,  // Mendelian error marker
  MT_N_TYPES
};
// Note: fully informative for both parents is impossible in biallelic genotype
// data. As we're using PLINK BED format data at present, we are using
// biallelic data.
// One can encode multiallelic variants using a series of biallelic markers. If
// this turns out to be a common case, it would be nice to let the user specify
// which of the biallelic markers actually correspond to a single multiallelic
// variant. In some settings you end up with a partly informative marker
// even though the series of markers combined is one fully informative (not just
// for one parent) marker. It may not make sense to go to the effort to code
// for this case since the number of such markers is often low, but it wouldn't
// necessarily be hard to do. At most four variants (for non-error markers)
// will be observed in the family data, and one can specifically handle this
// uncommon case in such a way that the standard PLINK representation isn't
// relied upon.

struct State;

// During backtracing, when multiple states lead to the same number of
// recombinations, we'd like to identify which inheritance vector values are
// uncertain. We do this by tracking the inheritance vector <iv> and ambiguous
// bits <ambig> that apply to each state at the previous or current maker
// (stored in <_curBTAmbigSet> and <_prevBTAmbigSet>) -- thus <stateIdx> is
// an index into the list of states in <_hmm> for the corresponding marker,
// whether that referred to by _cur* or _prev*.
// Storing these values allows us to diff (using bitwise xor) the <iv> values
// among all the states. Note that an advantage of storing these values separate
// from the states themselves is that multiple paths that lead to the same state
// but yield different propagated values are possible to store with this.
// An alternative, heavier-weight solution would be to make separate full States
// during back tracing, but the only differing values are those stored here, so
// this suffices. It's also immediately apparent when two paths coalesce on the
// same values since we use sets to store the objects.
struct BT_ambig_info {
  BT_ambig_info() { }
  BT_ambig_info(uint32_t s, uint64_t i, uint64_t a, uint8_t pp,
		uint8_t app) {
    stateIdx = s;
    iv = i;
    ambig = a;
    parPhase = pp;
    ambigParPhase = app;
  }
  uint32_t stateIdx;
  uint64_t iv;
  uint64_t ambig;
  uint8_t parPhase;
  uint8_t ambigParPhase;
};


class Phaser {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void init() {
      // TODO: optimization: set load factors for these next values?
      _stateHash.set_empty_key(NULL);
      _curBTAmbigSet = new BT_state_set();
      _prevBTAmbigSet = new BT_state_set();
      // TODO: would using pointers make this faster?
      // Note that this value is impossible as <State::ambig> will always have
      // pairs of bits set (or, for ambig1, the first order bit set in any pair)
      // instead of just the second bit as here:
      BT_ambig_info empty(UINT32_MAX, UINT64_MAX, 2, 1, 1);
      _curBTAmbigSet->set_empty_key(empty);
      _prevBTAmbigSet->set_empty_key(empty);

      int maxMarkers = 0;
      for(int c = 0; c < Marker::getNumChroms(); c++)
	if (Marker::getNumChromMarkers(c) > maxMarkers)
	  maxMarkers = Marker::getNumChromMarkers(c);
      _hmm.resize(maxMarkers);
      _hmmMarker.resize(maxMarkers);
      _genos.resize(maxMarkers);
      _ambigPrevLists.resize(maxMarkers);
    }
    static void run(NuclearFamily *theFam, int chrIdx, FILE *log);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void parBitsInit(int numChildren);
    static void getFamilyData(NuclearFamily *theFam, uint8_t missingPar,
			      int marker, uint8_t &parentData,
			      uint8_t &parentGenoTypes,
			      uint64_t childrenData[5], uint8_t &childGenoTypes,
			      int &numMissChildren);
    static int  getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes,
			      uint8_t &homParGeno);
    static void printMarkerType(int mt, FILE *log);
    static void makePartialStates(dynarray<State> &partialStates,
				  int markerTypes, uint8_t parentData,
				  uint8_t homParGeno,
				  const uint64_t childrenData[5]);
    static void makePartialFI1States(dynarray<State> &partialStates,
				     uint8_t parentData, uint8_t homParGeno,
				     const uint64_t childrenData[5]);
    static void makePartialPIStates(dynarray<State> &partialStates,
				    uint8_t parentData,
				    const uint64_t childrenData[5]);
    static void makeFullStates(const dynarray<State> &partialStates, int marker,
			       int firstMarker, const uint64_t childrenData[5],
			       uint8_t missingPar, int numChildren,
			       int numMissChildren);
    static uint8_t isIVambigPar(const State *state, uint8_t missingPar);
    static void mapPrevToFull(const State *prevState, uint8_t prevHMMIndex,
			      uint32_t prevIdx, const State &curPartial,
			      float minMaxRec[2],
			      const uint64_t childrenData[5],
			      uint8_t IVambigPar, uint8_t missingPar,
			      int numDataChildren, int numMarkersSincePrev,
			      bool &zeroRecombsThisPrev);
    static void checkPenalty(const State *prevState, const State &curPartial,
			     uint8_t isPI, uint8_t missingPar,
			     int numMarkersSincePrev, uint64_t &allButOneIV,
			     uint8_t &applyPenalty,
			     int16_t &numMarkersSinceNonHetPar,
			     int16_t &numMarkersSinceOneHetPar,
			     const uint64_t childrenData[5]);
    static void handlePI(const State *prevState, uint64_t &fullIV,
			 uint64_t &fullAmbig, uint64_t &recombs,
			 uint64_t parRecombs[2], uint64_t &propagateAmbig,
			 uint64_t &defaultPhaseHasRecomb,
			 uint64_t childPrevUnassigned[2],
			 uint64_t unambigHetRecombs[4],
			 const uint64_t childrenData[5]);
    static void calcHetChildPIRecombs(const uint64_t parRecombs[2],
				      const uint64_t unambigHets,
				      uint64_t unambigHetRecombs[4]);
    static void flipPIVals(uint64_t &fullIV, uint64_t &fullAmbig,
			   const uint64_t childrenData[5],
			   uint64_t propagateAmbig,
			   const uint64_t unambigHetRecombs[4],
			   const uint64_t childPrevUnassigned[2],
			   uint64_t defaultPhaseHasRecomb);
    static void fixRecombFromAmbig(uint64_t &fullIV, uint64_t &recombs,
				   const uint64_t parRecombs[2], uint8_t isPI,
				   uint64_t ambigOnlyPrev, uint8_t hetParent,
				   uint64_t &stdAmbigOnlyPrev,
				   uint64_t &ambig1PrevInfo,
				   uint64_t &ambig1Unassigned);
    static void updateStates(uint64_t fullIV, uint64_t fullAmbig,
			     uint64_t fullUnassigned, uint64_t ambig1Unassigned,
			     uint64_t recombs, uint64_t allButOneIV,
			     uint8_t applyPenalty,
			     const State *prevState, uint64_t stdAmbigOnlyPrev,
			     uint64_t ambig1PrevInfo, uint8_t hetParent,
			     uint8_t homParentGeno, uint8_t initParPhase,
			     uint8_t altPhaseType, uint8_t prevHMMIndex,
			     uint32_t prevIndex, uint8_t IVambigPar,
			     float minMaxRec[2], bool hetParentUndefined,
			     const uint64_t childrenData[5],int numDataChildren,
			     int16_t numMarkersSinceNonHetPar,
			     int16_t numMarkersSinceOneHetPar,
			     bool &zeroRecombsThisPrev);
    static inline uint8_t calcAmbigParHetBits(uint8_t hetPar1, uint8_t hetPar2,
					      uint8_t &hetParIsDiff);
    static State * lookupState(const uint64_t iv, const uint64_t ambig,
			       const uint64_t unassigned,
			       int &lowOrderChildBit);
    static int lowOrderUnambigUnassignedBit(const uint64_t allAmbig,
					    const uint64_t unassigned,
					    uint64_t &unambig,
					    uint64_t &ambigStd);
    static bool checkMinUpdate(uint64_t fullIV, uint64_t fullUnassigned,
			       uint64_t ambig1Unassigned, State *theState,
			       const State *prevState, uint8_t hetParent,
			       uint8_t homParentGeno, uint8_t curParPhase,
			       uint8_t altPhaseType, uint8_t ambigLocal,
			       uint8_t prevHMMIndex, uint32_t prevIndex,
			       uint8_t IVambigPar, float minMaxRec[2],
			       int16_t numMarkersSinceNonHetPar,
			       int16_t numMarkersSinceOneHetPar,
			       float totalRecombs, float totalLikehood,
			       size_t numRecombs, int lowOrderChildBit);
    static int8_t decideOptimalState(State *theState, const State *prevState,
				     uint8_t prevHMMIndex, float totalRecombs,
				     float totalLikelihood,
				     size_t newRecombs);
    static void updateAmbigPrev(State *theState, uint32_t prevIndex,
				bool newBestPrev);
    static void rmBadStatesCheckErrorFlag(dynarray<State*> &curStates,
					  float minMaxRec[2],
					  int numChildren);
    static void backtrace(NuclearFamily *theFam, uint8_t missingPar,
			  int chrFirstMarker, int chrLastMarker);
    static uint32_t findMinStates(dynarray<State*> &theStates);
    static void deleteStates(dynarray<State*> &theStates);
    static void propagateBackIV(uint64_t curIV, uint64_t curAmbig,
				State *prevState, uint64_t &newPrevIV,
				uint64_t &newPrevAmbig,
				uint8_t &newPrevParPhase,
				uint8_t &newPrevAmbigParPhase,
				uint8_t &numAmbig1Recombs);
    static void collectAmbigPrevIdxs(uint64_t curIV, uint64_t curAmbig,
				     uint8_t ambigPrev, uint32_t prevState,
				     const dynarray<State*> &prevStates,
				     BT_ambig_info &thePrevInfo,
				     uint8_t &numAmbig1Recombs);
    static void tracePrior(int hmmIndex, int stateIdx);
    static void printStates();
    static void printStates(int hmmIndex);


    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    static dynarray< dynarray<State*> > _hmm;
    static dynarray<int> _hmmMarker;
    static dynarray< dynarray<uint32_t> > _ambigPrevLists;

    // States at previous markers (up to CmdLineOpts::errorLength previous
    // informative sites) that (a) transition to the next marker with >= 1
    // recombination and therefore (b) have the potential to be a previous
    // state along a path that treats one or more (up to
    // CmdLineOpts::errorLength) markers as errors
    static std::deque< std::pair<int, uint32_t> > _prevErrorStates;

    // Container for parent and children's genotypes. Needed occasionally
    // during back tracing
    static dynarray< std::pair<uint8_t,uint64_t> > _genos;

    // Hash table (including all the book keeping stuff) to store states
    // to enable fast lookup of states a given previous state maps to
    typedef typename std::tuple<uint64_t,uint64_t,uint64_t>* iv_ambig;
    typedef typename std::tuple<uint64_t,uint64_t,uint64_t> iv_ambig_real;
    struct hashIVAmbig {
      size_t operator()(iv_ambig const key) const {
	// make a better hash function?
	return std::hash<uint64_t>{}(std::get<0>(*key)) +
	       std::hash<uint64_t>{}(std::get<1>(*key)) +
	       // unassigned bits aren't expected to be anything but 0 for most
	       // values so won't distinguish them.
	       // When they are non-zero, we also expect the other two values
	       // to distinguish the state and that it will only rarely be the
	       // case that a state with equivalent IV and ambig values will
	       // have differing unassigned bits.
	       std::get<2>(*key);
      }
    };
    struct eqIVAmbig {
      bool operator()(const iv_ambig k1, const iv_ambig k2) const {
	return k1 == k2 ||
	  (k1 && k2 && std::get<0>(*k1) == std::get<0>(*k2) &&
		       std::get<1>(*k1) == std::get<1>(*k2) &&
		       std::get<2>(*k1) == std::get<2>(*k2));
      }
    };
    typedef typename google::dense_hash_map<iv_ambig, State *, hashIVAmbig,
					    eqIVAmbig> state_ht;
    typedef typename state_ht::const_iterator state_ht_iter;
    static state_ht _stateHash;

    // for back tracing; see comment about struct BT_ambig_info declaration
    struct hashBTinfo {
      size_t operator()(BT_ambig_info const key) const {
	// make a better hash function?
	return std::hash<uint32_t>{}(key.stateIdx) +
	       std::hash<uint64_t>{}(key.iv) +
	       // during back tracing, most ambig bits will be 0 so won't
	       // really distinguish them, so we won't hash but use the raw
	       // value:
	       key.ambig +
	       std::hash<uint8_t>{}( (key.ambigParPhase << 2) | key.parPhase);
      }
    };
    struct eqBTinfo {
      bool operator()(const BT_ambig_info k1, const BT_ambig_info k2) const {
	return k1.stateIdx == k2.stateIdx && k1.iv == k2.iv &&
	       k1.ambig == k2.ambig && k1.parPhase == k2.parPhase &&
	       k1.ambigParPhase == k2.ambigParPhase;
      }
    };
    typedef typename google::dense_hash_set<BT_ambig_info,
					    hashBTinfo, eqBTinfo> BT_state_set;
    typedef typename BT_state_set::const_iterator state_set_iter;
    // For when there are ambiguous previous states
    static BT_state_set *_curBTAmbigSet;
    static BT_state_set *_prevBTAmbigSet;

    // What type of phasing are we doing?
    static PhaseMethod _phaseMethod;

    // Inheritance vector bits corresponding to parent 0 and 1 (indexed).
    // Index 2 corresponds to all inheritance vector bits set to 1.
    // Note these values are variable for each family depending on the number of
    // children in it.
    static uint64_t _parBits[3];
    // When changing the phase one or both parents, the inheritance vector
    // values corresponding to that parent are flipped. For the four different
    // possibilities, corresponding to the bits to be flipped for one child,
    // the following will contain the flip value that applies to all children
    // in the current family. (That is, the index value but repeated every 2
    // bits up to the 2*<numChildren>-th bit).
    static uint64_t _flips[4];
    // As above, but applies to ambiguous inheritance vector values. See the
    // comment in lookupState() which talks about the fact that there are
    // two distinct values only instead of four.
    static uint64_t _ambigFlips[4];

    // For tracking information about genetic distances
    static int _lastInformMarker;
    static int _curMarker;
};

struct State {
  // Inheritance vectors store two bits for each child, with bits 0,1 for child
  // 0; 1,2 for child 1, etc. Even bits indicate which haplotype (0 or 1)
  // parent 0 transmitted to the corresponding child, and odd bits indicate
  // the transmitted haplotype from parent 1.
  uint64_t iv;         // inheritance vector

  // Ambiguous bits: at partly informative markers, when a child is heterozygous
  // and recombines relative to the previous marker, which parent transmitted
  // a recombination -- and therefore what the iv value is -- is ambiguous.
  // This ambiguity needs to be tracked so that it can (a) be resolved using
  // information at later markers and (b) won't be treated as fixed with
  // respect to calculating numbers of recombinations at later markers.
  uint64_t ambig;      // ambiguous inheritance vector bits

  // TODO: comment on dual use of <unassigned>? In partial states it is a
  // mask for propagating from the previous marker. In full states it indicates
  // which bits shouldn't count for recombination

  // When a parent is homozygous, the corresponding inheritance vector bits
  // are unknown. At the beginning of the chromosome, some bits in <iv> are
  // invalid because there has not yet been a marker that is heterozygous for
  // one parent. This value indicates the bits that should not be considered
  // meaningful in <iv>.
  // Similarly, at the beginning of the chromosome, regardless of the
  // informativeness of the parent's genotype, if a child is missing data,
  // the corresponding <iv> values will not be assigned and so aren't considered
  // meaningful.
  // A more complex scenario revolves around ambig1 values. See the comments
  // in fixRecombFromAmbig() -- if an ambig1 assignment is made and a
  // subsequent fully informative marker has a recombination for the
  // corresponding child, this becomes equivalent to the child having one
  // parent unassigned so this value can be non-zero even after a partly
  // informative marker (where both parents are heterozygous)
  uint64_t unassigned; // iv bits that haven't been assigned up to this marker

  // index of previous state that leads to optimal phase here
  // See comment above <ambigPrev> below for the meaning of this value when
  // the previous state is ambiguous.
  // TODO: add checks to ensure the number of states does not grow beyond
  // the capacity here? Unlikely to be more than 2^32 states but could check
  uint32_t prevState;

  // relative HMM index -- many HMM indexes back -- for the previous state
  // If non-error, this will be 1
  uint8_t prevHMMIndex;

  // Ambiguous previous state? If this value is 1, then <prevState> stores
  // in every 6 bits starting from the lowest order bits, the indexes to
  // the ambiguous previous states.
  // This value is 2 whenever there are more than 32 / 6 == 5 ambiguous
  // previous states or if the index of any previous states is >=2^6 == 64,
  // the previous states are stored in a list inside <_ambigPrevLists> with
  // the value <prevState> being an index into that list.
  uint8_t  ambigPrev;  // fits in 2 bits

  // Minimum number of recombinations to reach this state
  float minRecomb;

  // For detecting cases where a parent without data has transmitted only one
  // haplotype to all children. Stores the number of markers (including
  // uninformative markers) since the last state that can reach this one that
  // was heterozygous for the parent in question.
  // This is signed and set to be negative when the penalty associated with
  // treating the IV as if all children had the same haplotype has been applied.
  int16_t numMarkersSinceNonHetPar;

  // For detecting cases where the children's IVs have been swapped from a
  // parent for which we have data to a parent without data. This manifests
  // as a series of markers that are either uninformative or where both parents
  // are heterozygous (since sites where the true heterozygous parent are
  // modeled as if they're heterozygous for the missing data parent).
  // Stores the number of markers since the last state that can reach this one
  // that was heterozygous for only one parent.
  // This is signed and set to be negative when the penalty applied to states
  // where both parents are heterozygous over a long stretch (which in practice
  // deletes that state path).
  int16_t numMarkersSinceOneHetPar;

  // We prefer to back trace to states that have recombinations immediately
  // before the current state as opposed to earlier. In general this decision
  // doesn't matter much since we track ambiguous <iv> values (see <ivFlippable>
  // in the backtrace() method).
  // Note: will not ever get more than 64 recombinations from the previous
  // marker since there are 64 bits in <iv>, so cannot exceed UINT8_MAX
  uint8_t  maxPrevRecomb;

  // Log likelihood of state
  float    maxLikelihood;

  // Which parent is heterozygous? Either 0, 1, or 2 for both, which corresponds
  // to MT_PI states
  uint8_t  hetParent;  // fits in 2 bits

  // Ambiguous as to which parent is heterozygous?
  // Each bit corresponds to a value of <hetParent>, which takes on values 0,1,2
  // If the corresponding bit is set to 1, it's possible for the state to have
  // the indicated <hetParent> value.
  uint8_t  ambigParHet;  // fits in 3 bits

  // Genotype of homozygous parent (may be missing)
  uint8_t  homParentGeno;  // fits in 2 bits

  // Four possible phase types for the parents: default heterozygous assignment
  // has allele 0 on haplotype 0 and allele 1 on haplotype 1. Can flip this
  // for each parent. Bit 0 here indicates (0 or 1) whether to flip parent 0
  // and bit 1 indicates whether to flip parent 1
  uint8_t  parentPhase;  // fits in 2 bits

  // Ambiguous parent phase at this marker? Because <hetParent> is subject to
  // change depending on the previous state (only when there is ambiguity),
  // we track three sets of values, stored sequentially. <hetParent> == 0 are
  // stored in bits 0 & 1, <hetParent> == 1 are in bits 2 & 3, and
  // <hetParent> == 2 are in bits 4-7. There are 2^H possible phase types where
  // H is the number of heterozygous parents. If the corresponding bit is set,
  // the phase type is valid.
  uint8_t  ambigParPhase;  // requires 8 bits

  // Arbitrary choice for parent of origin? Occurs at the beginning of the
  // chromosome when we don't have data for either parent. Also occurs within
  // and just after regions where the parents transmit IV values that are either
  // identical or exactly opposite.
  // In practice this means arbitrarily choosing one of the parents as
  // heterozygous at FI markers and one of two (equivalent but opposite)
  // phase assignments for the two parents at PI markers
  uint8_t  arbitraryPar;  // fits in 1 bit

  // Error state? Used in the context of <CmdLineOpts::max1MarkerRecomb>
  // A value of 1 means <this> is an error state.
  // A value of 2 means that some previous state leading to this one is an
  // error state. Tracking whether some state leading to this one is an error
  // state allows us to prefer non-error state paths when there are ambiguities
  uint8_t  error;  // fits in 2 bits
};

#endif // PHASER_H
