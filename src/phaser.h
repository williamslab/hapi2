// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/nuclearfamily.h>
#include <genetio/personbulk.h>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

#ifndef PHASER_H
#define PHASER_H

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

class Phaser {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void init() {
      // TODO: optimization: set load factors for these next values?
      _stateHash.set_empty_key(NULL);
      _curIdxSet = new state_idx_set();
      _prevIdxSet = new state_idx_set();
      _curIdxSet->set_empty_key(UINT16_MAX);
      _prevIdxSet->set_empty_key(UINT16_MAX);

      int maxMarkers = 0;
      for(int c = 0; c < Marker::getNumChroms(); c++)
	if (Marker::getNumChromMarkers(c) > maxMarkers)
	  maxMarkers = Marker::getNumChromMarkers(c);
      _hmm.resize(maxMarkers);
      _hmmMarker.resize(maxMarkers);
      _genos.resize(maxMarkers);
      _ambigPrevLists.resize(maxMarkers);
    }
    static void run(NuclearFamily *theFam, int chrIdx);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void parBitsInit(int numChildren);
    static void getFamilyData(NuclearFamily *theFam, int marker,
			      uint8_t &parentData, uint8_t &parentGenoTypes,
			      uint64_t childrenData[5],uint8_t &childGenoTypes);
    static int  getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes,
			      uint8_t &homParGeno);
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
			       const uint64_t childrenData[5],
			       NuclearFamily *theFam);
    static void mapPrevToFull(const State *prevState, int32_t prevIdx,
			      const State &curPartial,
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
			     uint64_t recombs, uint64_t prevUnassigned,
			     uint64_t stdAmbigOnlyPrev, uint64_t ambig1PrevInfo,
			     uint8_t hetParent, uint8_t homParentGeno,
			     uint8_t initParPhase, uint8_t altPhaseType,
			     int32_t prevIndex, uint16_t prevMinRecomb,
			     uint8_t prevError, bool hetParentUndefined,
			     const uint64_t childrenData[5]);
    static State * lookupState(const uint64_t iv, const uint64_t ambig,
			       const uint64_t unassigned);
    static void checkClearErrorFlag(dynarray<State*> &curStates);
    static void backtrace(NuclearFamily *theFam);
    static uint16_t findMinStates(dynarray<State*> &theStates);
    static void deleteStates(dynarray<State*> &theStates);


    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    static dynarray< dynarray<State*> > _hmm;
    static dynarray<int> _hmmMarker;
    static dynarray< dynarray<uint16_t> > _ambigPrevLists;

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
	return std::tr1::hash<uint64_t>{}(std::get<0>(*key)) +
	       std::tr1::hash<uint64_t>{}(std::get<1>(*key)) +
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

    struct eq_uint16 {
      bool operator()(const uint16_t k1, const uint16_t k2) const {
	return k1 == k2;
      }
    };
    typedef typename google::dense_hash_set<uint16_t, std::tr1::hash<uint16_t>,
					    eq_uint16> state_idx_set;
    typedef typename state_idx_set::const_iterator state_set_iter;
    // For when there are ambiguous previous states
    static state_idx_set *_curIdxSet;
    static state_idx_set *_prevIdxSet;

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
  uint64_t unassigned; // iv bits that haven't been assigned up to this marker

  // index of previous state that leads to optimal phase here
  // TODO: add checks to ensure the number of states does not grow beyond
  // the capacity here. Could enlarge if needed.
  // NOTE: if we need more than 2^16 slots for ambiguous previous states,
  // the following approaches suffice:
  // (1) have the field <ambigPrev> below represent a factor representing
  // multiples of 2^16 to add to the base <prevState> value (subtracted by
  // 1 so that <ambigPrev> == 1 is unshifted), thus expanding the allowable
  // range up to 2^(16+N-1), where N is the number of bits allotted to
  // <ambigPrev>
  // (2) Make this a uint16_t and compress the fields below which can fit
  // in very little space
  // (3) During a call to makeFullStates(), keep a list of ambigPrev values
  // to add to <_ambigPrevList>: some states are only ambiguous briefly and
  // wind up being unambiguous after examining all possible previous states
  uint16_t prevState;

  // Minimum number of recombinations to reach this state
  // TODO: add checks to ensure we never reach UINT16_MAX?
  uint16_t minRecomb;

  // Which parent is heterozygous? Either 0, 1, or 2 for both, which corresponds
  // to MT_PI states
  uint8_t  hetParent;     // Note: can fit in 2 bits

  // Genotype of homozygous parent (may be missing)
  uint8_t  homParentGeno; // Note: can fit in 2 bits

  // Four possible phase types for the parents: default heterozygous assignment
  // has allele 0 on haplotype 0 and allele 1 on haplotype 1. Can flip this
  // for each parent. Bit 0 here indicates (0 or 1) whether to flip parent 0
  // and bit 1 indicates whether to flip parent 1
  uint8_t  parentPhase;   // Note: can fit in 2 bits

  // TODO: optimization: more efficient to compress the above values instead?

  // Ambiguous parent phase at this marker? There are four possible parent
  // phase types, and each of the four bits in this value corresponds to one
  // type. If the corresponding bit is set, the phase type is valid.
  uint8_t  ambigParPhase : 4;

  // Ambiguous previous state? If non-zero, then <prevState> is an index
  // into <_ambigPrevLists> that gives a list of all the previous states
  // that map to this state with a total of <minRecomb>
  uint8_t  ambigPrev : 1;

  // Ambiguous as to which parent is heterozygous?
  uint8_t  ambigParHet : 1;

  // Error state? Used in the context of <CmdLineOpts::max1MarkerRecomb>
  // A value of 1 means <this> is an error state.
  // A value of 2 means that some previous state leading to this one is an
  // error state. Tracking whether some state leading to this one is an error
  // state allows us to prefer non-error state paths when there are ambiguities
  uint8_t  error : 2;
};

#endif // PHASER_H
