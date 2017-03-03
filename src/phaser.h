// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/nuclearfamily.h>
#include <genetio/personbulk.h>
#include <sparsehash/dense_hash_map>
#include <bitset>

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

// For reference, the following are the 2-bit values of all the genotypes
// in PLINK format data:
// 0 - homozygous for allele 0
// 1 - missing
// 2 - heterozygous
// 3 - homozygous for allele 1
enum Geno {
  G_HOM0 = 0,
  G_MISS = 1,
  G_HET  = 2,
  G_HOM1 = 3
};

struct State;

class Phaser {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void run(NuclearFamily *theFam, int chrIdx);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void init(int numChildren);
    static void getFamilyData(NuclearFamily *theFam, int marker,
			      uint8_t &parentData, uint8_t &parentGenoTypes,
			      uint64_t childrenData[5],uint8_t &childGenoTypes);
    static int  getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes);
    static void makePartialStates(dynarray<State> &partialStates,
				  int markerTypes, uint8_t parentData,
				  const uint64_t childrenData[5]);
    static void makePartialFI1States(dynarray<State> &partialStates,
				     uint8_t parentData,
				     const uint64_t childrenData[5]);
    static void makePartialPIStates(dynarray<State> &partialStates,
				    uint8_t parentData,
				    const uint64_t childrenData[5]);
    static void makeFullStates(const dynarray<State> &partialStates, int marker,
			       const uint64_t childrenData[5]);
    static void calcHetChildPIVals(const uint64_t recombs,
				   const uint64_t unambigHets,
				   uint64_t unambigHetRecombs[4]);
    static void flipPIVals(uint64_t &fullIV, uint64_t &fullAmbig,
			   const uint64_t childrenData[5],
			   uint64_t propagateAmbig,
			   const uint64_t unambigHetRecombs[4]);
    static void updateStates(uint64_t fullIV, uint64_t fullAmbig,
			     uint64_t fullUnassigned, uint64_t recombs,
			     uint8_t hetParent, uint8_t initParPhase,
			     uint8_t parPhaseFlip, uint16_t prevIndex,
			     uint16_t prevMinRecomb,
			     const uint64_t childrenData[5]);
    static State * lookupState(const uint64_t iv, const uint64_t ambig);
    static size_t popcount(uint64_t val) {
      // Currently using std::bitset, but the conversion is probably not free,
      // so optimize? TODO
      std::bitset<64> to_count(val);
      return to_count.count();
    }
    static void backtrace(NuclearFamily *theFam);
    static void findMinStates(int hmmIndex);


    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    static dynarray< dynarray<State*> > _hmm;
    static dynarray<int> _hmmMarker;
    static dynarray<uint16_t> _minStates;

    // Hash table (including all the book keeping stuff) to store states
    // to enable fast lookup of states a given previous state maps to
    typedef typename std::pair<uint64_t,uint64_t>* iv_ambig;
    typedef typename std::pair<uint64_t,uint64_t> iv_ambig_real;
    struct hashIVAmbig {
      size_t operator()(iv_ambig const key) const {
	// make a better hash function?
	return std::tr1::hash<uint64_t>{}(key->first) +
	       std::tr1::hash<uint64_t>{}(key->second);
      }
    };
    struct eqIVAmbig {
      bool operator()(const iv_ambig k1, const iv_ambig k2) const {
	return k1->first == k2->first && k1->second == k2->second;
      }
    };
    typedef typename google::dense_hash_map<iv_ambig, State *, hashIVAmbig,
					    eqIVAmbig> state_ht;
    typedef typename state_ht::const_iterator state_ht_iter;
    static state_ht _stateHash;

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
  uint16_t prevState;

  // Minimum number of recombinations to reach this state
  uint16_t minRecomb;

  // Which parent is heterozygous? Either 0, 1, or 2 for both, which corresponds
  // to MT_PI states
  uint8_t  hetParent  : 2;

  // Four possible phase types for the parents: default heterozygous assignment
  // has allele 0 on haplotype 0 and allele 1 on haplotype 1. Can flip this
  // for each parent. Bit 0 here indicates (0 or 1) whether to flip parent 0
  // and bit 1 indicates whether to flip parent 1
  uint8_t  parentPhase : 2;
};

#endif // PHASER_H
