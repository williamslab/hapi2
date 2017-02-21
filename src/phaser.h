// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/personbulk.h>

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

    static void run(PersonBulk::par_pair parents,
		    dynarray<PersonBulk*> &children, int chrIdx);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void init(int numChildren);
    static void getFamilyData(PersonBulk::par_pair parents,
			      dynarray<PersonBulk*> &children, int marker,
			      uint8_t &parentData, uint8_t &parentGenoTypes,
			      uint64_t childrenData[5],uint8_t &childGenoTypes);
    static int  getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes);
    static void makePartialStates(dynarray<State> &partialStates,
				  int markerTypes, uint8_t parentData,
				  uint64_t childrenData[5]);
    static void makePartialFI1States(dynarray<State> &partialStates,
				     uint8_t parentData,
				     uint64_t childrenData[5]);
    static void makePartialPIStates(dynarray<State> &partialStates,
				    uint8_t parentData,
				    uint64_t childrenData[5]);


    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    static dynarray< dynarray<State> > _hmm;
    static dynarray<int> _hmmMarker;

    // All inheritance vector bits set to 1 (variable value depending on the
    // number of children in the family)
    static uint64_t _allBitsSet;
    // Inheritance vector bits corresponding to parent 0 and 1 (indexed)
    static uint64_t _parBits[2];
};

struct State {
  // Inheritance vectors store two bits for each child, with bits 0,1 for child
  // 0; 1,2 for child 1, etc. Even bits indicate which haplotype (0 or 1)
  // parent 0 transmitted to the corresponding child, and odd bits indicate
  // the transmitted haplotype from parent 1.
  uint64_t iv;         // inheritance vector

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
  uint16_t prevState;
};

#endif // PHASER_H
