// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/personbulk.h>

#ifndef PHASER_H
#define PHASER_H

// Types of markers. We use bits set in an integer for markers that can
// have multiple types (occurs when one or both parents are missing data)
enum MarkerType {
  MT_FI_1, // Fully informative for one parent
  MT_PI,   // Partly informative (both parents heterozygous for same alleles)
  MT_UN,   // Uninformative (both parents homozygous)
  MT_AMBIG,// All children and observed parents are heterozygous or missing
  MT_ERROR,// Mendelian error marker
};
// Note: fully informative for both parents is impossible in biallelic genotype
// data. As we're using PLINK BED format data at present, we are using
// biallelic data.

// For reference, the following are the 2-bit values of all the genotypes
// in PLINK format data:
// 0 - homozygous for allele 0
// 1 - missing
// 2 - heterozygous
// 3 - homozygous for allele 1

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

    static int getMarkerType(uint8_t parentGenoTypes, uint8_t childGenoTypes);
};

#endif // PHASER_H
