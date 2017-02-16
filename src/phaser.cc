// HAPI 2: HAPlotype Inference
//
// This program is distributed under the terms of the GNU General Public License

#include "phaser.h"

void Phaser::run(PersonBulk::par_pair parents, dynarray<PersonBulk*> &children,
		 int chrIdx) {
  // Get first and last marker numbers for the chromosome to be phased:
  int firstMarker = Marker::getFirstMarkerNum(chrIdx);
  int lastMarker = Marker::getLastMarkerNum(chrIdx);

  int numChildren = children.length();

  // build states/phase each marker
  for(int m = firstMarker; m <= lastMarker; m++) {
    // Get data for dad (first 2 bits) and mom (bits 3-4)
    uint8_t dadData = parents->first->getBitGeno(m);
    uint8_t momData = parents->second->getBitGeno(m);
    uint8_t parentData = dadData + (momData << 2);

    // Which of the genotypes are present in the parents/children? There are
    // four possible genotypes, and the first four bits will have a value of 1
    // iff the corresponding genotype value is present in at least one child
    uint8_t parentGenoTypes = (1 << dadData) | (1 << momData);
    uint8_t childGenoTypes = 0;

    uint64_t childrenData = 0;
    for(int c = 0; c < numChildren; c++) {
      uint8_t curChildData = children[c]->getBitGeno(m);
      childrenData += curChildData << (c*2);
      childGenoTypes |= 1 << curChildData; // observed genotype <curChildData>
    }

    // Step 1: Determine marker type and check for Mendelian errors
    int mt = getMarkerType(parentGenoTypes, childGenoTypes);
    assert(mt > 0);

    if (mt & ((1 << MT_ERROR) | (1 << MT_AMBIG))) {
      // TODO: indicate that the marker is erroneous / ambiguous
    }
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
    // Parents are both homozygous for cases 1, 8, 9:
    case 1:  // only bit 0 set:    both homozygous for 0
      if (childGenoTypes & 12) // mask is 1100
	// Mendelian error -- invalid bits set: only hom for 0 and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;
      break;
    case 8:  // only bit 3 set:    both homozygous for 1
      if (childGenoTypes & 5) // mask is 0101
	// Mendelian error -- invalid bits set: only hom for 1 and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;
      break;
    case 9:  // bits 0 and 3 set:  one homozygous for 0, other 1
      if (childGenoTypes & 9) // mask is 1001
	// Mendelian error -- invalid bits set: only het and missing are
	// possible
	return 1 << MT_ERROR;
      else
	// both parents homzoygous: uninformative marker
	return 1 << MT_UN;

    //////////////////////////////////////////////////////////////////////////
    // One parent homozygous the other heterozygous for cases 5, 12
    case 5:  // bits 0 and 2 set:  one homozygous for 0, other heterozygous
      if (childGenoTypes & 8) // mask is 1000
	// Mendelian error -- child is homozygous for allele not present in
	// the homozygous parent
	return 1 << MT_ERROR;
      else
	// one parent heterozygous, other homozygous: informative for one parent
	return 1 << MT_FI_1;
    case 12: // bits 3 and 2 set:  one homozygous for 1, other heterozygous
      if (childGenoTypes & 1) // mask is 0001
	// Mendelian error -- child is homozygous for allele not present in
	// the homozygous parent
	return 1 << MT_ERROR;
      else
	// one parent heterozygous, other homozygous: informative for one parent
	return 1 << MT_FI_1;

    //////////////////////////////////////////////////////////////////////////
    // Both parents heterozygous for case 4
    case 4:  // only bit 2 set:    both heterozygous
      // both heterozygous for same alleles (by virtue of biallelic only data):
      // partly inforamtive
      //
      // No Mendelian errors possible in this case
      return 1 << MT_PI;

    //////////////////////////////////////////////////////////////////////////
    // Both parents missing for case 2
    case 2:  // only bit 1 set:    both missing
      // Further analysis done below using observed children genotypes
      missingType = 1;
      break;

    //////////////////////////////////////////////////////////////////////////
    // One parent missing, other non-missing for case 2
    // Further analysis done below using observed children genotypes
    // (including checks for Mendelian errors)
    case 3:  // bits 0 and 1 set: one parent homozgyous for 0, other missing
      missingType = 0;
      break;
    case 10: // bits 3 and 1 set: one parent homzoygous for 1, other missing
      missingType = 3;
      break;
    case 6:  // bits 2 and 1 set: one parent heterozygous, other missing
      missingType = 2;
      break;

    default:
      fprintf(stderr, "ERROR: got impossible parent genotype value %d\n",
	      parentGenoTypes);
      exit(5);
      break;
  }

  assert(missingType >= 0 && missingType <= 3);

  switch (childGenoTypes) {
    case 9:  // both homozygous types observed
    case 11: // above and missing type observed
    case 13: // both homozgyous types and heterozygous type observed
    case 15: // above and missing type observed
      // if children have both homozygous types, parents must both be
      // heterozygous: partly informative
      //
      // Ensure this accords with what we know about the parent genotypes:
      switch(missingType) {
	case 1: // both missing data -- no information to add -- consistent
	case 2: // known parent is heterozgyous -- consistent
	  return 1 << MT_PI;
	case 0: // one parent homozygous for 0
	case 3: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for
	  // allele that is not present in the parent. Since both homozygous
	  // types are present, this is violated.
	  return 1 << MT_ERROR;
      }
      break;
    case 1:  // only homozygous for 0 observed
    case 3:  // above and missing type observed
      // could in principle be anything; see if there is information on one
      // of the parents that constrains the type:
      switch (missingType) {
	case 1: // both missing data -- could be any type:
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case 2: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  return (1 << MT_UN) | (1 << MT_FI_1);
	case 3: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
      }
      break;
    case 8:  // only homozygous for 1 observed
    case 10: // above and missing type observed
      // could in principle be anything; see if there is information on one
      // of the parents that constrains the type:
      switch (missingType) {
	case 1: // both missing data -- could be any type:
	  return (1 << MT_UN) | (1 << MT_FI_1) | (1 << MT_PI);
	case 2: // one heterozygous parent
	  // can't be uninformative with one heterozygous parent:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case 3: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  return (1 << MT_UN) | (1 << MT_FI_1);
      }
      break;
    case 5:  // bits 0 and 2 set: homozygous and heterozygous types observed
    case 7:  // above and missing type observed
      // could be informative for one parent or partly informative; see if there
      // is information on one of the parents that constrains the type:
      switch (missingType) {
	case 1: // both missing data -- no information to add
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 2: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 0: // one parent homozygous for 0
	  // can't be partly informative with one homozygous parent:
	  return 1 << MT_FI_1;
	case 3: // one parent homozygous for 1
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
      }
      break;
    case 12: // bits 3 and 2 set: homozygous and heterozygous types observed
    case 14: // above and missing type observed
      // could be informative for one parent or partly informative; see if there
      // is information on one of the parents that constrains the type:
      switch (missingType) {
	case 1: // both missing data -- no information to add
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 2: // one heterozygous parent
	  // other parent could be either type, so no constraining:
	  return (1 << MT_FI_1) | (1 << MT_PI);
	case 0: // one parent homozygous for 0
	  // Mendelian error: can't get children that are homozygous for allele
	  // that is not present in the parent.
	  return 1 << MT_ERROR;
	case 3: // one parent homozygous for 1
	  // can't be partly informative with one homozygous parent:
	  return 1 << MT_FI_1;
      }
      break;
    case 4:  // only heterozygous types observed
    case 6:  // above and missing type observed
      // This is ambiguous unless we have a parent that is homozygous
      switch (missingType) {
	case 1: // both missing data -- no information to add
	  return 1 << MT_AMBIG;
	case 2: // one heterozygous parent
	  // no hints about phasing with all observed samples heterozygous
	  return 1 << MT_AMBIG;
	case 0: // one parent homozygous for 0
	case 3: // one parent homozygous for 1
	  // other parent must be heterozygous, so fully informative for that
	  // parent
	  return 1 << MT_FI_1;
      }
      break;
    case 2:  // all children missing data
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
