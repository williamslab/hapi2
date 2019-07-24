#include <stdint.h>

// For calculating canonical State::ambigParPhase values; need to update these
// bits, which are phase swap types, according to a <flipType> value -- the
// first index to swap01Phase and swap2Phase below. The <flipType> takes on
// values 0,1,2,3 corresponding, respectively, to inverting neither parent,
// inverting parent 0, inverting parent 1, and inverting both parents.

// The next value is for the first four bits of <State::ambigParPhase>, which
// indicate the possible phase swap types for parents 0 (bits 0 and 1) and
// 1 (bits 2 and 3).
uint8_t swap01Phase[4][16] = {  // swap type 0: no change:
				{  0,  1,  2,  3,
				   4,  5,  6,  7,
				   8,  9, 10, 11,
				  12, 13, 14, 15 },
				// swap type 1: invert parent 0
                                {  0,  2,  1,  3,   // parent 1 bits = 0
	                           4,  6,  5,  7,   // parent 1 bits = 1
			           8, 10,  9, 11,   // parent 1 bits = 2
			          12, 14, 13, 15 }, // parent 1 bits = 3
				// swap type 2: invert parent 1:
				{  0,  1,  2,  3,   // parent 1 bits = 0
				   8,  9, 10, 11,   // parent 1 bits = 1
				   4,  5,  6,  7,   // parent 1 bits = 2
				  12, 13, 14, 15 },
				// swap type 3: invert both parents:
				{  0,  2,  1,  3,   // only parent 0 bits here
				   8, 10,  9, 11,
				   4,  6,  5,  7,
				  12, 14, 13, 15 } };

// For the second set of four bits (bits 4-7) in State::ambigParPhase.
// This is complex: each bit in the value represents a phase swap type. Taking
// the log of the bits makes this transformation more intuitive, but we here we
// just hard code the conversion for efficiency. If we had taken the log, the
// comments below indicate the changes. (That is, the log2 of each bit
// successively gives the bit index, and the comments below reference these bit
// indices.)
uint8_t swap2Phase[4][16] = {
			      // swap type 0: no change:
			      {  0,  1,  2,  3,
			         4,  5,  6,  7,
			         8,  9, 10, 11,
			        12, 13, 14, 15 },
			      // swap type 1: invert parent 0
			      // 0 <-> 1 and 2 <-> 3  (00 <-> 01; 10 <-> 11)
			      {  0,    2,    1,    3,
			         8,  8+2,  8+1,  8+3,
			         4,  4+2,  4+1,  4+3,
			        12, 12+2, 12+1, 12+3 },
			      // swap type 2: invert parent 1
			      // 0 <-> 2 and 1 <-> 3  (00 <-> 10; 01 <-> 11)
			      {  0,   4,   8,   12,
			         1, 1+4, 1+8, 1+12,
			         2, 2+4, 2+8, 2+12,
			         3, 3+4, 3+8, 3+12 },
			      // swap type 3: invert both parents
			      // 0 <-> 3 and 1 <-> 2  (00 <-> 11; 01 <-> 10)
			      {  0,   8,   4,   12,
			         2, 2+8, 2+4, 2+12,
			         1, 1+8, 1+4, 1+12,
			         3, 3+8, 3+4, 3+12 } };
