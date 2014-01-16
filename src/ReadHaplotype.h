#ifndef _READHAPLOTYPE_H
#define _READHAPLOTYPE_H

#include <string>
#include <vector>
#include <deque>
#include "api/BamReader.h"
#include "ReadHaplotype.h"
#include "Variant.h"

class ReadHaplotype {
  // The ReadHaplotype class stores a reduced representation of a BAM
  // alignment, storing information about what variants they contain
  // and their positions. Each object contains as much information
  // about the read's haplotype as gatherable from the alignment. To
  // prevent excess object creation, a read's variants are stored as
  // pointers. Read level-variant information like quality is stored
  // as an accompanying private members.
  // Note: FIXME: this has to be changed, bad interface.
public:
  ReadHaplotype();
  ~ReadHaplotype();
  const std::vector<Variant*> GetReadVariants();
  const Variant* GetReadVariant(pos_t pos); // a read can only have a single variant per position
  
private:
  pos_t start; 
  pos_t end; // end of alignment, excluding soft clipped
  int32_t refid;
  std::vector< BamTools::CigarOp > cigar;
  std::string md;
  int nm;
  
  std::vector<Variant*> variants; // for all variants in read
};

#endif
