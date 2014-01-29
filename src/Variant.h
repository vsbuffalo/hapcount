#ifndef _VARIANT_H
#define _VARIANT_H

#include <string>
#include <memory>
#include "api/BamReader.h"

typedef long int pos_t;

class Variant; 

// For use in containers
typedef std::shared_ptr<Variant> VariantPtr;

enum class VariantType {
  SNP,
  Insertion,
  Deletion
};

class Variant {
  // Variants are any difference between reference and aligned read,
  // carried through the MD/CIGAR flags. These may be true genetic
  // variants or errors -- it's up to downstream code to figure this
  // out.
  // TODO
  // - add comparison operator for STL Set
public:
  Variant(VariantType type, int32_t refid, pos_t position, unsigned int length, 
	  double quality, std::string refseq, std::string altseq): 
    type(type),
    refid(refid),
    position(position),
    length(length),
    quality(quality),
    refseq(refseq),
    altseq(altseq) {};

  VariantType type;
  int32_t refid; // from BAM header
  pos_t position;
  unsigned int length;
  double quality;
  // For SNPs; altseq for insertions too
  std::string refseq; 
  std::string altseq;

  void print();
private: 
  Variant(const Variant&);
  Variant& operator=(const Variant&);
};

struct BedRegion {
  std::string seqname;
  pos_t start;
  pos_t end;  
  BedRegion(std::string sn, pos_t s, pos_t e)
    :seqname(sn), start(s), end(e) {};
};

#endif
