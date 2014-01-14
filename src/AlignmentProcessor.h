#ifndef _ALIGNMENTPROCESSOR_H
#define _ALIGNMENTPROCESSOR_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;


typedef long int pos_t;

enum VariantType {
  SNP,
  INSERTION,
  DELETION
};

class ReadHaplotype {
  // The ReadHaplotype class stores a reduced representation of a BAM
  // alignment, storing information about what variants they contain
  // and their positions. Each object contains as much information
  // about the read's haplotype as gatherable from the alignment.
public:
  ReadHaplotype(const BamAlignment& alignment);
  ~ReadHaplotype();
  const vector<Variant>& GetReadVariants(const Variant& variant);

private:
  pos_t start; 
  pos_t end; // end of alignment, excluding soft clipped
  int32_t refid;
  std::vector< BamTools::CigarOp > cigar;
  std::string md;
  int nm;

  vector<Variant> variants; // for all variants in read
};

class Variant {
  // Variants are any difference between reference and aligned read,
  // carried through the MD/CIGAR flags. These may be true genetic
  // variants or errors -- it's up to downstream code to figure this
  // out.
  // TODO
  // - add comparison operator for STL Set
public:
  Variant(const BamAlignment& alignment);
  ~Variant();
  VariantType type;
  int32_t refid; // from BAM header
  pos_t position;
  unsigned int length;
  double quality;
  // For SNPs; altseq for insertions too
  string refseq; 
  string altseq;
private: 
  //  static unsigned int count; // global variant counts
};

struct BedRegion {
  string seqname;
  pos_t start;
  pos_t end;  
  BedRegion(string sn, pos_t s, pos_t e)
    :seqname(sn), start(s), end(e) {};
};

class AlignmentProcessor {
public:
  string filename;
  
  AlignmentProcessor(const string& filename): filename(filename) {
    if (!reader.Open(filename)) {
      std::cerr << "count not open bamfile '" << filename << "'." << std::endl;
      std::abort();
    }
  };
  void blockReset(pos_t pos);
  void blockReset();
  pos_t processAlignment(const BamAlignment& alignment);
  bool isBlockEnd(const BamAlignment& alignment);
  int run();
  
private:
  const RefVector references;
  BamReader reader;

  int32_t curr_refid;
  pos_t last_aln_pos; // for ensuring everything is sorted

  // block-level info
  set<Variant>* block_variants;
  deque< vector<Variant> >* block_reads_variants;
  deque<BamAlignment>* block_reads; // for mapped reads
  pos_t block_start;
  pos_t block_end;
  pos_t block_last_variant_pos;

};

#endif
