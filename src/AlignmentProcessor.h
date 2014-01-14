#ifndef _ALIGNMENTPROCESSOR_H
#define _ALIGNMENTPROCESSOR_H

#include <string>
#include <vector>
#include <set>
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;


typedef unsigned long int pos_t;

enum VariantType {
  SNP,
  INSERTION,
  DELETION
};

class Variant {
public:
  VariantType type;
  int32_t refid; /* from BAM header */
  pos_t position;
  unsigned int length;
  double quality;
  /* For SNPs; altseq for insertions too */
  string refseq; 
  string altseq;
  
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
  void processAlignment(const BamAlignment& alignment);
  void isBlockEnd(const BamAlignment& alignment);
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
