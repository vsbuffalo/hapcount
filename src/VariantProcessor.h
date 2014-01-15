#ifndef _VARIANTPROCESSOR_H
#define _VARIANTPROCESSOR_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include "api/BamReader.h"

#include "Variant.h"

using namespace std;
using namespace BamTools;

class VariantProcessor {
public:
  string filename;
  
  VariantProcessor(const string& filename): filename(filename) {
    if (!reader.Open(filename)) {
      std::cerr << "count not open bamfile '" << filename << "'." << std::endl;
      std::abort();
    }
    last_aln_pos = 0; // TODO change when region added
  };
  void blockReset(); // for initialization
  void blockReset(pos_t pos);
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
