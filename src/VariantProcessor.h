#ifndef _VARIANTPROCESSOR_H
#define _VARIANTPROCESSOR_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include "api/BamReader.h"

#include "Variant.h"

class VariantProcessor {
public:
  std::string filename;
  
  VariantProcessor(const std::string& filename): filename(filename) {
    if (!reader.Open(filename)) {
      std::cerr << "count not open bamfile '" << filename << "'." << std::endl;
      std::abort();
    }
    last_aln_pos = 0; // TODO change when region added
  };
  void blockReset(); // for initialization
  void blockReset(pos_t pos);
  pos_t processAlignment(const BamTools::BamAlignment& alignment);
  bool isBlockEnd(const BamTools::BamAlignment& alignment);
  int run();
  
private:
  const BamTools::RefVector references;
  BamTools::BamReader reader;

  int32_t curr_refid;
  pos_t last_aln_pos; // for ensuring everything is sorted

  // block-level info
  std::set<Variant>* block_variants;
  std::deque< std::vector<Variant> >* block_reads_variants;
  std::deque<BamTools::BamAlignment>* block_reads; // for mapped reads
  pos_t block_start;
  pos_t block_end;
  pos_t block_last_variant_pos;

};

#endif
