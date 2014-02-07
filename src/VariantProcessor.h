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
  pos_t processMatchOrMismatch(const BamTools::BamAlignment& alignment,
			       std::vector<VariantPtr>& read_variants,
			       const uint32_t& op_length, const std::string& refseq,
			       const pos_t& refpos, const pos_t& readpos);
  pos_t processInsertion(const BamTools::BamAlignment& alignment,
			 std::vector<VariantPtr>& read_variants,
			 const uint32_t& op_length, const std::string& refseq,
			 const pos_t& refpos, const pos_t& readpos);
  pos_t processDeletion(const BamTools::BamAlignment& alignment,
			std::vector<VariantPtr>& read_variants,
			const uint32_t& op_length, const std::string& refseq,
			const pos_t& refpos, const pos_t& readpos);
  void processBlockAlignments();
  int run();

private:
  const BamTools::RefVector references;
  BamTools::BamReader reader;

  int32_t curr_refid;
  pos_t last_aln_pos; // for ensuring everything is sorted

  // block-level info
  std::set<VariantPtr> block_variants;
  std::deque< std::vector<VariantPtr> > block_read_variants;
  std::deque<BamTools::BamAlignment> block_alignments; // for mapped reads
  pos_t block_start;
  pos_t block_end;
  pos_t block_last_variant_pos;

  void addVariant();
};

#endif
