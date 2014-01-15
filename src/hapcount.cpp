#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <getopt.h>

#include "api/BamReader.h"
#include "VariantProcessor.h"
#include "Variant.h"

using namespace std;

int main(int argc, char *argv[]) {
  string bam_filename(argv[1]);
  VariantProcessor* ap = new VariantProcessor(bam_filename);

  ap->blockReset(0); // TODO add region support
  int test = ap->run();
  cerr << "there were " << test << " BAM entries processed" << endl;
  
  return 0;
}
