#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <getopt.h>

#include "api/BamReader.h"
#include "AlignmentProcessor.h"

using namespace std;

int main(int argc, char *argv[]) {
  string bam_filename(argv[1]);
  AlignmentProcessor* ap = new AlignmentProcessor(bam_filename);

  ap->blockReset();
  int test = ap->run();
  cerr << "there were " << test << " BAM entries processed" << endl;
  
  return 0;
}
