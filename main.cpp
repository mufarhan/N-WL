#include <iostream>
#include "NWL_WeisfeilerLehman.h"

using namespace std;

int main(int argc, char **argv) {

  NWL_WeisfeilerLehman nwl(argv[1], atoi(argv[2]), atoi(argv[3]));

  return 0;
}
