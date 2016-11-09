#include <fstream>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

int main() {
  const char* filename = "test.txt";
  std::ifstream inFile(filename);

  // Make sure the file stream is good
  if(!inFile) {
    cout << endl << "Failed to open file " << filename;
    return 1;
  }

  long n = 0;
  while(!inFile.eof()) {
    inFile >> n;
    cout << std::setw(10) << n;
  }
  cout << endl;
  return 0;
}

