#include "DNAAnalyzer.h"
#include <iostream>
#include <string>

using namespace std;

int main () {

	cout << "*** DNA Analyzer ***" << endl << endl;

	DNAAnalyzer cAnalyzer1 ("dnastrands.txt", ">C1");

	cout << "DNA Strings" << endl;
	cAnalyzer1.displayStrandSet (cout);

	cout << "Profile" << endl;
	cAnalyzer1.displayProfile (cout);

	cout << "Consensus" << endl;
	cAnalyzer1.displayConsensus (cout);

	return EXIT_SUCCESS;
}