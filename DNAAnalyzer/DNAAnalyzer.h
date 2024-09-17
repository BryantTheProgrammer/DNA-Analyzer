#pragma once

#include <iostream>
#include <fstream>
#include "DNAStrandSet.h"
#include "DNAStrand.h"

using namespace std;

class DNAAnalyzer {
public:
	enum class Bases : int {
		ADENINE_INDX = 0, CYTOSINE_INDX, GUANINE_INDX,
		THYMINE_INDX
	};

	DNAAnalyzer (const DNAStrandSet &rcDNAStrandSet, const string &rcID);
	DNAAnalyzer (const string &rcFileName, const string &rcID);

	void analyze (const string &rcID);
	void displayStrandSet (ostream &rcOutStream) const;
	void displayProfile (ostream &rcOutStream) const;
	void displayConsensus (ostream &rcOutStream) const;
	DNAStrand getConsensusStrand () const;
	static const char BASE_LOOKUP[];
private:
	DNAStrandSet mcDNAStrandSet;
	vector<vector<unsigned int>> mcProfile;
	DNAStrand mcConsensusStrand;
	void initProfile ();
	void calculateProfile ();
	DNAStrand calculateConsensus (const string &rcID) const;
};