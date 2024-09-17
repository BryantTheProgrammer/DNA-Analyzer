#pragma once
//*********************************************************************
// File name:		DNAStrandSet.h
// Author:			CS, Pacific University
// Date:
// Class:				CS 250
// Assignment:  02 DNAStrandSet Class
// Purpose:			Declare the interface for a set of DNAStrands.
//*********************************************************************
#include "DNAStrand.h"
#include <iostream>
#include <vector>

class DNAStrandSet {
public:
	DNAStrandSet ();
	DNAStrandSet (const DNAStrandSet &rcDNAStrandSet);

	void add (const DNAStrand &rcDNAStrand);
	unsigned int size () const;
	bool isIn (const DNAStrand &rcDNAStrand)const;
	DNAStrandSet setUnion (const DNAStrandSet &rcDNAStrandSet) const;
	DNAStrandSet setIntersection (const DNAStrandSet &rcDNAStrandSet) const;
	DNAStrand::FLAG read (istream &rcInStream);
	void write (ostream &rcOutStream) const;
	DNAStrand getStrand (unsigned whichStrand) const;
	bool equalSize () const;

private:
	vector<DNAStrand> mcDNAStrands;
};