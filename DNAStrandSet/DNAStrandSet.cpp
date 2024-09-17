/**********************************************************************
// File name:	 DNAStrandSet.cpp
// Author:		 Bryant Hayden
// Date:
// Class:			 CS 250
// Assignment: DNAStrandSet
// Purpose:		 Define each prototype for the DNAStrandSet class
**********************************************************************/

#include "DNAStrandSet.h"

/**********************************************************************
Function:    DNAStrandSet

Description: Provide an appropriate default constructor

Parameters:  none

Returned:    none
**********************************************************************/
DNAStrandSet::DNAStrandSet () {
	mcDNAStrands = {};
}
/**********************************************************************
Function:    DNAStrandSet

Description: Provide a constructor that accepts a DNAStrandSet and
						 initializes the calling object correctly.

Parameters:  rcDNAStrandSet

Returned:    none
**********************************************************************/
DNAStrandSet::DNAStrandSet (const DNAStrandSet &rcDNAStrandSet) {
	for (unsigned strand = 0; strand < rcDNAStrandSet.size (); strand++) {
		mcDNAStrands.push_back (rcDNAStrandSet.mcDNAStrands[strand]);
	}
}
/**********************************************************************
Function:    add

Description: Add a DNAStrandto the set. There is no limit to the number
						 of DNAStrand objects. Also, a set cannot contain a strand
						 with the same id

Parameters:  rcDNAStrand

Returned:    none
**********************************************************************/
void DNAStrandSet::add (const DNAStrand &rcDNAStrand) {
	if (!isIn (rcDNAStrand)) {
		mcDNAStrands.push_back (rcDNAStrand);
	}
}
/**********************************************************************
Function:    size

Description: Get the current number of items in the set

Parameters:  none

Returned:    size of set
**********************************************************************/
unsigned int DNAStrandSet::size () const {
	return mcDNAStrands.size ();//this size is a function of any vector
}
/**********************************************************************
Function:    isIn

Description: Determine if a particular DNAStrand is in the set

Parameters:  rcDNAStrand

Returned:    return true or false: if the DNAStrand is contained
**********************************************************************/
bool DNAStrandSet::isIn (const DNAStrand &rcDNAStrand) const {
	bool bContained = false;
	for (unsigned strand = 0; strand < mcDNAStrands.size () && !bContained;
		strand++) {
		if (mcDNAStrands[strand].isEqual (rcDNAStrand)) {
			bContained = true;
		}
	}
	return bContained;
}

/**********************************************************************
Function:    setUnion

Description: Provide the set union operation

Parameters:  rcDNAStrandSet

Returned:    DNAStrandSet
**********************************************************************/
DNAStrandSet DNAStrandSet::setUnion (const DNAStrandSet &rcDNAStrandSet)
const {
	DNAStrandSet cUnion (rcDNAStrandSet);
	for (unsigned strand = 0; strand < mcDNAStrands.size (); strand++) {
		if (!cUnion.isIn (mcDNAStrands[strand])) {
			cUnion.add (mcDNAStrands[strand]);
		}
	}
	return cUnion;
}

/**********************************************************************
Function:    setIntersection

Description: Provide the set intersection operation

Parameters:  rcDNAStrandSet

Returned:    intersection of both sets
**********************************************************************/
DNAStrandSet DNAStrandSet::setIntersection
(const DNAStrandSet &rcDNAStrandSet) const {
	DNAStrandSet cIntersection;
	for (unsigned strand = 0; strand < mcDNAStrands.size (); strand++) {
		if (rcDNAStrandSet.isIn (mcDNAStrands[strand])) {
			cIntersection.add (mcDNAStrands[strand]);
		}
	}
	return cIntersection;
}
/**********************************************************************
Function:    read

Description: Provide a read function to read an unknown number of
						 DNAStrand objects from the stream. Each DNAStrand is
						 separated by whitespace

Parameters:  rcInStream

Returned:    flag indicating any errors
**********************************************************************/
DNAStrand::FLAG DNAStrandSet::read (istream &rcInStream) {
	DNAStrand::FLAG flag = DNAStrand::FLAG::LEGAL;
	DNAStrand inStrand = { "","" };

	flag = inStrand.read (rcInStream);
	while (DNAStrand::FLAG::LEGAL == flag) {
		add (inStrand);
		inStrand = { "","" };
		flag = inStrand.read (rcInStream);
	}
	return flag;
}
/**********************************************************************
Function:    write

Description: Provide a write function to write each DNAStrandobject to
						 the stream, separated by a space (do not print a newline)

Parameters:  rcOutStream

Returned:    none
**********************************************************************/
void DNAStrandSet::write (ostream &rcOutStream) const {
	for (unsigned strand = 0; strand < mcDNAStrands.size (); strand++) {
		mcDNAStrands[strand].write (rcOutStream);
	}
}
/**********************************************************************
Function:    getStrand

Description: get a strand at provided location

Parameters:  which strand

Returned:    DNA Strand is passed back
**********************************************************************/
DNAStrand DNAStrandSet::getStrand (unsigned whichStrand) const {
	return mcDNAStrands[whichStrand];
}
/**********************************************************************
Function:    equalSize

Description: returns true if both strands are equal in size

Parameters:  none

Returned:    bEqual
**********************************************************************/
bool DNAStrandSet::equalSize () const {
	bool bEqual = true;
	for (unsigned strand = 0; strand < mcDNAStrands.size () && !bEqual;
		strand++) {
		if (mcDNAStrands[0].size () != mcDNAStrands[strand].size ()) {
			bEqual = false;
		}
	}
	return bEqual;
}
