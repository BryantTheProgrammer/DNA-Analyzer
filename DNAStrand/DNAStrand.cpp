//*********************************************************************
// File name:	 DNAStrand.cpp
// Author:		 Bryant Hayden
// Date:			 2/8/22-
// Class:			 CS-250-01
// Assignment: DNAStrand Class
// Purpose:		 Implementation of DNAstand.h
//*********************************************************************

#include "DNAStrand.h"

/**********************************************************************
Function:    DNAStrand

Description: Instantiate a DNAStrand using strings passed by reference

Parameters:  rcID - the Id used to construct the strand
						 rcrcDNAStrand - the DNAStrand used to construct the strand

Returned:    none
**********************************************************************/
DNAStrand::DNAStrand (const string& rcID,
	const string& rcDNAStrand)
{
	mID = rcID;
	mBases = rcDNAStrand;
}

/**********************************************************************
Function:    DNAStrand

Description: Instantiate a DNAStrand using preexisting DNAStrand by
						 Reference

Parameters:  rcDNAStrand - used to construct strand

Returned:    none
**********************************************************************/
DNAStrand::DNAStrand (const DNAStrand &rcDNAStrand) {
	mID = rcDNAStrand.mID;
	mBases = rcDNAStrand.mBases;
}
/**********************************************************************
Function:    getBase

Description: get a base at a particular location

Parameters:  whichBase - passed by value to give location in array

Returned:    char - base at location provided
**********************************************************************/
char DNAStrand::getBase (int whichBase) const {
	return mBases.at(whichBase);
}
/**********************************************************************
Function:    read

Description: read a strand from some input stream

Parameters:  rcInStream - some kind of istream

Returned:    flag - returns length of string. Or error code.
**********************************************************************/
 DNAStrand::FLAG DNAStrand::read (istream& rcInStream)
 {
	FLAG flag = LEGAL;
	//char inBase = ' ';
	//mBases = "";//needs to reset string since we are always appending it
	//rcInStream >> mID >> inBase;
	if (rcInStream >> mID >> mBases) {
		for (int base = 0; base < static_cast<int>(mBases.size ());
			base++) {
			if (ADENINE != mBases[base] && CYTOSINE != mBases[base] &&
				GUANINE != mBases[base] && THYMINE != mBases[base])
			{
				flag = ILLEGAL_BASE;
				break;
			}
		}
	}
	else if (mID.size() > NULL) {
		flag = ILLEGAL_EOF;
	}
	else {
		flag = LEGAL_EOF;
	}
	return flag;
}
/**********************************************************************
Function:    write

Description: Write the strand to some output stream

Parameters:  rcOutStream - some kind of ostream

Returned:    none
**********************************************************************/
void DNAStrand::write (ostream& rcOutStream) const
{
	rcOutStream << mID << ' ' << mBases << endl;
}
/**********************************************************************
Function:    writeBases

Description: Write the bases to some output stream

Parameters:  rcOutStream - some kind of ostream

Returned:    none
**********************************************************************/
void DNAStrand::writeBases (ostream& rcOutStream) const
{
	rcOutStream << mBases;
}
/**********************************************************************
Function:    reverse

Description: Return a reverse strand through the function name:i.e.
						 create a new strand that is the reverse of the original
						 strand and pass the reversed strand back. An ID is provided
						 as a parameter for the reversed strand that is passed back.

Parameters:  rcID - id for new string

Returned:    DNAStrand
**********************************************************************/
DNAStrand DNAStrand::reverse (const string& rcID) const
{
	DNAStrand cReverseStrand(rcID,mBases);
	int length = cReverseStrand.size();
	char temp;
	int mirrorBase;
	for (int base = 0; base < (length)/2; base++)
	{
		//create a index for the mirror of the middle of the strand
		mirrorBase = length - base-1;
		//preserve first char in temp
		temp = cReverseStrand.mBases.at (base);
		//set first equal to last
		cReverseStrand.mBases.at (base) =
			cReverseStrand.mBases.at (mirrorBase);
		//set last to persevered temp
		cReverseStrand.mBases.at (mirrorBase) = temp;
	}
	return cReverseStrand;
}
/**********************************************************************
Function:    complement

Description: Return a complement of a strand through the function
						 name:i.e. create a new strand that is the complement of the
						 original strand and pass the complemented strand back.An ID
						 is provided as a parameter for the reversed strand that is
						 passed back.

Parameters:  rcID - Id for new DNAStrand

Returned:    DNAStrand
**********************************************************************/
DNAStrand DNAStrand::complement (const string& rcID) const
{
	string complementStrand = mBases;
	for(int base = 0; base < static_cast<int>(mBases.length()); base++){
		complementStrand.at (base) =
			dnaBaseComplement(complementStrand.at (base));
	}
	DNAStrand cComplement(rcID, complementStrand);
	return cComplement;
}
/**********************************************************************
Function:    gcContent

Description: Calculate and return the GC content of a particular strand.
						 The GC content is the percentage of G's and C's in the
						 overall strand.

Parameters:  none

Returned:    gcContent
**********************************************************************/
double DNAStrand::gcContent () const
{
	const char C = 'C';
	const char G = 'G';
	double totalGAndC = 0.0;//to maintain double division
	int numberOfBases = static_cast<int>(mBases.size ());
	double gcContent = 0.0;
	for (int base = 0; base < numberOfBases; base++)
	{
		if (C == mBases.at(base) || G == mBases.at(base))
		{
			totalGAndC += 1;
		}
	}
	gcContent = totalGAndC / numberOfBases;
	return gcContent;
}
/**********************************************************************
Function:    hammingDistance

Description: Calculate and return the hamming distance of two strands:
						 The hamming distance between two DNA strands(s1, s2) of
						 equal length is the proportion of corresponding bases that
						 differ between the two strings.

Parameters:  rcDNAStrand

Returned:    hammingDistance
**********************************************************************/
double DNAStrand::hammingDistance (const DNAStrand& rcDNAStrand) const
{
	const double HAMMING_LENGTH_ERROR = -1.0;
	double hammingDistance;
	int numberOfBases = static_cast<int>(mBases.size ());
	double differenceCount = 0;

	for (int base = 0; base < numberOfBases; base++)
	{
		if (mBases.at (base) != rcDNAStrand.mBases.at (base)) {
			differenceCount++;
		}
	}
	hammingDistance = differenceCount / rcDNAStrand.mBases.size();

	if (mBases.size() != rcDNAStrand.size()) {
		hammingDistance = HAMMING_LENGTH_ERROR;
	}
	return hammingDistance;
}
/**********************************************************************
Function:    size

Description: Return the size of a strand(number of bases) through the
						 function name. The size is the number of bases in a
						 DNAStrand

Parameters:  none

Returned:    static_cast<int>(mBases.length ())
***********************************************************************/
unsigned int DNAStrand::size () const
{
	return static_cast<int>(mBases.length ());
}
/**********************************************************************
Function:    isEqual

Description: Return whether the IDs of two DNAStrands are the same.
						 Return the result through the function name. The function
						 accepts a single DNAStrand.

Parameters:  rcDNAStrand

Returned:    bISEqual
**********************************************************************/
bool DNAStrand::isEqual (const DNAStrand& rcDNAStrand) const
{
	return rcDNAStrand.mID == mID;
}
/**********************************************************************
Function:    dnaBaseComplement

Description: Return a complement of a strand through the function name:
						 i.e. create a new strand that is the complement of the
						 original strand and pass the completed strand back. An ID
						 is provided as a parameter for the reveres strand that is
						 passed back

Parameters:  base - base passed in to determine complement

Returned:    base
**********************************************************************/
char DNAStrand::dnaBaseComplement (char base) const
{

	switch (base)
	{
	case THYMINE:  base = ADENINE;
								 break;
	case ADENINE:  base = THYMINE;
								 break;
	case GUANINE:  base = CYTOSINE;
								 break;
	case CYTOSINE: base = GUANINE;
								 break;

	default:			 base = ILLEGAL_BASE;
	}

	return base;
}