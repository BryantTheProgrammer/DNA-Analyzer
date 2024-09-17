//############################################################################
//DNAStrand.h 
//############################################################################
#pragma once
//*********************************************************************
// File name:		DNAStrand.h
// Author:			CS, Pacific University
// Date:				2/7/22
// Class:				CS 250
// Assignment:  01 DNAStrand Class
// Purpose:			Declare the interface for a strand of DNA. A strand of
//							DNA contains a unique unknown number of chemical
//							bases adenine (A), guanine (G), cytosine (C), and
//							thymine (T). The ID will be a simplified FASTA string.
//*********************************************************************


#include <string>
#include <iostream>

using namespace std;

class DNAStrand {
public:
	enum FLAG : int {
		ILLEGAL_BASE = -2, ILLEGAL_EOF = -1, LEGAL = 0, LEGAL_EOF = 1
	};

	enum BASE : char {
		ADENINE = 'A', CYTOSINE = 'C', GUANINE = 'G',
		THYMINE = 'T', NO_BASE = ' '
	};

	DNAStrand (const string& rcID = "",//Why the [""]'s
		const string& rcDNAStrand = "");
	DNAStrand (const DNAStrand& rcDNAStrand);

	// string getID () const; is removed from the assignment
	char getBase (int whichBase) const;

	FLAG read (istream& rcInStream);
	void write (ostream& rcOutStream) const;
	void writeBases (ostream& rcOutStream) const;

	DNAStrand reverse (const string& rcID) const;
	DNAStrand complement (const string& rcID) const;
	double gcContent () const;
	double hammingDistance (const DNAStrand& rcDNAStrand) const;

	unsigned int size () const;
	bool isEqual (const DNAStrand& rcDNAStrand) const;

private:
	string mID;
	string mBases;

	char dnaBaseComplement (char base) const;
};   
  
 
//############################################################################
//DNAStrand.cpp 
//############################################################################
//**********************************************************************
// File name:	 DNAStrand.cpp
// Author:		 Bryant Hayden
// Date:			 2/8/22-
// Class:			 CS-250-01
// Assignment: DNAStrand Class
// Purpose:		 Implementation of DNAstand.h
//**********************************************************************

#include "DNAStrand.h"

/***********************************************************************
Function:    DNAStrand

Description: Instantiate a DNAStrand using strings passed by reference

Parameters:  rcID - the Id used to construct the strand
						 rcrcDNAStrand - the DNAStrand used to construct the strand

Returned:    none
***********************************************************************/
DNAStrand::DNAStrand (const string& rcID,
	const string& rcDNAStrand)
{
	mID = rcID;
	mBases = rcDNAStrand;
}

/***********************************************************************
Function:    DNAStrand

Description: Instantiate a DNAStrand using preexisting DNAStrand by
						 Reference

Parameters:  rcDNAStrand - used to construct strand

Returned:    none
***********************************************************************/
DNAStrand::DNAStrand (const DNAStrand& rcDNAStrand)
{
	mID = rcDNAStrand.mID;
	mBases = rcDNAStrand.mBases;
}
/***********************************************************************
Function:    getBase

Description: get a base at a particular location

Parameters:  whichBase - passed by value to give location in array

Returned:    char - base at location provided
***********************************************************************/
char DNAStrand::getBase (int whichBase) const {
	return mBases.at(whichBase);
}
/***********************************************************************
Function:    read

Description: read a strand from some input stream

Parameters:  rcInStream - some kind of istream

Returned:    flag - returns length of string. Or error code.
***********************************************************************/
 DNAStrand::FLAG DNAStrand::read (istream& rcInStream)
{
	FLAG flag = ILLEGAL_BASE;
	//char inBase = ' ';
	//mBases = "";//needs to reset string since we are always appending it
	//rcInStream >> mID >> inBase;
	if (rcInStream >> mID >> mBases) {
		for (int base = 0; base < mBases.size (); base++) {
			if (ADENINE == mBases[base] || CYTOSINE == mBases[base] ||
				GUANINE == mBases[base] || THYMINE == mBases[base])
			{
				flag = LEGAL;
			}
			else {
				flag = ILLEGAL_BASE;
				break;
			}
		}
	}
	else if (rcInStream >> mID) {
		flag = ILLEGAL_EOF;
	}
	else {
		flag = LEGAL_EOF;
	}
	return flag;
}
/***********************************************************************
Function:    write

Description: Write the strand to some output stream

Parameters:  rcOutStream - some kind of ostream

Returned:    none
***********************************************************************/
void DNAStrand::write (ostream& rcOutStream) const
{
	rcOutStream << mID << ' ' << mBases << endl;
}
/***********************************************************************
Function:    writeBases

Description: Write the bases to some output stream

Parameters:  rcOutStream - some kind of ostream

Returned:    none
***********************************************************************/
void DNAStrand::writeBases (ostream& rcOutStream) const
{
	rcOutStream << mBases;
}
/***********************************************************************
Function:    reverse

Description: Return a reverse strand through the function name:i.e.
						 create a new strand that is the reverse of the original
						 strand and pass the reversed strand back. An ID is provided
						 as a parameter for the reversed strand that is passed back.

Parameters:  rcID - id for new string

Returned:    DNAStrand
***********************************************************************/
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
/***********************************************************************
Function:    complement

Description: Return a complement of a strand through the function
						 name:i.e. create a new strand that is the complement of the
						 original strand and pass the complemented strand back.An ID
						 is provided as a parameter for the reversed strand that is
						 passed back.

Parameters:  rcID - Id for new DNAStrand

Returned:    DNAStrand
***********************************************************************/
DNAStrand DNAStrand::complement (const string& rcID) const
{
	string complementStrand = mBases;
	for(int base = 0; base < mBases.length(); base++){
		complementStrand.at (base) =
			dnaBaseComplement(complementStrand.at (base));
	}
	DNAStrand cComplement(rcID, complementStrand);
	return cComplement;
}
/***********************************************************************
Function:    gcContent

Description: Calculate and return the GC content of a particular strand.
						 The GC content is the percentage of G's and C's in the
						 overall strand.

Parameters:  none

Returned:    gcContent
***********************************************************************/
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
/***********************************************************************
Function:    hammingDistance

Description: Calculate and return the hamming distance of two strands:
						 The hamming distance between two DNA strands(s1, s2) of
						 equal length is the proportion of corresponding bases that
						 differ between the two strings.

Parameters:  rcDNAStrand

Returned:    hammingDistance
***********************************************************************/
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
/***********************************************************************
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
/***********************************************************************
Function:    isEqual

Description: Return whether the IDs of two DNAStrands are the same.
						 Return the result through the function name. The function
						 accepts a single DNAStrand.

Parameters:  rcDNAStrand

Returned:    bISEqual
***********************************************************************/
bool DNAStrand::isEqual (const DNAStrand& rcDNAStrand) const
{
	return rcDNAStrand.mID == mID;
}
/***********************************************************************
Function:    dnaBaseComplement

Description: Return a complement of a strand through the function name:
						 i.e. create a new strand that is the complement of the
						 original strand and pass the completed strand back. An ID
						 is provided as a parameter for the reveres strand that is
						 passed back

Parameters:  base - base passed in to determine complement

Returned:    base
***********************************************************************/
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
   
  
 
//############################################################################
//DNAStrandDriver.cpp 
//############################################################################
//**********************************************************************
// File name:	 DNAStrandDriver.cpp
// Author:		 Bryant Hayden
// Date:			 2/8/22-
// Class:			 CS-250-01
// Assignment: DNAStrand Class
// Purpose:		 driver to test DNA code
//**********************************************************************
#include "DNAStrand.h"
#include <iostream>
#include <iomanip>
#include <fstream>


using namespace std;

/***********************************************************************
Function:		  main

Description:	main driver used to test a DNAStrand object. This driver
							will accept a DNAStrand in FASTA format from the
							keyboard and then read in two strings from a file. The
							reverse complement, GC content, and hamming distance of
							two strands will be displayed.

Parameters:		none

Returned:			exit status
***********************************************************************/

int main () {
	const string INPUT_FILE_NAME = "dnastrands.txt";
	const int DECIMAL_PRECISION = 2;
	const int INSUFFICIENT_STRAND = 0;//might change
	const int MAX_STRAND_SIZE = 3;
	//Length are also being used to error check
	int flag1 = 0;//lengths are unccaery. Thats why we have size funtion
	int flag2 = 0;//implement flags not length for what is returned
	int flag3 = 0;
	string ID;
	string Bases;

	ifstream inFile;
	DNAStrand cS1 (">S1", "");
	DNAStrand cS2 ("", "");
	DNAStrand cS3 ("", "");

	cout << "*** DNA Strand Driver ***" << endl << endl;

	cout << "Enter DNA Strand (length = " << MAX_STRAND_SIZE << "): ";

	flag1 = cS1.read (cin);

	if (flag1 != DNAStrand::LEGAL || cS1.size() != MAX_STRAND_SIZE) {
		cout << endl << "Illegal Keyboard Input Format" << endl;
		exit (EXIT_FAILURE);
	}

		inFile.open (INPUT_FILE_NAME);
		if (inFile.fail ()) {
			cout << "ERROR Opening Input File";
			exit (EXIT_FAILURE);
		}
	//read lines 2&3 and save the returned value
	flag2 = cS2.read (inFile);
	flag3 = cS3.read (inFile);

	inFile.close ();

	if (flag2 != DNAStrand::FLAG::LEGAL ||
		(flag3 != DNAStrand::FLAG::LEGAL && flag3
			!= DNAStrand::FLAG::LEGAL_EOF)|| cS2.size() != MAX_STRAND_SIZE
		|| cS3.size() != MAX_STRAND_SIZE) {
		cout << endl << "Illegal Data File Format" << endl;
		inFile.close ();
		exit (EXIT_FAILURE);
	}
	cout << endl << "Original DNA Strands" << endl;
	cS1.write (cout);
	cS2.write (cout);
	cS3.write (cout);
	cout << endl;

	cout << "Reverse Complement of S1" << endl;
	cS1.reverse ("").complement(">S1").write (cout);

	cout << endl << "GC Content of S1" << endl;
	cout << fixed << setprecision(DECIMAL_PRECISION) << cS1.gcContent ()
		<< endl;

	cout << endl << "Hamming Distance (S1, S2)" << endl;
	cout << cS1.hammingDistance (cS2) << endl;

	return EXIT_SUCCESS;
}   
  
 
