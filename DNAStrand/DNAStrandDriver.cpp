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