/**********************************************************************
// File name:	 DNASetDriver.cpp
// Author:		 Bryant Hayden
// Date:
// Class:			 CS 250
// Assignment: DNAStrandSet
// Purpose:		 Driver for assignment
**********************************************************************/
#include "DNAStrandSet.h"
#include <iostream>
#include <iomanip>
#include <fstream>

void printMessage (const string& message, ostream& rcOutStream);

using namespace std;

int main () {
	const string INPUT1 = "dnastrands1.txt";
	const string INPUT2 = "dnastrands2.txt";

	ifstream inFile;

	DNAStrandSet cSet1;
	DNAStrandSet cSet2;

	DNAStrand::FLAG flag1;
	DNAStrand::FLAG flag2;

	//read 1
	inFile.open (INPUT1);

	if (inFile.fail ()) {
		cout << "ERROR OPENING INPUT FILE";
		exit (EXIT_FAILURE);
	}

	flag1 = cSet1.read (inFile);

	if (DNAStrand::FLAG::LEGAL > flag1) {
		cout << "Illegal DNAStrand Input" << endl;
		exit (EXIT_FAILURE);
	}
	inFile.close ();
	//end

	//Read 2
	inFile.open (INPUT2);

	if (inFile.fail ()) {
		cout << "FILE ERROR";
		exit (EXIT_FAILURE);
	}

	flag2 = cSet2.read (inFile);

	if (DNAStrand::FLAG::LEGAL > flag2) {
		cout << "Illegal DNAStrand Input" << endl;
		exit (EXIT_FAILURE);
	}
	inFile.close ();
	//end

	cout << "*** DNAStrandSet Analyzer ***" << endl;

	printMessage ("DNAStrandSet #1", cout);
	cSet1.write (cout);

	printMessage ("DNAStrandSet #2", cout);
	cSet2.write (cout);

	printMessage ("Set Union", cout);
	cSet1.setUnion (cSet2).write (cout);

	printMessage ("Set Intersection", cout);
	cSet1.setIntersection (cSet2).write (cout);

	return EXIT_SUCCESS;
}

void printMessage (const string& message, ostream& rcOutStream) {
	rcOutStream << endl << message << endl << setfill ('-')
		<< setw (message.length ()) << "" << endl;
}
