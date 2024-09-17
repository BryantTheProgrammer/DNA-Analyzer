#include "DNAAnalyzer.h"
#include <iomanip>

//we will be using 0,1,2,3 to represent A,C,G,T respectively
///IS THERE A WAY TO ACCESS THESE WITHOUT CREATING A LIBRARY
const char DNAAnalyzer::BASE_LOOKUP[] = { 'A','C','G','T' };
/**********************************************************************
Function:    DNAAnalyzer

Description: sets mcDNAStrandSet appropriately and calls analyze
						 passing in a strand id for the consensus strand

Parameters:  rcDNAStrandSet - set to be analyzed
						 rcID           - ID passed through function to be
															assigned in calculateConsensus

Returned:    none
**********************************************************************/
DNAAnalyzer::DNAAnalyzer (const DNAStrandSet &rcDNAStrandSet,
	const string &rcID) {
	mcDNAStrandSet = rcDNAStrandSet;
	analyze (rcID);
}
/**********************************************************************
Function:    DNAAnalyzer

Description: populates a DNAStrandSet by opening and reading the file
						 from rcFileName then call analyze as in 6.

Parameters:  rcFileName - name of file that contains a strand set
						 rcID       - ID passed through function to be
													assigned in calculateConsensus

Returned:    none
**********************************************************************/
DNAAnalyzer::DNAAnalyzer (const string &rcFileName, const string &rcID)
 {//{ on new line not done by accident/ran out of space
	ifstream inFile;

	inFile.open (rcFileName);
	if (inFile.fail ()) {
		cout << "Error: Cannot Open File: " << rcFileName << endl;
		exit (EXIT_FAILURE);
	}
	if (static_cast<int>(mcDNAStrandSet.read (inFile)) <
		DNAStrand::FLAG::LEGAL) {
		cout << "Error: ILLegal File Read: " << rcFileName << endl;
		exit (EXIT_FAILURE);
	}
	inFile.close ();
	analyze (rcID);
}
/**********************************************************************
Function:    analyze

Description: needs to make sure all DNAStrands are of equal length, or
						 an error is outputted to the display and the program is
						 terminated. If the strands are of equal length, then call
						 initProfile (), calculateProfile (), and
						 calculateConsensus () which assigns the mcConsensusStrand
						 to the consensus strand

Parameters:  rcID - ID passed through function to be assigned in
						 calculateConsensus

Returned:    none
**********************************************************************/
void DNAAnalyzer::analyze (const string &rcID) {
	//Check that the file is set is populated
	if (mcDNAStrandSet.size () <= 0) {
		cout << "No DNA Strands To Process" << endl;
		exit (EXIT_FAILURE);
	}
	//Size of 1st strand to be checked against all other strands
	const int STRAND_SIZE = mcDNAStrandSet.getStrand (0).size ();
	//iterate through strands(rows)
	for (int strand = 0;
		strand < static_cast<int>(mcDNAStrandSet.size ()); strand++) {
		//check that each strand matches in size
		if (STRAND_SIZE != mcDNAStrandSet.getStrand (strand).size ()) {
			cout << "Error: Non-Homologous Strands" << endl;
			exit (EXIT_FAILURE);
		}
	}
	initProfile ();
	calculateProfile ();
	mcConsensusStrand = calculateConsensus (rcID);
}
/**********************************************************************
Function:    displayStrandSet

Description: display the strand set

Parameters:  rcOutStream - stream to output strand set

Returned:    none
**********************************************************************/
void DNAAnalyzer::displayStrandSet (ostream &rcOutStream) const {
	mcDNAStrandSet.write (rcOutStream);
	rcOutStream << endl;
}
/**********************************************************************
Function:    displayProfile

Description: Display the profile

Parameters:  rcOutStream - oStream used to write Profile

Returned:    none
**********************************************************************/
void DNAAnalyzer::displayProfile (ostream &rcOutStream) const {
	for (unsigned int base = 0; base < mcProfile.size (); base++) {
		rcOutStream << BASE_LOOKUP[base];//start each row with base
		for (unsigned int location = 0;
			location < mcProfile.at (base).size (); location++) {
			rcOutStream << setw (3) << mcProfile.at (base).at (location);
		}
		rcOutStream << endl;
	}
	rcOutStream << endl;
}
/**********************************************************************
Function:    displayConsensus

Description: display the consensus strand

Parameters:  rcOutStream - Stream in which we write consensus strand

Returned:    none
**********************************************************************/
void DNAAnalyzer::displayConsensus (ostream &rcOutStream) const {
	mcConsensusStrand.write (cout);
}
/**********************************************************************
Function:    getConsensusStrand

Description: single line function, returns mcConsensusStrand

Parameters:  none

Returned:    mcConsensusStrand
**********************************************************************/
DNAStrand DNAAnalyzer::getConsensusStrand () const {
	return mcConsensusStrand;
}
/**********************************************************************
Function:    initProfile

Description: initialize matrix profile of 0's to display the
						 analysis

Parameters:  none

Returned:    none
**********************************************************************/
void DNAAnalyzer::initProfile () {
	//Vector of 0's
	//#of 0's equal to DNAStandsetSize of subzero
	vector <unsigned>templateVector;
	const int STRAND_SIZE = mcDNAStrandSet.getStrand (0).size ();
	for (int size = 0; size < STRAND_SIZE; size++) {
		templateVector.push_back (0);
	}
	for (int base = static_cast<int>(Bases::ADENINE_INDX);
		base <= static_cast<int>(Bases::THYMINE_INDX); base++) {
		mcProfile.push_back (templateVector);
	}
}
/**********************************************************************
Function:    calculateProfile

Description: Iterate through the strandSet by row by column. For each
						 base in the strandSet increment the appropriate counter in
						 mcProfile.

Parameters:  none

Returned:    none
**********************************************************************/
void DNAAnalyzer::calculateProfile () {
	const int BASE_COUNT = mcDNAStrandSet.getStrand (0).size ();
	const int SET_SIZE = mcDNAStrandSet.size ();
	Bases baseToIncrement = Bases::ADENINE_INDX;
	//iterate through strands of strandSet
	for (int strand = 0; strand < SET_SIZE; strand++) {
		//iterate through bases of a strand in strandSet
		for (int base = 0; base < BASE_COUNT; base++) {
			//Switch statement to take a base(char) and place it into an
			//baseIndex(int)
			switch (mcDNAStrandSet.getStrand (strand).getBase (base)) {
			case DNAStrand::BASE::ADENINE:
				baseToIncrement = Bases::ADENINE_INDX;
				break;
			case DNAStrand::BASE::CYTOSINE:
				baseToIncrement = Bases::CYTOSINE_INDX;
				break;
			case DNAStrand::BASE::GUANINE:
				baseToIncrement = Bases::GUANINE_INDX;
				break;
			case DNAStrand::BASE::THYMINE:
				baseToIncrement = Bases::THYMINE_INDX;
				break;
			default:cout << "calulateProfile Base Error";
				break;
			}
			//Increment count in BASE(row) at location(column)
			mcProfile.at (static_cast<int>(baseToIncrement)).at (base)++;
		}
	}
}
/**********************************************************************
Function:    calculateConsensus

Description: creates the consensus strand using the completed profile.
						 It is possible that the consensus strand is not unique.
						 For instance, in the first column of the above profile, if
						 instead of 3, 0, 1, 0 we had 2, 0, 2, 0,then either A or G
						 could be used in the consensus strand.I’m just looking for
						 a correct consensus strand so whichever base you use for a
						 tie is OK

Parameters:  rcID - the Id for new strand that we create

Returned:    the calculated consensus strand
**********************************************************************/
DNAStrand DNAAnalyzer::calculateConsensus (const string &rcID) const {
	string stringConsensus = "";
	int leadingBase = static_cast<int>(Bases::ADENINE_INDX);
	//iterate through the columns or location
	for (unsigned location = 0; location < mcProfile.at (0).size ();
		location++) {
		//iterate through each base count in a column to figure out the
		//most commonly appearing
		for (int base = static_cast<int>(Bases::ADENINE_INDX);
			base <= static_cast<int>(Bases::THYMINE_INDX);  base++) {
			if (mcProfile.at (base).at (location) >
				mcProfile.at (leadingBase).at (location)) {
				leadingBase = base;
			}
		}
		//After searching an entire column add the base to our consensus
		stringConsensus.push_back (BASE_LOOKUP[leadingBase]);
	}
	//Create a DNA strand with the rcID given and our calculated bases
	DNAStrand consensus (rcID, stringConsensus);
	return consensus;
}