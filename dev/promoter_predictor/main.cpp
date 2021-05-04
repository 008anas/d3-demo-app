#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/bind.hpp>
#include <ctime>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/assign.hpp>
#include <algorithm>



using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace boost::adaptors;


typedef vector < vector <double>* > table;
typedef map< string, table* > scoringMatrices;
typedef map< string, vector <double>* > mapOfScores;
typedef map< string, double > mapMat;

void loadScoringMatrices(scoringMatrices* sc, int numArgs, char* const argv[], vector <double>* maxvec){
	for(int i=4;i<numArgs;i++){
		ifstream mFile;
		mFile.open(argv[i]);
		string line, key="scoringMatrixNo"+lexical_cast<string>(i-3);
		
		table* tururu=new table();
		double maxs = 1.0;
		while(getline(mFile, line)){
			vector<string> sValues;
			split(sValues,line,is_any_of("\t "));
			if(sValues[0]!=">"){
				vector<double>* values=new vector<double>();
			
				BOOST_FOREACH(string tmp,sValues){
					values->push_back(lexical_cast<double>(tmp));
				}
				tururu->push_back(values);
				maxs *= *max_element(values->begin(), values->end());
			}
		}
		sc->insert(make_pair(key,tururu));
		maxvec->push_back(maxs);
	}
	
}


void loadEnergyMatrix(mapMat* emat, char* const argv[], int num){
	ifstream mFile;
	mFile.open(argv[num]);
	string line;
	string pairBases;
	double energyVal;
	while(getline(mFile, line)){
	//	cout << line << endl;
		vector<string> sValues;
		split(sValues,line,is_any_of("\t "));
		pairBases=sValues[0];
		energyVal=lexical_cast<double>(sValues[1]);
		emat->insert(make_pair(pairBases,energyVal));
		//cout << pairBases << endl;
	}
}


vector <double>* scanGenomeMotifs(scoringMatrices* sc, string keyVal, string sequence, int startval, double normValue){
	double sizeMotif = sc->at(keyVal)->size();
	vector <double>* scoreMotif = new vector <double>();
	for(int i=25-startval;i<(sequence.size()-sizeMotif);i++){
		double tmpScore=1;
		for(int j=0;j<sizeMotif;j++){
			char nucleotide = lexical_cast<char>(sequence.substr(i+j,1));
			switch(nucleotide){
				case 'A':
					tmpScore*=sc->at(keyVal)->at(j)->at(0);
					break;
				case 'C':
					tmpScore*=sc->at(keyVal)->at(j)->at(1);
					break;
				case 'G':
					tmpScore*=sc->at(keyVal)->at(j)->at(2);
					break;
				case 'T':
					tmpScore*=sc->at(keyVal)->at(j)->at(3);
					break;
			}
		}
		scoreMotif->push_back(-log2(tmpScore/normValue));
	}
	return scoreMotif;
}


vector <double>* scanDiNucleotides(string sequence, string dinucleotide, int posMin, int posMax){ //posMin and posMax relative to the tata start (i.e. -4,-1 for TG motif)
	vector <double>* scoreMotif = new vector <double>();
	int lenDin=dinucleotide.size();
	for(int x=24;x<(sequence.size()-14);x++){
		double tmpScore=0;
		for(int j=posMin;j<=posMax;j++){
			int A=x+j;
			if(sequence.substr(A, lenDin).compare(dinucleotide)==0){
				tmpScore=1;
				break;
			}
		}
		scoreMotif->push_back(tmpScore);
	}
	return scoreMotif;
}



vector <double>* scanGenomeEnergies(mapMat* em, string sequence, double GCcontent){
	vector <double>* scoreEnergy = new vector <double>();
	for(int i=25;i<(sequence.size()-15);i++){
		double tmpScore=1;
		for(int j=-25;j<15;j++){   // we calculate energy for 15nts, taken 2 by 2
			string dinucleotide = sequence.substr(i+j,2);
			tmpScore*=em->at(dinucleotide); 
		}
		scoreEnergy->push_back(-log2(tmpScore));
	}
	return scoreEnergy;
}

double countGCs(string sequence){
	double counts=0;
	for(int i=0;i<sequence.size();i++){
		if(sequence.substr(i,1).compare("C")==0 || sequence.substr(i,1).compare("G")==0){
			counts++;
		}
	}
	return counts;
}


vector <double>* scanSensitivePositions(string sequence, double GCcontent){
	vector <double>* sensitivePosBases = new vector <double>();
	for(int i=25;i<(sequence.size()-15);i++){
		double percentage=1;
		double numGCsAfterTATA;
		numGCsAfterTATA=countGCs(sequence.substr(i-5,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i-4,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i-3,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i+5,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i+6,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i+7,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		numGCsAfterTATA=countGCs(sequence.substr(i+8,2));
		percentage*=pow(0.2,numGCsAfterTATA);
		sensitivePosBases->push_back(-log2(percentage));
	}
	return sensitivePosBases;
}


vector <double>* scanminus45bases(string sequence, double GCcontent){
	vector <double>* minus45bases = new vector <double>();
	for(int h=25;h<=39;h++){
		double minus451=0;
		minus45bases->push_back(-log2(1-(minus451/15)));
	}
	for(int i=40;i<(sequence.size()-15);i++){
		int j=i-40;
		double minus452=countGCs(sequence.substr(j,15));
		if(minus452==15){
			minus452=14;
		}
		minus45bases->push_back(-log2(1-(minus452/15)));
	}
	
	return minus45bases;
}


string makeReverseStr(string sequence){
	string revComp=sequence;
	//cout << revComp.substr(0,5) << endl;
	reverse(revComp.begin(), revComp.end());
	for(int i=0;i<revComp.size();i++){
		char nucleotide = lexical_cast<char>(revComp.substr(i,1));
		switch(nucleotide){
			case 'A':
				revComp.replace(i,1,"T");
				break;
			case 'C':
				revComp.replace(i,1,"G");
				break;
			case 'G':
				revComp.replace(i,1,"C");
				break;
			case 'T':
				revComp.replace(i,1,"A");
				break;
		}
	}
	//cout << revComp.substr(0,5) << endl;
	return revComp;
}


int getPribnowStartVal(char* const argv[], int i){
	ifstream getPribnowStart;
	getPribnowStart.open(argv[i]);
	string linea;
	vector <string> splitted;
	int val;
	getline(getPribnowStart,linea);
	split(splitted, linea, is_any_of(" "));
	val=lexical_cast<int>(splitted[1]);
	getPribnowStart.close();
	return val; 
}

void normalizeData(vector <double>* dataset){
	double suma=0;
	for(int i=0;i<dataset->size();i++){
		suma+=dataset->at(i);
	}
	suma=suma/(dataset->size());
	for(int i=0;i<dataset->size();i++){
		dataset->at(i)=(dataset->at(i)/suma);
	}
}


int main (int argc, char * const argv[]) {
	if(argc<5){
		cout << "Calculates promoter scores according to different scoring matrices provided by the user" << endl;
		cout << "ARGS" << endl;
		cout << "1. Sequence file in FASTA format (only one sequence per file)" << endl;
		cout << "2. Energy Scoring Matrix for DNA stacking energy" << endl;
		cout << "3. Output file basename (no .txt)" << endl;
		cout << "4+. Scoring matrix/matrices (can give more than one), max motif length = 30 nt" << endl;
		cout << endl;
	}
	
	ifstream fastaSeq;
	string line,key;
	
	
	string TG="TG";
	string GC="GC";
	int pm=-4, pM=-1;
	
	fastaSeq.open(argv[1]);
	while(getline(fastaSeq, line)){
		if(line.substr(0,1).compare(">")!=0){
			cout << "Processing plus strand.. "<< endl;
			mapOfScores* resultToWrite = new mapOfScores();
	//		map <string, int> pribnowStart;
			
			double GCcontent = 100*countGCs(line)/line.size();
			
			mapMat* energies=new mapMat();
			loadEnergyMatrix(energies, argv, 2);
			
			scoringMatrices* scores=new scoringMatrices();
			vector <double>* maxvec = new vector <double>();  //this vector loads the max possible values for each motif so that we can normalize
			loadScoringMatrices(scores, argc, argv, maxvec);
			
			for(int i=4;i<argc;i++){
				int counter=0;
				int val = getPribnowStartVal(argv, i);
	//			pribnowStart[matID]=val;
				double maxval=lexical_cast<double>(maxvec->at(counter));
				string matID="scoringMatrixNo"+lexical_cast<string>(i-3);
				vector <double>* scoreForThisMotif = scanGenomeMotifs(scores, matID, line, val, maxval);
				resultToWrite->insert(make_pair(matID, scoreForThisMotif));
				cout << "Scoring matrix "<< i-3 << " processed" << endl;
				counter++;
			}
			// Now we calculate energies
			string energyID="Energy";
			string basePenaltiesID="Penalties";
			string beforeTG="TGx";
			string afterGC="xGC";
			string motif45ATs="ATs_-45";
			
			vector <double>* energyValuesForSeqs = scanGenomeEnergies(energies, line, GCcontent);
			vector <double>* basePenalties = scanSensitivePositions(line, GCcontent);
			vector <double>* TGprevious = scanDiNucleotides(line, TG, pm, pM);
			vector <double>* GCafter=scanDiNucleotides(line, GC, 10,13);
			vector <double>* motif45=scanminus45bases(line, GCcontent);
			
			normalizeData(energyValuesForSeqs);
			normalizeData(basePenalties);
			normalizeData(motif45);
			
			resultToWrite->insert(make_pair(energyID, energyValuesForSeqs));
			resultToWrite->insert(make_pair(basePenaltiesID, basePenalties));
			resultToWrite->insert(make_pair(beforeTG, TGprevious));
			resultToWrite->insert(make_pair(afterGC, GCafter));
			resultToWrite->insert(make_pair(motif45ATs, motif45));
			
			cout << "Energy matrix processed" <<endl;
			cout << "Base penalties applied" << endl;
			
			
			//write output file from start to end - 30 last nucleotides (max length of motifs = 30)
			string fileStreamPlus=lexical_cast<string>(argv[3])+"_plus.txt";
			ofstream outputFile;
			outputFile.open(fileStreamPlus.c_str());
			
			for(int j=0; j<(line.size()-50);j++){
				outputFile << j+25 << "\t"; 
				BOOST_FOREACH(string key, *resultToWrite | map_keys){
					outputFile << resultToWrite->at(key)->at(j) << "\t";

				}
				outputFile << "+" << endl; 
			}
			outputFile.close();
			
			mapOfScores* resultToWriteMinus = new mapOfScores();
			cout << "Processing minus strand.. "<< endl;
			
			string line2 = makeReverseStr(line);
			
			
			for(int i=4;i<argc;i++){
				int counter=0;
				int val = getPribnowStartVal(argv, i);
				double maxval=lexical_cast<double>(maxvec->at(counter));
				string matID="scoringMatrixNo"+lexical_cast<string>(i-3);
				vector <double>* scoreForThisMotif = scanGenomeMotifs(scores, matID, line2,val, maxval);
				resultToWriteMinus->insert(make_pair(matID, scoreForThisMotif));
				cout << "Scoring matrix "<< i-3 << " processed" << endl;
				counter++;
			}
			
			//string energyID2="Energy";
			//string basePenaltiesID2="Penalties";
			vector <double>* energyValuesForSeqs2 = scanGenomeEnergies(energies, line2, GCcontent);
			vector <double>* basePenalties2 = scanSensitivePositions(line2, GCcontent);
			vector <double>* TGprevious2 = scanDiNucleotides(line2, TG, pm, pM);
			vector <double>* GCafter2=scanDiNucleotides(line2, GC, 10,13);
			vector <double>* motif452=scanminus45bases(line2, GCcontent);
			
			normalizeData(energyValuesForSeqs2);
			normalizeData(basePenalties2);
			normalizeData(motif452);
			
			resultToWriteMinus->insert(make_pair(energyID, energyValuesForSeqs2));
			resultToWriteMinus->insert(make_pair(basePenaltiesID, basePenalties2));
			resultToWriteMinus->insert(make_pair(beforeTG, TGprevious2));
			resultToWriteMinus->insert(make_pair(afterGC, GCafter2));
			resultToWriteMinus->insert(make_pair(motif45ATs, motif452));
			
			cout << "Energy matrix processed" <<endl;
			cout << "Base penalties applied" << endl;
			
			//write output file from start to end - 30 last nucleotides (max length of motifs = 30) 
			string fileStreamMinus=lexical_cast<string>(argv[3])+"_minus.txt";
			ofstream outputFile2;
			outputFile2.open(fileStreamMinus.c_str());
			
			for(int j=0; j<(line2.size()-50);j++){
				outputFile2 << line2.size()-24-j << "\t"; 
				BOOST_FOREACH(string key, *resultToWriteMinus | map_keys){
					outputFile2 << resultToWriteMinus->at(key)->at(j) << "\t";
					
				}
				outputFile2 << "-" << endl; 
			}
			outputFile2.close();
			
			
		}
	}
	fastaSeq.close();  
	
    return 0;
}


