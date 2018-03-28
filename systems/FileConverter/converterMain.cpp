/*
 * converterMain.cpp
 *
 *  Created on: 16 Feb 2018
 *      Author: eabdelha
 */

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <set>
#include <stdlib.h>
#include <string>

using namespace std;

vector<string>* readLine(std::ifstream& file, int ignoreFirstColumn, int typeFrom, string delimiter) {

	std::string str;

	if(std::getline(file, str)) {

		//parse line
		vector<string>* v = new vector<string>();

		size_t pos = 0;
		std::string token;
		bool firstToken = true;
		while ((pos = str.find(delimiter)) != std::string::npos) {

			//parse consecutive tokens
		    token = str.substr(0, pos);
		    str.erase(0, pos + delimiter.length());

		    //ignore the first token
			if(firstToken && ignoreFirstColumn==1) {
				firstToken = false;
				continue;
			}
			firstToken = false;

			v->push_back(token);
		}
		v->push_back(str);

		return v;
	} else {
		return NULL;
	}
}

vector<string>* readLineSparse(std::ifstream& file, int ignoreFirstColumn, int typeFrom, unsigned int classIndexFrom) {

	std::string str;
	std::string delimiter = " ";

	if(std::getline(file, str)) {
		//parse line
		vector<string> temp;

		size_t pos = 0;
		std::string token;
		bool firstToken = true;
		while ((pos = str.find(delimiter)) != std::string::npos) {

			//parse consecutive tokens
		    token = str.substr(0, pos);
		    str.erase(0, pos + delimiter.length());

		    //ignore the first token
   			if(firstToken && ignoreFirstColumn==1) {
   				firstToken = false;
   				continue;
   			}
   			firstToken = false;

		    temp.push_back(token);
		}
		temp.push_back(str);

		vector<string>* v = new vector<string>();
		for(unsigned int i=0;i<temp.size();i++) {
			string token = temp[i];

			//get index and value
			int pos = token.find(":");
			unsigned int index = atoi(token.substr(0, pos).c_str());
			string value = token.substr(pos+1);
			while(v->size()<=index)
				v->push_back("0");
			(*v)[index] = value;
		}

		return v;
	} else {
		return NULL;
	}
}

vector<string>* readLineBin4(std::ifstream& file, int ignoreFirstColumn, int typeFrom, int nDim) {

	vector<string>* v = new vector<string>();

	bool firstCol = true;
	for(int i=0;i<nDim;i++) {
		float value = 0;
		file.read((char*)(&value), 4);
//		cout<<"size = "<<sizeof(float)<<endl;
//		cout<<value<<endl;
		if(file.eof())
			return NULL;
		if(firstCol && ignoreFirstColumn) {
			firstCol = false;
			continue;
		}

		std::ostringstream ss;
		ss << value;
		std::string s(ss.str());

		v->push_back(s);
	}

	if(v->size()==0) {
		delete v;
		return NULL;
	}

	return v;
}

vector<string>* readLineBin8(std::ifstream& file, int ignoreFirstColumn, int typeFrom, int nDim) {

	vector<string>* v = new vector<string>();

	bool firstCol = true;
	for(int i=0;i<nDim;i++) {
		double value;
		file.read((char*)(&value), 8);
		if(file.eof())
			return NULL;
		if(firstCol && ignoreFirstColumn) {
			firstCol = false;
			continue;
		}

		std::ostringstream ss;
		ss << value;
		std::string s(ss.str());

		v->push_back(s);
	}

	if(v->size()==0) {
		delete v;
		return NULL;
	}

	return v;
}

int main( int argc, char *argv[])
{
	/**
	 * File formats
	 * 0- csv file, each line contains values for each object with commas separated
	 * 1- csv-like file with the first line containing the number of objects then the number of dimensions. Followed by a series of lines, each represents an object
	 * 2- sparse list of items, each line has a class, then a series of index:value pairs. 0 values are ignored.
	 * 3- ARFF file format
	 * 4- Binary format; a header of 2 integers: 4-byte for #objects, the second 4-bytes is #coordinates. Followed by a series of 4-bytes float values that represent objects values.
	 * 5- Binary format; a header of 2 integers: 8-byte for #objects, the second 8-bytes is #coordinates. Followed by a series of 8-bytes float values that represent objects values.
	 */
	if (argc > 7) {
		char* inputFile = argv[1];
		char* outputFile = argv[2];
		int ignoreFirstColumn = 0;
		if(atoi(argv[3])==1)
			ignoreFirstColumn = 1;
		int typeFrom = atoi(argv[4]);
		int typeTo = atoi(argv[5]);

		int classIndexFrom;
		classIndexFrom = atoi(argv[6]);
		int classIndexTo;
		classIndexTo = atoi(argv[7]);

		int numObjects = 0;
		int numDimensions = 0;

		set<string> classValues;

		cout << "Input file: " << inputFile << std::endl;
		cout << "Output file: " << outputFile << std::endl;
		cout << "Ignore first column: "<<ignoreFirstColumn<<endl;
		cout << "Conversion from type = "<<typeFrom<<endl;
		cout << "Conversion to type = "<<typeTo<<endl;
		cout << "Class Index From = "<<classIndexFrom<<endl;
		cout << "Class Index To = "<<classIndexTo<<endl;

		//open input file
		std::ifstream file;

		if(typeFrom==4 || typeFrom==5) {
			file.open(inputFile, ios::binary | ios::in);
		}
		else
		{
			file.open(inputFile);
		}

		std::string str;
		//open output file
		ofstream tempFile;
		tempFile.open ("temp");

		//read header data
		switch(typeFrom) {
		case 0:
			break;
		case 1:
			{
			file>>numObjects;
			file>>numDimensions;
			string tempStr;
			std::getline(file, tempStr);
			}
			break;
		case 2:
			break;
		case 3:
			{
			std::string str;
			while(std::getline(file, str)) {
				if(str.compare("@DATA")==0) {
					break;
				}
			}
			}
			break;
		case 4:
		case 5:
			{
				file.read((char*)(&numObjects), 4);
				file.read((char*)(&numDimensions), 4);
				cout<<"binary read: "<<numObjects<<", "<<numDimensions<<endl<<flush;
			}
			break;
		}

		numObjects = 0;

		while(true) {
			vector<string>* v;
			string cls;
			if(typeFrom==2) {
				v = readLineSparse(file, ignoreFirstColumn, typeFrom, classIndexFrom);
			}
			else {
				if(typeFrom==0 || typeFrom==3)
					v = readLine(file, ignoreFirstColumn, typeFrom, ",");
				if(typeFrom==1)
					v = readLine(file, ignoreFirstColumn, typeFrom, " ");
				if(typeFrom==4) {
					v = readLineBin4(file, ignoreFirstColumn, typeFrom, numDimensions);
				}
				if(typeFrom==5) {
					v = readLineBin8(file, ignoreFirstColumn, typeFrom, numDimensions);
				}
			}
			if(v==NULL) break;

			if(typeFrom==2) {
				if(numDimensions<v->size())
					numDimensions=v->size();
			}

			if(numObjects==0)
				numDimensions = v->size();

			if(classIndexFrom>-1) {
				if(classIndexFrom>=v->size())
					cls = "0";
				else
					cls = (*v)[classIndexFrom];
				classValues.insert(cls);
			}

			for(unsigned int i=0;i<v->size();i++) {
				tempFile<<(*v)[i];
				if(i==(v->size()-1))
					tempFile<<"\n";
				else
					tempFile<<",";
			}

			delete v;
			numObjects++;
		}

		tempFile.close();
		file.close();

		cout<<"numObjects = "<<numObjects<<endl;
		cout<<"numDimensions = "<<numDimensions<<endl;

		//now read the temp file and output on the target file
		std::ifstream tempFile1("temp");

		//open the traget file
		ofstream outFile;
		outFile.open (outputFile);

		//write header file
		switch(typeTo) {
			case 0:
				break;
			case 1:
				outFile<<numObjects<<" "<<numDimensions<<endl;
				break;
			case 2:
				break;
			case 3:
				/*
				@RELATION name-rel1

				@ATTRIBUTE the NUMERIC
				...
				@ATTRIBUTE pos {P,N}

				@DATA
				 */
				outFile<<"@RELATION name-rel1"<<endl;
				outFile<<endl;
				for(int i=0;i<numDimensions;i++) {
					if(i==classIndexTo) {
						outFile<<"@ATTRIBUTE cls {";

						int t = 0;
						for(set<string>::iterator iter = classValues.begin(); iter!=classValues.end(); iter++) {
							if(t==0)
								t = 1;
							else
								outFile<<",";
							outFile<<(*iter);
						}

						outFile<<"}"<<endl;
					}
					else
						outFile<<"@ATTRIBUTE A"<<i<<" NUMERIC"<<endl;
				}
				outFile<<endl;
				outFile<<"@DATA"<<endl;
				break;
			case 4:
			case 5:
				outFile.write((char*)(&numObjects), 4);
				outFile.write((char*)(&numDimensions), 4);
				break;
		}

		while(true) {
			vector<string>* v = readLine(tempFile1, false, 0, ",");

			if(v==NULL) break;

			//if condition wether we do care about class
			if(classIndexFrom>-1) {
				string classValue = (*v)[classIndexFrom];
				(*v).erase(v->begin()+classIndexFrom);
				(*v).insert(v->begin()+classIndexTo, classValue);
			}

			bool valueWasWritten = false;
			unsigned int i = 0;
			for(;i<v->size();i++) {
				switch(typeTo) {
				case 0:
				case 3:
					outFile<<(*v)[i];
					if(i<v->size()-1)
						outFile<<",";
					break;
				case 1:
					outFile<<(*v)[i];
					if(i<v->size()-1)
						outFile<<" ";
					break;
				case 2:
					if(i==classIndexTo) {
						if(valueWasWritten) {
							outFile<<" ";
						}
						outFile<<(*v)[i];
					} else {
						if(((*v)[i]).compare("0")==0)
							continue;
						if(valueWasWritten) {
							outFile<<" ";
						}
						outFile<<(i-1)<<":"<<(*v)[i];
					}
					valueWasWritten = true;
					break;
				case 4:
					{
						float tempF = atof(((*v)[i]).c_str());
						outFile.write((char*)(&tempF), 4);
//						cout<<(*v)[i]<<":"<<tempF<<" | ";
					}
					break;
				case 5:
					double tempD = atof(((*v)[i]).c_str());
					outFile.write((char*)(&tempD), 8);
					break;
				}
			}
			if(typeTo==0 || typeTo==1 || typeTo==2 || typeTo==3)
				outFile<<"\n";
			delete v;
		}

		/*
		while(true) {
			vector<double>* v = readLine(file, ignoreFirstColumn, typeFrom);
			if(v==NULL) break;

			switch(type) {

			case 1:
				allValues.push_back(v);
				break;
			case 2:
				outFile<<(*v)[0];
				for(int i=1;i<v->size();i++) {
					if((*v)[i]==0)
						continue;
					outFile<<" "<<(i-1)<<":"<<(*v)[i];
				}
				outFile<<"\n";
				delete v;
				break;
			}
		}



		switch(type) {

		case 1:
			if(allValues.size()>0) {
				outFile<<allValues.size()<<" "<<allValues[0]->size()<<endl;
			} else {
				outFile<<"0 0"<<endl;
			}

			for(unsigned int i=0;i<allValues.size();i++)
			{
				vector<double>* currentV = allValues[i];
				for(unsigned int j=0;j<currentV->size();j++) {
					outFile<<(*currentV)[j];
					if(j<currentV->size()-1)
						outFile<<" ";
				}

				delete currentV;

				if(i<(allValues.size()-1))
					outFile<<endl;
			}
			break;
		}


*/

		cout<<"Converted file written on: "<<outputFile<<endl;
		outFile.close();
		tempFile.close();
	}
	else
	{
		cout<<"Error: You have to provide:\n1- input file\n2- output file\n3- a flag (1 to ignore first column and 0 otherwise)"<<endl;
		cout<<"4- from conversion type\n5- to conversion type\n";
		cout<<"\t0 -> csv\n";
		cout<<"\t1 -> csv preceded with a line containing #objects [space] #dimensions\n";
		cout<<"\t2 -> sparse list of items (feature number:feature value [space]) features with zero value are ignored\n";
		cout<<"\t3 -> ARFF file format\n";
		cout<<"\t4 -> Binary 4-bytes\n";
		cout<<"\t5 -> Binary 8-bytes\n";
		cout<<"6- Clas Index From\n7- Class Index To\n";
	}

	return 0;
}
