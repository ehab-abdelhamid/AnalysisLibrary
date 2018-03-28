/*
 * CSVReader.cpp
 *
 *  Created on: 27 Sep 2017
 *      Author: eabdelha
 */

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "CSVReader.h"

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

std::vector<Object*> CVReader::read(string filename, bool ignoreFirstLine)
{
	cout<<"Reading data ...";
	bool ignore = false;
	if(ignoreFirstLine)
		ignore = true;

	std::ifstream file(filename.c_str());

	vector<Object*> allObjects;
	CSVRow row;
	while(file >> row)
	{
		if(ignore) {
			ignore = false;
			continue;
		}
		Object* o = new Object();
		for(unsigned int i=0;i<row.size();i++) {
			o->addValue((float)row[i]);
		}
		//o->print(cout);
		allObjects.push_back(o);
	}

	cout<<" [#objects = "<<allObjects.size()<<"] ";
	cout<<"DONE."<<endl;
	return allObjects;
}

std::size_t CSVRow::size() const
{
	return m_data.size();
}
void CSVRow::readNextRow(std::istream& str)
{
	std::string         line;
	std::getline(str, line);

	std::stringstream   lineStream(line);
	std::string         cell;

	m_data.clear();
	while(std::getline(lineStream, cell, ','))
	{
		float cellF = atof(cell.c_str());
		m_data.push_back(cellF);
	}
	// This checks for a trailing comma with no data after it.
	if (!lineStream && cell.empty())
	{
		// If there was a trailing comma then add an empty element.
		m_data.push_back(-1);
	}
}
