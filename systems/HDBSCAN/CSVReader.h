/*
 * CSVReader.h
 *
 *  Created on: 27 Sep 2017
 *      Author: eabdelha
 */

#ifndef CSVREADER_H_
#define CSVREADER_H_

#include<vector>
#include<string>
#include "Object.h"

class CVReader {
	public:
		static std::vector<Object*> read(string filename, bool ignoreFirstLine);
};

class CSVRow
{
    public:
        float const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const;
        void readNextRow(std::istream& str);

    private:
        std::vector<float> m_data;
};

#endif /* CSVREADER_H_ */
