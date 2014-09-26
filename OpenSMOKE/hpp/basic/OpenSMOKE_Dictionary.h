/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
 *   alberto.cuoci@polimi.it   						                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef OPENSMOKE_DICTIONARY
#define OPENSMOKE_DICTIONARY

#include <vector>
#include "BzzMath.hpp"

class OpenSMOKE_1DMap
{
public:
    vector<string>  names;
};

class OpenSMOKE_2DMap
{
public:
	vector<double>  values;
    vector<string>  names;
};

class OpenSMOKE_Dictionary
{
public:

	OpenSMOKE_Dictionary();
    void SetupBase();
	void SetName(const std::string name);
    void Add(const std::string new_word, const char new_flag, const char new_type_id, const std::string new_comment);
    void Conflict(const std::string word_a, const std::string word_b);
    void Compulsory(const std::string word_a, const std::string word_b);
    void Compulsory(const std::string word_a, const std::string word_b, const std::string word_c);
    void Compulsory(const std::string word_a, const std::string word_b, const std::string word_c, const std::string word_d);
    void Compulsory(const std::string word_a, const std::string word_b, const std::string word_c, 
					const std::string word_d, const std::string word_e);
    void Compulsory(const std::string word_a, const std::string word_b, const std::string word_c, 
					const std::string word_d, const std::string word_e, const std::string word_f);
    void Lock();

    bool Return(const std::string new_word);
    bool Return(const std::string new_word, int  &value);
	bool Return(const std::string new_word, char &value);
    bool Return(const std::string new_word, double &value);
    bool Return(const std::string new_word, std::string &value);
    bool Return(const std::string new_word, double &number, std::string &unit);
	bool Return(const std::string new_word, vector<string> &names);
    bool Return(const std::string new_word, vector<double> &values, vector<string> &names);

    void ParseFile(const std::string fileName);
    void Help();

private:

	std::string name_object;

    bool iLocked;

    vector<string>  words;
    vector<char>    flag;
    vector<char>    type_id;
    vector<string>  comment;

    vector<string>  conflict_word_a;
    vector<string>  conflict_word_b;

    vector<string>  compulsory_word_a;
    vector<string>  compulsory_word_b;
    vector<string>  compulsory_word_c;
    vector<string>  compulsory_word_d;
    vector<string>  compulsory_word_e;
    vector<string>  compulsory_word_f;

    vector<string>  read_words;
    vector<int>     read_int;
    vector<char>    read_char;
    vector<double>  read_double;
    vector<string>  read_string;
    vector<char>    read_type_id;
    vector<OpenSMOKE_1DMap>  read_1DMap;
    vector<OpenSMOKE_2DMap>  read_2DMap;


    void ReadAdd(const std::string new_word);
    void ReadAdd(const std::string new_word, const int n);
	void ReadAdd(const std::string new_word, const char c);
    void ReadAdd(const std::string new_word, const double n);
    void ReadAdd(const std::string new_word, const std::string n);
    void ReadAdd(const std::string new_word, const double d, const std::string n);
	void ReadAdd(const std::string new_word, OpenSMOKE_1DMap _1DMap);
	void ReadAdd(const std::string new_word, OpenSMOKE_2DMap _2DMap);

    void CheckAdd(const std::string new_word);
    char Check(const std::string new_word);
    void CheckCompulsory(const vector<string> read_words);
    void CheckCompulsoryCouples(const vector<string> read_words);
    void CheckAll(const vector<string> read_words);
    void CheckConflicting(const vector<string> read_words);
    void CheckDouble(const vector<string> read_words);
    bool CheckForComment(const std::string read_word);
	bool CheckForKeyWord(const std::string read_word);

    void ErrorMessage(const std::string message);
    void WarningMessage(const std::string message);
};

#endif // OPENSMOKE_DICTIONARY


