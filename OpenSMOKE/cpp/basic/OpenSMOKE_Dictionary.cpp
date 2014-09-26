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

#include "basic/OpenSMOKE_Utilities.h"

#include "basic/OpenSMOKE_Dictionary.h"

void OpenSMOKE_Dictionary::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Dictionary"		<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Dictionary::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_Dictionary"		<< endl;
    cout << "Object:  " << name_object			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

OpenSMOKE_Dictionary::OpenSMOKE_Dictionary()
{
	name_object	= "[Name not assigned]";
}

void OpenSMOKE_Dictionary::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Dictionary::Lock()
{
    iLocked = true;
}

void OpenSMOKE_Dictionary::SetupBase()
{
    iLocked = false;

    Add("#NoVerbose",       'O', 'N', "");
    Add("#Nvideo",          'O', 'I', "Numer of steps for video output");
    Add("#Nfile",           'O', 'I', "Number of steps for file output");
    Add("#Nbackup",         'O', 'I', "Number of steps for backup");
    Add("#Help",            'O', 'N', "This help");
}

void OpenSMOKE_Dictionary::Add(const std::string new_word, const char new_flag, const char new_type_id, const std::string new_comment)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    CheckAdd(new_word);

    words.push_back(new_word);
    flag.push_back(new_flag);
    type_id.push_back(new_type_id);
    comment.push_back(new_comment);
}

void OpenSMOKE_Dictionary::Conflict(const std::string word_a, const std::string word_b)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    conflict_word_a.push_back(word_a);
    conflict_word_b.push_back(word_b);
}

void OpenSMOKE_Dictionary::Compulsory(const std::string word_a, const std::string word_b)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    compulsory_word_a.push_back(word_a);
    compulsory_word_b.push_back(word_b);
	compulsory_word_c.push_back("!");
    compulsory_word_d.push_back("!");
	compulsory_word_e.push_back("!");
    compulsory_word_f.push_back("!");
}

void OpenSMOKE_Dictionary::Compulsory(const std::string word_a, const std::string word_b, const std::string word_c)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    Check(word_c);
    compulsory_word_a.push_back(word_a);
    compulsory_word_b.push_back(word_b);
	compulsory_word_c.push_back(word_c);
    compulsory_word_d.push_back("!");
	compulsory_word_e.push_back("!");
    compulsory_word_f.push_back("!");
}

void OpenSMOKE_Dictionary::Compulsory(const std::string word_a, const std::string word_b, const std::string word_c, const std::string word_d)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    Check(word_c);
    Check(word_d);
    compulsory_word_a.push_back(word_a);
    compulsory_word_b.push_back(word_b);
	compulsory_word_c.push_back(word_c);
    compulsory_word_d.push_back(word_d);
	compulsory_word_e.push_back("!");
    compulsory_word_f.push_back("!");
}

void OpenSMOKE_Dictionary::Compulsory(const std::string word_a, const std::string word_b, const std::string word_c,
									  const std::string word_d, const std::string word_e)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    Check(word_c);
    Check(word_d);
    Check(word_e);
    compulsory_word_a.push_back(word_a);
    compulsory_word_b.push_back(word_b);
	compulsory_word_c.push_back(word_c);
    compulsory_word_d.push_back(word_d);
	compulsory_word_e.push_back(word_e);
    compulsory_word_f.push_back("!");
}

void OpenSMOKE_Dictionary::Compulsory(const std::string word_a, const std::string word_b, const std::string word_c,
									  const std::string word_d, const std::string word_e, const std::string word_f)
{
    if (iLocked == true)
        ErrorMessage("The dictionary is already locked");

    Check(word_a);
    Check(word_b);
    Check(word_c);
    Check(word_d);
    Check(word_e);
    Check(word_f);
    compulsory_word_a.push_back(word_a);
    compulsory_word_b.push_back(word_b);
	compulsory_word_c.push_back(word_c);
    compulsory_word_d.push_back(word_d);
	compulsory_word_e.push_back(word_e);
    compulsory_word_f.push_back(word_f);
}

char OpenSMOKE_Dictionary::Check(const std::string new_word)
{
    int n = words.size();
    for(int i=0; i<n; i++)
        if (words[i] == new_word)
            return type_id[i];
	
    ErrorMessage("The following option is not recognized: " + new_word);
    return 0;
}

void OpenSMOKE_Dictionary::CheckAdd(const std::string new_word)
{
    int n = words.size();
    for(int i=0; i<n; i++)
        if (words[i] == new_word)
            ErrorMessage("The following option is already in the dictionary: " + new_word);
}

void OpenSMOKE_Dictionary::CheckCompulsory(const vector<string> read_words)
{
    int nDictionary = words.size();
    int nRead       = read_words.size();

    for(int i=0; i<nDictionary; i++)
    {
        if(flag[i] == 'C')
        {
            int iFind = 0;
            for(int j=0; j<nRead; j++)
                if (words[i] == read_words[j])
                {
                    iFind = 1;
                    break;
                }

            if (iFind == 0)
                ErrorMessage("The following option is mandatory: " + words[i]);
        }
    }
}

void OpenSMOKE_Dictionary::CheckAll(const vector<string> read_words)
{
    int nDictionary = words.size();
    int nRead       = read_words.size();

    for(int j=0; j<nRead; j++)
    {
        int iFind = 0;
        for(int i=0; i<nDictionary; i++)
            if (words[i] == read_words[j])
            {
                iFind = 1;
                break;
            }

        if (iFind == 0)
			ErrorMessage("The following option is not recognized: " + read_words[j]);
    }
}

void OpenSMOKE_Dictionary::CheckConflicting(const vector<string> read_words)
{
    int nConflicting    = conflict_word_a.size();
    int nRead           = read_words.size();

    for(int j=0; j<nRead; j++)
        for(int i=0; i<nConflicting; i++)
        {
            if (conflict_word_a[i] == read_words[j])
                for(int jj=0; jj<nRead; jj++)
                {
                    if (conflict_word_b[i] == read_words[jj])
                        ErrorMessage("Conflicts between the following options: "
                                        + read_words[j] + " " + read_words[jj]);
                }
        }
}

void OpenSMOKE_Dictionary::CheckCompulsoryCouples(const vector<string> read_words)
{
    int nCompulsory     = compulsory_word_a.size();
    int nRead           = read_words.size();

	for(int i=0; i<nCompulsory; i++)
	{
		bool iFound = false;
		
		for(int j=0; j<nRead; j++)
            if (	compulsory_word_a[i] == read_words[j] || 
					compulsory_word_b[i] == read_words[j] ||
					compulsory_word_c[i] == read_words[j] ||
					compulsory_word_d[i] == read_words[j] ||
					compulsory_word_e[i] == read_words[j] ||
					compulsory_word_f[i] == read_words[j]	)	iFound = true;

		if (iFound == false)
		{
			if (compulsory_word_c[i] == "!")
				ErrorMessage("The following options are compulsory: " + compulsory_word_a[i] + " || " + compulsory_word_b[i]);
			else if (compulsory_word_d[i] == "!")
				ErrorMessage("The following options are compulsory: " + compulsory_word_a[i] + " || " + compulsory_word_b[i] + " || " + compulsory_word_c[i]);
			else if (compulsory_word_e[i] == "!")
				ErrorMessage("The following options are compulsory: " + compulsory_word_a[i] + " || " + compulsory_word_b[i] + " || " + compulsory_word_c[i] + " || " + compulsory_word_d[i]);
			else if (compulsory_word_f[i] == "!")
				ErrorMessage("The following options are compulsory: " + compulsory_word_a[i] + " || " + compulsory_word_b[i] + " || " + compulsory_word_c[i] + " || " + compulsory_word_d[i] + " || " + compulsory_word_e[i]);
			else
				ErrorMessage("The following options are compulsory: " + compulsory_word_a[i] + " || " + compulsory_word_b[i] + " || " + compulsory_word_c[i] + " || " + compulsory_word_d[i] + " || " + compulsory_word_e[i] + " || " + compulsory_word_f[i]);
		}
	}
}

void OpenSMOKE_Dictionary::CheckDouble(const vector<string> read_words)
{
    int nRead           = read_words.size();

    for(int j=0; j<nRead; j++)
        for(int i=0; i<nRead; i++)
            if(read_words[i] == read_words[j] && i!=j)
                ErrorMessage("The following options is used two times: " + read_words[j]);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word)
{
    read_words.push_back(new_word);
   	read_char.push_back(' ');
	read_int.push_back(0);
    read_double.push_back(0.);
    read_string.push_back("");
    read_type_id.push_back('N');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const int n)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(n);
    read_double.push_back(0.);
    read_string.push_back("");
    read_type_id.push_back('I');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const char c)
{
    read_words.push_back(new_word);
    read_char.push_back(c);
    read_int.push_back(0);
    read_double.push_back(0.);
    read_string.push_back("");
    read_type_id.push_back('H');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const double n)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(0);
    read_double.push_back(n);
    read_string.push_back("");
    read_type_id.push_back('D');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const std::string n)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(0);
    read_double.push_back(0.);
    read_string.push_back(n);
    read_type_id.push_back('S');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const double d, const std::string n)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(0);
    read_double.push_back(d);
    read_string.push_back(n);
    read_type_id.push_back('M');

	OpenSMOKE_1DMap _1DMap;
	OpenSMOKE_2DMap _2DMap;
	read_1DMap.push_back(_1DMap);
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, OpenSMOKE_1DMap _1DMap)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(0);
    read_double.push_back(0.);
    read_string.push_back("");
    read_type_id.push_back('V');
    read_1DMap.push_back(_1DMap);

	OpenSMOKE_2DMap _2DMap;
	read_2DMap.push_back(_2DMap);
}

void OpenSMOKE_Dictionary::ReadAdd(const std::string new_word, const OpenSMOKE_2DMap _2DMap)
{
    read_words.push_back(new_word);
    read_char.push_back(' ');
    read_int.push_back(0);
    read_double.push_back(0.);
    read_string.push_back("");
    read_type_id.push_back('L');
    read_2DMap.push_back(_2DMap);

	OpenSMOKE_1DMap _1DMap;
	read_1DMap.push_back(_1DMap);
}

bool OpenSMOKE_Dictionary::CheckForComment(const std::string read_word)
{
    if(read_word.at(0) == '/')
        if(read_word.at(1) == '/')
            return true;

    return false;
}

bool OpenSMOKE_Dictionary::CheckForKeyWord(const std::string read_word)
{
    if(read_word.at(0) == '#')
		return true;

    return false;
}


void OpenSMOKE_Dictionary::ParseFile(const std::string fileName)
{
    const int SIZE = 300;
    char comment[SIZE];

    char    flag;
    std::string  dummy;
    int     dummy_int;
    char    dummy_char;
    double  dummy_double;
    std::string  dummy_string;
	bool    iGoBack;

    ifstream fInput;
    openInputFileAndControl(fInput, fileName.c_str());
	
	iGoBack = false;
    do
    {
		if (iGoBack == false)
			fInput >> dummy;
		else dummy = dummy_string;

		iGoBack = false;

        if (dummy == "#Help")
            Help();

        if (dummy == "#END")
            break;
        else
        {
            if (!CheckForComment(dummy))
            {
                flag = Check(dummy);

                if(flag == 'N')
                {
                    ReadAdd(dummy);
                }

                if(flag == 'I')
                {
                    fInput >> dummy_int;
                    ReadAdd(dummy, dummy_int);
                }

                if(flag == 'H')
                {
                    fInput >> dummy_char;
                    ReadAdd(dummy, dummy_char);
                }

                if(flag == 'D')
                {
                    fInput >> dummy_double;
                    ReadAdd(dummy, dummy_double);
                }

                if(flag == 'S')
                {
                    fInput >> dummy_string;
                    ReadAdd(dummy, dummy_string);
                }

                if(flag == 'M')
                {
                    fInput >> dummy_double;
                    fInput >> dummy_string;
                    ReadAdd(dummy, dummy_double, dummy_string);
                }

				if(flag == 'L')
                {
					OpenSMOKE_2DMap _2DMap;

					for(;;)
					{
						fInput >> dummy_string;
						
						if (CheckForComment(dummy_string) == true)
						{
							fInput.getline(comment, SIZE);
							break;
						}

						if (CheckForKeyWord(dummy_string) == true)
						{
							
							iGoBack = true;
							break;
						}

						fInput >> dummy_double;
						_2DMap.names.push_back(dummy_string);
						_2DMap.values.push_back(dummy_double);
                    }
					
					ReadAdd(dummy, _2DMap);
                }

				if(flag == 'V')
                {
					OpenSMOKE_1DMap _1DMap;

					for(;;)
					{
						fInput >> dummy_string;
						
						if (CheckForComment(dummy_string) == true)
						{
							fInput.getline(comment, SIZE);
							break;
						}

						if (CheckForKeyWord(dummy_string) == true)
						{
							
							iGoBack = true;
							break;
						}

						_1DMap.names.push_back(dummy_string);
                    }
					 
					ReadAdd(dummy, _1DMap);
                }
            }
            else
                fInput.getline(comment, SIZE);
        }
    }
    while(!fInput.eof());

    fInput.close();

    CheckCompulsory(read_words);
	CheckCompulsoryCouples(read_words);
    CheckAll(read_words);
    CheckConflicting(read_words);
    CheckDouble(read_words);
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'N')
                ErrorMessage("The " + new_word + " is not a 'N' option!");
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, int &value)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'I')
                ErrorMessage("The " + new_word + " is not a 'I' option!");
            value = read_int[i];
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, char &value)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'H')
                ErrorMessage("The " + new_word + " is not a 'H' option!");
            value = read_char[i];
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, double &value)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'D')
                ErrorMessage("The " + new_word + " is not a 'D' option!");
            value = read_double[i];
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, std::string &value)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'S')
                ErrorMessage("The " + new_word + " is not a 'S' option!");
            value = read_string[i];
            return 1;
        }
    return 0;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, double &number, std::string &unit)
{
    int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'M')
                ErrorMessage("The " + new_word + " is not a 'M' option!");
            unit = read_string[i];
            number = read_double[i];
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, vector<double> &values, vector<string> &names)
{
	int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'L')
                ErrorMessage("The " + new_word + " is not a 'L' option!");
            values = read_2DMap[i].values;
            names  = read_2DMap[i].names;
            return true;
        }
    return false;
}

bool OpenSMOKE_Dictionary::Return(const std::string new_word, vector<string> &names)
{
	int n = read_words.size();
    for(int i=0; i<n; i++)
        if (read_words[i] == new_word)
        {
            if (read_type_id[i] != 'V')
                ErrorMessage("The " + new_word + " is not a 'V' option!");
            names  = read_1DMap[i].names;
            return true;
        }
    return false;
}

void OpenSMOKE_Dictionary::Help()
{
	int i;
    int n = words.size();

    cout << endl;
    cout << "OpenSMOKE Dictionary" << endl;
    cout << "------------------------------------------------------------------" << endl;

    cout << "Mandatory arguments:" << endl;
    for(i=0; i<n; i++)
        if(flag[i] == 'C')
        {
            if (type_id[i] == 'N')
                cout << "   " << words[i] << " \t\t\t" << comment[i] << endl;
            else if (type_id[i] == 'I')
                cout << "   " << words[i] << " [int]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'H')
                cout << "   " << words[i] << " [char]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'D')
                cout << "   " << words[i] << " [double]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'S')
                cout << "   " << words[i] << " [std::string]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'M')
                cout << "   " << words[i] << " [double] [std::string]\t" << comment[i] << endl;
        }
    cout << endl;

    cout << "Optional arguments:" << endl;
    for(i=0; i<n; i++)
        if(flag[i] == 'O')
        {
            if (type_id[i] == 'N')
                cout << "   " << words[i] << " \t\t\t" << comment[i] << endl;
            else if (type_id[i] == 'I')
                cout << "   " << words[i] << " [int]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'H')
                cout << "   " << words[i] << " [char]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'D')
                cout << "   " << words[i] << " [double]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'S')
                cout << "   " << words[i] << " [std::string]\t\t" << comment[i] << endl;
            else if (type_id[i] == 'M')
                cout << "   " << words[i] << " [double] [std::string]\t" << comment[i] << endl;
        }
    cout << endl;

    cout << "Press enter to continue..." << endl;
    getchar();
    exit(-1);
}
