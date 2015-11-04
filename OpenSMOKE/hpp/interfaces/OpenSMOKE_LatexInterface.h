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

#ifndef OPENSMOKE_LATEXINTERFACE
#define OPENSMOKE_LATEXINTERFACE

#include "BzzMath.hpp"
using namespace std;

class OpenSMOKE_LatexInterface
{
public:

        OpenSMOKE_LatexInterface();
        OpenSMOKE_LatexInterface(std::string _fileName);
        void include_figure(std::string fileName, std::string caption);
        void include_table(std::string *titlesUp, std::string *titlesLeft, BzzMatrix &matrix);
        void setup();
        void new_section(std::string section_name);
        void add(std::string message);
        void new_page();
        void close();
        void create_pdf();

private:
        std::string fileName;
        ofstream latex_file;
        static const std::string newline;
};

#endif // OPENSMOKE_LATEXINTERFACE


