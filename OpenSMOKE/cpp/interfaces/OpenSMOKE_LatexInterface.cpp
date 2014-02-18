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

#include "interfaces/OpenSMOKE_LatexInterface.h"

const string OpenSMOKE_LatexInterface::newline = "\n";

OpenSMOKE_LatexInterface::OpenSMOKE_LatexInterface()
{
    fileName = "report";
    string dummy = fileName + ".tex";
    latex_file.open(dummy.c_str());
}

OpenSMOKE_LatexInterface::OpenSMOKE_LatexInterface(string _fileName)
{
    fileName = _fileName;
    string dummy = fileName + ".tex";
    latex_file.open(dummy.c_str());
}

void OpenSMOKE_LatexInterface::include_figure(string fileName, string caption)
{
    string message;

    message  = newline;
    message += "   \\begin{figure}"                         + newline;
    message += "      \\centering"                          + newline;
    message += "      \\includegraphics{" + fileName + "}"  + newline;
    message += "      \\caption{" + caption + "}"           + newline;
    message += "      \\label{fig:" + fileName + "}"        + newline;
    message += "   \\end{figure}"                           + newline;
    message += newline;

    latex_file << message;
}

void OpenSMOKE_LatexInterface::setup()
{
    string head_file;

    head_file  = "\\documentclass[a4paper,12pt]{article}" + newline;
    head_file += "\\usepackage{graphicx}" + newline+newline;
    head_file += "\\begin{document}" + newline + newline + newline;

    latex_file << head_file;
}

void OpenSMOKE_LatexInterface::new_section(string section_name)
{
    string section;

    section  = newline;
    section += "   \\section{" + section_name + "}" + newline;
    section += "   \\label{" + section_name + "}" + newline;
    section += newline;

    latex_file << section;
}

void OpenSMOKE_LatexInterface::add(string message)
{
    latex_file << newline;
    latex_file << message;
    latex_file << newline;
}

void OpenSMOKE_LatexInterface::new_page()
{
    latex_file << newline;
    latex_file << "   \\newpage";
    latex_file << newline;
}

void OpenSMOKE_LatexInterface::close()
{
    latex_file << newline;
    latex_file << "\\end{document}";
    latex_file << newline;
    latex_file.close();
}

void OpenSMOKE_LatexInterface::create_pdf()
{
    string latex_creation   = "latex "  + fileName;
    string pdf_creation     = "dvipdf " + fileName + ".dvi " + fileName + ".pdf";

    system(latex_creation.c_str());
    system(pdf_creation.c_str());
}


void OpenSMOKE_LatexInterface::include_table(string *titlesUp, string *titlesLeft, BzzMatrix &matrix)
{
    int i, j;

    string message;
    for(j=1;j<=matrix.Columns();j++)
        message += "c";

    latex_file << "   \\begin{tabular}{|c|" + message + "|}" << endl;
    latex_file << "      \\hline" << endl;

    latex_file << "  ";
    for(j=1;j<=matrix.Columns();j++)
        latex_file << " & " << titlesUp[j];
    latex_file << "\\\\" << endl;

    latex_file << "\\hline" << endl;

    for(i=1;i<=matrix.Rows();i++)
    {
        latex_file << titlesLeft[i];
        for(j=1;j<=matrix.Columns();j++)
            latex_file << " & " << matrix[i][j] ;
        latex_file << "\\\\" << endl;
    }

    latex_file << "      \\hline" << endl;
    latex_file << "   \\end{tabular}" << endl;
}
