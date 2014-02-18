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

#ifndef OPENSMOKE_GNUPLOTINTERFACE
#define OPENSMOKE_GNUPLOTINTERFACE

#include "BzzMath.hpp"

class OpenSMOKE_GnuPlotInterface
{
public:

    OpenSMOKE_GnuPlotInterface();
    OpenSMOKE_GnuPlotInterface(string folderName);

    void setPlot(string _title, string _xlabel, string _ylabel);

    void setKey(string _label_1);
    void setKey(string _label_1, string _label_2);
    void setKey(string _label_1, string _label_2, string _label_3);
    void setKey(string _label_1, string _label_2, string _label_3, string _label_4);
    void setKey(string _label_1, string _label_2, string _label_3,
                string _label_4, string _label_5);
    void setKey(string _label_1, string _label_2, string _label_3,
                string _label_4, string _label_5, string _label_6);

    void setKind( char _kind_1);
    void setKind( char _kind_1, char _kind_2);
    void setKind( char _kind_1, char _kind_2, char _kind_3);
    void setKind( char _kind_1, char _kind_2, char _kind_3, char _kind_4);
    void setKind( char _kind_1, char _kind_2, char _kind_3,
                  char _kind_4, char _kind_5);
    void setKind( char _kind_1, char _kind_2, char _kind_3,
                  char _kind_4, char _kind_5, char _kind_6);

    void setLogScaleX();
    void setLogScaleY();

    void plot(string fileFigure, string fileData, int x1, int y1);
    void plot(string fileFigure, string fileData, int x1, int y1, int x2, int y2);
    void plot(string fileFigure, string fileData, int x1, int y1, int x2, int y2, int x3, int y3);
    void plot(string fileFigure, string fileData, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4);
    void plot(string fileFigure, string fileData, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, 
					          int x5, int y5);
    void plot(string fileFigure, string fileData, int x1, int y1, int x2, int y2, 
						  int x3, int y3, int x4, int y4, int x5, int y5, int x6, int y6);

    void plot(string fileFigure, string fileData, BzzVectorInt &x, BzzVectorInt &y);

    void plot(string fileFigure, string fileData_1, string fileData_2, int x, int y);

    void plot(  string fileFigure, string fileData_1, string fileData_2, string fileData_3,
                int x, int y);

private:

        void ResetValues();
        void plot(string message_a, string message_b);
        string i2s(int number);

        string  folderFigure;
        string *label;
        char   *kind;
        string title;
        string xlabel;
        string ylabel;

        int iSetPlot;
        int iSetLabel;
        int iSetKind;
        int iSetLogScaleX;
        int iSetLogScaleY;
        int iSetRangeX;
        int iSetRangeY;

};

#endif // OPENSMOKE_GNUPLOTINTERFACE
