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

#include "interfaces/OpenSMOKE_GnuPlotWrap.h"
#include <sstream>

#include "interfaces/OpenSMOKE_GnuPlotInterface.h"


OpenSMOKE_GnuPlotInterface::OpenSMOKE_GnuPlotInterface()
{
    folderFigure = "GnuPlot";

    std::string message = "mkdir " + folderFigure;

    label	= new std::string[7];
    kind	= new char[7];
    ResetValues();
}

OpenSMOKE_GnuPlotInterface::OpenSMOKE_GnuPlotInterface(std::string folderName)
{
    folderFigure = folderName;

    std::string message = "mkdir " + folderFigure;

    label = new std::string[7];
    kind = new char[7];
    ResetValues();
}

void OpenSMOKE_GnuPlotInterface::setPlot(std::string _title, std::string _xlabel, std::string _ylabel)
{
    title = _title;
    xlabel = _xlabel;
    ylabel = _ylabel;

    iSetPlot = 1;
}

void OpenSMOKE_GnuPlotInterface::setKey(std::string _label_1)
{
    label[1] = _label_1;
    iSetLabel = 1;
}

void OpenSMOKE_GnuPlotInterface::setKey(std::string _label_1, std::string _label_2)
{
    label[1] = _label_1;    label[2] = _label_2;
    iSetLabel = 2;
}

void OpenSMOKE_GnuPlotInterface::setKey(std::string _label_1, std::string _label_2, std::string _label_3)
{
    label[1] = _label_1;    label[2] = _label_2;    label[3] = _label_3;
    iSetLabel = 3;
}

void OpenSMOKE_GnuPlotInterface::setKey(std::string _label_1, std::string _label_2, std::string _label_3,
                          std::string _label_4)
{
    label[1] = _label_1;    label[2] = _label_2;    label[3] = _label_3;
    label[4] = _label_4;
    iSetLabel = 4;
}

void OpenSMOKE_GnuPlotInterface::setKey(std::string _label_1, std::string _label_2, std::string _label_3,
                          std::string _label_4, std::string _label_5)
{
    label[1] = _label_1;    label[2] = _label_2;    label[3] = _label_3;
    label[4] = _label_4;    label[5] = _label_5;
    iSetLabel = 5;
}

void OpenSMOKE_GnuPlotInterface::setKey(  std::string _label_1, std::string _label_2, std::string _label_3,
                            std::string _label_4, std::string _label_5, std::string _label_6)
{
    label[1] = _label_1;    label[2] = _label_2;    label[3] = _label_3;
    label[4] = _label_4;    label[5] = _label_5;    label[6] = _label_6;
    iSetLabel = 6;
}

void OpenSMOKE_GnuPlotInterface::setKind(char _kind_1)
{
    kind[1]  = _kind_1;
    iSetKind = 1;
}

void OpenSMOKE_GnuPlotInterface::setKind( char _kind_1, char _kind_2)
{
    kind[1]  = _kind_1;    kind[2] = _kind_2;
    iSetKind = 2;
}

void OpenSMOKE_GnuPlotInterface::setKind( char _kind_1, char _kind_2, char _kind_3)
{
    kind[1]  = _kind_1;    kind[2] = _kind_2;    kind[3] = _kind_3;
    iSetKind = 3;
}

void OpenSMOKE_GnuPlotInterface::setKind( char _kind_1, char _kind_2, char _kind_3, char _kind_4)
{
    kind[1]  = _kind_1;    kind[2] = _kind_2;    kind[3] = _kind_3;
    kind[4]  = _kind_4;
    iSetKind = 4;
}

void OpenSMOKE_GnuPlotInterface::setKind( char _kind_1, char _kind_2, char _kind_3,
                            char _kind_4, char _kind_5)
{
    kind[1]  = _kind_1;    kind[2] = _kind_2;    kind[3] = _kind_3;
    kind[4]  = _kind_4;    kind[5] = _kind_5;
    iSetKind = 5;
}

void OpenSMOKE_GnuPlotInterface::setKind( char _kind_1, char _kind_2, char _kind_3,
                            char _kind_4, char _kind_5, char _kind_6)
{
    kind[1]  = _kind_1;    kind[2] = _kind_2;    kind[3] = _kind_3;
    kind[4]  = _kind_4;    kind[5] = _kind_5;    kind[6] = _kind_6;
    iSetKind = 6;
}

void OpenSMOKE_GnuPlotInterface::setLogScaleX()
{
    iSetLogScaleX = 1;
}

void OpenSMOKE_GnuPlotInterface::setLogScaleY()
{
    iSetLogScaleY = 1;
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1)
{
    BzzVectorInt x(1, x1);
    BzzVectorInt y(1, y1);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1, int x2, int y2)
{
    BzzVectorInt x(2, x1, x2);
    BzzVectorInt y(2, y1, y2);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1, int x2, int y2, int x3, int y3)
{
    BzzVectorInt x(3, x1, x2, x3);
    BzzVectorInt y(3, y1, y2, y3);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4)
{
    BzzVectorInt x(4, x1, x2, x3, x4);
    BzzVectorInt y(4, y1, y2, y3, y4);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int x5, int y5)
{
    BzzVectorInt x(5, x1, x2, x3, x4, x5);
    BzzVectorInt y(5, y1, y2, y3, y4, y5);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, int x1, int y1, int x2, int y2, 
							int x3, int y3, int x4, int y4, int x5, int y5, int x6, int y6)
{
    BzzVectorInt x(6, x1, x2, x3, x4, x5, x6);
    BzzVectorInt y(6, y1, y2, y3, y4, y5, y6);
    plot(fileFigure, fileData, x, y);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData, BzzVectorInt &x, BzzVectorInt &y)
{
    int i;

    std::string message_a = "set output \"";
    message_a += folderFigure + "/" + fileFigure + ".eps\"";

    std::string message_b = "p ";
    for (i=1;i<=x.Size();i++)
    {
        message_b += " '" + fileData + "' " + " u ";
        message_b += i2s(x[i]) + ":" + i2s(y[i]);
        if(iSetKind > 0)
            if      (kind[i] == 'l') message_b +=  " w l lw 5";
            else if (kind[i] == 'p') message_b +=  " w p ps 3";
        if(iSetLabel > 0) message_b += " ti \"" + label[i] + "\"";

        if (i!=x.Size())    message_b += " , ";
    }

    plot(message_a, message_b);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure, std::string fileData_1, std::string fileData_2, int x, int y)
{
    std::string message_a = "set output \"";
    message_a += folderFigure + "/" + fileFigure + ".eps\"";

    std::string message_b = "p ";
    message_b += " '" + fileData_1 + "' " + " u ";
    message_b += i2s(x) + ":" + i2s(y);
    if(iSetKind > 0)
        if      (kind[1] == 'l') message_b +=  " w l lw 5";
        else if (kind[1] == 'p') message_b +=  " w p ps 3";
    if(iSetLabel > 0) message_b += " ti \"" + label[1] + "\"";
    message_b += " , ";

    message_b += " '" + fileData_2 + "' " + " u ";
    message_b += i2s(x) + ":" + i2s(y);
    if(iSetKind > 0)
        if      (kind[2] == 'l') message_b +=  " w l lw 5";
        else if (kind[2] == 'p') message_b +=  " w p ps 3";
    if(iSetLabel > 0) message_b += " ti \"" + label[2] + "\"";

    plot(message_a, message_b);
}

void OpenSMOKE_GnuPlotInterface::plot(std::string fileFigure,
                        std::string fileData_1, std::string fileData_2, std::string fileData_3,
                        int x, int y)
{
    std::string message_a = "set output '";
    message_a += folderFigure + "/" + fileFigure + ".eps'";

    std::string message_b = "p ";
    message_b += " '" + fileData_1 + "' " + " u ";
    message_b += i2s(x) + ":" + i2s(y);
    if(iSetKind > 0)
        if      (kind[1] == 'l') message_b +=  " w l lw 5";
        else if (kind[1] == 'p') message_b +=  " w p ps 3";
    if(iSetLabel > 0) message_b += " ti '" + label[1] + "'";
    message_b += " , ";

    message_b += " '" + fileData_2 + "' " + " u ";
    message_b += i2s(x) + ":" + i2s(y);
    if(iSetKind > 0)
        if      (kind[2] == 'l') message_b +=  " w l lw 5";
        else if (kind[2] == 'p') message_b +=  " w p ps 3";
    if(iSetLabel > 0) message_b += " ti '" + label[2] + "'";
    message_b += " , ";

    message_b += " '" + fileData_3 + "' " + " u ";
    message_b += i2s(x) + ":" + i2s(y);
    if(iSetKind > 0)
        if      (kind[3] == 'l') message_b +=  " w l lw 5";
        else if (kind[3] == 'p') message_b +=  " w p ps 3";
    if(iSetLabel > 0) message_b += " ti '" + label[3] + "'";


    plot(message_a, message_b);
}


std::string OpenSMOKE_GnuPlotInterface::i2s(int number)
{
    stringstream a;
    a << number;
    return a.str();
}

void OpenSMOKE_GnuPlotInterface::plot(std::string message_a, std::string message_b)
{
    Gnuplot g = Gnuplot("");
    g.cmd("set terminal postscript eps enhanced color solid");

    g.cmd("set grid");
    if (iSetPlot==1)
    {
        g.set_title(title);
        g.set_xlabel(xlabel).set_ylabel(ylabel);
    }

    if (iSetLogScaleX == 1) g.cmd("set logscale x");
    if (iSetLogScaleY == 1) g.cmd("set logscale y");

    if (iSetRangeX == 1)    g.cmd("");
    if (iSetRangeY == 1)    g.cmd("");

    g.cmd(message_a);
    g.cmd(message_b);

    ResetValues();
}

void OpenSMOKE_GnuPlotInterface::ResetValues()
{
    iSetPlot        = 0;
    iSetLabel       = 0;
    iSetLogScaleX   = 0;
    iSetLogScaleY   = 0;
    iSetRangeX      = 0;
    iSetRangeY      = 0;
}
