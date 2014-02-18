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


#include "basic/OpenSMOKE_Constants.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "basic/OpenSMOKE_Dictionary.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Grid1D.h"

#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_Kinetics.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "engine/OpenSMOKE_GlobalKinetics.h"


#include "idealreactors/OpenSMOKE_GasStream.h"
#include "idealreactors/OpenSMOKE_UD_Profile.h"
#include "idealreactors/OpenSMOKE_0DReactor.h"

#include "idealreactors/pfr/OpenSMOKE_PFR_Geometry.h"
//#include "idealreactors/pfr/OpenSMOKE_PFR.h"
//#include "idealreactors/batch/OpenSMOKE_Batch.h"

#include "addons/OpenSMOKE_2EModel.h"
#include "addons/OpenSMOKE_SootManager.h"

#include "preprocessing/OpenSMOKE_PreProcessorIdealGas.h"
#include "preprocessing/OpenSMOKE_PreProcessorKinetics.h"
#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"

#include "interfaces/OpenSMOKE_GnuPlotInterface.h"
#include "interfaces/OpenSMOKE_GnuPlotWrap.h"
#include "interfaces/OpenSMOKE_LatexInterface.h"
#include "interfaces/SimpleGlob.h"
#include "interfaces/SimpleOpt.h"

#include "qmom/OpenSMOKE_QMOM_Module.h"

