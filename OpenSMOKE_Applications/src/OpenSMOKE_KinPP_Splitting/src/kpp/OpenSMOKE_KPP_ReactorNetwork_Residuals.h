/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
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

#ifndef OpenSMOKE_KPP_ReactorNetwork_Residuals_H
#define OpenSMOKE_KPP_ReactorNetwork_Residuals_H

#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"

class OpenSMOKE_KPP_ReactorNetwork;

class OpenSMOKE_KPP_ReactorNetwork_Residuals
{

public:

	OpenSMOKE_KPP_ReactorNetwork_Residuals(OpenSMOKE_KPP_ReactorNetwork& network);

	void SetNumberOfTimeSteps(const int nTimeStep_);
	void Analysis();	
	void WriteResidualsOnVideo();

	inline const BzzVector& residuals()			const { return residuals_; }
	
	inline const BzzVector& speciesResiduals()	const { return speciesResiduals_; }
	inline const BzzVector& speciesNormInf()	const { return speciesNormInf_; }
	inline const BzzVector& speciesNorm1()		const { return speciesNorm1_; }
	inline const BzzVector& speciesNorm2()		const { return speciesNorm2_; }

	inline const BzzVector& reactorResiduals()	const { return reactorResiduals_; }
	inline const BzzVector& reactorNormInf()	const { return reactorNormInf_; }
	inline const BzzVector& reactorNorm1()		const { return reactorNorm1_; }
	inline const BzzVector& reactorNorm2()		const { return reactorNorm2_; }

	inline double normInf() const	{ return normInf_; }
	inline double norm1()	const	{ return norm1_; }
	inline double norm2()	const	{ return norm2_; }

	inline double RatioNormInf()	const { return ratioNormInf->MeanRatios(); }
	inline double RatioNorm1()		const { return ratioNorm1->MeanRatios(); }

private:

	OpenSMOKE_KPP_ReactorNetwork& network_;
	
	BzzVector residuals_;

	double normInf_;
	double norm1_;
	double norm2_;

	BzzVector speciesResiduals_;
	BzzVector speciesNormInf_;
	BzzVector speciesNorm1_;
	BzzVector speciesNorm2_;
	
	BzzVector reactorResiduals_;
	BzzVector reactorNormInf_;
	BzzVector reactorNorm1_;
	BzzVector reactorNorm2_;

	BzzVector historyNormInf_;
	BzzVector historyNorm1_;
	BzzVector historyNorm2_;

	int nTimeSteps_;
	int currentTimeStep_;
	int nFrequencyRefinedAnalysis_;

	double maxSpecies;
	double minSpecies;

	ofstream fResiduals;
	ofstream fResidualsSpecies;
	ofstream fResidualsReactorStatistics;
	ofstream fResidualsSpeciesStatistics;

	void AnalysisRefined();

	void WriteResidualsLabels();
	void WriteResidualsSpeciesLabels();
	void WriteResidualsStatisticsLabels();

	void WriteResidualsOnFile();
	void WriteResidualsSpeciesOnFile();

	void Statistics();

	BzzVector ranges_;

	KPP_Network_Status currentStatus_;
	int localIteration_;
	double startingNormInf_;
	double startingNorm1_;
	double startingNorm2_;

	RingVector<double>* ratioNormInf;
	RingVector<double>* ratioNorm1;
        
        static const int nRanges_;
        static const double endingRange_;
        static const double startingRange_;
        
private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};

#endif // OpenSMOKE_KPP_ReactorNetwork_Residuals_H