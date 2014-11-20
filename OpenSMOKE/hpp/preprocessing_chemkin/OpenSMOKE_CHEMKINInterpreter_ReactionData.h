/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_REACTIONDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_REACTIONDATA_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_CHEMKINInterpreter_KineticsData;
class OpenSMOKE_CHEMKINInterpreter_ThermoData;
class OpenSMOKE_CHEMKINInterpreter_UnitsData;
class OpenSMOKE_CHEMKINInterpreter_TransportData;

class OpenSMOKE_CHEMKINInterpreter_ReactionData
{
public:

	OpenSMOKE_CHEMKINInterpreter_ReactionData();
	void Setup(const int _iReaction, const int _iLine, OpenSMOKE_CHEMKINInterpreter_KineticsData *_ptKinetics, OpenSMOKE_CHEMKINInterpreter_UnitsData  *_ptUnits);
	void Setup(const double _energy_factor, const double _frequency_factor);
	void ReverseRateFromReaction(const int _iReactionBis, OpenSMOKE_CHEMKINInterpreter_ReactionData &r);

	virtual ~OpenSMOKE_CHEMKINInterpreter_ReactionData();

	std::string reaction_string_clean;
	std::string reaction_string_complete;

	double A;
	double Beta;
	double E;

	double energy_factor;
	double frequency_factor;

	int nReactants;
	vector<string>	nameDirect;
	BzzVectorInt	indexDirect;
	BzzVector		nuDirect;
	int nProducts;
	vector<string>	nameInverse;
	BzzVectorInt	indexInverse;
	BzzVector		nuInverse;

	int nGlobalReactants;
	vector<string>	nameGlobalDirect;
	BzzVectorInt	indexGlobalDirect;
	BzzVector		lambdaGlobalDirect;

	int nGlobalProducts;
	vector<string>	nameGlobalInverse;
	BzzVectorInt	indexGlobalInverse;
	BzzVector		lambdaGlobalInverse;

	bool iReversible;
	bool iGlobal;
	bool iThermodynamicConsistency;

	BzzVectorInt	indexDirectOrdered;
	BzzVectorInt	indexInverseOrdered;
	BzzVectorInt	indexDirectOrderedSequence;
	BzzVectorInt	indexInverseOrderedSequence;

	// Pressure dependent reactions
	kindOfPressureDependence kindPressureDependence;
	BzzVector troe;
	BzzVector sri;

	// High-pressure limits (Chemically Activated Bimolecular Reactions)
	bool iHighPressure;
	BzzVector high_pressure;

	// Low-pressure limits (Fall-off Reactions)
	bool iLowPressure;
	BzzVector low_pressure;

	// Third-body
	bool			iThirdBody;
	BzzVector		efficiencies_coefficients;
	BzzVectorInt	efficiencies_indices;

	void Chebishev(vector<string> instructions);
	void CollisionalEfficiency(vector<string> instructions);
	void Duplicate(vector<string> instructions);
	void EnergyLossParameter(vector<string> instructions);
	void Fit1(vector<string> instructions);
	void ForwardReactionKineticParameter(vector<string> instructions);
	void HighPressureLimit(vector<string> instructions);
	void JanevLangerReactionRate(vector<string> instructions);
	void LowPressureLimit(vector<string> instructions);
	void LandauTellerReaction(vector<string> instructions);
	void MomentumTransferCollisionFrequency(vector<string> instructions);
	void ChebishevPressureLimits(vector<string> instructions);
	void PressureLogarithmicInterpolation(vector<string> instructions);
	void ReverseRate(vector<string> instructions);
	void LandauTellerReverseReaction(vector<string> instructions);
	void BackwardReactionKineticParameter(vector<string> instructions);
	void ChebishevTemperatureLimits(vector<string> instructions);
	void SRIFallOffReaction(vector<string> instructions);
	void SpeciesTemperatureDependence(vector<string> instructions);
	void TROEFallOffReaction(vector<string> instructions);
	void Units(vector<string> instructions);
	void UserRateSubRoutine(vector<string> instructions);
	void CollisionCrossSection(vector<string> instructions);
	void ThirdBodyEfficiencies(vector<string> instructions);
	void ThirdBodySingle();
	void ThermodynamicConsistency(vector<string> instructions);
	void TAR(vector<string> instructions);

	void Check();
	void CheckStoichiometry(OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, BzzVectorInt &indices, vector<string> &element_names, vector<string> &isotope_names);
	void CompactStoichiometry();
	void ReactionString();
	void Complete();
	void PrintOnFile(ofstream &fOutput);
	void PrintOnBinaryFile(BzzSave &outputFile, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter_TransportData &transport);
	void SummaryOnFile(ofstream &fOutput);
	void WriteWarningMessageOnFile(const std::string message);



	OpenSMOKE_CHEMKINInterpreter_KineticsData	*ptKinetics;
	OpenSMOKE_CHEMKINInterpreter_UnitsData		*ptUnits;
	int iLine;
	int iReaction;
	int iReactionBis;
	int iReverseReaction;

	double sumNu;
	double lambdaTotal;
	double conversion_frequency_factor;

	bool iDuplicate;

	bool iLandauTeller;
	BzzVector landau_teller;

	bool iJanevLanger;
	BzzVector janev_langer;

	bool iPowerSeries;
	BzzVector power_series;

	bool iPressureLogarithmic;
	BzzMatrix pressure_logarithmic;

	bool iChebishevPolynomials;
	BzzVector chebishev_vector;
	BzzVector chebishev_pressure;
	BzzVector chebishev_temperature;

	bool iReverseRate;
	BzzVector reverse_rate;

	bool iReverseRateLandauTeller;
	BzzVector reverse_rate_landau_teller;

	bool iConversionEnergy;
	bool iConversionFrequencyFactor;

	bool iCollisionEfficiency;
	double kStar_collision_frequency;
	double dAB;
	double WAB;

	bool iTAR;
	BzzVector TAR_series;

	bool IsReversibleReaction();
	bool IsThirdBodyReaction();
	bool IsLandauTellerReaction();

	std::string thirdbodysingle;

private:
	
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	std::string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_REACTIONDATA_H)
