/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_DROPLETMICROGRAVITY_DAE_OBJECTS_H)
#define OPENSMOKE_DROPLETMICROGRAVITY_DAE_OBJECTS_H

#include "BzzMath.hpp"

class OpenSMOKE_Droplet;
class OpenSMOKE_DropletMicrogravity;

class OpenSMOKE_DropletMicrogravity_MyDaeSystemEigenValue : public BzzDaeSystemObject
{
public:
	void assignDroplet(OpenSMOKE_DropletMicrogravity *droplet);

	OpenSMOKE_DropletMicrogravity *ptDroplet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_DropletMicrogravity_MyDaeSystemUnsteadyBatch : public BzzDaeSystemObject
{
public:
	void assignDroplet(OpenSMOKE_DropletMicrogravity *droplet);

	OpenSMOKE_DropletMicrogravity *ptDroplet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Droplet_MyDaeSystemNoMomentum : public BzzDaeSystemObject
{
public:
	void assignDroplet(OpenSMOKE_Droplet *droplet);

	OpenSMOKE_Droplet *ptDroplet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Droplet_MyDaeSystemOnlyTemperature : public BzzDaeSystemObject
{
public:
	void assignDroplet(OpenSMOKE_Droplet *droplet);

	OpenSMOKE_Droplet *ptDroplet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

#endif // !defined(OPENSMOKE_DROPLETMICROGRAVITY_DAE_OBJECTS_H)


/*

		else
		{
			double time0=0.;
			double TEND = 0.;
			double deltat=data->deltaTime;
			setInitialConditions();
			Summary();	
			PrepareSystem();
			daeSystemEigenValue.assignDroplet(this);
			BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &daeSystemEigenValue, dimensionBlock);	
		//	DAESystemSolution(&o, data->tEnd);  
			DAESystemSolution(&o, time0+deltat);  

			ofstream fDynamic;
			openOutputFileAndControl(fDynamic, "Dynamic.out");
			fDynamic.setf(ios::scientific);

			ofstream fMap;
			openOutputFileAndControl(fMap, "Map.out");
			fMap.setf(ios::scientific);

			//iDamping = false;
			for(;;)
			{
				int iTMax, iOHMax;
				double TMax = T.Max(&iTMax);
				double OHMax = OmegaGas.GetColumn(mix->recognize_species("OH")).Max(&iOHMax);

				cout << "Time: " << time0 << endl;
				fDynamic << setw(16) << left << time0;
				fDynamic << setw(16) << left << diameterDroplet*1000.;												// [mm]
				fDynamic << setw(16) << left << BzzPow2(diameterDroplet/data->diameterDroplet0);					// [-]
				fDynamic << setw(16) << left << mdot*1e3;																// [g/s]
				fDynamic << setw(16) << left << 4./Constants::pi*m[1]/LiquidDensity()/diameterDroplet * 1e6;			// dD2/dt[mm2/s]
				fDynamic << setw(16) << left << 4.*m[1]/LiquidDensity()/diameterDroplet * 1e6;							// dA/dt[mm2/s]
				fDynamic << setw(16) << left << 2./Constants::pi*m[1]/LiquidDensity()/BzzPow2(diameterDroplet) * 1e3;	// dD/dt[mm/s]
				fDynamic << setw(16) << left << massDroplet*1000.;														// [g]
				fDynamic << setw(16) << left << massDroplet/data->massDroplet0;											// [-]
				fDynamic << setw(16) << left << 2.*grid.x[iTMax]/diameterDroplet;									// [-]
				fDynamic << setw(16) << left << TMax;																	// [K]
				fDynamic << setw(16) << left << 2.*grid.x[iOHMax]/diameterDroplet;									// [-]
				fDynamic << setw(16) << left << OHMax;																	// [-]
				fDynamic << endl;

				int jH2O = mix->recognize_species("H2O");
				int jCO2 = mix->recognize_species("CO2");
				int jOH = mix->recognize_species("OH");
				for(int i=1;i<=N;i++)
				{
					fMap << setw(16) << left << time0;
					fMap << setw(16) << left << grid.x[i];
					fMap << setw(16) << left << grid.x[i]/(data->diameterDroplet0/2.);
					fMap << setw(16) << left << T[i];
					for(int j=1;j<=NCLiquid;j++)
						fMap << setw(16) << left << OmegaGas[i][data->jFuel[j]];
					fMap << setw(16) << left << OmegaGas[i][data->jOxidizer];
					fMap << setw(16) << left << OmegaGas[i][jH2O];
					fMap << setw(16) << left << OmegaGas[i][jCO2];
					fMap << setw(16) << left << OmegaGas[i][jOH];
					fMap << endl;
				}
				fMap << endl;

				cout << "Before" << endl;
				cout << " Mass:     " << massDroplet << endl;
				cout << " Diameter: " << diameterDroplet << endl;
				cout << " M:        " << mdot << " " << m[1] << endl;
			
				for(int j=1;j<=NCLiquid;j++)
				{
					rhoDroplet[j]		= data->liquid_species[j].rho(TDroplet[1]);
					dHVaporization[j]	= data->liquid_species[j].Hv(TDroplet[1]);
					pvDroplet[j]		= data->liquid_species[j].Pv(TDroplet[1]);
					cpDroplet[j]		= data->liquid_species[j].Cp(TDroplet[1]);
				}
				massDroplet    -= deltat*mdot;
	//			data->TDroplet += 6./rhoDroplet/data->diameterDroplet/cpDroplet * 
		//						  (-mdot/Constants::pi/BzzPow2(data->diameterDroplet)*dHVaporization + lambda[1]*dT_over_dr[1]);
				diameterDroplet = pow(6.*massDroplet/Constants::pi/LiquidDensity(),1./3.); 
				grid.Construct(N, data->ratioRadii*(data->diameterDroplet0/2.)-diameterDroplet/2., data->stretchingFactor, diameterDroplet/2.);
	
				cout << "After" << endl;
				cout << " Mass:     " << massDroplet << endl;
				cout << " Diameter: " << diameterDroplet << endl;
				cout << " M:        " << mdot << " " << m[1] << endl;

				time0 += deltat;
				o(xFirstGuess, time0, inDerAlg, &daeSystemEigenValue, dimensionBlock);	
				DAESystemSolutionAgain(&o, time0+deltat); 
			
			
				if (diameterDroplet/data->diameterDroplet0 < 0.10)
					break;

				deltat *= data->ratioDeltaTime;; 
			}

			fDynamic.close();
			fMap.close();

			PrintFinalSolution();
			PrintAdditionalFinalSolution();
		}
*/

/*

void OpenSMOKE_DropletMicrogravity::DAESystemSolutionAgain(BzzDaeSparseObject *o, double tEnd)
{	
	o->StepPrint(DAE_Print);
	o->SetMinimumConstraints(xMin);
	o->SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) 100*MachEps()
	o->SetTollRel(data->relTolerances);
	o->SetTollAbs(data->absTolerances);

	double timeStart = BzzGetCpuTime();
	bzzStop = 0;  
	xFirstGuess = (*o)(tEnd, tEnd);

	cout << endl;
	cout << "Number of steps: "					<< o->GetNumStep() << endl;
	cout << "Number of function for Jacobian: " << o->GetNumFunctionForJacobian() << endl;
	cout << "Numerical Jacobians: "				<< o->GetNumNumericalJacobian() << endl;
	
	cout << "Time DAE solution: "				<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
}

*/