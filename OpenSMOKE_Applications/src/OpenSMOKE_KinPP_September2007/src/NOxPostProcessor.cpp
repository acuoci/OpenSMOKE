
#include "BzzMath.hpp"

#include "myBzzCSTRNetwork.hpp"
#include <stdio.h>
#include <iostream>

using namespace std;


int main(void)
{

	cout << endl;
	cout << "----------------------------------------------------------" << endl;
	cout << "-                                                        -" << endl;
	cout << "-                    NOX Post Processor                  -" << endl;
	cout << "-                                                        -" << endl;
	cout << "-          Politecnico di Milano - Maggio 2007           -" << endl;
	cout << "-                                                        -" << endl;
	cout << "----------------------------------------------------------" << endl;
	cout << endl;


	BzzPrint("\n\n\nCSTR Network");
	{
		bzzFileOut = stdout; 
		RemoveWarningWindow();
		int cicloCluster,cicloDiffusion,iaia,relaxation;
		double startUser = BzzGetUserTime();
		double startKernel = BzzGetKernelTime();
		double start = BzzGetCpuTime();
		
		myBzzCSTRNetwork cstr;
		cstr.SetMemoTemperature();

		cout << "--------------------------------------------" << endl;
		cout << "Diffusion Relaxation:                       " << endl;
		cout << "--------------------------------------------" << endl;
		cout << " 0 = Without Relaxation					 " << endl;
		cout << " 1 = With Relaxation                        " << endl;
		cout << " 2 = Start again from previous results      " << endl;
		cout << endl;
		cout << " Enter your choice: ";
		cin >> relaxation;
		cout << endl << endl;
		
		cout << "--------------------------------------------" << endl;
		cout << "Network Clustering                          " << endl;
		cout << "--------------------------------------------" << endl;
		cout << "0  = Original number of reactors            " << endl;
		cout << ">0 = Network Clustering                     " << endl;
		cout << endl;
		cout << " Enter your choice: ";
		cin  >> cicloCluster;
		cout << endl << endl;

		// No Relaxation
		if(relaxation == 0)
		{
			cicloDiffusion = 1;

			cout << "-----------------------------------------------------" << endl;
			cout << "Initial Data                                         " << endl;
			cout << "-----------------------------------------------------" << endl;
			cout << "0  = Data from fluent results (file FirstGuess.bzz)  " << endl;
			cout << "1  = Data from previous iteration (file mass.tmp)    " << endl;
			cout << endl;
			cout << " Enter your choice: ";
			cin  >> iaia;
			cout << endl << endl;

			
			cstr("Input/CFDNetwork.bzz","Input/FirstGuess.bzz", cicloCluster, cicloDiffusion, iaia,relaxation);
			cstr.OutputPrint(0);
			cstr.Save('*',"Temp/mass.tmp");
		}

		// Start from previous solution
		else if(relaxation == 2)
		{
			cicloDiffusion = 1;
			iaia = 1;
			cstr("Input/CFDNetwork.bzz","Input/FirstGuess.bzz", cicloCluster, cicloDiffusion, iaia, relaxation);
			cstr.OutputPrint(0);
			cstr.Save('*',"Temp/mass.tmp");
		}
		
		// Relaxation
		else if(relaxation == 1)
		{
			cicloDiffusion = 4;
			
			cout << "-----------------------------------------------------" << endl;
			cout << "Initial Data                                         " << endl;
			cout << "-----------------------------------------------------" << endl;
			cout << "0  = Data from fluent results (file FirstGuess.bzz)  " << endl;
			cout << "1  = Data from previous iteration (file mass.tmp)    " << endl;
			cout << endl;
			cout << " Enter your choice: ";
			cin  >> iaia;
			cout << endl << endl;
			
			cstr("Input/CFDNetwork.bzz","Input/FirstGuess.bzz", cicloCluster, cicloDiffusion, iaia, relaxation);
			cstr.OutputPrint(0);
			cstr.Save('*',"Temp/mass.tmp");
			iaia = 1;
			for(cicloDiffusion = 3;cicloDiffusion > 0;cicloDiffusion--)
			{
				printf("\ncicloDiffusion %d",cicloDiffusion);
				cstr("Input/CFDNetwork.bzz","Input/FirstGuess.bzz", cicloCluster, cicloDiffusion, iaia, relaxation);
				cstr.OutputPrint(0);
				cstr.Save('*',"Temp/mass.tmp");
			}
		}

		::BzzPrint("\nCpu Seconds for complete solution: %e",BzzGetCpuTime() - start);
		::BzzPrint("\nUser Seconds for complete solution: %e",BzzGetUserTime() - startUser);
		::BzzPrint("\nKernel Seconds for complete solution: %e",BzzGetKernelTime() - startKernel);
		BzzPause();

	}

	return 0;
}