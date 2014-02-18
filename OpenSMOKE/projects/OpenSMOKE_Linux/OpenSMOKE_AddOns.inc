# COMPILER
PATH_OBJECTS 	= OpenSMOKE_AddOns/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/addons
PATH_HPP 	= ../../hpp/addons


# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_2EModel.o \
	$(PATH_OBJECTS)/OpenSMOKE_ElementFluxAnalysis.o \
	$(PATH_OBJECTS)/OpenSMOKE_InverseKinetics.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Flame1D.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_General.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_PFR.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_CSTR.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ShockTube.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Batch.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ICEM_MultiZone.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.o \
	$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.o \
	$(PATH_OBJECTS)/OpenSMOKE_RateOfProductionAnalysis.o \
	$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis.o \
	$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.o \
	$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.o \
	$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.o \
	$(PATH_OBJECTS)/OpenSMOKE_SootManager.o \
	$(PATH_OBJECTS)/OpenSMOKE_PolimiSoot.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_AddOns_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_AddOns_$(STRING_NAME).a $(OBJS)


# CHEMISTRY

#1
$(PATH_OBJECTS)/OpenSMOKE_2EModel.o : $(PATH_CPP)/OpenSMOKE_2EModel.cpp $(PATH_HPP)/OpenSMOKE_2EModel.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_2EModel.cpp -o $(PATH_OBJECTS)/OpenSMOKE_2EModel.o

#2
$(PATH_OBJECTS)/OpenSMOKE_ElementFluxAnalysis.o : $(PATH_CPP)/OpenSMOKE_ElementFluxAnalysis.cpp $(PATH_HPP)/OpenSMOKE_ElementFluxAnalysis.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ElementFluxAnalysis.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ElementFluxAnalysis.o

#3
$(PATH_OBJECTS)/OpenSMOKE_InverseKinetics.o : $(PATH_CPP)/OpenSMOKE_InverseKinetics.cpp $(PATH_HPP)/OpenSMOKE_InverseKinetics.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_InverseKinetics.cpp -o $(PATH_OBJECTS)/OpenSMOKE_InverseKinetics.o

#4
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor.o : $(PATH_CPP)/OpenSMOKE_PostProcessor.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor.o

#5
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Flame1D.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_Flame1D.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_Flame1D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_Flame1D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Flame1D.o

#6
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_General.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_General.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_General.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_General.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_General.o

#7
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_PFR.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_PFR.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_PFR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_PFR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_PFR.o

#8
$(PATH_OBJECTS)/OpenSMOKE_RateOfProductionAnalysis.o : $(PATH_CPP)/OpenSMOKE_RateOfProductionAnalysis.cpp $(PATH_HPP)/OpenSMOKE_RateOfProductionAnalysis.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_RateOfProductionAnalysis.cpp -o $(PATH_OBJECTS)/OpenSMOKE_RateOfProductionAnalysis.o

#9
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.o

#10
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.o

#11
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.o	

#12
$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis.o : $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis.cpp $(PATH_HPP)/OpenSMOKE_SensitivityAnalysis.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis.o

#13
$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.o : $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.cpp $(PATH_HPP)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.o

#14
$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.o : $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.cpp $(PATH_HPP)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.o

#15
$(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.o : $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.cpp $(PATH_HPP)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.o

#16
$(PATH_OBJECTS)/OpenSMOKE_SootManager.o : $(PATH_CPP)/OpenSMOKE_SootManager.cpp $(PATH_HPP)/OpenSMOKE_SootManager.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SootManager.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SootManager.o

#17


#18
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_CSTR.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_CSTR.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_CSTR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_CSTR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_CSTR.o

#19
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ShockTube.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_ShockTube.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_ShockTube.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_ShockTube.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ShockTube.o	

#20
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Batch.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_Batch.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_Batch.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_Batch.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_Batch.o

#21
$(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ICEM_MultiZone.o : $(PATH_CPP)/OpenSMOKE_PostProcessor_ICEM_MultiZone.cpp $(PATH_HPP)/OpenSMOKE_PostProcessor_ICEM_MultiZone.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PostProcessor_ICEM_MultiZone.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PostProcessor_ICEM_MultiZone.o

#22
$(PATH_OBJECTS)/OpenSMOKE_PolimiSoot.o : $(PATH_CPP)/OpenSMOKE_PolimiSoot.cpp $(PATH_HPP)/OpenSMOKE_PolimiSoot.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PolimiSoot.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PolimiSoot.o
		 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_AddOns_$(STRING_NAME).a

