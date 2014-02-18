#COMPILER
PATH_OBJECTS 	= OpenSMOKE_SymbolicKinetics/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/symbolickinetics
PATH_HPP 	= ../../hpp/symbolickinetics

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_SymbolicKinetics.o \
	$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.o \
	$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_Jacobian.o \
	$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_ReactionRates.o \
	$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.o \
	$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_Jacobian.o \
	$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_ReactionRates.o \
	$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.o \
	$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_Jacobian.o \
	$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.o \
	$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_Jacobian.o \
	$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_ReactionRates.o \
	$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.o \
	$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_Jacobian.o \
	$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_Jacobian.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_Jacobian.o \
	$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.o \
	$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_Jacobian.o \
	$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.o \
	$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_Jacobian.o \
	$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.o \
	$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_Jacobian.o \
	$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_ReactionRates.o \
	$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.o \
	$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_Jacobian.o \
	$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_ReactionRates.o


# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_SymbolicKinetics_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_SymbolicKinetics_$(STRING_NAME).a $(OBJS)


# MAKE

$(PATH_OBJECTS)/OpenSMOKE_SymbolicKinetics.o : $(PATH_CPP)/OpenSMOKE_SymbolicKinetics.cpp $(PATH_HPP)/OpenSMOKE_SymbolicKinetics.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_SymbolicKinetics.cpp -o $(PATH_OBJECTS)/OpenSMOKE_SymbolicKinetics.o

#GRI30
$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.o : $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.cpp $(PATH_HPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.h
		 	 $(CCP) $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.cpp -o $(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.o

$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_Jacobian.o : $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_Jacobian.cpp $(PATH_HPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.h
		 	 $(CCP) $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_Jacobian.cpp -o $(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_Jacobian.o

$(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_ReactionRates.o : $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_ReactionRates.cpp $(PATH_HPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30.h
		 	 $(CCP) $(PATH_CPP)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_ReactionRates.cpp -o $(PATH_OBJECTS)/gri30/OpenSMOKE_SymbolicKinetics_GRI30_ReactionRates.o

#GRI12
$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.o : $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.cpp $(PATH_HPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.h
		 	 $(CCP) $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.cpp -o $(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.o

$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_Jacobian.o : $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_Jacobian.cpp $(PATH_HPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.h
		 	 $(CCP) $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_Jacobian.cpp -o $(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_Jacobian.o

$(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_ReactionRates.o : $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_ReactionRates.cpp $(PATH_HPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12.h
		 	 $(CCP) $(PATH_CPP)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_ReactionRates.cpp -o $(PATH_OBJECTS)/gri12/OpenSMOKE_SymbolicKinetics_GRI12_ReactionRates.o

#FLUENT_DRM22_POLIMI_NOX
$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.o : $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.cpp $(PATH_HPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h
		 	 $(CCP) $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.cpp -o $(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.o

$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_Jacobian.o : $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_Jacobian.cpp $(PATH_HPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h
		 	 $(CCP) $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_Jacobian.cpp -o $(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_Jacobian.o

$(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_ReactionRates.o : $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_ReactionRates.cpp $(PATH_HPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h
		 	 $(CCP) $(PATH_CPP)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_ReactionRates.cpp -o $(PATH_OBJECTS)/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX_ReactionRates.o


#Polimi_H2CONOX_1101
$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.o : $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.cpp $(PATH_HPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.cpp -o $(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.o

$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_Jacobian.o : $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_Jacobian.cpp $(PATH_HPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_Jacobian.o

$(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_ReactionRates.o : $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_ReactionRates.cpp $(PATH_HPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101_ReactionRates.o


#Polimi_H2CONOX_NoThermal_1101
$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.o : $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.cpp $(PATH_HPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.o

$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_Jacobian.o : $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_Jacobian.cpp $(PATH_HPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_Jacobian.o

$(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_ReactionRates.o : $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_ReactionRates.cpp $(PATH_HPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101_ReactionRates.o


#Polimi_H2CONOX_NoNNH_1101
$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.o : $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.cpp $(PATH_HPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.o

$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_Jacobian.o : $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_Jacobian.cpp $(PATH_HPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_Jacobian.o

$(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_ReactionRates.o : $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_ReactionRates.cpp $(PATH_HPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101_ReactionRates.o


#Polimi_H2CONOX_NoN2O_1101
$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.o : $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.cpp $(PATH_HPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.cpp -o $(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.o

$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_Jacobian.o : $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_Jacobian.cpp $(PATH_HPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_Jacobian.o

$(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_ReactionRates.o : $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_ReactionRates.cpp $(PATH_HPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101_ReactionRates.o



#SanDiego_AVIO
$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.o : $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.cpp $(PATH_HPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.h
		 	 $(CCP) $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.cpp -o $(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.o

$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_Jacobian.o : $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_Jacobian.cpp $(PATH_HPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.h
		 	 $(CCP) $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_Jacobian.cpp -o $(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_Jacobian.o

$(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_ReactionRates.o : $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_ReactionRates.cpp $(PATH_HPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.h
		 	 $(CCP) $(PATH_CPP)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_ReactionRates.cpp -o $(PATH_OBJECTS)/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO_ReactionRates.o

#Polimi_c1c3htnox_1101
$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.o : $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.cpp $(PATH_HPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.o

$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_Jacobian.o : $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_Jacobian.cpp $(PATH_HPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_Jacobian.o

$(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_ReactionRates.o : $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_ReactionRates.cpp $(PATH_HPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101_ReactionRates.o


#Polimi_nc7_avio
$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.o : $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.cpp $(PATH_HPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.h
		 	 $(CCP) $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.cpp -o $(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.o

$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_Jacobian.o : $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_Jacobian.cpp $(PATH_HPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.h
		 	 $(CCP) $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_Jacobian.o

$(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_ReactionRates.o : $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_ReactionRates.cpp $(PATH_HPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.h
		 	 $(CCP) $(PATH_CPP)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO_ReactionRates.o

#Polimi_c1c3htnox_avio_0702
$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.o : $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.cpp $(PATH_HPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.o

$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_Jacobian.o : $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_Jacobian.cpp $(PATH_HPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_Jacobian.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_Jacobian.o

$(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_ReactionRates.o : $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_ReactionRates.cpp $(PATH_HPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.h
		 	 $(CCP) $(PATH_CPP)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_ReactionRates.cpp -o $(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702_ReactionRates.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o $(PATH_OBJECTS)/sandiego_avio/*.o $(PATH_OBJECTS)/gri30/*.o $(PATH_OBJECTS)/gri12/*.o $(PATH_OBJECTS)/fluent_drm22_polimi_nox/*.o $(PATH_OBJECTS)/polimi_h2conox_1101/*.o $(PATH_OBJECTS)/polimi_h2conox_nothermal_1101/*.o $(PATH_OBJECTS)/polimi_h2conox_nonnh_1101/*.o $(PATH_OBJECTS)/polimi_h2conox_non2o_1101/*.o $(PATH_OBJECTS)/polimi_c1c3htnox_1101/*.o $(PATH_OBJECTS)/polimi_nc7_avio/*.o $(PATH_OBJECTS)/polimi_c1c3htnox_avio_0702/*.o  ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_SymbolicKinetics_$(STRING_NAME).a

