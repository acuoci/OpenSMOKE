make -f OpenSMOKE_AddOns_GNU clean
make -f OpenSMOKE_AddOns_GNU

make -f OpenSMOKE_Basic_GNU clean
make -f OpenSMOKE_Basic_GNU

make -f OpenSMOKE_BatchReactor_GNU clean
make -f OpenSMOKE_BatchReactor_GNU

make -f OpenSMOKE_CHEMKINInterpreter_GNU clean
make -f OpenSMOKE_CHEMKINInterpreter_GNU

make -f OpenSMOKE_Distributions_GNU clean
make -f OpenSMOKE_Distributions_GNU

make -f OpenSMOKE_Engine_GNU clean
make -f OpenSMOKE_Engine_GNU

make -f OpenSMOKE_Droplet_GNU clean
make -f OpenSMOKE_Droplet_GNU

make -f OpenSMOKE_Flame1D_GNU clean
make -f OpenSMOKE_Flame1D_GNU
 
make -f OpenSMOKE_Flamelet_GNU clean
make -f OpenSMOKE_Flamelet_GNU

make -f OpenSMOKE_KinPP_GNU clean
make -f OpenSMOKE_KinPP_GNU

make -f OpenSMOKE_IdealReactors_GNU clean
make -f OpenSMOKE_IdealReactors_GNU

make -f OpenSMOKE_Interfaces_GNU clean
make -f OpenSMOKE_Interfaces_GNU

make -f OpenSMOKE_MixManager_GNU clean
make -f OpenSMOKE_MixManager_GNU

make -f OpenSMOKE_PlugFlowReactor_GNU clean
make -f OpenSMOKE_PlugFlowReactor_GNU

make -f OpenSMOKE_CSTR_GNU clean
make -f OpenSMOKE_CSTR_GNU

make -f OpenSMOKE_PreProcessing_GNU clean
make -f OpenSMOKE_PreProcessing_GNU

make -f OpenSMOKE_QMOM_GNU clean
make -f OpenSMOKE_QMOM_GNU

make -f OpenSMOKE_ShockTube_GNU clean
make -f OpenSMOKE_ShockTube_GNU

make -f OpenSMOKE_LiquidProperties_GNU clean
make -f OpenSMOKE_LiquidProperties_GNU

make -f OpenSMOKE_ICEM_GNU clean
make -f OpenSMOKE_ICEM_GNU

make -f OpenSMOKE_SurfaceChemistry_GNU clean
make -f OpenSMOKE_SurfaceChemistry_GNU

make -f OpenSMOKE_SymbolicKinetics_GNU clean
make -f OpenSMOKE_SymbolicKinetics_GNU

make -f OpenSMOKE_SurfaceChemistry_GNU clean
make -f OpenSMOKE_SurfaceChemistry_GNU 

ar rc ../../lib/linux/gnu/libOpenSMOKE_PostProcessor_GNU.a  OpenSMOKE_AddOns/gnu/*.o OpenSMOKE_Flame1D/gnu/*.o OpenSMOKE_Basic/gnu/*.o OpenSMOKE_Engine/gnu/*.o OpenSMOKE_Engine/gnu/*.o
