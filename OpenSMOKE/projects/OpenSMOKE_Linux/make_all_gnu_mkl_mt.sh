make -f OpenSMOKE_AddOns_GNU_MKL_MT clean
make -f OpenSMOKE_AddOns_GNU_MKL_MT

make -f OpenSMOKE_Basic_GNU_MKL_MT clean
make -f OpenSMOKE_Basic_GNU_MKL_MT

make -f OpenSMOKE_BatchReactor_GNU_MKL_MT clean
make -f OpenSMOKE_BatchReactor_GNU_MKL_MT

make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL_MT clean
make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL_MT

make -f OpenSMOKE_Distributions_GNU_MKL_MT clean
make -f OpenSMOKE_Distributions_GNU_MKL_MT

make -f OpenSMOKE_Engine_GNU_MKL_MT clean
make -f OpenSMOKE_Engine_GNU_MKL_MT

make -f OpenSMOKE_Droplet_GNU_MKL_MT clean
make -f OpenSMOKE_Droplet_GNU_MKL_MT

make -f OpenSMOKE_Flame1D_GNU_MKL_MT clean
make -f OpenSMOKE_Flame1D_GNU_MKL_MT
 
make -f OpenSMOKE_Flamelet_GNU_MKL_MT clean
make -f OpenSMOKE_Flamelet_GNU_MKL_MT

make -f OpenSMOKE_KinPP_GNU_MKL_MT clean
make -f OpenSMOKE_KinPP_GNU_MKL_MT

make -f OpenSMOKE_IdealReactors_GNU_MKL_MT clean
make -f OpenSMOKE_IdealReactors_GNU_MKL_MT

make -f OpenSMOKE_Interfaces_GNU_MKL_MT clean
make -f OpenSMOKE_Interfaces_GNU_MKL_MT

make -f OpenSMOKE_MixManager_GNU_MKL_MT clean
make -f OpenSMOKE_MixManager_GNU_MKL_MT

make -f OpenSMOKE_PlugFlowReactor_GNU_MKL_MT clean
make -f OpenSMOKE_PlugFlowReactor_GNU_MKL_MT

make -f OpenSMOKE_CSTR_GNU_MKL_MT clean
make -f OpenSMOKE_CSTR_GNU_MKL_MT

make -f OpenSMOKE_PreProcessing_GNU_MKL_MT clean
make -f OpenSMOKE_PreProcessing_GNU_MKL_MT

make -f OpenSMOKE_QMOM_GNU_MKL_MT clean
make -f OpenSMOKE_QMOM_GNU_MKL_MT

make -f OpenSMOKE_ShockTube_GNU_MKL_MT clean
make -f OpenSMOKE_ShockTube_GNU_MKL_MT

make -f OpenSMOKE_LiquidProperties_GNU_MKL_MT clean
make -f OpenSMOKE_LiquidProperties_GNU_MKL_MT

make -f OpenSMOKE_ICEM_GNU_MKL_MT clean
make -f OpenSMOKE_ICEM_GNU_MKL_MT

make -f OpenSMOKE_SurfaceChemistry_GNU_MKL_MT clean
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL_MT

make -f OpenSMOKE_SurfaceChemistry_GNU_MKL_MT clean
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL_MT 

make -f OpenSMOKE_SymbolicKinetics_GNU_MKL_MT clean
make -f OpenSMOKE_SymbolicKinetics_GNU_MKL_MT

ar rc ../../lib/linux/GNU_MKL_MT/libOpenSMOKE_PostProcessor_GNU_MKL_MT.a  OpenSMOKE_AddOns/GNU_MKL_MT/*.o OpenSMOKE_Flame1D/GNU_MKL_MT/*.o OpenSMOKE_Basic/GNU_MKL_MT/*.o OpenSMOKE_Engine/GNU_MKL_MT/*.o OpenSMOKE_Engine/GNU_MKL_MT/*.o
