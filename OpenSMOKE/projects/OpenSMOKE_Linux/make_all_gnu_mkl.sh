make -f OpenSMOKE_AddOns_GNU_MKL clean
make -f OpenSMOKE_AddOns_GNU_MKL

make -f OpenSMOKE_Basic_GNU_MKL clean
make -f OpenSMOKE_Basic_GNU_MKL

make -f OpenSMOKE_BatchReactor_GNU_MKL clean
make -f OpenSMOKE_BatchReactor_GNU_MKL

make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL clean
make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL

make -f OpenSMOKE_Distributions_GNU_MKL clean
make -f OpenSMOKE_Distributions_GNU_MKL

make -f OpenSMOKE_Engine_GNU_MKL clean
make -f OpenSMOKE_Engine_GNU_MKL

make -f OpenSMOKE_Droplet_GNU_MKL clean
make -f OpenSMOKE_Droplet_GNU_MKL

make -f OpenSMOKE_Flame1D_GNU_MKL clean
make -f OpenSMOKE_Flame1D_GNU_MKL
 
make -f OpenSMOKE_Flamelet_GNU_MKL clean
make -f OpenSMOKE_Flamelet_GNU_MKL

make -f OpenSMOKE_KinPP_GNU_MKL clean
make -f OpenSMOKE_KinPP_GNU_MKL

make -f OpenSMOKE_IdealReactors_GNU_MKL clean
make -f OpenSMOKE_IdealReactors_GNU_MKL

make -f OpenSMOKE_Interfaces_GNU_MKL clean
make -f OpenSMOKE_Interfaces_GNU_MKL

make -f OpenSMOKE_MixManager_GNU_MKL clean
make -f OpenSMOKE_MixManager_GNU_MKL

make -f OpenSMOKE_PlugFlowReactor_GNU_MKL clean
make -f OpenSMOKE_PlugFlowReactor_GNU_MKL

make -f OpenSMOKE_CSTR_GNU_MKL clean
make -f OpenSMOKE_CSTR_GNU_MKL

make -f OpenSMOKE_PreProcessing_GNU_MKL clean
make -f OpenSMOKE_PreProcessing_GNU_MKL

make -f OpenSMOKE_QMOM_GNU_MKL clean
make -f OpenSMOKE_QMOM_GNU_MKL

make -f OpenSMOKE_ShockTube_GNU_MKL clean
make -f OpenSMOKE_ShockTube_GNU_MKL

make -f OpenSMOKE_LiquidProperties_GNU_MKL clean
make -f OpenSMOKE_LiquidProperties_GNU_MKL

make -f OpenSMOKE_ICEM_GNU_MKL clean
make -f OpenSMOKE_ICEM_GNU_MKL

make -f OpenSMOKE_SurfaceChemistry_GNU_MKL clean
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL

make -f OpenSMOKE_SurfaceChemistry_GNU_MKL clean
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL 

make -f OpenSMOKE_SymbolicKinetics_GNU_MKL clean
make -f OpenSMOKE_SymbolicKinetics_GNU_MKL

ar rc ../../lib/linux/gnu_mkl/libOpenSMOKE_PostProcessor_GNU_MKL.a  OpenSMOKE_AddOns/gnu_mkl/*.o OpenSMOKE_Flame1D/gnu_mkl/*.o OpenSMOKE_Basic/gnu_mkl/*.o OpenSMOKE_Engine/gnu_mkl/*.o OpenSMOKE_Engine/gnu_mkl/*.o
