mkdir ../../lib
mkdir ../../lib/linux
mkdir ../../lib/linux/gnu

mkdir OpenSMOKE_AddOns
mkdir OpenSMOKE_AddOns/gnu
make -f OpenSMOKE_AddOns_GNU clean
make -f OpenSMOKE_AddOns_GNU

mkdir OpenSMOKE_Basic
mkdir OpenSMOKE_Basic/gnu
make -f OpenSMOKE_Basic_GNU clean
make -f OpenSMOKE_Basic_GNU

mkdir OpenSMOKE_BatchReactor
mkdir OpenSMOKE_BatchReactor/gnu
make -f OpenSMOKE_BatchReactor_GNU clean
make -f OpenSMOKE_BatchReactor_GNU

mkdir OpenSMOKE_CHEMKINInterpreter
mkdir OpenSMOKE_CHEMKINInterpreter/gnu
make -f OpenSMOKE_CHEMKINInterpreter_GNU clean
make -f OpenSMOKE_CHEMKINInterpreter_GNU

mkdir OpenSMOKE_Distributions
mkdir OpenSMOKE_Distributions/gnu
make -f OpenSMOKE_Distributions_GNU clean
make -f OpenSMOKE_Distributions_GNU

mkdir OpenSMOKE_Engine
mkdir OpenSMOKE_Engine/gnu
make -f OpenSMOKE_Engine_GNU clean
make -f OpenSMOKE_Engine_GNU

mkdir OpenSMOKE_Droplet
mkdir OpenSMOKE_Droplet/gnu
make -f OpenSMOKE_Droplet_GNU clean
make -f OpenSMOKE_Droplet_GNU

mkdir OpenSMOKE_Flame1D
mkdir OpenSMOKE_Flame1D/gnu
make -f OpenSMOKE_Flame1D_GNU clean
make -f OpenSMOKE_Flame1D_GNU
 
mkdir OpenSMOKE_Flamelet
mkdir OpenSMOKE_Flamelet/gnu
make -f OpenSMOKE_Flamelet_GNU clean
make -f OpenSMOKE_Flamelet_GNU

mkdir OpenSMOKE_KinPP
mkdir OpenSMOKE_KinPP/gnu
make -f OpenSMOKE_KinPP_GNU clean
make -f OpenSMOKE_KinPP_GNU

mkdir OpenSMOKE_IdealReactors
mkdir OpenSMOKE_IdealReactors/gnu
make -f OpenSMOKE_IdealReactors_GNU clean
make -f OpenSMOKE_IdealReactors_GNU

mkdir OpenSMOKE_Interfaces
mkdir OpenSMOKE_Interfaces/gnu
make -f OpenSMOKE_Interfaces_GNU clean
make -f OpenSMOKE_Interfaces_GNU

mkdir OpenSMOKE_MixManager
mkdir OpenSMOKE_MixManager/gnu
make -f OpenSMOKE_MixManager_GNU clean
make -f OpenSMOKE_MixManager_GNU

mkdir OpenSMOKE_PlugFlowReactor
mkdir OpenSMOKE_PlugFlowReactor/gnu
make -f OpenSMOKE_PlugFlowReactor_GNU clean
make -f OpenSMOKE_PlugFlowReactor_GNU

mkdir OpenSMOKE_CSTR
mkdir OpenSMOKE_CSTR/gnu
make -f OpenSMOKE_CSTR_GNU clean
make -f OpenSMOKE_CSTR_GNU

mkdir OpenSMOKE_PreProcessing
mkdir OpenSMOKE_PreProcessing/gnu
make -f OpenSMOKE_PreProcessing_GNU clean
make -f OpenSMOKE_PreProcessing_GNU

mkdir OpenSMOKE_QMOM
mkdir OpenSMOKE_QMOM/gnu
make -f OpenSMOKE_QMOM_GNU clean
make -f OpenSMOKE_QMOM_GNU

mkdir OpenSMOKE_ShockTube
mkdir OpenSMOKE_ShockTube/gnu
make -f OpenSMOKE_ShockTube_GNU clean
make -f OpenSMOKE_ShockTube_GNU

mkdir OpenSMOKE_LiquidProperties
mkdir OpenSMOKE_LiquidProperties/gnu
make -f OpenSMOKE_LiquidProperties_GNU clean
make -f OpenSMOKE_LiquidProperties_GNU

mkdir OpenSMOKE_ICEM
mkdir OpenSMOKE_ICEM/gnu
make -f OpenSMOKE_ICEM_GNU clean
make -f OpenSMOKE_ICEM_GNU

mkdir OpenSMOKE_SurfaceChemistry
mkdir OpenSMOKE_SurfaceChemistry/gnu
make -f OpenSMOKE_SurfaceChemistry_GNU clean
make -f OpenSMOKE_SurfaceChemistry_GNU

mkdir OpenSMOKE_SymbolicKinetics
mkdir OpenSMOKE_SymbolicKinetics/gnu
mkdir OpenSMOKE_SymbolicKinetics/gnu/sandiego_avio
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_nc7_avio
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_h2conox_nothermal_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_h2conox_nonnh_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_h2conox_non2o_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_h2conox_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_c1c3htnox_avio_0702
mkdir OpenSMOKE_SymbolicKinetics/gnu/polimi_c1c3htnox_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu/gri30
mkdir OpenSMOKE_SymbolicKinetics/gnu/gri12
mkdir OpenSMOKE_SymbolicKinetics/gnu/fluent_drm22_polimi_nox
make -f OpenSMOKE_SymbolicKinetics_GNU clean
make -f OpenSMOKE_SymbolicKinetics_GNU


ar rc ../../lib/linux/gnu/libOpenSMOKE_PostProcessor_GNU.a  OpenSMOKE_AddOns/gnu/*.o OpenSMOKE_Flame1D/gnu/*.o OpenSMOKE_Basic/gnu/*.o OpenSMOKE_Engine/gnu/*.o OpenSMOKE_Engine/gnu/*.o
