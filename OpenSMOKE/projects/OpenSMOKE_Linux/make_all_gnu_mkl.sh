mkdir ../../lib
mkdir ../../lib/linux
mkdir ../../lib/linux/gnu_mkl

mkdir OpenSMOKE_AddOns
mkdir OpenSMOKE_AddOns/gnu_mkl
make -f OpenSMOKE_AddOns_GNU_MKL clean
make -f OpenSMOKE_AddOns_GNU_MKL

mkdir OpenSMOKE_Basic
mkdir OpenSMOKE_Basic/gnu_mkl
make -f OpenSMOKE_Basic_GNU_MKL clean
make -f OpenSMOKE_Basic_GNU_MKL

mkdir OpenSMOKE_BatchReactor
mkdir OpenSMOKE_BatchReactor/gnu_mkl
make -f OpenSMOKE_BatchReactor_GNU_MKL clean
make -f OpenSMOKE_BatchReactor_GNU_MKL

mkdir OpenSMOKE_CHEMKINInterpreter
mkdir OpenSMOKE_CHEMKINInterpreter/gnu_mkl
make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL clean
make -f OpenSMOKE_CHEMKINInterpreter_GNU_MKL

mkdir OpenSMOKE_Distributions
mkdir OpenSMOKE_Distributions/gnu_mkl
make -f OpenSMOKE_Distributions_GNU_MKL clean
make -f OpenSMOKE_Distributions_GNU_MKL

mkdir OpenSMOKE_Engine
mkdir OpenSMOKE_Engine/gnu_mkl
make -f OpenSMOKE_Engine_GNU_MKL clean
make -f OpenSMOKE_Engine_GNU_MKL

mkdir OpenSMOKE_Droplet
mkdir OpenSMOKE_Droplet/gnu_mkl
make -f OpenSMOKE_Droplet_GNU_MKL clean
make -f OpenSMOKE_Droplet_GNU_MKL

mkdir OpenSMOKE_Flame1D
mkdir OpenSMOKE_Flame1D/gnu_mkl
make -f OpenSMOKE_Flame1D_GNU_MKL clean
make -f OpenSMOKE_Flame1D_GNU_MKL
 
mkdir OpenSMOKE_Flamelet
mkdir OpenSMOKE_Flamelet/gnu_mkl
make -f OpenSMOKE_Flamelet_GNU_MKL clean
make -f OpenSMOKE_Flamelet_GNU_MKL

mkdir OpenSMOKE_KinPP
mkdir OpenSMOKE_KinPP/gnu_mkl
make -f OpenSMOKE_KinPP_GNU_MKL clean
make -f OpenSMOKE_KinPP_GNU_MKL

mkdir OpenSMOKE_IdealReactors
mkdir OpenSMOKE_IdealReactors/gnu_mkl
make -f OpenSMOKE_IdealReactors_GNU_MKL clean
make -f OpenSMOKE_IdealReactors_GNU_MKL

mkdir OpenSMOKE_Interfaces
mkdir OpenSMOKE_Interfaces/gnu_mkl
make -f OpenSMOKE_Interfaces_GNU_MKL clean
make -f OpenSMOKE_Interfaces_GNU_MKL

mkdir OpenSMOKE_MixManager
mkdir OpenSMOKE_MixManager/gnu_mkl
make -f OpenSMOKE_MixManager_GNU_MKL clean
make -f OpenSMOKE_MixManager_GNU_MKL

mkdir OpenSMOKE_PlugFlowReactor
mkdir OpenSMOKE_PlugFlowReactor/gnu_mkl
make -f OpenSMOKE_PlugFlowReactor_GNU_MKL clean
make -f OpenSMOKE_PlugFlowReactor_GNU_MKL

mkdir OpenSMOKE_CSTR
mkdir OpenSMOKE_CSTR/gnu_mkl
make -f OpenSMOKE_CSTR_GNU_MKL clean
make -f OpenSMOKE_CSTR_GNU_MKL

mkdir OpenSMOKE_PreProcessing
mkdir OpenSMOKE_PreProcessing/gnu_mkl
make -f OpenSMOKE_PreProcessing_GNU_MKL clean
make -f OpenSMOKE_PreProcessing_GNU_MKL

mkdir OpenSMOKE_QMOM
mkdir OpenSMOKE_QMOM/gnu_mkl
make -f OpenSMOKE_QMOM_GNU_MKL clean
make -f OpenSMOKE_QMOM_GNU_MKL

mkdir OpenSMOKE_ShockTube
mkdir OpenSMOKE_ShockTube/gnu_mkl
make -f OpenSMOKE_ShockTube_GNU_MKL clean
make -f OpenSMOKE_ShockTube_GNU_MKL

mkdir OpenSMOKE_LiquidProperties
mkdir OpenSMOKE_LiquidProperties/gnu_mkl
make -f OpenSMOKE_LiquidProperties_GNU_MKL clean
make -f OpenSMOKE_LiquidProperties_GNU_MKL

mkdir OpenSMOKE_ICEM
mkdir OpenSMOKE_ICEM/gnu_mkl
make -f OpenSMOKE_ICEM_GNU_MKL clean
make -f OpenSMOKE_ICEM_GNU_MKL

mkdir OpenSMOKE_SurfaceChemistry
mkdir OpenSMOKE_SurfaceChemistry/gnu_mkl
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL clean
make -f OpenSMOKE_SurfaceChemistry_GNU_MKL

mkdir OpenSMOKE_SymbolicKinetics
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/sandiego_avio
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_nc7_avio
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_h2conox_nothermal_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_h2conox_nonnh_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_h2conox_non2o_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_h2conox_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_c1c3htnox_avio_0702
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/polimi_c1c3htnox_1101
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/gri30
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/gri12
mkdir OpenSMOKE_SymbolicKinetics/gnu_mkl/fluent_drm22_polimi_nox
make -f OpenSMOKE_SymbolicKinetics_GNU_MKL clean
make -f OpenSMOKE_SymbolicKinetics_GNU_MKL


ar rc ../../lib/linux/gnu_mkl/libOpenSMOKE_PostProcessor_GNU_MKL.a  OpenSMOKE_AddOns/gnu_mkl/*.o OpenSMOKE_Flame1D/gnu_mkl/*.o OpenSMOKE_Basic/gnu_mkl/*.o OpenSMOKE_Engine/gnu_mkl/*.o OpenSMOKE_Engine/gnu_mkl/*.o
