#!/bin/bash

#python morpholopy.py -d  /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe1_ZdepBoost_04Zsun \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu2_Fe05_ZdepBoost_04Zsun \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu2_Fe1_ZdepBoost_04Zsun \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_04Zsun \
#                      -n nu1p6Fe1Zdep nu2Fe05Zdep nu2Fe1Zdep nu1p6Fe05Zdep \
#                      -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                      -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                      -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ComparisonZdepDTDnuFe

python morpholopy.py -d  /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu2_Fe05_ZdepBoost_04Zsun \
                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu2_Fe05O07Mg07 \
                      -n nu2Fe05ZdepBoost nu1p6Fe1O1Mg1 nu2Fe05O07Mg07 \
                      -c halo_0036.properties halo_0036.properties halo_0036.properties \
                      -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
                      -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/L25N188ComparisonABCModels

#python3 morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin10hi40_DTDExp2nu1p6_vBF \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin10hi40_DTDExp2nu1p6_Fe05O07Mg07 \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                         /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05O07Mg07 \
#                      -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/MinSNMassTest \
#                      -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                      -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                      -n SN10to40Fe1Mg1O1 SN10to40Fe05Mg07O07 SN8to40Fe1Mg1O1 SN8to40Fe05Mg07O07 \
#                      -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L050N376_NewYields_SNIILin8hi40_DTDExp2nu1p6_ZdepB \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L050N376_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ConvergenceTest2 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n L50DTDtau2nu1p6Zdep L50DTDtau2nu1p6vBF \
#                     -m 5e9

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L050N376_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ConvergenceTest1 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n L50DTDtau2nu1p6 L25DTDtau2nu1p6 \
#                     -m 1e9

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L050N376_NewYields_SNIILin8hi40_DTDExp2nu1p6_ZdepB \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_04Zsun \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ConvergenceTest2 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n L50DTDtau2nu1p6Zdep L25DTDtau2nu1p6Zdep \
#                     -m 1e9

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L012N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L012N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBFO07Mg07 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L012N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_ZdepB \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/L12Comparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe1O1Mg1 DTDtau2nu1p6Fe05O07Mg07 DTDtau2nu1p6Fe05Zdep

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_08Zsun \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_04Zsun \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/MetalDepBoostFactorsComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe05 DTDtau2nu1p6Fe05O05Mg05 DTDtau2nu1p6Fe05Zdep08Zsun DTDtau2nu1p6Fe05Zdep04Zsun \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/L025SNII8to40Tau2nu1p6vBF \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNII8to40Tau2nu1p6Fe05vBF SNII8to40Tau2nu1p6Fe1vBF \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe025 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe075 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe1 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestConstBoostFactorTau2nu1p6 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe025 DTDtau2nu1p6Fe05 DTDtau2nu1p6Fe075 DTDtau2nu1p6Fe1 \
#                     -m 7e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05Zhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalDepBoostFactorZTau2nu2 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu2Fe05Z0p4Zsun DTDtau2nu2Fe05Z0p8Zsun DTDtau2nu2Fe05Z1p5Zsun DTDtau2nu2Fe05 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05M10to40 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05M20to40 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05M30to40 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMassDepBoostFactorFe05 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu1p6Fe05M10to40 DTDExptau2nu1p6Fe05M20to40 DTDExptau2nu1p6Fe05M30to40 DTDExptau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe4Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe4Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetallicityDepBoostFactor2 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu1p6Fe4Zhigh1 DTDExptau2nu1p6Fe4Zhigh2 DTDExptau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05FeZhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05FeZhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05FeZhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalDepBoostFactorZhigh1 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu2Fe05Fe4Z0p4Zsun DTDtau2nu2Fe05Fe4Z0p6Zsun DTDtau2nu2Fe05Fe4Z0p8Zsun DTDtau2nu2Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2Zhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalDepBoostFactorZhigh3 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe2Z0p4Zsun DTDtau2nu1p6Fe2Z0p6Zsun DTDtau2nu1p6Fe2Z0p8Zsun DTDtau2nu1p6 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05O05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_O15Mg15 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/DTDExpTau2nu2Comparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu2 DTDtau2nu2Fe05Si05C15 DTDExptau2nu2Fe05O05 DTDExptau2nu2Mg15O15 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_O05Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1_O05Fe05Si05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau4nu1_O05Fe05Si05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLawBeta05nu1O05Fe05 DTDExptau1p5nu1p6O05Fe05 DTDExptau2nu1O05Fe05 DTDExptau4nu1O05Fe05 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau1nu2 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_new_plots/L025NYSampleA \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNII8to40Tau2nu1p6 SNII10to40Tau2nu1p6 SNII8to40Tau2nu2 SNII8to40Tau1nu2 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_new_plots/L025NYSampleD \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n PowerLawBeta08 PowerLawBeta08Fe05 Tau1nu1p5Fe05O05Mg05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi100_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIMassComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNIIlo8hi40 SNIIlo10hi40 SNIIlo8hi100 SNIIlo10hi100 \
#                     -m 5e8
#
#
#

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau4nu1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau4nu1_O05Fe05Si05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDtau4nu1 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau4nu1 DTDtau4nu1O05Fe05Si05 \
#                     -m 5e8
#
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L012N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -n L025SNIILinlo8hi100Tau2Nu1p6 L012SNIILinlo8hi100Tau2Nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/Convergence_NewYields_SNIILinlo8hi100 \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L012N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L012N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -n SNIILinlo8hi40Tau2Nu1p6 SNIILinlo8hi100Tau2Nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L012N188_NewYields_SNIILin_Mass_Range_comparison \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L012N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties -n L12SNIIlo8hi40DTDtau2nu1p6 L25SNIIlo8hi40DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/NewYieldsSNIILinlo8hi40_DTDtau2nu1p6 \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n SNIILinlo8hi40DTDPowerLawBeta05nu1 SNIILinlo8hi40DTDPowerLawBeta05nu1Fe05 SNIILinlo8hi40DTDtau1p5nu1p6Fe05O05Mg05 SNIILinlo8hi40DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L25NYSNIILinlo8hi40DTDComparison \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -n SNIILinlo8hi100Tau2Nu1p6 SNIILinlo8hi100Tau2Nu1p6O05Mg05Si05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_OMgSi_boost \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05_Fe15C15 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05_Fe2C2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n O05Mg05Si05 O05Mg05Si05Fe15C15 O05Mg05Si06Fe2C2 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_new_plots/L025NYSampleE \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_O05Mg05Si05 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -n SNIILinlo8hi100Tau2Nu1p6Fe13 SNIILinlo8hi100Tau2Nu1p6Fe16 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/Test_boost \
#                     -m 1e9

#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_Fe15 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6_Fe2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n SNIILinlo8hi100Tau2Nu1p6 SNIILinlo8hi100Tau2Nu1p6Fe15 SNIILinlo8hi100Tau2Nu1p6Fe2 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/Test_boost \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NewYieldsSNIILinlo8hi100Tau2Nu1p2 NewYieldsSNIILinlo8hi100Tau2Nu1p4 NewYieldsSNIILinlo8hi100Tau2Nu1p6 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/abundance_plots/L025N188_NewYields_SNIILinlo8hi100_DTDnu_comparison \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau1p5nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau3nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau4nu1p4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NewYieldsSNIILinlo8hi100Tau1p5Nu1p4 NewYieldsSNIILinlo8hi100Tau2Nu1p4 NewYieldsSNIILinlo8hi100Tau3Nu1p4 NewYieldsSNIILinlo8hi100Tau4Nu1p4 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/abundance_plots/L025N188_NewYields_SNIILinlo8hi100_DTDtau_comparison \
#                     -m 1e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NewYieldsSNIILinlo8hi100Tau2Nu1p2 NewYieldsSNIILinlo8hi100Tau2Nu1p4 NewYieldsSNIILinlo8hi100Tau2Nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L25N188NY/SNIILin8100DTDnuComparison \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau1p5nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau3nu1p4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau4nu1p4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NewYieldsSNIILinlo8hi100Tau1p5Nu1p4 NewYieldsSNIILinlo8hi100Tau2Nu1p4 NewYieldsSNIILinlo8hi100Tau3Nu1p4 NewYieldsSNIILinlo8hi100Tau4Nu1p4 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L25N188NY/SNIILin8100DTDtauComparison \
#                     -m 5e8

#python morpholopy.py -d /cosma6/data/dp004/Eagle/wmfw23/snap_bkp/newchem/volumes/chem_6s/L06N094_SNIIlin_lo_8_hi_40_SNIaNu_1p6_tau_3p5/proc/ \
#                        /cosma6/data/dp004/Eagle/wmfw23/snap_bkp/newchem/volumes/chem_6s/L06N094_SNIIlin_lo_8_hi_70_SNIaNu_1p6_tau_3p5/proc/ \
#                        /cosma6/data/dp004/Eagle/wmfw23/snap_bkp/newchem/volumes/chem_6s/L06N094_SNIIlin_lo_8_hi_100_SNIaNu_1p6_tau_3p5/proc/ \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/Morpholopy_plots/JamesProduction/L06_SNII_limit_comparison \
#                     -c velociraptor_0034.properties.0 velociraptor_0034.properties.0 velociraptor_0034.properties.0 \
#                     -s colibre_0034.hdf5 colibre_0034.hdf5 colibre_0034.hdf5 \
#                     -n KobayashiYieldslo8hi40 KobayashiYieldslo8hi70 KobayashiYieldslo8hi100

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoostTest_newSNIIform \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/Morpholopy_plots/L012NewSNIIFormComparison \
#                     -c halo_0036.properties.0 halo_0036.properties.0 halo_0036.properties.0 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n NomotoYieldsO05Fe1 NomotoYieldsO05Fe05 KobayashiNewFormO05Fe05

#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/prediction2_full_noZ \
#                     /cosma7/data/dp004/dc-chai1/HAWK/prediction1 \
#                     -s colibre_0023.hdf5 colibre_0023.hdf5 -c halo_0023.properties halo_0023.properties \
#                      -n wave2ColibreDTDFe1Mg4NoZDependence wave1ColibreDTDFe1Mg4 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/BestWaveRuns/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_Test_Kobayashi_yields_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_Test_Kobayashi_yields_2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                    -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n ThielemannSNIapowerlaw TestKobayashiSNIaEAGLEDTD TestKobayashiSNIapowerlaw \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/TestKobayashi/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_01_COLIBRE/L012N188_SNIIBoost_4 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/NomotoTestBoostFactors \
#                     -c halo_0036.properties.0 halo_0036.properties.0 halo_0036.properties.0 halo_0036.properties.0 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNIaEagleDTDO05Fe1 SNIaEagleDTDO05Fe05 SNIaPowerLawO2Mg2Fe05 SNIaPowerLawO15Mg15Fe05 \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_2 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_3 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n SNIIYieldsMmin13Mmax40 SNIIYieldsMmin10Mmax40 SNIIYieldsMmin7Mmax40 SNIIYieldsMmin13Mmax100 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/NomotoYieldsMassLimitTest \
#                     -m 1e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0033.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0033.properties halo_0036.properties \
#                     -n SNIINewYieldsSNIapowerlaw05 SNIINewYieldsSNIapowerlaw15 SNIINewYieldsSNIapowerlaw08 SNIINewYieldsSNIapowerlaw12 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/NomotoYieldsSNIaPowerLaw

#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/NomotoYieldsSNIaPowerLaw/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_NewYields_SNIaEAGLEDTD \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_NewYields_SNIapowerlaw_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_NewYields_SNIapowerlaw_2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NewYieldsSNIaEAGLEDTD NewYieldsSNIaPowerLaw05 NewYieldsSNIaPowerLaw15 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/L12TestSNIa/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_Test_C_Ba_Eu_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_Test_C_Ba_Eu_2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n TestCBaEu1 TestCBaEu2 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/TestCBaEu/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_SNIa_power_law_2 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n SNIINewYieldsSNIapowerlaw05 SNIINewYieldsSNIapowerlaw15 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/NomotoYieldsSNIaPowerLaw/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_2 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_3 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/New_yields_Model_limit_Mass_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n SNIIYieldsMmin13Mmax40 SNIIYieldsMmin10Mmax40 SNIIYieldsMmin7Mmax40 SNIIYieldsMmin13Mmax100 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/NomotoYieldsMassLimitTest/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_Model_1 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_Model_2 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_Model_3 \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L006N094_New_yields_Model_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n NomotoYieldsMmin10Mmax40Fe15 NomotoYieldsMmin10Mmax40Fe2 NomotoYieldsMmin10Mmax40SNIa1Gyr NomotoYieldsMmin10Mmax40SNIa05Gyr \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/NomotoYields/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_1 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_2 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_3 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n L12SNIatau1Gyr L12SNIatau2Gyr L12SNIatau3Gyr L12SNIatau4Gyr \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n L12newAGB L12oldyields \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#

#python morpholopy.py -d ../../ \
#                     -s colibre_0034.hdf5 -c halo_0034.properties -n test \
#                     -o ../../Morpholopy_Plots/

#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties \
#                     -n L12newAGB \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n L12newAGB L12oldyields \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d ../../ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n DTDExptau2nu1p6Fe05Z0p4Zsun \
#                     -o ../../Morpholopy_Plots/ \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/123b2_121b2_SNII_Mg_2_Fe_1_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 123b2121b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/116b2_115b2_Fe1_new_cooling_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 116b2115b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/106b3_104b2_norm_0p3_new_cooling_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 106b3104b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/122b2_121b2_SNIa_best_fit_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 122b2121b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_boost_1/ \
#                     -s colibre_0034.hdf5 -c halo_0034.properties -n SNIaFe05Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_boost_2/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIaFe025Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_11_COLIBRE/EAGLE_SNIa_diffusion_1/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIaFe1Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_2/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe1Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_5/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe05Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_6/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe025Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa2GyrFe05Mg2L12N188 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_025_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa2GyrFe025Mg2L12N188 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/


