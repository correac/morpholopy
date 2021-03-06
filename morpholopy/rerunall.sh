#!/bin/bash

# Runs MorpholoPy on a bunch of simulations:

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01_1SDT_ESN3e51"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01_1SDT_ESN3e51/output_0823"
snap="0823"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L006N376_FKIN00_NOEOS_VKICK000SLOPE00NORM00"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/XMAS2020_L006N376_FKIN00_NOEOS_VKICK000SLOPE00NORM00/output_2729"
snap="2729"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/AGN5_L006188_00_UVB_dust1_CR0_G0_shield1"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/AGN5_L006188_00_UVB_dust1_CR0_G0_shield1/output_2729"
snap="2729"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/38_AGNdT_3e8_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/38_AGNdT_3e8_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01_1SDT_ESN3e51"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01_1SDT_ESN3e51/output_0623"
snap="0623"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/AGN5_L006188_00_UVB_dust1_CR0_G0_shield1"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/AGN5_L006188_00_UVB_dust1_CR0_G0_shield1/output_2729"
snap="2729"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/01_reference_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/15_alphavir_0p5_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/15_alphavir_0p5_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/39_only_UVB_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/39_only_UVB_L006N0188/output_016"
snap="16"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/03_SNE_4e51_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/03_SNE_4e51_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/05_delay_40Myr_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/05_delay_40Myr_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/09_steepEOS_1e3_pressurelaw_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/09_steepEOS_1e3_pressurelaw_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/11_fkin_0p1_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/11_fkin_0p1_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/33_BHboost_densdep_beta0p5_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/33_BHboost_densdep_beta0p5_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/36_SF_efficiency_0p001_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/36_SF_efficiency_0p001_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/35_SF_efficiency_0p003_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/35_SF_efficiency_0p003_L006N0188/output_023/"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/29_SNE_densdep_2_4_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/29_SNE_densdep_2_4_L006N0188/output_023"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/31_SNE_metdep_2_4_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/31_SNE_metdep_2_4_L006N0188/output_023/"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/37_same_gr_softening_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/37_same_gr_softening_L006N0188/output_023/"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output

folder="/cosma7/data/dp004/dc-chai1/HAWK/40_alphavir_0p25_L006N0188"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/40_alphavir_0p25_L006N0188/output_023/"
snap="23"
python morpholopy.py -d=$folder -n=$snap -o=$output


