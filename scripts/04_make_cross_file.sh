cd ~/projects/NLB_MDH
# make cross file for RILs
Rscript scripts/convert_rqtl2_to_rqtl.R data/RIL.csv analyses/IBM_NLB_BLUPS_RILs.csv data/ibm302map.csv analyses/RIL_cross.csv
# make cross file for B73 crosses
Rscript scripts/convert_rqtl2_to_rqtl.R data/B73_hybrids.csv analyses/IBM_NLB_BLUPS_B73BC.csv data/ibm302map.csv analyses/B73_cross.csv  
# make cross file for Mo17 crosses
Rscript scripts/convert_rqtl2_to_rqtl.R data/Mo17_hybrids.csv analyses/IBM_NLB_BLUPS_Mo17BC.csv data/ibm302map.csv analyses/Mo17_cross.csv
# make cross file for old RIL data
Rscript scripts/convert_rqtl2_to_rqtl.R data/RIL.csv data/old_nlb_data.csv data/ibm302map.csv analyses/old_nlb_data_cross.csv
# make cross file for env BLUPs
Rscript scripts/convert_rqtl2_to_rqtl.R data/RIL.csv analyses/line_blups_envs.csv data/ibm302map.csv analyses/env_blups_cross.csv   