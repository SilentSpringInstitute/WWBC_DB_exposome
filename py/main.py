from os import path
from prepForFrag import prepForFrag
from WWBC_database import WWBC_database, matchChemicals
from toolbox import pathFolder
from filterAnnotation import filterAnnotation
from mapAfterFrag import mapAfterFrag
from checkPoint import checkDB

# define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_DATA = PR_ROOT + "data/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")

# SET OF CRITERIA #
###################
# may need to be place as in input file

#minMW = 100
#maxMW = 1000
#lipinski_violation = 3
#list_chemicals_metabolite = ["Drug_UCSF_PXYS", "Drug_most comon and haz", "Disinfectant", "FRs", "PFAS", "MC", "MGDev", "ERactive_bin", "E2Up_bin", "P4Up_bin", "pesticidemammarytumors_bin", "nitroPAH_bin", "pellizzari_bin", "phthalate"]


#  RUN - prepare database  #
############################
#c_db = WWBC_database.WWWBC_database(PR_DATA + "WWBC_MS_database_4.7.21.csv", minMW, maxMW, lipinski_violation, list_chemicals_metabolite, PR_RESULTS, "WWBC_MS_database_6.30.21")
#c_db.main()


## check composition of the DB ###
##################################
## to test
p_DB = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/check/chem_in_DB.txt"

l_check = ["DTXSID1022556", "DTXSID9022601","DTXSID9020116","DTXSID3020910","DTXSID3046742","DTXSID4024983","DTXSID1022845","DTXSID5020364","DTXSID0020365","DTXSID3022877","DTXSID5020364","DTXSID7022883","DTXSID70227388","DTXSID2021735","DTXSID8021480","DTXSID1020562","DTXSID5023035","DTXSID9023049","DTXSID3020625","DTXSID4039657","DTXSID2020634","DTXSID9044299","DTXSID6020648","DTXSID8041032","DTXSID4034150","DTXSID4025371","DTXSID6025438","DTXSID7020760","DTXSID9037664","DTXSID6020804","DTXSID4020822","DTXSID2020892","DTXSID00858942","DTXSID4041070","DTXSID0045703","DTXSID9023413","DTXSID7023431","DTXSID8020541","DTXSID8023135","DTXSID4021185","DTXSID8023557","DTXSID5040910","DTXSID6034186","DTXSID5046354","DTXSID1034187","DTXSID3023631","DTXSID8022371","DTXSID8023688","DTXSID8048288","DTXSID30154863","DTXSID1032278","DTXSID5023742","DTXSID9046023"]
k_search = "DTXSID"

c_check = checkDB.checkDB(p_DB)
c_check.searchInDB(k_search, l_check, p_out)
here


# PREP FOR FRAGMENTATION #
##########################

# Filter anotation #
####################

#from checkPoint import checkTable
#stophere

## nurse ##
# pos
#p_nurse_matched_pos = PR_DATA + "result_NTA_tofilter/List_matched_no_filter_pos_N_FB_07.21_rar.csv"
#p_nurse_matched_pos = PR_DATA + "result_NTA_tofilter/List_matched_no_filter_pos_N_FB_07.21_rar.csv"

# neg
#p_nurse_matched_neg = PR_DATA + "result_NTA_tofilter/List_matched_no_filter_neg_FB_07.21.csv"

#### FF ####
############
# neg
#p_ff_matched_neg = PR_DATA + "result_NTA_tofilter/20210716Neg_features_FF_limited_filtered.csv"

#pos
#p_ff_matched_pos = PR_DATA + "result_NTA_tofilter/20210716POS_features_FF_limited_filtered.csv"

##
# filter setting for analysis
#p_filter = PR_DATA + "result_NTA_tofilter/filter_criteria.txt"

## NEG mode
#pr_out = pathFolder.createFolder(PR_RESULTS + "forFragNeg/")
#l_p_annotation = [p_nurse_matched_neg, p_ff_matched_neg]

#c_prepForFrag = prepForFrag.prepForFrag(l_p_annotation, p_filter, "NEG", pr_out)
#c_prepForFrag.filteredAnnotation()
#c_prepForFrag.mergeAnnotationForFrag()

# POS mode
#pr_out = pathFolder.createFolder(PR_RESULTS + "forFragPos/")
#l_p_annotation = [p_nurse_matched_pos, p_ff_matched_pos]

#c_prepForFrag = prepForFrag.prepForFrag(l_p_annotation, p_filter, "POS", pr_out)
#c_prepForFrag.filteredAnnotation()
#c_prepForFrag.mergeAnnotationForFrag()


## nurse ##
###########

# neg
#p_nurse_matched_neg = PR_DATA + "result_NTA_tofilter/List_matched_no_filter_neg_FB_07.21.csv"
#pr_out = pathFolder.createFolder(PR_RESULTS + "filter_annotation/")

#c_filteria = filterAnnotation.filterAnnotation(p_nurse_matched_neg, p_filter, pr_out)
#c_filteria.loadCriteria()
#c_filteria.filterByCriteriaA()
#c_filteria.removeDuplicateBasedOnCriteriaD()



########### Rematch after frag ################
###################################################
pr_out_frag = pathFolder.createFolder(PR_RESULTS + "map_frag_10-8-21/")

p_pos_frag = PR_DATA + "frag_results_10-8-21/node_attributes_table_neg.tsv"
p_neg_frag = PR_DATA + "frag_results_10-8-21/node_attributes_table_pos1.tsv"

p_combine_for_frag_neg = PR_RESULTS + "forFragNeg/NEG_filter_features_3MW_30deltaRT_NurseFFCombined.csv"
p_combine_for_frag_pos = PR_RESULTS + "forFragPos/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"

c_map = mapAfterFrag.mapAfterFrag(p_pos_frag, p_combine_for_frag_pos, "POS", pr_out_frag)
c_map.map()

c_map = mapAfterFrag.mapAfterFrag(p_neg_frag, p_combine_for_frag_neg, "NEG", pr_out_frag)
c_map.map()


STOPHEREINMAIN
from checkPoint import checkTable

# Extract substructure from matches #
#####################################
pr_out = pathFolder.createFolder(PR_RESULTS + "match_chem/")
d_substructures = {"nitro benzene":["[O-][N+](=O)C1=CC=CC=C1"], "nitro":["[O-][N+](=O)C"], "sulfate":["S(=O)(=O)([O-])C"], "aniline":["Nc1ccccc1"], "guanidium":["NC(N)=N"]}

p_neg_match = PR_DATA + "match_chemicals/20210501list_matched_NEGfeaturesFF_Priority.xlsx"
c_match = matchChemicals.matchChemicals(p_neg_match, "NEG_match",pr_out)
c_match.loadChem("20210501list_matched_NEGfeature")
c_match.searchSubstructure(d_substructures)


p_pos_match = PR_DATA + "match_chemicals/20210501list_matched_POSfeatures_FF_prioritize.xlsx"
c_match = matchChemicals.matchChemicals(p_pos_match, "POS_match",pr_out)
c_match.loadChem("20210501list_matched_POSfeature")
c_match.searchSubstructure(d_substructures)
