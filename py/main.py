from os import path
from prepForFrag import prepForFrag
from WWBC_database import WWBC_database, matchChemicals
from toolbox import pathFolder
from filterAnnotation import filterAnnotation
from mapAfterFrag import mapAfterFrag


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
