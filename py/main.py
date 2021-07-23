from os import path
from WWBC_database import WWBC_database, matchChemicals
from toolbox import pathFolder
from filterAnnotation import filterAnnotation



# define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
PR_DATA = PR_ROOT + "data/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")

# SET OF CRITERIA #
###################
# may need to be place as in input file

minMW = 100
maxMW = 1000
lipinski_violation = 3
list_chemicals_metabolite = ["Drug_UCSF_PXYS", "Drug_most comon and haz", "Disinfectant", "FRs", "PFAS", "MC", "MGDev", "ERactive_bin", "E2Up_bin", "P4Up_bin", "pesticidemammarytumors_bin", "nitroPAH_bin", "pellizzari_bin", "phthalate"]


#  RUN - prepare database  #
############################
#c_db = WWBC_database.WWWBC_database(PR_DATA + "WWBC_MS_database_4.7.21.csv", minMW, maxMW, lipinski_violation, list_chemicals_metabolite, PR_RESULTS, "WWBC_MS_database_6.30.21")
#c_db.main()


# Filter anotation #
####################

p_filter = PR_DATA + "result_NTA_tofilterforseg/filter_criteria.txt"
p_matched_neg = PR_DATA + "result_NTA_tofilterforseg/List_matched_no_filter_neg_FB_07.21.csv"
pr_out = pathFolder.createFolder(PR_RESULTS + "filter_annotation/")

c_filteria = filterAnnotation.filterAnnotation(p_matched_neg, p_filter, pr_out)
c_filteria.loadCriteria()
c_filteria.filterByCriteria()

from checkPoint import checkTable

STOPHEREINMAIN
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
