from os import path
from prepForFrag import prepForFrag
from WWBC_database import WWBC_database, matchChemicals
from toolbox import pathFolder, toolbox
from filterAnnotation import filterAnnotation
from mapAfterFrag import mapAfterFrag
from checkPoint import checkDB, checkTable
from confirmFrag import confirmFrag

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
## added mapping file for drug in nurse / survey and 
#l_d_fileToMap = [{"pfile":PR_DATA + "buildDB/UCSF_Haz_Drugs Inventory_FINAL 2_2018_wDTSXIDs.csv", "headertoMap": "DTXSID (CompTox)", "addCol": "UCSF_haz_drug"}, {"pfile":PR_DATA + "buildDB/nurses_otherdrugs_rar_10272021.csv", "headertoMap":"DTSXID (from wwbc database plus other sources)", "addCol": "Survey_drug"}]
#c_db = WWBC_database.WWBC_database(PR_DATA + "buildDB/WWBC_MS_database_10.28.21.csv", minMW, maxMW, lipinski_violation, list_chemicals_metabolite, PR_RESULTS, l_d_fileToMap, "WWBC_MS_database_10.28.21")
#c_db.main()


# MAP NURSE RESULTS ON MODE FROM MS
#p_filin = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/nurse_drug_selectionRR_10-29-21/all_considered_for_target_N_10_29_21.csv"
#p_filout = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/nurse_drug_selectionRR_10-29-21/all_considered_for_target_N_10_29_21_mappedOnMode.csv"

#p_pos = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/nurse_drug_selectionRR_10-29-21/List_matched_no_filter_neg_FB_06.21_rar.csv"
#p_neg = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/nurse_drug_selectionRR_10-29-21/List_matched_no_filter_pos_FB_06.21_rar.csv"


#d_tomap = [{"pfile":p_pos, "headertoMap": "DB_name", "addCol": "POS mode"},{"pfile":p_neg, "headertoMap": "DB_name", "addCol": "NEG mode"} ]
#toolbox.manualMapping(p_filin, d_tomap, p_filout)



## check composition of the DB ###
##################################
## to test
#p_DB = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
#p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/check/chem_in_DB.txt"

#l_check = ["DTXSID1022556", "DTXSID9022601","DTXSID9020116","DTXSID3020910","DTXSID3046742","DTXSID4024983","DTXSID1022845","DTXSID5020364","DTXSID0020365","DTXSID3022877","DTXSID5020364","DTXSID7022883","DTXSID70227388","DTXSID2021735","DTXSID8021480","DTXSID1020562","DTXSID5023035","DTXSID9023049","DTXSID3020625","DTXSID4039657","DTXSID2020634","DTXSID9044299","DTXSID6020648","DTXSID8041032","DTXSID4034150","DTXSID4025371","DTXSID6025438","DTXSID7020760","DTXSID9037664","DTXSID6020804","DTXSID4020822","DTXSID2020892","DTXSID00858942","DTXSID4041070","DTXSID0045703","DTXSID9023413","DTXSID7023431","DTXSID8020541","DTXSID8023135","DTXSID4021185","DTXSID8023557","DTXSID5040910","DTXSID6034186","DTXSID5046354","DTXSID1034187","DTXSID3023631","DTXSID8022371","DTXSID8023688","DTXSID8048288","DTXSID30154863","DTXSID1032278","DTXSID5023742","DTXSID9046023"]
#k_search = "DTXSID"

#c_check = checkDB.checkDB(p_DB)
#c_check.searchInDB(k_search, l_check, p_out)


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



####################################
## nurse ##
###########

# neg
#p_nurse_matched_neg = PR_DATA + "result_NTA_tofilter/List_matched_no_filter_neg_FB_07.21.csv"
#pr_out = pathFolder.createFolder(PR_RESULTS + "filter_annotation/")

#c_filteria = filterAnnotation.filterAnnotation(p_nurse_matched_neg, p_filter, pr_out)
#c_filteria.loadCriteria()
#c_filteria.filterByCriteriaA()
#c_filteria.removeDuplicateBasedOnCriteriaD()
###############################################


########### Rematch after frag ################
###################################################
#pr_out_frag = pathFolder.createFolder(PR_RESULTS + "map_frag_10-28-21/")

#p_pos_frag = PR_DATA + "frag_results_10-8-21/node_attributes_table_neg.tsv"
#p_neg_frag = PR_DATA + "frag_results_10-8-21/node_attributes_table_pos1.tsv"

#p_combine_for_frag_neg = PR_RESULTS + "forFragNeg/NEG_filter_features_3MW_30deltaRT_NurseFFCombined.csv"
#p_combine_for_frag_pos = PR_RESULTS + "forFragPos/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"

#c_map = mapAfterFrag.mapAfterFrag(p_pos_frag, p_combine_for_frag_pos, "POS", pr_out_frag)
#p_pos = c_map.map()

#c_map = mapAfterFrag.mapAfterFrag(p_neg_frag, p_combine_for_frag_neg, "NEG", pr_out_frag)
#p_neg = c_map.map()

# check duplicate
#d_pos = toolbox.loadMatrix(p_pos)
#d_neg = toolbox.loadMatrix(p_neg)

#print(list(set(list(d_pos.keys())) & set(list(d_neg.keys()))))

# Extract substructure from matches #
#####################################
#pr_out = pathFolder.createFolder(PR_RESULTS + "match_chem/")
#d_substructures = {"nitro benzene":["[O-][N+](=O)C1=CC=CC=C1"], "nitro":["[O-][N+](=O)C"], "sulfate":["S(=O)(=O)([O-])C"], "aniline":["Nc1ccccc1"], "guanidium":["NC(N)=N"]}

#p_neg_match = PR_DATA + "match_chemicals/20210501list_matched_NEGfeaturesFF_Priority.xlsx"
#c_match = matchChemicals.matchChemicals(p_neg_match, "NEG_match",pr_out)
#c_match.loadChem("20210501list_matched_NEGfeature")
#c_match.searchSubstructure(d_substructures)


#p_pos_match = PR_DATA + "match_chemicals/20210501list_matched_POSfeatures_FF_prioritize.xlsx"
#c_match = matchChemicals.matchChemicals(p_pos_match, "POS_match",pr_out)
#c_match.loadChem("20210501list_matched_POSfeature")
#c_match.searchSubstructure(d_substructures)


#check merge for frag
# check comfirmation list
pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/all_considered_for_target_N_11.8.21_LH ranking_rar_ETF merge_rar2_AB_ETF_rar_LB.csv"
pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/N_FFConfirmationAndSemitarget_list.csv"
p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/log_diff.txt"
p_pos_mode = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"
p_neg_mode = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/NEG_filter_features_3MW_30deltaRT_NurseFFCombined.csv"

#cdiff = checkTable.checkTable(pfile1, pfile2, p_out)
#cdiff.loadFile()
#cdiff.compareCol()
#cdiff.diffValueInCol("DTXSID", ["X"] )#["name_original", "ID", "SMILES", "SMILES_cleaned", "Molweight_cleaned", "ID_raw"])
#cdiff.checkDuplicate(pfile=2, col_header="DTXSID", l_merge_col = ["featureid.x", "featureid.y", "parent.mass", "RTMean", "MW", "RT", "revised.priority.based.on.likelihood.to.work.on.lc.ms...1st.20", "revised.priority.based.on.likelihood.to.work.on.lc.ms...2st.20", "choose.for.targeted.long.and.older.", "POS.mode",
#                                                                                 "NEG.mode", "neg.pos.from.JT", "LH.rank", "LH.notes", "DTSC.and.other.chemists..any.concern.with.compatability.with.your.methods..chemical.stability.in.serum", 
#                                                                                 "DTSC...has.standards.", "Elissia.ID.supplier", "Link.to.order", "QC.EPA", "In.volatile.EPA.list", "In.pubchem", "Polarizability", "LC.MS.compatible..Elissia.",
#                                                                                 "X.1", "X.2", "mode", "name", "group", "mean_visit1", "mean_visit2", "mean_visit3", "mean_N", "mean_OW"], p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/log_dup.csv")

#cdiff.checkDupAndTrackFeatureOrigin(pfile=2, col_header="DTXSID", l_merge_col = ["featureid.x", "featureid.y", "parent.mass", "RTMean", "MW", "RT", "revised.priority.based.on.likelihood.to.work.on.lc.ms...1st.20", "revised.priority.based.on.likelihood.to.work.on.lc.ms...2st.20", "choose.for.targeted.long.and.older.", "POS.mode",
#                                                                                 "NEG.mode", "neg.pos.from.JT", "LH.rank", "LH.notes", "DTSC.and.other.chemists..any.concern.with.compatability.with.your.methods..chemical.stability.in.serum", 
#                                                                                 "DTSC...has.standards.", "Elissia.ID.supplier", "Link.to.order", "QC.EPA", "In.volatile.EPA.list", "In.pubchem", "Polarizability", "LC.MS.compatible..Elissia.",
#                                                                                 "X.1", "X.2", "mode", "name", "group", "mean_visit1", "mean_visit2", "mean_visit3", "mean_N", "mean_OW"], l_features_col=["featureid.x", "featureid.y"], d_origin = {"POS":toolbox.loadMatrix(p_pos_mode), "NEG": toolbox.loadMatrix(p_neg_mode)}, p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/log_dup_mode.csv")


#cdiff.checkIfInfragList(pfile1, check_col = "DTXSID", flag_col= "choose for targeted 1st 20", d_frag = {"POS":toolbox.loadMatrix(p_pos_mode, sep = ","), "NEG": toolbox.loadMatrix(p_neg_mode, sep = ",")}, p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/log_mode.csv")




## Step 5 => confirm from frag
pr_data_confirm = PR_DATA + "toConformFrag_12-21-21/"
p_chem_frag = pr_data_confirm + "all_considered_for_target_N_11.8.21_LH_ranking_rar_ETF_merge_rar2_AB_ETF_rar_LB.csv"#List of target chemicals (only need lines 1-42 - 1 n column H or I): https://berkeley.box.com/s/wsqe11x7de9quea15un2ta5n4bmsyval
p_nurse_pos = pr_data_confirm + "List_matched_no_filter_pos_FB_06.21_rar.csv"#Nurse pos features:https://berkeley.box.com/s/r1b5z7vs9nt8alf6s9yyw2huz1n0n9yv
p_nurse_neg = pr_data_confirm + "List_matched_no_filter_neg_FB_06.21_rar.csv"#nurse negative features: https://berkeley.box.com/s/a0h4qxl8arbn7wlel6jc92d26v3pd2hr
p_FF_neg = pr_data_confirm +  "20210716Neg_features_FF_limited.csv"#FF pos features https://berkeley.box.com/s/ot52cqpdc1di4nc01oxg9coqvwftr078
p_FF_pos = pr_data_confirm +  "20210716POS_features_FF_limited.csv"#FF neg features  https://berkeley.box.com/s/gbd7lvpl5m7mp7i3r5e8tia4y01by8pc
p_node_neg = pr_data_confirm +  "node_attributes_table_neg.tsv"
p_node_pos = pr_data_confirm +  "node_attributes_table_pos1.tsv"


#revised priority based on likelihood to work on lc/ms - 1st 20	


l_select_target_chem = ["revised priority based on likelihood to work on lc/ms - 1st 20", "revised priority based on likelihood to work on lc/ms - 2st 20"]
l_col_target_chem = ["revised priority based on likelihood to work on lc/ms - 1st 20", "revised priority based on likelihood to work on lc/ms - 2st 20", "CASRN","name_original","formula_original","Elissia.ID.supplier", "Link.to.order", "QC EPA", "In volatile EPA list", "In pubchem", "Polarizability", "LC-MS compatible (Elissia)", "ionization", "comments"]
d_col_study = {"nurse":["mz","time", "mean_N", "mean_OW", "mean_FB"], "FF":["mz.y","rt.y", "mean_allvisit", "mean_FB"]}

#l_select_target_chem, l_col_target_chem, l_col_analysis, p_out
p_data_confirm = pr_data_confirm + "mapping_out_1-5-22.csv"

c_confirm_frag = confirmFrag.confirmFrag(p_chem_frag, p_nurse_pos, p_nurse_neg, p_FF_pos, p_FF_neg, p_node_neg, p_node_pos)
c_confirm_frag.loadDataSet()
c_confirm_frag.mapFeatureToChem(l_select_target_chem, l_col_target_chem, p_data_confirm)

