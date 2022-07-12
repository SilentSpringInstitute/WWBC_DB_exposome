from os import path
from toolbox import pathFolder, toolbox
import CompDesc

PR_ROOT = path.abspath("../../") + "/"
PR_DATA = PR_ROOT + "data/"
PR_RESULTS = pathFolder.createFolder(PR_ROOT + "results/")






def fragMatch(p_frag_output, p_db, p_original_db, digit_mass, mode, pr_out, delta_mass = 0):

    l_d_frag_output = toolbox.loadMatrixToList(p_frag_output)
    l_d_db = toolbox.loadMatrixToList(p_db, sep = ",")

    ### Load the original DB for SMILES access
    d_original_db = toolbox.loadMatrix(p_original_db)

    # frag results
    # format MW with digit and split dtx
    for d_frag_output in l_d_frag_output:
        d_frag_output["parent.mass"] = round(abs(float(d_frag_output["parent.mass"])) + delta_mass, digit_mass)
        d_frag_output["RTMean"] = round(float(d_frag_output["RTMean"]), 1)
        d_frag_output["ConsensusID"] = d_frag_output["ConsensusID"].split(",")

        # split SMILES
        d_frag_output["ConsensusSMILES"] = d_frag_output["ConsensusSMILES"].split(",")


    # database
    # define feature id, split DTXSID, format MW
    for d_db in l_d_db:
        d_db["featureid"] = ["%s-%s"%(round(float(d_db["MW"]), 3), round(float(d_db["RT"]), 1))]
        
        # format MW
        d_db["MW"] = round(float(d_db["MW"]), digit_mass)

        # format list dtxsid
        d_db["DB_name"] = d_db["DB_name"].split(";")

        # format list name
        d_db["name"] = d_db["name"].split(";")

        d_db["SMILES"] = []
        # map on DB
        for ID in d_original_db.keys():
            if d_original_db[ID]["DB_name"] in d_db["DB_name"]:
                d_db["SMILES"].append(d_original_db[ID]["SMILES_cleaned"])

        d_db["ESI"] = d_db["ESI"].split(";")

    d_merge_mass = {}
    # merge DB on frag results
    for d_frag_output in l_d_frag_output:
        MW_frag = d_frag_output["parent.mass"]
        if d_frag_output["ConsensusID"] == ['']:
            continue

        for d_db in l_d_db:
            MW_db = d_db["MW"]

            if MW_frag == MW_db:
                #print(d_db["featureid"])
                #print(MW_frag, MW_db)
                #print(d_frag_output["ConsensusID"])
                #print(d_db["DB_name"])

                try: d_merge_mass[MW_db]
                except:
                    d_merge_mass[MW_db] = {}
                    d_merge_mass[MW_db]["DB"] = []
                    d_merge_mass[MW_db]["frag_result"] = []
                
                if not d_db in d_merge_mass[MW_db]["DB"]:
                    d_merge_mass[MW_db]["DB"].append(d_db) 
                
                if not d_frag_output in d_merge_mass[MW_db]["frag_result"]:
                    d_merge_mass[MW_db]["frag_result"].append(d_frag_output) 
                
    
    # confirm match with DTXID
    #DB_name_new==ConsensusID -> confirmed
    #DB_name_new == NA -> unknown

    matched = 0
    for MW_matched in d_merge_mass.keys():
        l_dtx_DB = [d_db["DB_name"] for d_db in d_merge_mass[MW_matched]["DB"]]
        l_dtx_DB = list(set([item for sublist in l_dtx_DB for item in sublist if item != ""]))

        l_dtx_consensus = [d_frag["ConsensusID"] for d_frag in d_merge_mass[MW_matched]["frag_result"]]
        l_dtx_consensus = list(set([item for sublist in l_dtx_consensus for item in sublist if item != ""]))
        
        l_db_name = [d_db["name"] for d_db in d_merge_mass[MW_matched]["DB"]]
        l_db_name = list(set([item.replace("--", ",") for sublist in l_db_name for item in sublist if item != ""]))

        l_feature_id =  [d_db["featureid"] for d_db in d_merge_mass[MW_matched]["DB"]]
        l_feature_id = list(set([item.replace("--", ",") for sublist in l_feature_id for item in sublist if item != ""])) 

        l_smi_DB = [d_db["SMILES"] for d_db in d_merge_mass[MW_matched]["DB"]]
        l_smi_DB = list(set([item for sublist in l_smi_DB for item in sublist if item != ""]))

        l_smi_consensus = [d_db["ConsensusSMILES"] for d_db in d_merge_mass[MW_matched]["frag_result"]]
        l_smi_consensus = list(set([item for sublist in l_smi_consensus for item in sublist if item != ""]))

        l_ESI = [d_db["ESI"] for d_db in d_merge_mass[MW_matched]["DB"]]
        l_ESI = list(set([item for sublist in l_ESI for item in sublist if item != ""]))

        d_merge_mass[MW_matched]["l_dtx_consensus"] = l_dtx_consensus
        d_merge_mass[MW_matched]["l_dtx_DB"] = l_dtx_DB
        d_merge_mass[MW_matched]["l_name_DB"] = l_db_name
        d_merge_mass[MW_matched]["l_feature_id"] = l_feature_id
        d_merge_mass[MW_matched]["l_smi_db"] = l_smi_DB
        d_merge_mass[MW_matched]["l_smi_consensus"] = l_smi_consensus
        d_merge_mass[MW_matched]["ESI"] = l_ESI

        # find intersept
        l_inter = list(set(l_dtx_DB) & set(l_dtx_consensus))
        l_inter = list(set(l_inter))

        if l_dtx_DB == [] and l_dtx_consensus != []:
            d_merge_mass[MW_matched]["match"] = "Unknown - no chemical in DB"
            d_merge_mass[MW_matched]["list match"] = []
        elif  l_dtx_consensus == []:
            d_merge_mass[MW_matched]["match"] = "No match in frag"
            d_merge_mass[MW_matched]["list match"] = []
        elif l_inter == []:
            # need to check 
            # check SMILES matching
            flag_match = 0
            for smi_db in d_merge_mass[MW_matched]["l_smi_db"]:
                
                if flag_match == 1.0:
                    break
                
                c_smidb = CompDesc.CompDesc(smi_db, "")
                c_smidb.prepChem()
                c_smidb.computeFP("MACCS")
                
                for smi_consensus in d_merge_mass[MW_matched]["l_smi_consensus"]:
                    c_smi_consensus = CompDesc.CompDesc(smi_consensus, "")
                    c_smi_consensus.prepChem()
                    c_smi_consensus.computeFP("MACCS")
                    
                    tanimoto_match = c_smidb.computeSimilarityFP(c_smi_consensus, "MACCS", "Tanimoto")
                    print(smi_db, smi_consensus)
                    print(tanimoto_match)
                    if tanimoto_match == 1.0:
                        d_merge_mass[MW_matched]["match"] = "Matched with SMILES"
                        d_merge_mass[MW_matched]["list match"] = [smi_db, smi_consensus]
                        matched = matched + 1
                        flag_match = 1
                        break

            if flag_match == 0.0:
                d_merge_mass[MW_matched]["match"] = "No match intercept with DTXSID or SMILES"
                d_merge_mass[MW_matched]["list match"] = []
        
        else:
            d_merge_mass[MW_matched]["match"] = "Matched with DTXSID"
            d_merge_mass[MW_matched]["list match"] = l_inter
            matched = matched + len(l_inter)

    print("Mode %s\nFor %s mass digit\n%s matches\n"%(mode, digit_mass, matched))



    p_filout = "%s%s_Digit-%s_deltaMass-%s.csv"%(pr_out, mode, digit_mass, delta_mass)

    filout = open(p_filout, "w")
    filout.write("MW_matched\tDTXID in DB\tName in DB\tConsensus DTXID\tfeatureid\tDTXSID matched\tMatch\tESI\n")
    for MW_matched in d_merge_mass.keys():
        filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(MW_matched, ";".join(d_merge_mass[MW_matched]["l_dtx_DB"]), ";".join(d_merge_mass[MW_matched]["l_name_DB"]), ";".join(d_merge_mass[MW_matched]["l_dtx_consensus"]), ";".join(d_merge_mass[MW_matched]["l_feature_id"]),  ";".join(d_merge_mass[MW_matched]["list match"]), d_merge_mass[MW_matched]["match"], ";".join(d_merge_mass[MW_matched]["ESI"])))
    filout.close()


def checkResult(p_jess, p_alex):
    """
    Check on the feature id
    """


    l_d_jess = toolbox.loadMatrixToList(p_jess, sep = ",")
    l_feature_id_jess = []
    for d_jess in l_d_jess:
        if d_jess["confirmed_match"] == "1":
            l_feature_id_jess.append(d_jess["featureid"])
    
    l_d_alex = toolbox.loadMatrixToList(p_alex)
    l_feature_id_alex = []
    for d_alex in l_d_alex:
        if d_alex["Match"] == "Matched with DTXSID":
            l_feature_id_alex = l_feature_id_alex + d_alex["featureid"].split(";")

    l_feature_id_alex = list(set(l_feature_id_alex))
    l_feature_id_jess = list(set(l_feature_id_jess))

    l_inter = list(set(l_feature_id_alex) & set(l_feature_id_jess))

    print("Alex input")
    print(len(l_feature_id_alex))
    print(l_feature_id_alex)

    print("Jess input")
    print(len(l_feature_id_jess))
    print(l_feature_id_jess)

    print("Intercept")
    print(len(l_inter))
    print(l_inter)

    print("Unique in Alex")
    print([x for x in l_feature_id_alex if x not in l_feature_id_jess])

    print("Unique in Jess")
    print([x for x in l_feature_id_jess if x not in l_feature_id_alex])

# frag match positive node
pr_out = pathFolder.createFolder(PR_RESULTS + "fragMatch/")
p_output_frag_pos = PR_DATA + "resultFrag10-2021/node_attributes_table_pos1.tsv" #~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_pos1.tsv
p_db_frag_pos = PR_DATA + "resultFrag10-2021/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv" #"~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation list to DTSC/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"
p_original_DB = PR_RESULTS + "WWBC_MS_database_6.30.21_prepForAnotation.csv"

fragMatch(p_output_frag_pos, p_db_frag_pos, p_original_DB, 1, "Positive", pr_out)
fragMatch(p_output_frag_pos, p_db_frag_pos, p_original_DB, 2, "Positive", pr_out)
fragMatch(p_output_frag_pos, p_db_frag_pos, p_original_DB, 3, "Positive", pr_out)

p_output_frag_neg = PR_DATA + "resultFrag10-2021/node_attributes_table_neg.tsv" #~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_pos1.tsv
p_db_frag_neg = PR_DATA + "resultFrag10-2021/NEG_filter_features_3MW_30deltaRT_NurseFFCombined.csv" #"~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation list to DTSC/POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"

fragMatch(p_output_frag_neg, p_db_frag_neg, p_original_DB, 1, "Negative", pr_out, 2.0)
fragMatch(p_output_frag_neg, p_db_frag_neg, p_original_DB, 2, "Negative", pr_out, 2.0)
fragMatch(p_output_frag_neg, p_db_frag_neg, p_original_DB, 3, "Negative", pr_out, 2.0)


# comparison output
#p_alex = pr_out + "Negative_Digit-3_deltaMass-2.0.csv"
#p_jess = pr_out + "NEG_matched_confirmed_3digits.csv"

#checkResult(p_jess, p_alex)


#p_alex = pr_out + "Negative_Digit-2_deltaMass-2.0.csv"
#p_jess = pr_out + "NEG_matched_confirmed_2digits.csv"

#checkResult(p_jess, p_alex)