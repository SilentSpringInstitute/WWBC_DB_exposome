from os import path
from toolbox import toolbox




PR_ROOT = path.abspath("../../") + "/"
PR_DATA = PR_ROOT + "data/"
pr_in = PR_DATA + "comparison_2rounds/"


def comparisonRound():

    p_neg_2020 = pr_in + "list_DTXSID_fragmentation_nurses_FF_neg.csv"
    p_neg_2021 = pr_in + "NEG_filter_features_3MW_30deltaRT_NurseFFCombined.csv"

    p_pos_2020 = pr_in + "list_DTXSID_fragmentation_nurses_FF_pos.csv"
    p_pos_2021 = pr_in + "POS_filter_features_3MW_30deltaRT_NurseFFCombined.csv"


    l_d_neg_2020 = toolbox.loadMatrixToList(p_neg_2020, sep = ",")
    l_d_neg_2021 = toolbox.loadMatrixToList(p_neg_2021, sep = ",")


    l_d_pos_2020 = toolbox.loadMatrixToList(p_pos_2020, sep = ",")
    l_d_pos_2021 = toolbox.loadMatrixToList(p_pos_2021, sep = ",")


    l_mass_neg_2020 = []
    for d_neg_2020 in l_d_neg_2020:
        mass = float(d_neg_2020["mz"])
        mass = round(mass, 2)
        l_mass_neg_2020.append(mass)


    l_mass_neg_2021 = []
    for d_neg_2021 in l_d_neg_2021:
        mass = float(d_neg_2021["MW"])
        mass = round(mass, 2)
        l_mass_neg_2021.append(mass)

    # remove duplicate
    l_mass_neg_2021 = list(set(l_mass_neg_2021))
    l_mass_neg_2020 = list(set(l_mass_neg_2020))

    l_inter_neg = list(set(l_mass_neg_2021) & set(l_mass_neg_2020))

    print("===negative mode===")
    print("Nb mass in 2020 round: ", len(l_mass_neg_2020))
    print("Nb mass in 2021 round: ", len(l_mass_neg_2021))

    print("Intersection: ", len(l_inter_neg))
    print(l_inter_neg)


    l_mass_pos_2020 = []
    for d_pos_2020 in l_d_pos_2020:
        mass = float(d_pos_2020["mz"])
        mass = round(mass, 2)
        l_mass_pos_2020.append(mass)


    l_mass_pos_2021 = []
    for d_pos_2021 in l_d_pos_2021:
        mass = float(d_pos_2021["MW"])
        mass = round(mass, 2)
        l_mass_pos_2021.append(mass)

    # remove duplicate
    l_mass_pos_2021 = list(set(l_mass_pos_2021))
    l_mass_pos_2020 = list(set(l_mass_pos_2020))

    l_inter_pos = list(set(l_mass_pos_2021) & set(l_mass_pos_2020))

    print("===positive mode===")
    print("Nb mass in 2020 round: ", len(l_mass_pos_2020))
    print("Nb mass in 2021 round: ", len(l_mass_pos_2021))

    print("Intersection: ", len(l_inter_pos))
    print(l_inter_pos)



    # check inter

    p_tsv_neg =  pr_in + "node_attributes_table_neg.tsv"
    p_tsv_pos =  pr_in + "node_attributes_table_pos1.tsv"

    l_d_tsv_neg = toolbox.loadMatrixToList(p_tsv_neg, sep = "\t")
    l_d_tsv_pos = toolbox.loadMatrixToList(p_tsv_pos, sep = "\t")

    l_mass_tsv_neg = []
    for d_tsv_neg in l_d_tsv_neg:
        mass = float(d_tsv_neg["parent.mass"])
        l_mass_tsv_neg.append(round(mass, 2))

    l_mass_tsv_neg = list(set(l_mass_tsv_neg))
    l_inter_tsv_neg = list(set(l_inter_neg) & set(l_mass_tsv_neg))

    print("=====TSV file neg=====")
    print("Nb unique parents mass:", len(l_mass_tsv_neg))
    print("Nb intercept: ", len(l_inter_tsv_neg))
    print(l_inter_tsv_neg)



    l_mass_tsv_pos = []
    for d_tsv_pos in l_d_tsv_pos:
        mass = float(d_tsv_pos["parent.mass"])
        l_mass_tsv_pos.append(round(mass, 2))

    l_mass_tsv_pos = list(set(l_mass_tsv_pos))
    l_inter_tsv_pos = list(set(l_inter_pos) & set(l_mass_tsv_pos))

    print("=====TSV file pos=====")
    print("Nb unique parents mass:", len(l_mass_tsv_pos))
    print("Nb intercept: ", len(l_inter_tsv_pos))
    print(l_inter_tsv_pos)



def merge20and21(p_tsv_20, p_tsv_21, n_dgit, p_filout):

    l_d_tsv_21 =  toolbox.loadMatrixToList(p_tsv_21, sep = "\t")
    l_d_tsv_20 = toolbox.loadMatrixToList(p_tsv_20, sep = "\t")

    #merge on the parent.mass with n digit   

    for d_tsv_20 in l_d_tsv_20:
        d_tsv_20["parent.mass"] = abs(round(float(d_tsv_20["parent.mass"]), n_dgit))
    
    l_k20 = list(d_tsv_20.keys())

    for d_tsv_21 in l_d_tsv_21:
        d_tsv_21["parent.mass"] = abs(round(float(d_tsv_21["parent.mass"]), n_dgit))

    l_k21 = list(d_tsv_21.keys())

    l_k = list(set(l_k20) & set(l_k21))

    filout = open(p_filout, "w")
    filout.write("%s\n"%("\t".join([str(k) + "20\t" + str(k) + "21" for k in l_k])))
    for d_tsv_20 in l_d_tsv_20:
        mass_parents = d_tsv_20["parent.mass"]
        print("20", mass_parents)

        for d_tsv_21 in l_d_tsv_21:
            print("21", d_tsv_21["parent.mass"])
            if mass_parents == d_tsv_21["parent.mass"]:
                filout.write("%s\n"%("\t".join([str(d_tsv_20[k]) + "\t" + str(d_tsv_21[k]) for k in l_k])))
    
    filout.close()



    return 


#### run main ####
##################


#comparisonRound()


## POSITIVE ##
p_tsv_21_pos = pr_in + "node_attributes_table_pos1.tsv" #"~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_pos1.tsv
p_tsv_20_pos = pr_in + "node_attributes_table_positive.tsv" #"~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/DoNotUse_2020Fragmentation/node_attributes_table_positive.tsv"

#merge20and21(p_tsv_20_pos, p_tsv_21_pos, 1, pr_in + "mergePos20-21_1digit.csv")
#merge20and21(p_tsv_20_pos, p_tsv_21_pos, 2, pr_in + "mergePos20-21_2digit.csv")


## NEGATIVE ##
p_tsv_21_pos = pr_in + "node_attributes_table_neg.tsv" #~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_neg.tsv
p_tsv_20_pos = pr_in + "node_attributes_table_negative.tsv" #"~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/DoNotUse_2020Fragmentation/node_attributes_table_negative.tsv"

merge20and21(p_tsv_20_pos, p_tsv_21_pos, 1, pr_in + "mergeNEG20-21_1digit.csv")
merge20and21(p_tsv_20_pos, p_tsv_21_pos, 2, pr_in + "mergeNEG20-21_2digit.csv")


