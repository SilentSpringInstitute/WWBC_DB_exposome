from toolbox import toolbox
from re import search


class checkTable:
    def __init__(self, p_table1, p_table2, p_log):
        self.p_table1 = p_table1
        self.p_table2 = p_table2
        
        # prepare the report
        self.f_log = open(p_log, "w")
        self.f_log.write("File1 => %s\n"%(self.p_table1))
        self.f_log.write("File2 => %s\n"%(self.p_table2))

    def loadFile(self):

        self.d_table1 = toolbox.loadMatrix(self.p_table1, sep = ",")
        self.d_table2 = toolbox.loadMatrix(self.p_table2, sep = ",")

    def compareCol(self):

        self.f_log.write("\n=====Compare column names=====\n")

        l_col_table1 = list(self.d_table1[list(self.d_table1.keys())[0]].keys())
        l_col_table2 = list(self.d_table2[list(self.d_table2.keys())[0]].keys())

        l_inter_col = list(set(l_col_table1) & set(l_col_table2))
        print("INTERCOL")
        print(l_inter_col)

        self.l_inter_col = l_inter_col

        self.f_log.write("File 1 => %s rows\n"%(len(list(self.d_table1.keys()))))
        self.f_log.write("File 2 => %s rows\n"%(len(list(self.d_table2.keys()))))
        self.f_log.write("File 1 => %s col\n"%(len(l_col_table1)))
        self.f_log.write("File 2 => %s col\n"%(len(l_col_table2)))
        self.f_log.write("Intersection => %s col\n"%(len(l_inter_col)))
        self.f_log.write("Col intersection => %s\n"%(",".join(l_inter_col)))
        self.f_log.write("\n")


    def diffValueInCol(self, join_table, l_skip_col):

        self.f_log.write("\n=====Compare value in intersection col=====\nTable join with col: %s\n"%(join_table))
        
        count_similar = 0
        count_diff = 0
        lw = []
        for k1 in self.d_table1.keys():
            join_k1 = self.d_table1[k1][join_table]
            for k2 in self.d_table2.keys():
                join_k2 = self.d_table2[k2][join_table]

                if join_k1 == join_k2:
                    flag=0
                    l_col_diff = []
                    for intercol in self.l_inter_col:
                        if self.d_table1[k1][intercol] == "-": self.d_table1[k1][intercol] = ""

                        try: self.d_table1[k1][intercol] = str(int(float(self.d_table1[k1][intercol])))
                        except:pass
                        try: self.d_table2[k2][intercol] = str(int(float(self.d_table2[k2][intercol])))
                        except:pass
                        if self.d_table2[k2][intercol] == "-": self.d_table2[k2][intercol] = ""
                        if self.d_table1[k1][intercol] != self.d_table2[k2][intercol] and not intercol in l_skip_col :
                            if not search(self.d_table1[k1][intercol], self.d_table2[k2][intercol]):
                                flag = 1
                                l_col_diff.append(intercol)
                            
                            
                    if flag == 0:
                        count_similar = count_similar + 1
                    else:
                        count_diff = count_diff + 1
                        lw.append("<<Diff %s\n%s\n%s => %s\n%s => %s\nCol diff: %s\n>>\n"%(join_k1, ";".join(l_col_diff), ";".join([self.d_table1[k1][k] for k in l_col_diff]), self.p_table1, ";".join([self.d_table2[k2][k] for k in l_col_diff]), self.p_table2, ";".join(l_col_diff)))


        self.f_log.write("Similar rows: %s\n"%(count_similar))
        self.f_log.write("Different rows: %s\n\n"%(count_diff))

        self.f_log.write("\n".join(lw))

        self.f_log.close()

    def checkDuplicate(self, pfile, col_header, l_merge_col, p_out):
        
        if pfile == 2:
            d_filin = self.d_table2
        else:
            d_filin = self.d_table1
        
        d_out = {}
        for chem in d_filin.keys():
            header = d_filin[chem][col_header]
            if not header in list(d_out.keys()):
                d_out[header] = {}
                for merge_col in l_merge_col:
                    d_out[header][merge_col] = []
                    d_out[header]["count"] = 0

            for merge_col in l_merge_col:
                l_val_col = d_filin[chem][merge_col].split(",")
                for val_col in l_val_col:
                    if not val_col in  d_out[header][merge_col]:
                        d_out[header][merge_col].append(val_col)
            d_out[header]["count"] = d_out[header]["count"] + 1
        
        filout = open(p_out, "w")
        filout.write("%s\t%s\tCount\n"%(col_header, "\t".join(["merge_" + merge_col for merge_col in l_merge_col])))
        for h in d_out.keys():
            filout.write("%s\t%s\t%s\n"%(h,"\t".join([",".join(d_out[h][merge_col]) for merge_col in l_merge_col]), d_out[h]["count"]))
        filout.close()
    

    def checkDupAndTrackFeatureOrigin(self, pfile, col_header, l_merge_col, l_features_col, d_origin, p_out):
        
        
        if pfile == 2:
            d_filin = self.d_table2
        else:
            d_filin = self.d_table1
        
        d_out = {}
        for chem in d_filin.keys():
            header = d_filin[chem][col_header]
            if not header in list(d_out.keys()):
                d_out[header] = {}
                for merge_col in l_merge_col:
                    d_out[header][merge_col] = []
                    d_out[header]["count"] = 0
                    d_out[header]["feature"] = []

            for merge_col in l_merge_col:
                l_val_col = d_filin[chem][merge_col].split(",")
                for val_col in l_val_col:
                    if not val_col in  d_out[header][merge_col]:
                        d_out[header][merge_col].append(val_col)
            d_out[header]["count"] = d_out[header]["count"] + 1

            #check feature 
            for mode in d_origin.keys():
                l_features_mode = d_origin[mode]["featureid"].split(";")
                for feature_col in l_features_col:
                    l_features_tocheck = d_filin[chem][feature_col].split(",")
                    
                    #interception
                    l_feature_inter = list(set(l_features_mode) & set(l_features_tocheck))
                    
                    if len(l_feature_inter) > 0:
                        d_out[header]["feature"].append(mode)
        
        filout = open(p_out, "w")
        filout.write("%s\t%s\tFeature_origin\tCount\n"%(col_header, "\t".join(["merge_" + merge_col for merge_col in l_merge_col])))
        for h in d_out.keys():
            filout.write("%s\t%s\t%s\n"%(h,"\t".join([",".join(d_out[h][merge_col]) for merge_col in l_merge_col]), ",".join(d_out[h]["feature"]), d_out[h]["count"]))
        filout.close()


    def checkIfInfragList(self, pfile1, check_col, flag_col, d_frag, p_out):
        
        d_out = {}
        d_filin = toolbox.loadMatrix(pfile1, sep = ",")
        for chem in d_filin.keys():
            if d_filin[chem][flag_col] != "1":
                continue
            DTXSID = d_filin[chem][check_col]
            d_out[DTXSID] = []
            
            for mode in d_frag.keys():
                for chem_frag in d_frag[mode].keys():
                    l_DTXSID_frag = d_frag[mode][chem_frag]["DB_name"].split(";")
                    print(l_DTXSID_frag)
                    if DTXSID in  l_DTXSID_frag:
                        d_out[DTXSID].append(mode)
        
        filout = open(p_out, "w")
        filout.write("DTXSID\tMode\n")
        for chem in d_out.keys():
            filout.write("%s\t%s\n"%(chem, ",".join(d_out[chem])))
        filout.close()
                    





pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"

# check positive nurse
#pfile2 = "./../../data/result_NTA/20210627FragList_FF_nursesPOS_cleaned.csv"
#p_out = "./../../results/diff_table/log_positive_nurse.txt"

# check parent compounds
#pfile2 = "./../../data/result_NTA/ssi_db_withParentCompunds.csv"
#p_out = "./../../results/diff_table/log_ssi_parentCompound.txt"


# check database from Jessica
#pfile2 = "./../../data/result_NTA/WWBC_MS_database_4.7.21_prepForAnotation_JAT.csv"
#p_out = "./../../results/diff_table/log_DB_jessica.txt"

# check database from Jessica
#pfile2 = "./../../data/result_NTA/2021.06.27ssi_db_withParentCompunds.csv"
#p_out = "./../../results/diff_table/log_2021.06.27ssi.txt"

# vincent
#pfile1 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_06.21_rar.csv"
#pfile2 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_06.21.csv"
#p_out = "./../../results/diff_table/log_rar_vincent.txt"

#vincent rerun 07-21
#pfile2 = "./../../data/result_NTA/List_matched_no_filter_pos_FB_07.21.csv"
#p_out = "./../../results/diff_table/log_FB_pos_vincent_07-21.txt"


# check diff database
#pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
#pfile2 = "./../../results/WWBC_MS_database_4.7.21_prepForAnotation.csv"
#p_out = "./../../results/diff_table/DB_diff.txt"


# chech jessica
#pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
#pfile2 = "./../../data/result_NTA/2021.07.09ssi_db_withParentCompunds.csv"
#p_out = "./../../results/diff_table/log_2021.07.09ssi_db_withParentCompunds.txt"

#vincent rerun 07-21
#pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
#pfile2 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_07.21.csv"
#p_out = "./../../results/diff_table/log_FB_neg_vincent_07-21.txt"


# concatene POS_NEG_FF
#pfile2 = "./../../data/result_NTA/20210713_FullFragList_FF_nursesPOS_sheet2.csv"
#p_out = "./../../results/diff_table/log_20210713_FullFragList_FF_nursesPOS_sheet2.txt"


# check after filtering file
#pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/exposome/data/result_NTA_tofilterforseg/list_filtered_jessica_dup.csv"
#pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/exposome/results/filter_annotation/filtered_annotation.csv"
#p_out = "./../../results/diff_table/log_diff_annotation.txt"


# check after all filtering
#pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/result_NTA/20210723_FF_Nurse_FragList_NEG.csv"
#pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/forFragNeg/NEG_filter_features_3MW_30deltaRT.csv"
#p_out = "./../../results/diff_table/log_diff_annotation.txt"


# check comfirmation list
#pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/all_considered_for_target_N_11.8.21_LH ranking_rar_ETF merge_rar2_AB_ETF_rar_LB.csv"
#pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/N_FFConfirmationAndSemitarget_list.csv"
#p_out = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/merge_check_12-15-21/log_diff.csv"


#cdiff = checkTable(pfile1, pfile2, p_out)
#cdiff.loadFile()
#cdiff.compareCol()
#cdiff.diffValueInCol("mz", ["name_original", "ID", "SMILES", "SMILES_cleaned", "Molweight_cleaned", "ID_raw"])






