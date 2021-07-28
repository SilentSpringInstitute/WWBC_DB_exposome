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

        self.d_table1 = toolbox.loadMatrix(self.p_table1, sep = "\t")
        self.d_table2 = toolbox.loadMatrix(self.p_table2, sep = ",")

    def compareCol(self):

        self.f_log.write("\n=====Compare column names=====\n")

        l_col_table1 = list(self.d_table1[list(self.d_table1.keys())[0]].keys())
        l_col_table2 = list(self.d_table2[list(self.d_table2.keys())[0]].keys())

        l_inter_col = list(set(l_col_table1) & set(l_col_table2))

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


    

pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"

# check positive nurse
#pfile2 = "./../../data/result_NTA/20210627FragList_FF_nursesPOS_cleaned.csv"
#p_out = "./../../results/diff_table/log_positive_nurse.txt"

# check parent compounds
pfile2 = "./../../data/result_NTA/ssi_db_withParentCompunds.csv"
p_out = "./../../results/diff_table/log_ssi_parentCompound.txt"


# check database from Jessica
pfile2 = "./../../data/result_NTA/WWBC_MS_database_4.7.21_prepForAnotation_JAT.csv"
p_out = "./../../results/diff_table/log_DB_jessica.txt"

# check database from Jessica
pfile2 = "./../../data/result_NTA/2021.06.27ssi_db_withParentCompunds.csv"
p_out = "./../../results/diff_table/log_2021.06.27ssi.txt"

# vincent
#pfile1 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_06.21_rar.csv"
#pfile2 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_06.21.csv"
#p_out = "./../../results/diff_table/log_rar_vincent.txt"

#vincent rerun 07-21
pfile2 = "./../../data/result_NTA/List_matched_no_filter_pos_FB_07.21.csv"
p_out = "./../../results/diff_table/log_FB_pos_vincent_07-21.txt"


# check diff database
pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
pfile2 = "./../../results/WWBC_MS_database_4.7.21_prepForAnotation.csv"
p_out = "./../../results/diff_table/DB_diff.txt"


# chech jessica
pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
pfile2 = "./../../data/result_NTA/2021.07.09ssi_db_withParentCompunds.csv"
p_out = "./../../results/diff_table/log_2021.07.09ssi_db_withParentCompunds.txt"

#vincent rerun 07-21
pfile1 = "./../../results/WWBC_MS_database_6.30.21_prepForAnotation.csv"
pfile2 = "./../../data/result_NTA/List_matched_no_filter_neg_FB_07.21.csv"
p_out = "./../../results/diff_table/log_FB_neg_vincent_07-21.txt"


# concatene POS_NEG_FF
pfile2 = "./../../data/result_NTA/20210713_FullFragList_FF_nursesPOS_sheet2.csv"
p_out = "./../../results/diff_table/log_20210713_FullFragList_FF_nursesPOS_sheet2.txt"


# check after filtering file
pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/exposome/data/result_NTA_tofilterforseg/list_filtered_jessica_dup.csv"
pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/exposome/results/filter_annotation/filtered_annotation.csv"
p_out = "./../../results/diff_table/log_diff_annotation.txt"


# check after all filtering
pfile2 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/data/result_NTA/20210723_FF_Nurse_FragList_NEG.csv"
pfile1 = "/mnt/c/Users/AlexandreBorrel/research/SSI/NTA/results/forFragNeg/NEG_filter_features_3MW_30deltaRT.csv"
p_out = "./../../results/diff_table/log_diff_annotation.txt"


cdiff = checkTable(pfile1, pfile2, p_out)
cdiff.loadFile()
cdiff.compareCol()
cdiff.diffValueInCol("mz", ["name_original", "ID", "SMILES", "SMILES_cleaned", "Molweight_cleaned", "ID_raw"])






