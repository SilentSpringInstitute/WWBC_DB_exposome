import toolbox



class checkTable:
    def __init__(self, p_table1, p_table2, pr_out):
        self.p_table1 = p_table1
        self.p_table2 = p_table2
        self.pr_out = pr_out
        
        # prepare the report
        self.f_log = open(pr_out + "log.txt", "w")
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
                    for intercol in self.l_inter_col:
                        if self.d_table1[k1][intercol] != self.d_table2[k2][intercol] and not intercol in l_skip_col:
                            flag = 1
                            break
                    
                    if flag == 0:
                        count_similar = count_similar + 1
                    else:
                        count_diff = count_diff + 1
                        lw.append("<<Diff %s\n%s\n%s => %s\n%s => %s\n>>\n"%(join_k1, ",".join(self.l_inter_col), ",".join([self.d_table1[k1][k] for k in self.l_inter_col]), self.p_table1, ",".join([self.d_table2[k2][k] for k in self.l_inter_col]), self.p_table2))


        self.f_log.write("Similar rows: %s\n"%(count_similar))
        self.f_log.write("Different rows: %s\n\n"%(count_diff))

        self.f_log.write("\n".join(lw))

        self.f_log.close()


    
pfile1 = "./../../results/WWBC_MS_database_4.7.21_prepForAnotation.csv"
pfile2 = "./../../data/result_NTA/20210627FragList_FF_nursesPOS_cleaned.csv"
prout = "./../../results/diff_table/"

cdiff = checkTable(pfile1, pfile2, prout)
cdiff.loadFile()
cdiff.compareCol()
cdiff.diffValueInCol("DB_name", ["name_original", "ID"])

