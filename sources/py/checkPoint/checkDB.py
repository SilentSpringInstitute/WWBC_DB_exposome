from toolbox import toolbox





class checkDB:
    def __init__(self, p_DB):
        self.p_DB = p_DB
        self.d_DB = toolbox.loadMatrix(self.p_DB)
    
        self.l_col = list(self.d_DB[list(self.d_DB.keys())[0]].keys())


    def searchInDB(self, k_search, l_vals, p_out):

        filout = open(p_out, "w")
        filout.write("Search in DB: %s\n"%(self.p_DB))
        filout.write("Search in col: %s\n\n"%(k_search))

        if not k_search in self.l_col:
            filout.write("ERROR: %s is not a valid col\n"%(k_search))


        for val in l_vals:
            flag_in = 0
            for chem in self.d_DB.keys():
                if self.d_DB[chem][k_search] == val:
                    filout.write("%s in DB\n"%(val))
                    flag_in = 1
                    break
            if flag_in == 0:
                filout.write("%s not found\n"%(val))
        
        filout.close()




