from os import name
from re import search
import numpy

from toolbox import toolbox
from copy import deepcopy


class filterAnnotation:
    def __init__(self, p_matched, p_filter_criteria, pr_out):

        self.p_matched = p_matched
        self.p_filter_criteria = p_filter_criteria
        self.pr_out = pr_out


    def loadCriteria(self):

        d_filteria = {}
        f_criteria = open(self.p_filter_criteria, "r")
        l_criteria = f_criteria.readlines()
        f_criteria.close()


        for line_criteria in l_criteria:
            if search("^Criteria", line_criteria):
                l_filters = line_criteria.strip().split(" ")
                d_filteria[l_filters[0]] = " ".join(l_filters[1:])
        
        self.d_criteria = d_filteria


    def filterByCriteria(self):

        
        self.d_prefilter = toolbox.loadMatrix(self.p_matched, sep = ",")
        self.l_colname = list(list(self.d_prefilter.values())[0].keys())


        l_filter_by_criteria = []
        for criteria_filteria in self.d_criteria.keys():

            # define list to work on 
            self.l_work = list(self.d_prefilter.keys())
            l_criteria = self.d_criteria[criteria_filteria].split(" AND ")
            print(l_criteria)
            
            l_l_and = []
            for criterion in l_criteria:
                l_instances = criterion.split(" OR ")
                
                l_l_or = []
                for instance in l_instances:
                    l_elem = instance.split(" ")
                    if len(l_elem) == 3:
                        l_l_or.append(self.filterTable(l_elem[0], l_elem[1], l_elem[2]))
                
                    elif "AVG" in l_elem:
                        value_avg = self.computeAVG(l_elem[l_elem.index("AVG") + 1])
                        print(value_avg)
                        del l_elem[l_elem.index("AVG") + 1]
                        l_elem[l_elem.index("AVG")] = value_avg 
                        if len(l_elem) == 3:
                            l_l_or.append(self.filterTable(l_elem[0], l_elem[1], l_elem[2]))
                              
                # take union of OR and redefine l_work
                if len(l_l_or) == 1:
                    l_l_and.append(l_l_or[0])
                    self.l_work = l_l_or[0]
                else:
                    l_l_and.append(list(set().union(*l_l_or)))
                    self.l_work = list(set().union(*l_l_or))


            # if it is a main criteria apply directly to matrix
            if search("CriteriaM", criteria_filteria):
                if len(l_l_and) == 1:
                    l_chem = l_l_and[0] 
                else:
                    l_chem = list(set(l_l_and[0]).intersection(*l_l_and))
                
                # reduce the matrix
                for chem in list(self.d_prefilter.keys()):
                    if not chem in l_chem:
                        del self.d_prefilter[chem]

            else:
                if len(l_l_and) == 1:
                    l_filter_by_criteria.append(l_l_and[0])
                    
                else:
                    l_filter_by_criteria.append(list(set(l_l_and[0]).intersection(*l_l_and)))
        
            




        # union criteria
        print([len(l) for l in l_filter_by_criteria])
        #print(l_filter_by_criteria)
        
        if len(l_filter_by_criteria) == 1:
            l_filter_all_criteria = l_filter_by_criteria[0]     
        else:
            l_filter_all_criteria = list(set().union(*l_filter_by_criteria))




        p_filout = self.pr_out + "filtered_annotation.csv"
        filout = open(p_filout, "w")

        filout.write("%s\n"%("\t".join(self.l_colname)))
        for ID in l_filter_all_criteria:
            filout.write("%s\t%s\n"%(ID, "\t".join(["%s"%(self.d_prefilter[ID][col]) for col in self.l_colname[1:]])))
        filout.close()



    
    def filterTable(self, name_col, condition, value):

        # test first if name_col is in the table
        if search("/", name_col):
            l_function_col_div = name_col.split("/")
                
        elif not name_col in self.l_colname:
            print("ERROR colname: %s"%(name_col))
            return
        


        l_work = deepcopy(self.l_work)

        i = 0
        imax = len(l_work)
        while i < imax:
            if condition == "==":
                if not name_col in self.l_colname and "l_function_col_div" in locals():
                    if float(self.d_prefilter[l_work[i]][l_function_col_div[0]])/float(self.d_prefilter[l_work[i]][l_function_col_div[1]]) != float(value):
                        del l_work[i]
                        imax = imax - 1
                        continue                    
                
                else:
                    if self.d_prefilter[l_work[i]][name_col] != value:
                        del l_work[i]
                        imax = imax - 1
                        continue
            
            else:
                # need to transform in value
                value_f = float(value)
                
                # need to check if equation
                if "l_function_col_div" in locals():
                    try:value_table = float(self.d_prefilter[l_work[i]][l_function_col_div[0]])/float(self.d_prefilter[l_work[i]][l_function_col_div[1]])
                    except:
                        del l_work[i]
                        imax = imax - 1
                        continue  
                    #print(value_table, i)
                else:
                    try:
                        value_table = float(self.d_prefilter[l_work[i]][name_col])
                    except:
                        del l_work[i]
                        imax = imax - 1
                        continue

                # check value
                if condition == ">=":
                    if value_table < value_f:
                        del l_work[i]
                        imax = imax - 1
                        continue
                
                elif condition == "<=": 
                    if value_table > value_f:
                        del l_work[i]
                        imax = imax - 1
                        continue

                
                elif condition == ">":
                    if value_table <= value_f:
                        del l_work[i]
                        imax = imax - 1
                        continue
                
                elif condition == "<": 
                    if value_table >= value_f:
                        del l_work[i]
                        imax = imax - 1
                        continue
            i = i + 1
        
        return l_work


    def computeAVG(self, col_name):

        l_val = []
        for chem in self.d_prefilter.keys():
            try: l_val.append(float(self.d_prefilter[chem][col_name]))
            except: pass
        
        return numpy.average(l_val)
