from re import search
from filterAnnotation import filterAnnotation

class prepForFrag:
    def __init__(self, l_p_annotated_files, p_filter_criteria, name_file, pr_out):
        self.l_p_annotated_filtered_files = l_p_annotated_files
        self.pr_out = pr_out
        self.name_file = name_file
        self.p_filter_criteria = p_filter_criteria
        self.verbose = 0

    def filteredAnnotation(self):
        """
        Apply anotation filters on several annotation files organized in a list
        """


        d_d_filtered_annotation = {}
        # load class criteria before the loop to avoid memory duplicate
        self.c_filteria = filterAnnotation.filterAnnotation("", self.p_filter_criteria, self.pr_out)
        self.c_filteria.loadCriteria()
        for p_annotated_file in self.l_p_annotated_filtered_files:
            self.c_filteria.p_matched = p_annotated_file
            self.c_filteria.filterByCriteriaA() # apply criteria on the mapping
            d_d_filtered_annotation[p_annotated_file.split("/")[-1][0:-4]] = self.c_filteria.d_filtered
        self.d_d_filtered_annotation = d_d_filtered_annotation   
        
        if self.verbose == 1:
            print("Filtered annotation")
            for k_filtered in self.d_d_filtered_annotation.keys():
                print( k_filtered, " filtered:", len(list(d_d_filtered_annotation[k_filtered].keys())))

    def mergeAnnotationForFrag(self):
        """
        Merge annotation for fragmentation based on features and lost of the chemicals mapping
        arg: none
            work on the self.d_d_filtered_annotation => annotation need to be computed before
        return:
            self.d_annotationMerged and write in a file for exportation in the self.pr_out folder
        """

        # step 1 => extract criteria round digit for MW and rt delta
        for criteria in self.c_filteria.d_criteria.keys():
            if search("^CriteriaD", criteria):
                if search ("mz", self.c_filteria.d_criteria[criteria]):
                    mw_digits = int(self.c_filteria.d_criteria[criteria].split("mz => ")[-1].split(" ")[0])
                elif search("time", self.c_filteria.d_criteria[criteria]):
                    rt_delta = int(self.c_filteria.d_criteria[criteria].split("time => ")[-1].split(" ")[0])
                elif search("intensity", self.c_filteria.d_criteria[criteria]):
                    intensity_min = float(self.c_filteria.d_criteria[criteria].strip().split("intensity >= ")[-1].split(" ")[0])


        # step2 group by MW
        d_MW = {}
        for k_filtered in self.d_d_filtered_annotation.keys():
            print(k_filtered)
            for ID in self.d_d_filtered_annotation[k_filtered].keys():
                ##### HERE format table for FF to have the same structure
                # => need to change in early pipeline step ==> colname are not consistant 
                if not "mz" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    self.d_d_filtered_annotation[k_filtered][ID]["mz"] = self.d_d_filtered_annotation[k_filtered][ID]["mz.x"]
                if not "time" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    self.d_d_filtered_annotation[k_filtered][ID]["time"] = self.d_d_filtered_annotation[k_filtered][ID]["rt.x"]
                if not "mean_N_OW" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    self.d_d_filtered_annotation[k_filtered][ID]["mean_N_OW"] = self.d_d_filtered_annotation[k_filtered][ID]["mean_visit1"]
                if not "featureid" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    self.d_d_filtered_annotation[k_filtered][ID]["featureid"] = str(self.d_d_filtered_annotation[k_filtered][ID]["mz"]) + "-" + str(self.d_d_filtered_annotation[k_filtered][ID]["time"])
                
                if not "ESI" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    self.d_d_filtered_annotation[k_filtered][ID]["ESI"] = "FF"
                else:
                    self.d_d_filtered_annotation[k_filtered][ID]["ESI"] = "N"


                # check intensity
                # for nurse
                if "mean_N_OW" in list(self.d_d_filtered_annotation[k_filtered][ID].keys()):
                    try:intensity_val = float(self.d_d_filtered_annotation[k_filtered][ID]["mean_N_OW"])
                    except: continue # case where NA and not able to convert in float
                    if intensity_val <= intensity_min:
                        continue

                MW = str(round(float(self.d_d_filtered_annotation[k_filtered][ID]["mz"]), mw_digits))
                if not MW in list(d_MW.keys()):
                    d_MW[MW] = []
                d_add = {"mz":[self.d_d_filtered_annotation[k_filtered][ID]["mz"]], "featureid":[self.d_d_filtered_annotation[k_filtered][ID]["featureid"]], "time":[self.d_d_filtered_annotation[k_filtered][ID]["time"]], "formula_cleaned": [self.d_d_filtered_annotation[k_filtered][ID]["formula_cleaned"]], "DB_name": [self.d_d_filtered_annotation[k_filtered][ID]["DB_name"]], "ESI":  [self.d_d_filtered_annotation[k_filtered][ID]["ESI"]]}
                d_MW[MW].append(d_add)


        # step 3 with the RT and the dofference 30s 
        ## only one turn of RT
        for MW in d_MW.keys():
            i = 0
            imax = len(d_MW[MW])
            while i < imax:
                RT_i = float(d_MW[MW][i]["time"][0])
                j = i + 1
                while j < imax:
                    RT_j = float(d_MW[MW][j]["time"][0])

                    if abs(RT_i - RT_j) <= rt_delta:
                        for k in d_MW[MW][i].keys():
                            d_MW[MW][i][k] = d_MW[MW][i][k] + d_MW[MW][j][k]
                        imax = imax - 1
                        del d_MW[MW][j]
                        continue
                    else:
                        j = j + 1
                i = i + 1
        
        if self.verbose == 1:
            print("Check RT and MW filter:")
            for MW in d_MW.keys():
                print(MW, " time: ", ",".join([str(d_MW_t["time"]) for d_MW_t in d_MW[MW]]))
                print(MW, " DBname: ", ",".join([str(d_MW_t["DB_name"]) for d_MW_t in d_MW[MW]]))

        # write results
        p_filout = self.pr_out + "%s_filter_features_%sMW_%sdeltaRT.csv"%(self.name_file, mw_digits, rt_delta)
        filout = open(p_filout, "w")
        l_col = list(d_MW[list(d_MW.keys())[0]][0].keys())
        filout.write("MW\tRT\t" + "\t".join(l_col) + "\n")
        l_MW = list(d_MW.keys())
        l_MW.sort()
        for MW in l_MW:
            for d_RT in d_MW[MW]:
                filout.write("%s\t%s\t%s\n"%(MW, d_RT["time"][0], "\t".join([';'.join(d_RT[k]) for k in l_col])))
        filout.close()



        return 