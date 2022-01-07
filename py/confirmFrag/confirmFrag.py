from toolbox import toolbox






class confirmFrag:
    def __init__(self, p_chem_target, p_pos_nurse, p_neg_nurse, p_pos_FF, p_neg_FF, p_neg_node, p_pos_node):
        
        self.p_chem_target = p_chem_target
        self.p_pos_nurse = p_pos_nurse
        self.p_neg_nurse = p_neg_nurse
        self.p_pos_FF = p_pos_FF
        self.p_neg_FF = p_neg_FF
        self.p_neg_node = p_neg_node
        self.p_pos_node = p_pos_node
        
        pass
    
    def loadDataSet(self):
        
        self.d_chem_target = toolbox.loadMatrix(self.p_chem_target, sep = ",")
        self.d_pos_nurse = toolbox.loadMatrix(self.p_pos_nurse, sep = ",")
        self.d_neg_nurse = toolbox.loadMatrix(self.p_neg_nurse, sep = ",")
        self.d_pos_FF = toolbox.loadMatrix(self.p_pos_FF, sep = ",")
        self.d_neg_FF = toolbox.loadMatrix(self.p_neg_FF, sep = ",")
        self.d_node_pos = toolbox.loadMatrix(self.p_pos_node, sep = "\t")
        self.d_node_neg = toolbox.loadMatrix(self.p_neg_node, sep = "\t")
        
    def mapFeatureToChem(self, l_select_target_chem, l_col_target_chem, p_out):
        
        d_out = {}
        for chem_target in self.d_chem_target.keys():
            
            
            flag = 0
            for select_target_chem in l_select_target_chem:
                if self.d_chem_target[chem_target][select_target_chem] != "":
                    flag=1
            if flag == 0:
                continue
            
            
            
            DTXSID = self.d_chem_target[chem_target]["DTXSID"]
            d_out[DTXSID] = {}
            
            d_out[DTXSID]["target_chem"] = {}
            for col_target_chem in l_col_target_chem:
                d_out[DTXSID]["target_chem"][col_target_chem] = self.d_chem_target[chem_target][col_target_chem]

            #pos nurse
            d_out[DTXSID]["NURSE_POS"] = {}
            for id_pos_nurse in self.d_pos_nurse.keys():
                if self.d_pos_nurse[id_pos_nurse]["DB_name"] == DTXSID:
                    feature = self.d_pos_nurse[id_pos_nurse]["featureid"]
                    if not feature in list(d_out[DTXSID]["NURSE_POS"].keys()):
                        d_out[DTXSID]["NURSE_POS"][feature] = {}
                        d_out[DTXSID]["NURSE_POS"][feature]["RT"] = self.d_pos_nurse[id_pos_nurse]["time"]
                        d_out[DTXSID]["NURSE_POS"][feature]["MZ"] = self.d_pos_nurse[id_pos_nurse]["mz"]
                        d_out[DTXSID]["NURSE_POS"][feature]["NURSE_mean_FB"] = self.d_pos_nurse[id_pos_nurse]["mean_FB"]
                        d_out[DTXSID]["NURSE_POS"][feature]["NURSE_mean_OW"] = self.d_pos_nurse[id_pos_nurse]["mean_OW"]
                        d_out[DTXSID]["NURSE_POS"][feature]["NURSE_mean_N"] = self.d_pos_nurse[id_pos_nurse]["mean_N"]
                        d_out[DTXSID]["NURSE_POS"][feature]["FF_mean_visit1"] = ""
                        d_out[DTXSID]["NURSE_POS"][feature]["FF_mean_visit2"] = ""
                        d_out[DTXSID]["NURSE_POS"][feature]["FF_mean_visit3"] = ""
                        d_out[DTXSID]["NURSE_POS"][feature]["FF_mean_FB"] = ""
                        d_out[DTXSID]["NURSE_POS"][feature]["RT_mean"] = []
                        d_out[DTXSID]["NURSE_POS"][feature]["parent_mass"] = []
                        
            #neg nurse
            d_out[DTXSID]["NURSE_NEG"] = {}
            for id_neg_nurse in self.d_neg_nurse.keys():
                if self.d_neg_nurse[id_neg_nurse]["DB_name"] == DTXSID:
                    feature = self.d_neg_nurse[id_neg_nurse]["featureid"]
                    if not feature in list(d_out[DTXSID]["NURSE_NEG"].keys()):
                        d_out[DTXSID]["NURSE_NEG"][feature] = {}
                        d_out[DTXSID]["NURSE_NEG"][feature]["RT"] = self.d_neg_nurse[id_neg_nurse]["time"]
                        d_out[DTXSID]["NURSE_NEG"][feature]["MZ"] = self.d_neg_nurse[id_neg_nurse]["mz"]
                        d_out[DTXSID]["NURSE_NEG"][feature]["NURSE_mean_FB"] = self.d_neg_nurse[id_neg_nurse]["mean_FB"]
                        d_out[DTXSID]["NURSE_NEG"][feature]["NURSE_mean_OW"] = self.d_neg_nurse[id_neg_nurse]["mean_OW"]
                        d_out[DTXSID]["NURSE_NEG"][feature]["NURSE_mean_N"] = self.d_neg_nurse[id_neg_nurse]["mean_N"]
                        d_out[DTXSID]["NURSE_NEG"][feature]["FF_mean_visit1"] = ""
                        d_out[DTXSID]["NURSE_NEG"][feature]["FF_mean_visit2"] = ""
                        d_out[DTXSID]["NURSE_NEG"][feature]["FF_mean_visit3"] = ""
                        d_out[DTXSID]["NURSE_NEG"][feature]["FF_mean_FB"] = ""
                        d_out[DTXSID]["NURSE_NEG"][feature]["RT_mean"] = []
                        d_out[DTXSID]["NURSE_NEG"][feature]["parent_mass"] = []
            
            #pos FF
            d_out[DTXSID]["FF_POS"] = {}
            for id_pos_FF in self.d_pos_FF.keys():
                if self.d_pos_FF[id_pos_FF]["DB_name"] == DTXSID:
                    feature = self.d_pos_FF[id_pos_FF]["featureid"]
                    if not feature in list(d_out[DTXSID]["FF_POS"].keys()):
                        d_out[DTXSID]["FF_POS"][feature] = {}
                        d_out[DTXSID]["FF_POS"][feature]["RT"] = self.d_pos_FF[id_pos_FF]["rt"]
                        d_out[DTXSID]["FF_POS"][feature]["MZ"] = self.d_pos_FF[id_pos_FF]["mz"]
                        d_out[DTXSID]["FF_POS"][feature]["FF_mean_FB"] = self.d_pos_FF[id_pos_FF]["mean_FB"]
                        d_out[DTXSID]["FF_POS"][feature]["NURSE_mean_OW"] = ""
                        d_out[DTXSID]["FF_POS"][feature]["NURSE_mean_FB"] = ""
                        d_out[DTXSID]["FF_POS"][feature]["NURSE_mean_N"] = ""
                        d_out[DTXSID]["FF_POS"][feature]["FF_mean_visit1"] = self.d_pos_FF[id_pos_FF]["mean_visit1"]
                        d_out[DTXSID]["FF_POS"][feature]["FF_mean_visit2"] = self.d_pos_FF[id_pos_FF]["mean_visit2"]
                        d_out[DTXSID]["FF_POS"][feature]["FF_mean_visit3"] = self.d_pos_FF[id_pos_FF]["mean_visit3"]
                        d_out[DTXSID]["FF_POS"][feature]["RT_mean"] = []
                        d_out[DTXSID]["FF_POS"][feature]["parent_mass"] = []
            #neg FF
            d_out[DTXSID]["FF_NEG"] = {}
            for id_neg_FF in self.d_neg_FF.keys():
                if self.d_neg_FF[id_neg_FF]["DB_name"] == DTXSID:
                    feature = self.d_neg_FF[id_neg_FF]["featureid"]
                    if not feature in list(d_out[DTXSID]["FF_NEG"].keys()):
                        d_out[DTXSID]["FF_NEG"][feature] = {}
                        d_out[DTXSID]["FF_NEG"][feature]["RT"] = self.d_neg_FF[id_neg_FF]["rt.x"]
                        d_out[DTXSID]["FF_NEG"][feature]["MZ"] = self.d_neg_FF[id_neg_FF]["mz.x"]
                        d_out[DTXSID]["FF_NEG"][feature]["FF_mean_FB"] = self.d_neg_FF[id_neg_FF]["mean_FB"]
                        d_out[DTXSID]["FF_NEG"][feature]["NURSE_mean_OW"] = ""
                        d_out[DTXSID]["FF_NEG"][feature]["NURSE_mean_FB"] = ""
                        d_out[DTXSID]["FF_NEG"][feature]["NURSE_mean_N"] = ""
                        d_out[DTXSID]["FF_NEG"][feature]["FF_mean_visit1"] = self.d_neg_FF[id_neg_FF]["mean_visit1"]
                        d_out[DTXSID]["FF_NEG"][feature]["FF_mean_visit2"] = self.d_neg_FF[id_neg_FF]["mean_visit2"]
                        d_out[DTXSID]["FF_NEG"][feature]["FF_mean_visit3"] = self.d_neg_FF[id_neg_FF]["mean_visit3"]
                        d_out[DTXSID]["FF_NEG"][feature]["RT_mean"] = []
                        d_out[DTXSID]["FF_NEG"][feature]["parent_mass"] = [] 
        
        ### add RT_mean and parent mass from node file
        l_neg = ["FF_NEG", "NURSE_NEG"]
        l_pos = ["FF_POS", "NURSE_POS"]
        
        for DTXSID in d_out.keys():
            # neg mode
            for neg_node in l_neg:
                for k in self.d_node_neg.keys():
                    l_consensus_id = self.d_node_neg[k]["ConsensusID"].split(",")
                    if DTXSID in l_consensus_id:
                        RT_mean = self.d_node_neg[k]["RTMean"]
                        mass_parent = self.d_node_neg[k]["parent.mass"]
                        for feature in d_out[DTXSID][neg_node].keys():
                            if not RT_mean in d_out[DTXSID][neg_node][feature]["RT_mean"]:
                                d_out[DTXSID][neg_node][feature]["RT_mean"].append(RT_mean)
                            if not mass_parent in d_out[DTXSID][neg_node][feature]["parent_mass"]:
                                d_out[DTXSID][neg_node][feature]["parent_mass"].append(mass_parent)
                            
            for pos_node in l_pos:
                for k in self.d_node_pos.keys():
                    l_consensus_id = self.d_node_pos[k]["ConsensusID"].split(",")
                    if DTXSID in l_consensus_id:
                        RT_mean = self.d_node_pos[k]["RTMean"]
                        mass_parent = self.d_node_pos[k]["parent.mass"]
                        for feature in d_out[DTXSID][pos_node].keys():
                            if not RT_mean in d_out[DTXSID][pos_node][feature]["RT_mean"]:
                                d_out[DTXSID][pos_node][feature]["RT_mean"].append(RT_mean)
                            if not mass_parent in d_out[DTXSID][pos_node][feature]["parent_mass"]:
                                d_out[DTXSID][pos_node][feature]["parent_mass"].append(mass_parent)
        
        
        
        # write in a file by feature ID
        filout = open(p_out, "w")
        filout.write("featureid\tDTXSID\tNURSE_NEG\tNURSE_POS\tFF_NEG\tFF_POS\tMZ\tRT\tNURSE_mean_FB\tNURSE_mean_N\tNURSE_mean_OW\tFF_mean_FB\tFF_mean_visit1\tFF_mean_visit2\tFF_mean_visit3\tRT_mean\tparent_mass\t" + "\t".join(l_col_target_chem) + "\n")
        
        l_mode = ["NURSE_POS", "NURSE_NEG", "FF_POS", "FF_NEG"]
        
        for DTXSID in d_out.keys():
            flag_mode = 0
            for mode in l_mode:
                if d_out[DTXSID][mode] == {}:
                    continue
                else:
                    flag_mode = 1

                for featureid in d_out[DTXSID][mode].keys():
                    if mode == "NURSE_NEG":
                        filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(featureid, DTXSID, 1, 0, 0, 0, d_out[DTXSID][mode][featureid]["MZ"], d_out[DTXSID][mode][featureid]["RT"],
                                                                                     d_out[DTXSID][mode][featureid]["NURSE_mean_FB"], d_out[DTXSID][mode][featureid]["NURSE_mean_N"], d_out[DTXSID][mode][featureid]["NURSE_mean_OW"], d_out[DTXSID][mode][featureid]["FF_mean_FB"],
                                                                                     d_out[DTXSID][mode][featureid]["FF_mean_visit1"], d_out[DTXSID][mode][featureid]["FF_mean_visit2"], d_out[DTXSID][mode][featureid]["FF_mean_visit3"],
                                                                                     ",".join(d_out[DTXSID][mode][featureid]["RT_mean"]), ",".join(d_out[DTXSID][mode][featureid]["parent_mass"]),
                                                                                     "\t".join([d_out[DTXSID]["target_chem"][col_chem_target] for col_chem_target in l_col_target_chem])))
                    elif mode == "NURSE_POS":
                        filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(featureid, DTXSID, 0, 1, 0, 0, d_out[DTXSID][mode][featureid]["MZ"], d_out[DTXSID][mode][featureid]["RT"],
                                                                                     d_out[DTXSID][mode][featureid]["NURSE_mean_FB"], d_out[DTXSID][mode][featureid]["NURSE_mean_N"], d_out[DTXSID][mode][featureid]["NURSE_mean_OW"], d_out[DTXSID][mode][featureid]["FF_mean_FB"],
                                                                                     d_out[DTXSID][mode][featureid]["FF_mean_visit1"], d_out[DTXSID][mode][featureid]["FF_mean_visit2"], d_out[DTXSID][mode][featureid]["FF_mean_visit3"],
                                                                                     ",".join(d_out[DTXSID][mode][featureid]["RT_mean"]), ",".join(d_out[DTXSID][mode][featureid]["parent_mass"]),
                                                                                     "\t".join([d_out[DTXSID]["target_chem"][col_chem_target] for col_chem_target in l_col_target_chem])))
                    elif mode == "FF_POS":
                        filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(featureid, DTXSID, 0, 0, 0, 1, d_out[DTXSID][mode][featureid]["MZ"], d_out[DTXSID][mode][featureid]["RT"],
                                                                                     d_out[DTXSID][mode][featureid]["NURSE_mean_FB"], d_out[DTXSID][mode][featureid]["NURSE_mean_N"], d_out[DTXSID][mode][featureid]["NURSE_mean_OW"], d_out[DTXSID][mode][featureid]["FF_mean_FB"],
                                                                                     d_out[DTXSID][mode][featureid]["FF_mean_visit1"], d_out[DTXSID][mode][featureid]["FF_mean_visit2"], d_out[DTXSID][mode][featureid]["FF_mean_visit3"],
                                                                                     ",".join(d_out[DTXSID][mode][featureid]["RT_mean"]), ",".join(d_out[DTXSID][mode][featureid]["parent_mass"]),
                                                                                     "\t".join([d_out[DTXSID]["target_chem"][col_chem_target] for col_chem_target in l_col_target_chem])))
                    elif mode == "FF_NEG":
                        filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(featureid, DTXSID, 0, 0, 1, 0, d_out[DTXSID][mode][featureid]["MZ"], d_out[DTXSID][mode][featureid]["RT"],
                                                                                     d_out[DTXSID][mode][featureid]["NURSE_mean_FB"], d_out[DTXSID][mode][featureid]["NURSE_mean_N"], d_out[DTXSID][mode][featureid]["NURSE_mean_OW"], d_out[DTXSID][mode][featureid]["FF_mean_FB"],
                                                                                     d_out[DTXSID][mode][featureid]["FF_mean_visit1"], d_out[DTXSID][mode][featureid]["FF_mean_visit2"], d_out[DTXSID][mode][featureid]["FF_mean_visit3"],
                                                                                     ",".join(d_out[DTXSID][mode][featureid]["RT_mean"]), ",".join(d_out[DTXSID][mode][featureid]["parent_mass"]),
                                                                                     "\t".join([d_out[DTXSID]["target_chem"][col_chem_target] for col_chem_target in l_col_target_chem])))
            
            print(DTXSID, flag_mode)
            if flag_mode == 0:
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("NA", DTXSID, 0, 0, 0, 0, "NA", "NA",
                                                                                     "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "\t".join([d_out[DTXSID]["target_chem"][col_chem_target] for col_chem_target in l_col_target_chem])))
                
        
        filout.close()
        