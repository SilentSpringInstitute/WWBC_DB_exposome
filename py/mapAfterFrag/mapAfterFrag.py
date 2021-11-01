from toolbox import toolbox

class mapAfterFrag:
    def __init__(self, p_result_frag, p_input_frag, mode, pr_out):

        self.p_result_frag = p_result_frag
        self.p_input_frag = p_input_frag
        self.pr_out = pr_out
        self.mode = mode

        self.l_input_for_frag = toolbox.loadMatrixToList(self.p_input_frag)
        self.d_frag_results = toolbox.loadMatrix(self.p_result_frag)


    def map(self):

        p_filout = self.pr_out + "mapping_" + self.mode + ".csv"
        filout = open(p_filout, "w")
        filout.write("ConsensusID\tFoundInFragInput\tname\tESI\n")

        l_in = []
        for cluster_index in self.d_frag_results.keys():
            l_consensus_ID = self.d_frag_results[cluster_index]["ConsensusID"].split(",")
            for consensus_ID in l_consensus_ID:
                if consensus_ID in l_in:
                    continue
                flag = 0
                for d_input_frag in self.l_input_for_frag:
                    l_db_name = d_input_frag["DB_name"].split(";")
                    for db_name in l_db_name:
                        if db_name == consensus_ID:
                            filout.write("%s\t1\t%s\t%s\n"%(consensus_ID, d_input_frag["name"], d_input_frag["ESI"])) 
                            flag = 1
                            l_in.append(consensus_ID)
                    if flag == 1 :
                        break
                if flag == 0:
                        filout.write("%s\t0\t%s\t%s\n"%(consensus_ID, "NA", "NA"))
        filout.close()

        return p_filout