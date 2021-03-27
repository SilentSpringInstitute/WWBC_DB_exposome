from os import path
import runExternal
import toolbox

class Biotransformer:
    def __init__(self, smi_in, inchikey, pr_out):
        self.smi_in = smi_in
        self.inchikey = inchikey
        self.pr_out = pr_out

    def predBiotransformer(self, type_transfomation):

        p_filout = "%s%s_%s.csv"%(self.pr_out, self.inchikey, type_transfomation)
        if not path.exists(p_filout):
            runExternal.BioTransformer(self.smi_in, type_transfomation, p_filout)

            if not path.exists(p_filout):
                filout = open(p_filout, "w")
                filout.write("Log\nNO TRANSFORMATION")
                filout.close()

       