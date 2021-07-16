import toolbox
from rdkit import Chem

import CompDesc

class matchChemicals:

    def __init__(self, p_chem, name, pr_out):
        self.p_chem = p_chem
        self.name = name
        self.pr_out = pr_out

    def loadChem(self, name_sheet):
        d_chem = toolbox.loadExcelSheet(self.p_chem, name_sheet, "DTXSID")
        self.d_chem = d_chem

    def searchSubstructure(self, d_substructure):

        p_filout = self.pr_out + self.name  + ".csv"

        d_out = {}
        for dtxsid in self.d_chem.keys():
            try:c_chem = Chem.MolFromSmiles(self.d_chem[dtxsid]["SMILES_cleaned"])
            except:continue
            d_out[dtxsid] = {}
            for sub_name in d_substructure.keys():
                d_out[dtxsid][sub_name] = "0"
                for sub_smiles in d_substructure[sub_name]:
                    c_sub = Chem.MolFromSmiles(sub_smiles)
                    
                    # do the search
                    in_chem = c_chem.HasSubstructMatch(c_sub)
                    if in_chem == True:
                        d_out[dtxsid][sub_name] = "1"


        filout = open(p_filout, "w")
        filout.write("DTXSID\tCASRN\tname_original\t%s\n"%("\t".join(list(d_substructure.keys()))))
        for dtxsid in d_out.keys():
            filout.write("%s\t%s\t%s\t%s\n"%(dtxsid, self.d_chem[dtxsid]["CASRN"], self.d_chem[dtxsid]["name_original"], "\t".join([d_out[dtxsid][sub] for sub in d_substructure.keys()])))
        filout.close()
