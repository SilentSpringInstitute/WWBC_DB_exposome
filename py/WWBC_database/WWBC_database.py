from os import name, path
from random import shuffle
from Desc1D2D import molproperty
from rdkit import Chem
from copy import deepcopy
from re import search

import toolbox
import CompDesc
import pathFolder
import smi2name
import Biotransformer
import prepSMILES


class WWWBC_database:
    def __init__(self, p_database, minMW, maxMW, lipinski_violation, list_chemicals_metabolite, pr_out, name_DB = ""):

        self.p_database = p_database
        if name_DB == "":
            name_DB = p_database.split("/")[-1][:-4]
        self.name_DB = name_DB
        self.pr_out = pr_out
        self.minMW = minMW
        self.maxMW = maxMW
        self.lipinski_violation = lipinski_violation
        self.list_chemicals_metabolite = list_chemicals_metabolite
        self.pr_desc = pathFolder.createFolder(pr_out + "DESC/")
    
    def main(self):
        self.loadData()
        self.cleanStructure()
        self.processMW(self.minMW, self.maxMW)

        # compute biotransformation by chemicals
        pr_biotransformation = pathFolder.createFolder(self.pr_out + "biotransformer_output/")
        #self.computeBiotransformation("phaseII", pr_biotransformation)

        # organize by type source
        self.organizeBiotransformationByChemSources(pr_biotransformation,  self.list_chemicals_metabolite)
        self.processMetaboliteFilter(self.minMW, self.maxMW, self.lipinski_violation)
        self.fusionDBAndMetabolite()

    def loadData(self):

        d_data = toolbox.loadMatrix(self.p_database, ",")
        #self.d_data ={k: d_data[k] for k in list(d_data)[:100]}
        self.d_data =  d_data      

    def cleanStructure(self):

        self.pr_prepDB = pathFolder.createFolder(self.pr_out + "prepDB/")

        p_filout = self.pr_prepDB + self.name_DB + "_clean_structure.csv"
        p_flog = self.pr_prepDB + self.name_DB + "_clean_structure.log"

        if path.exists(p_filout):
            d_prep = toolbox.loadMatrix(p_filout)
            self.d_prep = d_prep
        else:
            l_chem = list(self.d_data.keys())
            i = 0
            imax = len(l_chem)
            self.d_prep = {}
            while i < imax:
                if i % 100 == 0:
                    print(i)
                
                # search for phthalate
                if search("phthalate", self.d_data[l_chem[i]]["Compound"].lower()):
                    self.d_data[l_chem[i]]["phthalate"] = "1"
                else:
                    self.d_data[l_chem[i]]["phthalate"] = "0"


                # need to remove duplicate in SMILES
                l_SMILES = self.d_data[l_chem[i]]["SMILES"].split(", ")
                l_formula = self.d_data[l_chem[i]]["Formula"].split(", ")
                
                if len(l_SMILES) > 1 or len(l_formula) > 1:
                    l_SMILES = list(set(l_SMILES))
                    l_formula = list(set(l_formula))
                    if len(l_SMILES) == 1 or len(l_formula) == 1: # case of mixture duplicate
                        self.d_data[l_chem[i]]["Formula"] = self.d_data[l_chem[i]]["Formula"].split(", ")[0]
                        self.d_data[l_chem[i]]["CASRN"] = self.d_data[l_chem[i]]["CASRN"].split(", ")[0]
                        self.d_data[l_chem[i]]["Compound"] = self.d_data[l_chem[i]]["Compound"].split(", ")[0]
                        self.d_data[l_chem[i]]["SMILES"] = self.d_data[l_chem[i]]["SMILES"].split(", ")[0]
                        self.d_data[l_chem[i]]["MONOISOTOPIC_MASS"] = self.d_data[l_chem[i]]["MONOISOTOPIC_MASS"].split(", ")[0]

                    else:# case duplicate the entry
                        j = 1
                        jmax = len(l_SMILES)
                        while j < jmax:
                            k_add = str(imax + 1)
                            imax = imax + 1
                            l_chem.append(k_add)
                            self.d_data[k_add] = deepcopy(self.d_data[l_chem[i]])
                            self.d_data[k_add]["ID"] = k_add
                            try:self.d_data[k_add]["Formula"] = list(set(self.d_data[k_add]["Formula"].split(", ")))[j]
                            except:self.d_data[k_add]["Formula"] = list(set(self.d_data[k_add]["Formula"].split(", ")))[0]
                            try:self.d_data[k_add]["CASRN"] = list(set(self.d_data[k_add]["CASRN"].split(", ")))[j]
                            except:self.d_data[k_add]["CASRN"] = list(set(self.d_data[k_add]["CASRN"].split(", ")))[0]
                            try:self.d_data[k_add]["Compound"] = list(set(self.d_data[k_add]["Compound"].split(", ")))[j]
                            except:self.d_data[k_add]["Compound"] = list(set(self.d_data[k_add]["Compound"].split(", ")))[0]
                            try:self.d_data[k_add]["SMILES"] = list(set(self.d_data[k_add]["SMILES"].split(", ")))[j]
                            except:self.d_data[k_add]["SMILES"] = list(set(self.d_data[k_add]["SMILES"].split(", ")))[0]
                            try:self.d_data[k_add]["MONOISOTOPIC_MASS"] = list(set(self.d_data[k_add]["MONOISOTOPIC_MASS"].split(", ")))[j]
                            except:self.d_data[k_add]["MONOISOTOPIC_MASS"] = list(set(self.d_data[k_add]["MONOISOTOPIC_MASS"].split(", ")))[0]
                            self.d_data[k_add]["DTXSID"] = self.d_data[k_add]["DTXSID"] + "_%s"%(j + 1)
                            j = j + 1
                        
                        self.d_data[l_chem[i]]["Formula"] = list(set(self.d_data[k_add]["Formula"].split(", ")))[0]
                        self.d_data[l_chem[i]]["CASRN"] = list(set(self.d_data[l_chem[i]]["CASRN"].split(", ")))[0]
                        self.d_data[l_chem[i]]["Compound"] = list(set(self.d_data[l_chem[i]]["Compound"].split(", ")))[0]
                        self.d_data[l_chem[i]]["SMILES"] = list(set(self.d_data[l_chem[i]]["SMILES"].split(", ")))[0]
                        self.d_data[l_chem[i]]["MONOISOTOPIC_MASS"] = list(set(self.d_data[l_chem[i]]["MONOISOTOPIC_MASS"].split(", ")))[0]
                        self.d_data[l_chem[i]]["DTXSID"] =  self.d_data[l_chem[i]]["DTXSID"] + "_1"       
                            
                # change the smiles prep based on the new RDKIT update ==> https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
                c_chem = CompDesc.CompDesc(self.d_data[l_chem[i]]["SMILES"], self.pr_desc)
                c_chem.prepChem()
                if not l_chem[i] in list(self.d_prep.keys()):
                    self.d_prep[l_chem[i]] = {}
                
                if c_chem.err == 0:
                    # Apply new prep phase
                    smi_mol_cleaned = prepSMILES.prepSMILES(c_chem.mol)
                    c_chem.smi = smi_mol_cleaned[0]
                    c_chem.mol = smi_mol_cleaned[1]
                    self.d_prep[l_chem[i]]["SMILES_cleaned"] = c_chem.smi
                    self.d_prep[l_chem[i]]["formula"] = Chem.rdMolDescriptors.CalcMolFormula(c_chem.mol)
                    self.d_prep[l_chem[i]]["name_cleaned"] = smi2name.pubchempySmiles2name(c_chem.smi)
                    
                    # compute only MW not all other descriptor
                    try:MW = molproperty.getExactMolWt(c_chem.mol)
                    except: c_chem.err = 1
                    #c_chem.computeAll2D()
                    if c_chem.err == 0:
                        self.d_prep[l_chem[i]]["Molweight_cleaned"] = MW
                    else:
                        self.d_prep[l_chem[i]]["Molweight_cleaned"] = ""
                else:
                    self.d_prep[l_chem[i]]["SMILES_cleaned"] = ""
                    self.d_prep[l_chem[i]]["Molweight_cleaned"] = ""
                    self.d_prep[l_chem[i]]["name_cleaned"] = ""
                    self.d_prep[l_chem[i]]["formula"] = ""
                i = i + 1 

            l_prop_data = list(self.d_data[list(self.d_data.keys())[0]].keys())
            l_prop_data.remove("DTXSID")
            l_prop_data.remove("Compound")
            l_prop_data.remove("Formula")
            l_prop_data.remove("CASRN")
            l_prop_data.remove("SMILES")
            l_prop_data.remove("ID")


            filout = open(p_filout, "w")
            filout.write("ID\tDTXSID\tCASRN\tname_original\tformula\tSMILES\tSMILES_cleaned\tname_cleaned\tformula_cleaned\tMolweight_cleaned\t%s\n"%("\t".join(l_prop_data)))
            for chem in self.d_prep.keys():
                # reformate CAS properly
                if search("/", self.d_data[chem]["CASRN"]):
                    self.d_data[chem]["CASRN"] = toolbox.searchCASRNFromDTXSID(self.d_data[chem]["DTXSID"])
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.d_data[chem]["ID"], self.d_data[chem]["DTXSID"],  self.d_data[chem]["CASRN"], self.d_data[chem]["Compound"], self.d_data[chem]["Formula"], self.d_data[chem]["SMILES"], self.d_prep[chem]["SMILES_cleaned"], self.d_prep[chem]["name_cleaned"], self.d_prep[chem]["formula"], self.d_prep[chem]["Molweight_cleaned"], "\t".join(str(self.d_data[chem][k_prop]) for k_prop in l_prop_data)))
            filout.close()

            filog = open(p_flog, "w")
            filog.write("ID\tDTXSID\tCASRN\tname_original\tformula\tSMILES\tSMILES_cleaned\tname_cleaned\tformula_cleaned\tMolweight_cleaned\t%s\n"%("\t".join(l_prop_data)))
            for chem in self.d_prep.keys():
                if self.d_prep[chem]["Molweight_cleaned"] == "":
                    filog.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.d_data[chem]["ID"], self.d_data[chem]["DTXSID"],  self.d_data[chem]["CASRN"], self.d_data[chem]["Compound"], self.d_data[chem]["Formula"], self.d_data[chem]["SMILES"], self.d_prep[chem]["SMILES_cleaned"], self.d_prep[chem]["name_cleaned"], self.d_prep[chem]["formula"], self.d_prep[chem]["Molweight_cleaned"], "\t".join(str(self.d_data[chem][k_prop]) for k_prop in l_prop_data)))
            filog.close()

    def computeBiotransformation(self, type_biotransformation, pr_biotransformer):

        #p_filout = self.pr_out + "biotransformation_phaseII.csv"
        #filout = open(p_filout, "w")

        d_transformation = {}
        l_chem =  list(self.d_prep.keys())
        shuffle(l_chem)
        i = 0
        for chem in l_chem:
            if i % 100 == 0:
                print(i)
            i = i + 1
            if self.d_prep[chem]["SMILES_cleaned"] == "":
                continue
            c_chem = CompDesc.CompDesc(self.d_prep[chem]["SMILES_cleaned"], "")
            c_chem.smi = self.d_prep[chem]["SMILES_cleaned"]
            c_chem.mol = Chem.MolFromSmiles(self.d_prep[chem]["SMILES_cleaned"])
            c_chem.generateInchiKey()
            if c_chem.err == 1:
                continue
            c_biotransformation = Biotransformer.Biotransformer(self.d_prep[chem]["SMILES_cleaned"], c_chem.inchikey, pr_biotransformer)
            d_transformation[chem] = c_biotransformation.predBiotransformer(type_biotransformation)
        
    def organizeBiotransformationByChemSources(self, pr_biotransformation, l_chem_classes):

        pr_out = pathFolder.createFolder(self.pr_out + "prepPhaseII/")
        self.pr_metabolite = pr_out
        p_filout = pr_out + "PhaseII_metabo_" + self.name_DB + ".csv"
        if path.exists(p_filout):
            self.d_metabolite = toolbox.loadMatrix(p_filout)
            return 

        d_out = {}
        l_chem = list(self.d_prep.keys())
        i = 0
        imax = len(l_chem)
        while i < imax:
            if i % 100 == 0:
                print(i)

            for chem_class in l_chem_classes:
                if self.d_prep[l_chem[i]][chem_class] == "1":
                    SMILES_clean = self.d_prep[l_chem[i]]["SMILES_cleaned"]
                    if SMILES_clean == "":
                        i = i + 1 
                        continue
                    c_chem = CompDesc.CompDesc(SMILES_clean, "")
                    c_chem.smi = SMILES_clean
                    c_chem.mol = Chem.MolFromSmiles(SMILES_clean)
                    c_chem.generateInchiKey()
                    if c_chem.err == 1:
                        i = i + 1
                        continue
                    inchikey = c_chem.inchikey
                    d_biotransformation = toolbox.loadMatrix(pr_biotransformation + inchikey + "_phaseII.csv", ",")
                    if "NO TRANSFORMATION" in list(d_biotransformation.keys()):
                        i = i + 1
                        continue
                    else:
                        for inch in d_biotransformation.keys():
                            try:SMILES_metabo = d_biotransformation[inch]["SMILES"]
                            except: continue
                            # clean SMILES from metabo
                            c_chem = CompDesc.CompDesc(SMILES_metabo, self.pr_desc)
                            c_chem.prepChem()
                            if c_chem.err == 1:
                                continue
                            else:
                                SMILES_metabo_cleaned = c_chem.smi
                                DTXSID = toolbox.searchDTXIDFromSMILES(SMILES_metabo_cleaned)
                            if not SMILES_metabo_cleaned in list(d_out.keys()):
                                d_out[SMILES_metabo_cleaned] = {}
                                d_out[SMILES_metabo_cleaned] = d_biotransformation[inch]
                                d_out[SMILES_metabo_cleaned]["name_cleaned"] = smi2name.pubchempySmiles2name(SMILES_metabo_cleaned)
                                d_out[SMILES_metabo_cleaned]["DTXSID"] = DTXSID
                                d_out[SMILES_metabo_cleaned]["Precursors"] = []
                                d_out[SMILES_metabo_cleaned]["List Reaction"] = []
                            if not self.d_data[l_chem[i]]["DTXSID"] in d_out[SMILES_metabo_cleaned]["Precursors"]:
                                d_out[SMILES_metabo_cleaned]["Precursors"].append(self.d_data[l_chem[i]]["DTXSID"])
                                d_out[SMILES_metabo_cleaned]["List Reaction"].append(d_biotransformation[inch]["Reaction ID"])
            
            i = i + 1

        
        filout = open(p_filout, "w")
        filout.write("name\tDTXSID\tSMILES_metabolite\tname_cleaned\tPUBCHEM_CID\tMolecular formula\tMajor Isotope Mass\tALogP\tLipinski_Violations\tMetabolite ID\tReaction\tReaction ID\tEnzyme(s)\tPrecursors\n")
        for chem_metabo in d_out.keys():
            chem_name = "__".join(["%s-%s"%(d_out[chem_metabo]["Precursors"][i], d_out[chem_metabo]["List Reaction"][i]) for i in range(0, len(d_out[chem_metabo]["Precursors"]))])
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chem_name, d_out[chem_metabo]["DTXSID"], chem_metabo, d_out[chem_metabo]["name_cleaned"], d_out[chem_metabo]["PUBCHEM_CID"], d_out[chem_metabo]["Molecular formula"], d_out[chem_metabo]["Major Isotope Mass"], d_out[chem_metabo]["ALogP"], d_out[chem_metabo]["Lipinski_Violations"], d_out[chem_metabo]["Metabolite ID"], d_out[chem_metabo]["Reaction"], d_out[chem_metabo]["Reaction ID"], d_out[chem_metabo]["Enzyme(s)"], ",".join(d_out[chem_metabo]["Precursors"])))
        filout.close()

        self.f_metabolite = p_filout
        self.d_metabolite = toolbox.loadMatrix(p_filout)

    def processMetaboliteFilter(self, minMW, maxMW, maxLipinskiViolation):
        
        p_filout =  "%sPhaseII_metabo_MW-%s-%s_Lip%s_%s.csv"%(self.pr_metabolite, minMW, maxMW, maxLipinskiViolation, self.name_DB)
        p_flog = "%sPhaseII_metabo_MW-%s-%s_Lip%s_%s.log"%(self.pr_metabolite, minMW, maxMW, maxLipinskiViolation, self.name_DB)

        if path.exists(p_filout) and path.exists(p_flog):
            self.d_metabolite_cleaned = toolbox.loadMatrix(p_filout)
            return 

        l_out = []
        l_in = []
        l_chem = list(self.d_metabolite.keys())
        i = 0
        imax = len(l_chem)
        while i < imax:
            MW = self.d_metabolite[l_chem[i]]["Major Isotope Mass"]
            lipinski = int(self.d_metabolite[l_chem[i]]["Lipinski_Violations"])
            try:
                MW = float(MW)
            except:
                l_out.append(l_chem[i])
                i = i + 1 
                continue

            if MW >= minMW and MW <= maxMW:
                if lipinski <= maxLipinskiViolation:
                    l_in.append(l_chem[i])
                else:
                    l_out.append(l_chem[i])
            else:
                l_out.append(l_chem[i])

            i = i + 1
        
        l_prop_data = list(self.d_metabolite[list(self.d_metabolite.keys())[0]].keys())
        l_prop_data.remove("name")
        l_prop_data.remove("SMILES_metabolite")

        filout = open(p_filout, "w")
        filout.write("name\tSMILES_metabolite\t%s\n"%("\t".join(l_prop_data)))
        for chem in l_in:
            filout.write("%s\t%s\t%s\n"%(chem, self.d_metabolite[chem]["SMILES_metabolite"], "\t".join(str(self.d_metabolite[chem][k_prop]) for k_prop in l_prop_data)))
        filout.close()

        flog = open(p_flog, "w")
        flog.write("name\tSMILES_metabolite\t%s\n"%("\t".join(l_prop_data)))
        for chem in l_out:
            flog.write("%s\t%s\t%s\n"%(chem, self.d_metabolite[chem]["SMILES_metabolite"], "\t".join(str(self.d_metabolite[chem][k_prop]) for k_prop in l_prop_data)))
        flog.close()    

        self.d_metabolite_cleaned = toolbox.loadMatrix(p_filout)
            
    def processMW(self, minMW, maxMW):
        """
        Add column DBname
        """
        
        p_filout = self.pr_prepDB + self.name_DB + "_filterMW.csv"
        p_flog = self.pr_prepDB + self.name_DB + "_filterMW.log"

        if path.exists(p_filout) and path.exists(p_flog):
            self.d_DB_cleaned = toolbox.loadMatrix(p_filout)
            return 

        l_out = []
        l_in = []
        l_chem = list(self.d_data.keys())
        i = 0
        imax = len(l_chem)
        while i < imax:
            MW = self.d_prep[l_chem[i]]["Molweight_cleaned"]
            if MW == "":
                MW = self.d_data[l_chem[i]]["MONOISOTOPIC_MASS"]
            try:
                MW = float(MW)
            except:
                l_out.append(l_chem[i])
                i = i + 1 
                continue

            if MW >= minMW and MW <= maxMW:
                l_in.append(l_chem[i])
            else:
                l_out.append(l_chem[i])

            i = i + 1
        
        l_prop_data = list(self.d_data[list(self.d_data.keys())[0]].keys())
        l_prop_data.remove("DTXSID")
        l_prop_data.remove("Compound")
        l_prop_data.remove("Formula")
        l_prop_data.remove("CASRN")
        l_prop_data.remove("SMILES")
        l_prop_data.remove("ID")



        filout = open(p_filout, "w")
        filout.write("ID\tDTXSID\tCASRN\tname_original\tformula\tSMILES\tSMILES_cleaned\tname_cleaned\tformula_cleaned\tMolweight_cleaned\t%s\n"%("\t".join(l_prop_data)))
        for chem in l_in:
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.d_data[chem]["ID"], self.d_data[chem]["DTXSID"],  self.d_data[chem]["CASRN"], self.d_data[chem]["Compound"].replace(" ", ","), self.d_data[chem]["Formula"], self.d_data[chem]["SMILES"], self.d_prep[chem]["SMILES_cleaned"], self.d_prep[chem]["name_cleaned"], self.d_prep[chem]["formula"], self.d_prep[chem]["Molweight_cleaned"], "\t".join(str(self.d_data[chem][k_prop]) for k_prop in l_prop_data)))
        filout.close()

        flog = open(p_flog, "w")
        flog.write("ID\tDTXSID\tCASRN\tname_original\tformula\tSMILES\tSMILES_cleaned\tname_cleaned\tformula_cleaned\tMolweight_cleaned\t%s\n"%("\t".join(l_prop_data)))
        for chem in l_out:
            flog.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.d_data[chem]["ID"], self.d_data[chem]["DTXSID"],  self.d_data[chem]["CASRN"], self.d_data[chem]["Compound"].replace(" ", ","), self.d_data[chem]["Formula"], self.d_data[chem]["SMILES"], self.d_prep[chem]["SMILES_cleaned"], self.d_prep[chem]["name_cleaned"], self.d_prep[chem]["formula"], self.d_prep[chem]["Molweight_cleaned"], "\t".join(str(self.d_data[chem][k_prop]) for k_prop in l_prop_data)))
        flog.close()        

        self.d_DB_cleaned = toolbox.loadMatrix(p_filout)

    def fusionDBAndMetabolite(self):

        l_chem_DB = list(self.d_DB_cleaned.keys())
        l_chem_metabo = list(self.d_metabolite_cleaned.keys())


        l_prop_data = list(self.d_DB_cleaned[list(self.d_DB_cleaned.keys())[0]].keys())
        l_prop_data.remove("DTXSID")
        l_prop_data.remove("CASRN")
        l_prop_data.remove("name_original")
        l_prop_data.remove("formula")
        l_prop_data.remove("SMILES_cleaned")
        l_prop_data.remove("name_cleaned")
        l_prop_data.remove("SMILES")
        l_prop_data.remove("formula_cleaned")
        l_prop_data.remove("Molweight_cleaned")
        l_prop_data.remove("ID")



        p_filout = self.pr_out + self.name_DB + "_prepForAnotation.csv"
        filout = open(p_filout, "w")
        filout.write("ID\tDB_name\tDTXSID\tCASRN\tname_original\tformula_original\tSMILES_original\tSMILES_cleaned\tname_cleaned\tformula_cleaned\tMolweight_cleaned\tparent_id\t%s\n"%("\t".join(l_prop_data)))

        i = 0
        ID = 1
        imax = len(l_chem_DB)

        while i < imax:
            filout.write("%s\t%s\t%s\t'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ID, self.d_DB_cleaned[l_chem_DB[i]]["DTXSID"], self.d_DB_cleaned[l_chem_DB[i]]["DTXSID"], self.d_DB_cleaned[l_chem_DB[i]]["CASRN"], self.d_DB_cleaned[l_chem_DB[i]]["name_original"], self.d_DB_cleaned[l_chem_DB[i]]["formula"], self.d_DB_cleaned[l_chem_DB[i]]["SMILES"], self.d_DB_cleaned[l_chem_DB[i]]["SMILES_cleaned"], self.d_DB_cleaned[l_chem_DB[i]]["name_cleaned"], self.d_DB_cleaned[l_chem_DB[i]]["formula"], self.d_DB_cleaned[l_chem_DB[i]]["Molweight_cleaned"], "", "\t".join(str(self.d_DB_cleaned[l_chem_DB[i]][k_prop]) for k_prop in l_prop_data)))
            i = i + 1
            ID = ID + 1

        i = 0
        imax = len(l_chem_metabo)
        while i < imax:
            if self.d_metabolite_cleaned[l_chem_metabo[i]]["name_cleaned"] == "None":
                self.d_metabolite_cleaned[l_chem_metabo[i]]["name_cleaned"] = ""
            filout.write("%s\t%s\t%s\t'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ID, self.d_metabolite_cleaned[l_chem_metabo[i]]["name"], "MTX%s"%(i) , "", "", self.d_metabolite_cleaned[l_chem_metabo[i]]["Molecular formula"], self.d_metabolite_cleaned[l_chem_metabo[i]]["SMILES_metabolite"], self.d_metabolite_cleaned[l_chem_metabo[i]]["SMILES_metabolite"], self.d_metabolite_cleaned[l_chem_metabo[i]]["name_cleaned"], self.d_metabolite_cleaned[l_chem_metabo[i]]["Molecular formula"], self.d_metabolite_cleaned[l_chem_metabo[i]]["Major Isotope Mass"],self.d_metabolite_cleaned[l_chem_metabo[i]]["Precursors"], "\t".join("metabolite" if k_prop == "source" else  "" for k_prop in l_prop_data)))
            i = i + 1
            ID = ID + 1


        filout.close()