from os import path
import WWBC_database
import pathFolder


# define folder #
#################
PR_ROOT = path.abspath("../../") + "/"
P_DATA = PR_ROOT + "data/"
pr_out = pathFolder.createFolder(PR_ROOT + "results/")

# SET OF CRITERIA #
###################

minMW = 100
maxMW = 1000
lipinski_violation = 3
list_chemicals_metabolite = ["Drug_UCSF_PXYS", "Drug_most comon and haz", "Disinfectant", "FRs", "PFAS", "MC", "MGDev", "ERactive_bin", "E2Up_bin", "P4Up_bin", "pesticidemammarytumors_bin", "nitroPAH_bin", "pellizzari_bin", "phthalate"]


#  RUN  #
#########
c_db = WWBC_database.WWWBC_database(P_DATA + "WWBC_MS_database_4.7.21.csv", minMW, maxMW, lipinski_violation, list_chemicals_metabolite, pr_out)
c_db.main()