from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
uncharger = rdMolStandardize.Uncharger()

### need to be integrated in the COMPDESC lib

def prepSMILES(m):

    #m = Chem.MolFromSmiles(smi_in,sanitize=False)
    #m.UpdatePropertyCache(strict=False)
    #Chem.SanitizeMol(m,sanitizeOps=(Chem.SANITIZE_ALL^Chem.SANITIZE_CLEANUP^Chem.SANITIZE_PROPERTIES))
    #cm = rdMolStandardize.Normalize(m)
    cm = uncharger.uncharge(m)
    smi_out = Chem.MolToSmiles(cm)

    return [smi_out, cm]
