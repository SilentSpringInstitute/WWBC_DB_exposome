import urllib.request
import pubchempy


def CIRconvert(smi):
    try:
        url ="https://cactus.nci.nih.gov/chemical/structure/" + smi+"/iupac_name" 
        ans = urllib.request.urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Name Not Available'

#smiles  = 'CCCCC(C)CC'
#print(smiles, CIRconvert(smiles))


def pubchempySmiles2name(smiles):

    try:
        compounds = pubchempy.get_compounds(smiles, namespace='smiles')
        match = compounds[0]
        return match.iupac_name
    except:
        return "-"
