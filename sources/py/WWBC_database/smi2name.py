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
    
    try: compounds = pubchempy.get_compounds(smiles, namespace='smiles')
    except: return "None"
    try:match = compounds[0]
    except: return "Noone"
    name = match.iupac_name
    if name == None:
        return "None"
    else:
        return name
