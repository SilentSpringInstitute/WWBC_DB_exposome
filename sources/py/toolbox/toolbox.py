import toolbox
from .DBrequest import DBrequest
import pandas
from re import search


def loadMatrixToList(pmatrixIn, sep = "\t", header_out = False):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    l_out = []
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["ID"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "ID"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        j = 0
        if len(lvalues) != len(lheaders):
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
        jmax = len(lheaders)
        dtemp = {}
        while j < jmax:
            try:dtemp[lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        l_out.append(dtemp)
        i += 1
    if header_out == True:
        return [l_out, lheaders]
    else:
        return l_out

def loadMatrix(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1) -1):
        lheaders = ["X"] + lheaders


    i = 0
    while i < len(lheaders):
        if lheaders[i] == "":
            lheaders[i] = "X"
        i += 1

    i = 1
    imax = len(llinesMat)
    while i < imax:
        if llinesMat[1][0] == "\"" and i + 1 < imax and llinesMat[i + 1][0] != "\"":
            llinesMat[i] = llinesMat[i].strip() + " " +  llinesMat[i+1]

        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        if len(lvalues) != len(lheaders):
            print("====")
            print(pmatrixIn)
            print("ERROR - line: ", i)
            print(lvalues)
            print(lheaders)
            print("====")

        jmax = len(lheaders)
        while j < jmax:
            try:dout[kin][lheaders[j]] = lvalues[j].replace("--", ",")
            except:pass
            j += 1
        i += 1

    return dout

def formatLine(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + "--"
        else:
            linenew = linenew + linein[i]

        i += 1

    linenew = linenew.replace('\"', "")
    return linenew

def writeMatrix(ddesc, pdescAct, sep = "\t"):


    filout = open(pdescAct, "w")
    lheader = list(ddesc[list(ddesc.keys())[0]].keys())

    # put header in first

    if "CAS" in lheader:
        del lheader[lheader.index("CAS")]
        lheader = ["CAS"] + lheader
    elif "CASID" in lheader:
        del lheader[lheader.index("CASID")]
        lheader = ["CASID"] + lheader
    else:
        lheader = ["CASID"] + lheader
        for casID in list(ddesc.keys()):
            ddesc[casID]["CASID"] = casID


    filout.write(sep.join(lheader) + "\n")
    for casID in list(ddesc.keys()):
        lval = [str(ddesc[casID][i]) for i in lheader]
        filout.write(sep.join(lval) + "\n")
    filout.close()

def formatLineDataset(linein):

    linein = linein.replace("\n", "")
    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == '"' and flagchar == 0:
            flagchar = 1
        elif linein[i] == '"' and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace('\"', "")
    return linenew


def mergeDict(l_dict):

    dout = {}

    lk = list(l_dict[0].keys())

    for k in lk:
        lval = [l_dict[i][k] for i in range(0,len(l_dict))]
        #identic
        if lval.count(lval[0]) == len(lval):
            dout[k] = lval[0]
        else:
            # remove duplicate
            lval = list(dict.fromkeys(lval))
            dout[k] = "----".join(lval)

    return dout

def searchDTXIDFromSMILES(SMILES_cleaned):
    cDB = DBrequest()
    extract_DTXSID = cDB.execCMD("SELECT dsstox_id FROM chemicals WHERE smiles_clean='%s'"%(SMILES_cleaned))

    if extract_DTXSID == []:
        return ""
    else:
        DTXSID = extract_DTXSID[0][0]
        if DTXSID == None:
            return ""
        else:
            return DTXSID

def searchCASRNFromDTXSID(dtxsid):

    cDB = DBrequest()
    extract_CASRN = cDB.execCMD("SELECT casn FROM chemicals WHERE dsstox_id='%s'"%(dtxsid))

    if extract_CASRN == []:
        return ""
    else:
        CASRN = extract_CASRN[0][0]
        if CASRN == None:
            return ""
        else:
            return CASRN  

def loadExcelSheet(p_excel, name_sheet, k_head):
    """
    TO DO: Add check duplicate rownames
    """

    d_out = {}
    # load MC list
    data_frame = pandas.read_excel(p_excel, name_sheet, engine='openpyxl')
    data_size = data_frame.shape
    nb_row = data_size[0]
    nb_col = data_size[1]
    l_col = data_frame.columns

    i = 0
    while i < nb_row:
        rowname = data_frame.iloc[i][k_head]
        if not rowname in list(d_out.keys()):
            d_out[rowname] = {}
        j = 0
        while j < nb_col:
            colname = l_col[j]
            if search("Unnamed:", str(colname)):
                j = j + 1
                continue
            else:
                d_out[rowname][colname] = data_frame.iloc[i][l_col[j]]
            j = j + 1
        i = i + 1
    
    return d_out


def manualMapping(p_filToMap, l_dfile_to_map, p_filout):
    """
    function used to map new function in a previous table
    [{"pfile":"", "headertoMap": "", "addCol": ""}]
    """
    
    l_d_filin = loadMatrixToList(p_filToMap, header_out=True, sep = ",")
    l_header = l_d_filin[1]
    l_d_filin = l_d_filin[0]

    l_header_add = []
    for dfile_to_map in l_dfile_to_map:
        l_header_add.append(dfile_to_map["addCol"])
        l_dtomap = loadMatrixToList(dfile_to_map["pfile"], sep = ",")
        for d_filin in l_d_filin:
            print(d_filin)
            for d_tomap in l_dtomap:
                if d_filin["DTXSID"] == d_tomap[dfile_to_map["headertoMap"]]:
                    d_filin[dfile_to_map["addCol"]] = "1"
                    break
            try: d_filin[dfile_to_map["addCol"]]
            except: d_filin[dfile_to_map["addCol"]] = ""
    
    # rewrite the dataset
    filout = open(p_filout, "w")
    filout.write("\t".join(l_header) + "\t" + "\t".join(l_header_add) + "\n")
    for d_filin in l_d_filin:
        filout.write("%s\t%s\n"%("\t".join([d_filin[h] for h in l_header]), "\t".join([d_filin[h] for h in l_header_add])))
    filout.close()



