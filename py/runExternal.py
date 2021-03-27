from os import name, path, getcwd, system, chdir


if name == "nt":
    PR_BIOTRANSFORMER = "C://Users/aborr/research/PFAS_project/Silent_Spring/BioTransformerJar/biotransformerjar/"
else:
    PR_BIOTRANSFORMER = "/mnt/c/Users/aborr/research/PFAS_project/Silent_Spring/BioTransformerJar/biotransformerjar/"


##############
# run biotransformer tool

def BioTransformer(smi, btType, p_out, nsteps=1):

    p_out = path.abspath(p_out)
    workdir = getcwd()
    chdir(PR_BIOTRANSFORMER)
    if name == "nt":
        cmd = "java.exe -jar ./biotransformer-1.1.5.jar -k pred -b %s -ismi \"%s\" -ocsv %s -s %s"%(btType, smi, p_out, nsteps)
    else:
        cmd = "java -jar ./biotransformer-1.1.5.jar -k pred -b %s -ismi \"%s\" -ocsv %s -s %s"%(btType, smi, p_out, nsteps)
    print("********")
    print(cmd)
    print("***********")
    system(cmd)
    chdir(workdir)