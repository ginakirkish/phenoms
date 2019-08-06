import argparse
import pbr
import os
from glob import glob
from pbr.base import _get_output
import json, uuid
from nipype.utils.filemanip import load_json
heuristic = load_json(os.path.join(os.path.split(pbr.__file__)[0], "heuristic.json"))["filetype_mapper"]
from subprocess import Popen, PIPE, check_call, check_output
import pandas as pd
import shutil

base_dir = '/data/henry12/phenoms/Basel/'


def calc_first(seg, num):
    num1 = num - .5
    num2 = num + .5
    cmd = ["fslstats",seg,"-l", "{}".format(num1), "-u","{}".format(num2), "-V"]
    proc = Popen(cmd, stdout=PIPE)
    area = [l.decode("utf-8").split() for l in proc.stdout.readlines()[:]][0][0]
    print(seg, area)
    return(area)

def get_first_values(first_path):
    L_thal,L_caud,L_put,L_pall,L_hipp, L_amy, L_acc,  R_thal, R_caud, R_put,R_pall, R_hipp, R_amy, R_acc,BS = '','','','','','','','','','','','','','',''
    for files in os.listdir(first_path):
        if files.endswith("firstseg.nii.gz") or files.endswith("firstsegs.nii.gz") and not "_N4_N4" in files:
            print(files)
            seg = first_path + files
            L_thal = calc_first(seg, int(10))
            L_caud = calc_first(seg, 11)
            L_put = calc_first(seg, 12)
            L_pall = calc_first(seg, 13)
            L_hipp = calc_first(seg, 17)
            L_amy = calc_first(seg, 18)
            L_acc = calc_first(seg,26)
            R_thal = calc_first(seg, 49)
            R_caud = calc_first(seg, 50)
            R_put = calc_first(seg, 51)
            R_pall = calc_first(seg, 52)
            R_hipp = calc_first(seg, 53)
            R_amy = calc_first(seg, 54)
            R_acc = calc_first(seg,58)
            BS = calc_first(seg, 16)
            print(L_thal,L_caud,L_put,L_pall,L_hipp, L_amy, L_acc,  R_thal, R_caud, R_put,R_pall, R_hipp, R_amy, R_acc,BS)

    return [L_thal,L_caud,L_put,L_pall,L_hipp, L_amy, L_acc,  R_thal, R_caud, R_put,R_pall, R_hipp, R_amy, R_acc,BS]


def get_sienax(report):
    VS, PG, VCSF, GM,WM,BV = "","","","","",""
    with open(report, "r") as f:
        lines = [line.strip() for line in f.readlines()]
        for line in lines:
            if not len(line) >= 1:
                continue
            if line.startswith("VSCALING"):
                VS = line.split()[1]
            elif line.startswith("pgrey"):
                PG = line.split()[2]
            elif line.startswith("vcsf"):
                VCSF = line.split()[2]
            elif line.startswith("GREY"):
                GM = line.split()[2]
            elif line.startswith("WHITE"):
                WM = line.split()[2]
            elif line.startswith("BRAIN"):
                BV = line.split()[2]

    return [ VS, PG, VCSF, GM, WM, BV]

def get_t1(ID):
    t1 = ""
    subject_path = base_dir +'/'+ ID
    if os.path.exists(subject_path):
        print(subject_path)
        for f in os.listdir(subject_path):
            new_path = subject_path +'/'+ f
            for files in os.listdir(new_path):
                if files.endswith(".nii.gz") and files.startswith("o") and not files.endswith("N4.nii.gz"):
                    #print(files)
                    print("**************",t1)
    return t1

def write_csv(c,out):
    t1 = ""
    df = pd.read_csv("{}".format(c))
    for idx in range(len(df)):

        id = df.loc[idx, 'subjects']
        ID = str(id).split("_")[0]
        tp = str(id).split("_")[1]
        
        print(ID, "*****************")
        df.loc[idx, 'subjects'] = ID
        df.loc[idx, 'timepoint'] = tp
        #t1 = get_t1(ID)

        first_path =  glob('/data/henry12/phenoms/Basel/first_output/{}*/'.format(id))
        print(first_path)

        if len(first_path)>0:
            first_path = first_path[0]
            try:
                t1 = first_path.split(ID)[-1].split("_")[:-1]
                df.loc[idx, "T1"] = t1
                print(t1, "T1")
            except:
                pass
            print("FIRST PATH", first_path)
            F = get_first_values(first_path)
            df.loc[idx, 'Left-Thalamus-Proper new'] = F[0]
            df.loc[idx, 'Left-Caudate new'] = F[1]
            df.loc[idx, 'Left-Putamen new'] = F[2]
            df.loc[idx, 'Left-Pallidum new'] = F[3]
            df.loc[idx, 'Left-Hippocampus new'] = F[4]
            df.loc[idx, 'Left-Amygdala new'] = F[5]
            df.loc[idx, 'Left-Accumbens new'] = F[6]
            df.loc[idx, 'Right-Thalamus-Proper new'] = F[7]
            df.loc[idx, 'Right-Caudate new'] = F[8]
            df.loc[idx, 'Right-Putamen new'] = F[9]
            df.loc[idx, 'Right-Pallidum new'] = F[10]
            df.loc[idx, 'Right-Hippocampus new'] = F[11]
            df.loc[idx, 'Right-Amygdala new'] = F[12]
            df.loc[idx, 'Right-Accumbens new'] = F[13]
            df.loc[idx, 'Brain Stem new'] = F[14]
        sienax_path = glob(base_dir + '/sienax_output/{}*'.format(id))
        if len(sienax_path) >0:
            report = sienax_path[0] + '/report.sienax'
            SX = get_sienax(report)
            df.loc[idx, 'V Scale'] = SX[0]
            df.loc[idx, 'pGM'] = SX[1]
            df.loc[idx, 'CSF'] = SX[2]
            df.loc[idx, 'GM'] = SX[3]
            df.loc[idx, 'WM'] = SX[4]
            df.loc[idx, "BV"] = SX[5]


    df.to_csv('{}'.format(out))





if __name__ == '__main__':
    parser = argparse.ArgumentParser('This code allows you to grab the mse, scanner, birthdate, sex  given a csv (with date and msid)')
    parser.add_argument('-i', help = 'csv containing the msid and date')
    parser.add_argument('-o', help = 'output csv for siena data')
    parser.add_argument
    args = parser.parse_args()
    c = args.i
    out = args.o
    print(c, out)
    write_csv(c, out)