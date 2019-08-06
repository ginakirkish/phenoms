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

base_dir = '/data/henry12/phenoms/OSLO/'

def rename(final_path):
    print(final_path, final_path.replace(")",""))
    new_name = final_path.replace(")","")
    os.rename(final_path, new_name)

def get_brain(t1_path):
    BM = glob(t1_path + 'r*brain*.nii.gz')
    if len(BM) >= 1:
        for b in BM:
            if not "mask" in b:
                BM = b
    else:
        BM = ""
    return BM


def run_sienax(t1, odir):

    if not os.path.exists(odir + '/I_brain.nii.gz'):
        if os.path.exists(t1):
            #cmd = ["sienax_optibet", t1, "-r", "-d", "-o", t1_out ]
            print("sienax_optibet", t1, "-r", "-d", "-o", odir )
        """try:
            Popen(cmd).wait()
        except:
            pass"""
    #else:
        #print(t1_out, "SIENAX FILE EXISTS")


def run_first(t1,odir):

    if not os.path.exists(odir):
        os.mkdir(odir)
    #if len(glob('{}/*first*'.format(odir))) == 0:


        fname = t1.split('/')[-1].split('.nii')[0]
        oname = os.path.join(odir, fname)

        omat = oname+'_to_std_sub'

        flirt_job = ['first_flirt', t1, omat]
        mat_name = omat+'.mat'

        if not os.path.exists(mat_name):
            print(flirt_job)
            check_call(flirt_job)

        imsfirst = []
        imscorr = []
        for s in ['L_Accu', 'L_Amyg', 'L_Caud', 'L_Hipp', 'L_Puta', 'L_Thal', 'R_Accu', 'L_Pall',
                  'R_Amyg', 'R_Caud', 'R_Hipp', 'R_Pall', 'R_Puta', 'R_Thal', 'BrStem']:

            modelN = '336'
            bcorr = '1'
            intref = '0'
            if 'Accu' in s:
                nmodes = '50'
            elif 'Amyg' in s:
                nmodes = '50'
                intref = '1'
            elif 'Caud' in s:
                nmodes = '30'
                intref= '1'
            elif 'Hipp' in s:
                nmodes = '30'
                intref = '1'
            elif 'Late' in s:
                nmodes = '40'
                intref = '1'
            elif 'Pall' in s:
                nmodes = '40'
                bcorr = '0'
            elif 'Puta' in s:
                nmodes = '40'
                bcorr = '0'
            elif 'Thal' in s:
                nmodes = '40'
                bcorr = '0'
            elif 'Stem' in s:
                nomdes = '40'
            else:
                raise ValueError('The structure {} is not in the structure list'.format(s))

            imfirst = '{}-{}_first'.format(oname, s)
            imsfirst.append(imfirst)

            FSLDIR = '/netopt/fsl5'

            if intref == '0':
                if bcorr == '1':
                    cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', t1, '-t', mat_name,
                           '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/{}_bin.bmv'.format(FSLDIR, modelN, s)]
                else:
                     cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', t1, '-t', mat_name,
                           '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/05mm/{}_05mm.bmv'.format(FSLDIR, modelN, s)]
            else:
                cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', t1, '-t', mat_name,
                       '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/intref_thal/{}.bmv'.format(FSLDIR, modelN, s), '-intref', '{}/data/first/models_{}_bin/05mm/{}_Thal_05mm.bmv'.format(FSLDIR, modelN, s.split('_')[0])]
            print('\nRunning Segmentation {} with: {}'.format(s,cmd))
            check_call(cmd)

            imcorr = '{}-{}_corr'.format(oname, s)
            imscorr.append(imcorr)

            btype = 'fast'
            if bcorr != '1': btype='none'
            cmd = ['{}/bin/first_boundary_corr'.format(FSLDIR), '-s', imfirst, '-o', imcorr,
                   '-i', t1, '-b', btype]
            print(cmd)
            check_call(cmd)

        cmd = ['{}/bin/fslmerge'.format(FSLDIR), '-t', '{}_all_{}_firstsegs'.format(oname, btype)]


        for imCORR in imscorr:
            if os.path.exists(imCORR + '.nii.gz'):
                cmd.append(imCORR+'.nii.gz')
                print(cmd)

        print(cmd)
        check_call(cmd)

        cmd = ['{}/bin/fslmerge'.format(FSLDIR), '-t', '{}_all_{}_origsegs'.format(oname, btype)]
        for imname in imsfirst:
            if os.path.exists(imname + '.nii.gz'):
                cmd.append(imname+'.nii.gz')
        print(cmd)
        check_call(cmd)

        cmd =  ['{}/bin/first_mult_bcorr'.format(FSLDIR),'-i', t1, '-u', '{}_all_{}_origsegs'.format(oname, btype),
                '-c', '{}_all_{}_firstsegs'.format(oname, btype), '-o', '{}_all_{}_firstsegs'.format(oname, btype)]
        print('{}/bin/first_mult_bcorr'.format(FSLDIR),'-i', t1, '-u', '{}_all_{}_origsegs'.format(oname, btype),
                '-c', '{}_all_{}_firstsegs'.format(oname, btype), '-o', '{}_all_{}_firstsegs'.format(oname, btype))
        check_call(cmd)


def bias_corr(t1):
    N4 = t1.replace(".nii","_N4.nii")
    if not "N4" in t1:
        if not os.path.exists(N4):
            #print("N4BiasFieldCorrection", "-d", "3", "-i", t1,"-o",N4 )
            cmd = ["N4BiasFieldCorrection", "-d", "3", "-i",t1,"-o",N4]
            Popen(cmd).wait()
    return N4




def convert_dcm2nii(E):
    cmd = ["dcm2nii", E ]
    print(cmd)
    Popen(cmd).wait()

def run_all(t1, sub):
    try:
        N4 = bias_corr(t1)
        odir = '/data/henry12/phenoms/OSLO/first_output/' + sub
        run_first(N4,odir)
        #odir = '/data/henry12/phenoms/OSLO/sienax_output/' + sub
        #run_sienax(N4,odir)
    except:
        pass


"""RSI_1 = base_dir + '/Oslo_MS_RSI_1/'
for annon in os.listdir(RSI_1):
    x = RSI_1 +'/'+ annon
    for E in os.listdir(RSI_1 + '/'+ annon):
        o = RSI_1 + '/'+ annon + '/'+ E
        #rename(o)
        for sq in os.listdir(o):
            sq_path = o +'/'+ sq
            #convert_dcm2nii(sq_path)
        if os.path.exists(o + '/FSPGR_SAG_TI450_3/'):
            for t1 in os.listdir(o + '/FSPGR_SAG_TI450_3/'):
                if t1.startswith("o"):
                    t1 = o + '/FSPGR_SAG_TI450_3/'+ t1
                    #print(t1)
                    sub = '/Oslo_MS_RSI_1/'+ t1.split('/')[8]

                    run_all(t1, sub)"""


"""print("^^^^^^^^^^^^")
RSI_2 =  base_dir + '/Oslo_MS_RSI_2/'
for mr in os.listdir(RSI_2):

    x = RSI_2 +'/'+ mr
    print(x)
    if os.path.exists(x + '/Sag_T1_BRAVO_1iso_3/'):
        o = x + '/Sag_T1_BRAVO_1iso_3/'
        for t1 in os.listdir(o ):
            if t1.startswith("o"):
                t1 = o + t1
                print("T1", t1)
                sub = '/Oslo_MS_RSI_2/'+ t1.split('/')[8]
                print("sub",sub)

                run_all(t1, sub)"""


MSX_1 = base_dir + '/Oslo_MSX_1/'
for ms in os.listdir(MSX_1):
    x = MSX_1 +'/'+ ms
    #print(x)
    for t1 in os.listdir(x):
        if "MPRAGE" in t1:
            t1_path = x +'/'+ t1
            print(t1_path)
            #convert_dcm2nii(t1_path)
            t1 = glob(t1_path + '/*.nii.gz')
            if len(t1)>=1:
                t1 = t1[0]
                t1_orient = t1.replace(".nii","_reorient.nii")
                sub = '/Oslo_MSX_1/'+ t1.split('/')[8]

                if not os.path.exists(t1_orient):
                    cmd = ["fslreorient2std",t1, t1_orient]
                    print(cmd)
                    Popen(cmd).wait()
                print("SUB", sub)
                run_all(t1_orient, sub)










