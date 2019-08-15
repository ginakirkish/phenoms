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

base_dir = '/data/henry12/phenoms/Leuven'

def rename(final_path):
    print(final_path, final_path.replace(" ","-"))
    os.rename(final_path, final_path.replace(" ","-"))

def get_brain(t1_path):
    BM = glob(t1_path + 'r*brain*.nii.gz')
    if len(BM) >= 1:
        for b in BM:
            if not "mask" in b:
                BM = b
    else:
        BM = ""
    return BM


def run_sienax(t1,sub):
    #sub = t1.split('_')[4]
    t1_out =  base_dir+ '/sienax_output/' + sub
    if not os.path.exists(t1_out + '/I_brain.nii.gz'):
        if os.path.exists(t1):
            cmd = ["sienax_optibet", t1, "-r", "-d", "-o", t1_out ]
            print("sienax_optibet", t1, "-r", "-d", "-o", t1_out )
        """try:
            Popen(cmd).wait()
        except:
            pass"""
    else:
        print(t1_out, "SIENAX FILE EXISTS")


def run_first(t1,sub):

    #sub = t1.split('_')[3:4]
    odir =  base_dir + '/first_output/' + sub
    print("output", odir)
    if not os.path.exists(odir): os.mkdir(odir)
    #if len(glob('{}/*first*'.format(odir))) == 0:
    x = 0
    if x == 0:
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

def run_reorient(t1):
    cmd = ["fslreorient2std", t1, t1.replace(".nii","reorient.nii")]
    print(cmd)
    Popen(cmd).wait()
    t1 = t1.replace(".nii","reorient.nii")
    return t1

def run_all(t1,sub):
    reorient = run_reorient(t1)
    N4 = bias_corr(reorient)
    #run_sienax(N4, sub)
    try:
        run_first(N4, sub)
    except:
        pass


def convert_dcm2nii(E):
    cmd = ["dcm2nii", E]
    Popen(cmd).wait()

subjects = glob('{}/files*/*3DT1*'.format(base_dir))
for t1 in subjects:
    #print(s)
    #for files in os.listdir(s):
    if t1.endswith(".nii"):
        #print(files)

        sub = t1.split('/')[5].split('_')[4]
        print(t1, sub)
        run_all(t1, sub)